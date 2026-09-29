#include <ParseLQCData/MesonPipeline.h>
#include <MesonAnalysis/MesonExtractor.h>
#include <Statistics/Resampling.h>
#include <IOdata/DirectoryScanner.h>
#include <IOdata/FileReader.h>

#include <algorithm>
#include <charconv>
#include <future>
#include <map>
#include <regex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace lqcd::meson {

namespace {

struct ChannelMeta {
    std::string direction;
    std::string start_tag;
};

[[nodiscard]] size_t infer_num_lines_from_directory(const std::filesystem::path& input_dir) {
    static const std::regex lattice_size_pattern(R"((?:^|[^[:alnum:]])L?([0-9]+)(?:x|X|T|t)([0-9]+))");
    for (auto path = input_dir; !path.empty(); path = path.parent_path()) {
        const auto name = path.filename().string();
        std::smatch match;
        if (std::regex_search(name, match, lattice_size_pattern)) {
            size_t spatial_extent = 0;
            const auto value = match[1].str();
            const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), spatial_extent);
            if (error == std::errc{} && end == value.data() + value.size() && spatial_extent > 0) {
                return spatial_extent;
            }
        }
        const auto parent = path.parent_path();
        if (parent == path) break;
    }
    throw std::invalid_argument("Cannot infer spatial extent (Ns) from input directory: " + input_dir.string());
}

[[nodiscard]] MesonAnalysisResult summarize_channel_matrix(
    const std::vector<double>& combined_matrix,
    const size_t num_lines,
    const size_t num_cfgs,
    const size_t binsize) {
    MesonAnalysisResult result;
    if (combined_matrix.empty() || num_cfgs == 0 || num_lines == 0) {
        return result;
    }

    const auto binned = lqcd::stats::compute_block_binning(combined_matrix, num_lines, num_cfgs, binsize);
    const size_t new_cols = num_cfgs - (num_cfgs % binsize);
    const size_t num_bins = (binsize <= 1) ? num_cfgs : (new_cols / binsize);
    result.n_bins = num_bins;

    const auto jk = lqcd::stats::compute_jackknife(binned, num_lines, num_bins);
    result.folded_jk_matrix = lqcd::stats::compute_folding(jk, num_lines, num_bins);

    result.means.resize(num_lines, 0.0);
    result.errors.resize(num_lines, 0.0);
    lqcd::stats::compute_row_mean_and_jackknife_error(
        result.folded_jk_matrix, num_lines, num_bins,
        result.means, result.errors
    );

    return result;
}

[[nodiscard]] std::vector<ChannelMeta> build_channel_metas(
    std::span<const MesonChannelMapping> channels,
    bool is_single_source) {
    std::vector<ChannelMeta> metas;
    metas.reserve(channels.size());
    for (const auto& ch : channels) {
        if (is_single_source) {
            metas.push_back({
                ch.dir,
                "--- " + ch.type + " to " + ch.type + " --- spatial:" + ch.dir + " --- 0/0/0/0"
            });
        } else {
            metas.push_back({
                ch.dir,
                "--- " + ch.type + " to " + ch.type + " --- spatial:" + ch.dir + " ---"
            });
        }
    }
    return metas;
}

} // namespace

[[nodiscard]] MesonAnalysisResult run_meson_pipeline(
    const std::filesystem::path& input_dir,
    std::span<const MesonChannelMapping> channels,
    const size_t binsize,
    const size_t num_lines,
    size_t thread_count,
    const bool is_single_source) {

    MesonAnalysisResult result;
    if (channels.empty()) {
        return result;
    }

    const size_t resolved_num_lines = num_lines == 0 ? infer_num_lines_from_directory(input_dir) : num_lines;

    // 1. 自然排序扫描所有强子关联函数文件
    const auto files = iodata::scan_natural_sorted(input_dir, "test1_lhadrons_", "_mesons_multi_src");
    const size_t num_cfgs = files.size();
    if (num_cfgs == 0) {
        return result;
    }
    result.n_raw_cfgs = num_cfgs;

    if (thread_count == 0) {
        thread_count = std::max(1u, std::thread::hardware_concurrency());
    }

    // 矩阵数据：num_lines 行 x N 列 (num_lines x num_cfgs)
    std::vector<double> combined_matrix(resolved_num_lines * num_cfgs, 0.0);

    const auto metas = build_channel_metas(channels, is_single_source);

    // 2. 多线程并发读取并处理各个配置文件的信道
    const size_t chunk_size = (num_cfgs + thread_count - 1) / thread_count;
    std::vector<std::future<void>> futures;
    futures.reserve(thread_count);

    for (size_t t = 0; t < thread_count; ++t) {
        const size_t start_idx = t * chunk_size;
        const size_t end_idx = std::min(start_idx + chunk_size, num_cfgs);
        if (start_idx >= end_idx) continue;

        futures.push_back(std::async(std::launch::async, [&, start_idx, end_idx]() {
            for (size_t cfg_i = start_idx; cfg_i < end_idx; ++cfg_i) {
                // 每个文件只需单次读取到内存
                const std::string content = iodata::read_file_to_string(files[cfg_i]);
                if (content.empty()) continue;

                // 累加该文件所有信道的结果
                for (const auto& meta : metas) {
                    std::vector<double> blk;
                    if (is_single_source) {
                        blk = extract_single_file_singlesrc_block(content, meta.start_tag, resolved_num_lines);
                    } else {
                        blk = extract_single_file_averaged_block(content, meta.direction, meta.start_tag, resolved_num_lines);
                    }
                    if (blk.size() == resolved_num_lines) {
                        for (size_t r = 0; r < resolved_num_lines; ++r) {
                            combined_matrix[r * num_cfgs + cfg_i] += blk[r];
                        }
                    }
                }
            }
        }));
    }

    for (auto& f : futures) {
        f.get();
    }

    // 各信道平均
    const double inv_num_channels = 1.0 / static_cast<double>(channels.size());
    for (size_t i = 0; i < combined_matrix.size(); ++i) {
        combined_matrix[i] *= inv_num_channels;
    }

    return summarize_channel_matrix(combined_matrix, resolved_num_lines, num_cfgs, binsize);
}

[[nodiscard]] std::map<std::string, MesonAnalysisResult> run_meson_pipeline_batch(
    const std::filesystem::path& input_dir,
    const std::map<std::string, std::vector<MesonChannelMapping>>& channel_defs,
    const size_t binsize,
    const size_t num_lines,
    size_t thread_count,
    const bool is_single_source) {

    std::map<std::string, MesonAnalysisResult> results;
    if (channel_defs.empty()) {
        return results;
    }

    const size_t resolved_num_lines = num_lines == 0 ? infer_num_lines_from_directory(input_dir) : num_lines;

    const auto files = iodata::scan_natural_sorted(input_dir, "test1_lhadrons_", "_mesons_multi_src");
    const size_t num_cfgs = files.size();
    if (num_cfgs == 0) {
        return results;
    }

    if (thread_count == 0) {
        thread_count = std::max(1u, std::thread::hardware_concurrency());
    }

    std::map<std::string, std::vector<double>> matrices;
    for (const auto& [ch, mappings] : channel_defs) {
        if (!mappings.empty()) {
            matrices[ch] = std::vector<double>(resolved_num_lines * num_cfgs, 0.0);
        }
    }

    const size_t chunk_size = (num_cfgs + thread_count - 1) / thread_count;
    std::vector<std::future<void>> futures;
    futures.reserve(thread_count);

    for (size_t t = 0; t < thread_count; ++t) {
        const size_t start_idx = t * chunk_size;
        const size_t end_idx = std::min(start_idx + chunk_size, num_cfgs);
        if (start_idx >= end_idx) continue;

        futures.push_back(std::async(std::launch::async, [&, start_idx, end_idx]() {
            for (size_t cfg_i = start_idx; cfg_i < end_idx; ++cfg_i) {
                const std::string content = iodata::read_file_to_string(files[cfg_i]);
                if (content.empty()) continue;

                for (const auto& [ch, mappings] : channel_defs) {
                    auto it = matrices.find(ch);
                    if (it == matrices.end() || mappings.empty()) continue;

                    const auto metas = build_channel_metas(mappings, is_single_source);
                    for (const auto& meta : metas) {
                        std::vector<double> blk;
                        if (is_single_source) {
                            blk = extract_single_file_singlesrc_block(content, meta.start_tag, resolved_num_lines);
                        } else {
                            blk = extract_single_file_averaged_block(content, meta.direction, meta.start_tag, resolved_num_lines);
                        }
                        if (blk.size() != resolved_num_lines) continue;
                        for (size_t r = 0; r < resolved_num_lines; ++r) {
                            it->second[r * num_cfgs + cfg_i] += blk[r];
                        }
                    }
                }
            }
        }));
    }

    for (auto& f : futures) {
        f.get();
    }

    for (const auto& [ch, mappings] : channel_defs) {
        auto it = matrices.find(ch);
        if (it == matrices.end()) continue;
        auto& matrix = it->second;
        const double inv_num_channels = 1.0 / static_cast<double>(mappings.size());
        for (auto& value : matrix) {
            value *= inv_num_channels;
        }
        results[ch] = summarize_channel_matrix(matrix, resolved_num_lines, num_cfgs, binsize);
        results[ch].n_raw_cfgs = num_cfgs;
    }

    return results;
}

} // namespace lqcd::meson

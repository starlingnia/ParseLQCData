#include <ParseLQCData/MesonAnalysis/MesonPipeline.h>
#include <ParseLQCData/MesonAnalysis/MesonExtractor.h>
#include <ParseLQCData/Statistics/Resampling.h>
#include <IOdata/DirectoryScanner.h>
#include <IOdata/FileReader.h>

#include <algorithm>
#include <future>
#include <thread>
#include <vector>

namespace lqcd::meson {

[[nodiscard]] MesonAnalysisResult run_meson_pipeline(
    const std::filesystem::path& input_dir,
    std::span<const MesonChannelMapping> channels,
    const size_t binsize,
    const size_t num_lines,
    size_t thread_count) {

    MesonAnalysisResult result;
    if (channels.empty()) {
        return result;
    }

    // 1. 自然排序扫描所有多源强子关联函数文件
    const auto files = iodata::scan_natural_sorted(input_dir, "test1_lhadrons_", "_mesons_multi_src");
    const size_t num_cfgs = files.size();
    if (num_cfgs == 0) {
        return result;
    }
    result.n_raw_cfgs = num_cfgs;

    if (thread_count == 0) {
        thread_count = std::max(1u, std::thread::hardware_concurrency());
    }

    // 矩阵数据：48 行 x N 列 (num_lines x num_cfgs)
    std::vector<double> combined_matrix(num_lines * num_cfgs, 0.0);

    // 预先构造各 channel 对应的 start_tag 与 direction
    struct ChannelMeta {
        std::string direction;
        std::string start_tag;
    };
    std::vector<ChannelMeta> metas;
    metas.reserve(channels.size());
    for (const auto& ch : channels) {
        metas.push_back({
            ch.dir,
            "--- " + ch.type + " to " + ch.type + " --- spatial:" + ch.dir + " ---"
        });
    }

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
                    const auto blk = extract_single_file_averaged_block(content, meta.direction, meta.start_tag, num_lines);
                    for (size_t r = 0; r < num_lines; ++r) {
                        combined_matrix[r * num_cfgs + cfg_i] += blk[r];
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

    // 3. Binning 分块平滑
    const auto binned = lqcd::stats::compute_block_binning(combined_matrix, num_lines, num_cfgs, binsize);
    const size_t new_cols = num_cfgs - (num_cfgs % binsize);
    const size_t num_bins = (binsize <= 1) ? num_cfgs : (new_cols / binsize);
    result.n_bins = num_bins;

    // 4. Jackknife 留一重采样
    const auto jk = lqcd::stats::compute_jackknife(binned, num_lines, num_bins);

    // 5. 时间反演对称折叠
    result.folded_jk_matrix = lqcd::stats::compute_folding(jk, num_lines, num_bins);

    // 6. 最终均值与 Jackknife 统计误差计算
    result.means.resize(num_lines, 0.0);
    result.errors.resize(num_lines, 0.0);
    lqcd::stats::compute_row_mean_and_jackknife_error(
        result.folded_jk_matrix, num_lines, num_bins,
        result.means, result.errors
    );

    return result;
}

} // namespace lqcd::meson

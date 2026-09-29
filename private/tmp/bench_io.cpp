// 临时性能基准：量化 ParseLQCData C++ IO / 解析路径的真实瓶颈
// 用法: bench_io <input_dir> [threads]
#include <IOdata/DirectoryScanner.h>
#include <IOdata/FileReader.h>
#include <IOdata/FastParser.h>
#include <MesonAnalysis/MesonExtractor.h>
#include <ParseLQCData/MesonPipeline.h>
#include <Statistics/Resampling.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <future>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using Clock = std::chrono::steady_clock;
static double secs(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
}

using Meta = std::pair<std::string, std::string>; // (type, dir)
static const std::vector<std::pair<std::string, std::vector<Meta>>> CHANNELS = {
    {"AV",  {{"AVector2","DIRX"},{"AVector3","DIRX"},{"AVector3","DIRY"},{"AVector1","DIRY"},{"AVector1","DIRZ"},{"AVector2","DIRZ"}}},
    {"S",   {{"TVector4","DIRZ"},{"TVector4","DIRX"},{"TVector4","DIRY"}}},
    {"Tt",  {{"TVector1","DIRX"},{"TVector3","DIRZ"},{"TVector2","DIRY"}}},
    {"PS",  {{"TAVector4","DIRZ"},{"TAVector4","DIRX"},{"TAVector4","DIRY"}}},
    {"Xt",  {{"TAVector1","DIRX"},{"TAVector3","DIRZ"},{"TAVector2","DIRY"}}},
    {"Vec", {{"Vector2","DIRX"},{"Vector3","DIRX"},{"Vector3","DIRY"},{"Vector1","DIRY"},{"Vector1","DIRZ"},{"Vector2","DIRZ"}}},
};

// ---------- 单遍解析器：一次线性扫描同时抽取全部信道 ----------
struct FlatMeta {
    size_t channel;
    std::string type;
    std::string dir;
    int axis; // 0=x 1=y 2=z
};

static int dir_axis(const std::string& d) {
    if (d == "DIRX") return 0;
    if (d == "DIRY") return 1;
    return 2;
}

// 返回 [channel][meta] 的 48 维移位平均向量
static std::vector<std::vector<double>> parse_all_metas_one_pass(
    std::string_view content,
    const std::vector<FlatMeta>& metas,
    size_t num_lines) {

    std::vector<std::vector<double>> sums(metas.size(), std::vector<double>(num_lines, 0.0));
    std::vector<size_t> counts(metas.size(), 0);

    std::vector<double> blk(num_lines, 0.0);
    long cur = -1;
    int cur_shift = 0;
    size_t cur_fill = 0;

    iodata::for_each_line(content, [&](std::string_view line) {
        if (!line.empty() && line[0] == '-') {
            cur = -1;
            cur_fill = 0;
            if (line.size() < 5) return;
            const size_t sp = line.find(' ', 4);
            if (sp == std::string_view::npos) return;
            const std::string_view type = line.substr(4, sp - 4);
            const size_t spos = line.find("spatial:");
            if (spos == std::string_view::npos) return;
            const std::string_view dir = line.substr(spos + 8, 4);
            for (size_t m = 0; m < metas.size(); ++m) {
                if (metas[m].type.size() == type.size() && metas[m].dir.size() == dir.size() &&
                    std::string_view(metas[m].type) == type && std::string_view(metas[m].dir) == dir) {
                    cur = static_cast<long>(m);
                    break;
                }
            }
            if (cur < 0) return;
            // 解析末尾 x/y/z/t
            const size_t last = line.rfind("--- ");
            if (last == std::string_view::npos) { cur = -1; return; }
            std::string_view coord = line.substr(last + 4);
            int nums[4] = {0,0,0,0};
            size_t idx = 0, start = 0;
            while (idx < 4 && start < coord.size()) {
                const size_t slash = coord.find('/', start);
                const std::string_view tok = (slash == std::string_view::npos) ? coord.substr(start)
                                                                              : coord.substr(start, slash - start);
                int v = 0;
                if (iodata::parse_int(tok, v)) nums[idx++] = v; else break;
                if (slash == std::string_view::npos) break;
                start = slash + 1;
            }
            cur_shift = nums[metas[static_cast<size_t>(cur)].axis];
            return;
        }
        if (cur < 0 || cur_fill >= num_lines) return;

        size_t token_idx = 0;
        double real_val = 0.0;
        bool got = false;
        iodata::for_each_token(line, [&](std::string_view token) {
            if (token_idx == 1) { if (iodata::parse_double(token, real_val)) got = true; }
            ++token_idx;
        });
        if (!got) return;
        blk[cur_fill++] = real_val;
        if (cur_fill == num_lines) {
            auto& s = sums[static_cast<size_t>(cur)];
            const int n = static_cast<int>(num_lines);
            const int sh = ((cur_shift % n) + n) % n;
            for (int i = 0; i < n; ++i) {
                s[static_cast<size_t>(i)] += blk[static_cast<size_t>((i + sh) % n)];
            }
            ++counts[static_cast<size_t>(cur)];
            cur = -1;
            cur_fill = 0;
        }
    });

    for (size_t m = 0; m < metas.size(); ++m) {
        if (counts[m] == 0) continue;
        const double inv = 1.0 / static_cast<double>(counts[m]);
        for (auto& v : sums[m]) v *= inv;
    }
    return sums;
}

static void run_stats(const std::vector<double>& matrix, size_t num_lines, size_t num_cfgs,
                      size_t binsize, std::vector<double>& means, std::vector<double>& errors) {
    const auto binned = lqcd::stats::compute_block_binning(matrix, num_lines, num_cfgs, binsize);
    const size_t new_cols = num_cfgs - (num_cfgs % binsize);
    const size_t num_bins = (binsize <= 1) ? num_cfgs : (new_cols / binsize);
    const auto jk = lqcd::stats::compute_jackknife(binned, num_lines, num_bins);
    const auto folded = lqcd::stats::compute_folding(jk, num_lines, num_bins);
    means.assign(num_lines, 0.0);
    errors.assign(num_lines, 0.0);
    lqcd::stats::compute_row_mean_and_jackknife_error(folded, num_lines, num_bins, means, errors);
}

int main(int argc, char** argv) {
    if (argc < 2) { std::printf("usage: %s <input_dir> [threads]\n", argv[0]); return 1; }
    const std::filesystem::path dir = argv[1];
    const size_t T = (argc > 2) ? static_cast<size_t>(std::stoul(argv[2]))
                                : std::max(1u, std::thread::hardware_concurrency());
    const size_t num_lines = 48;
    const size_t binsize = 4;

    auto t0 = Clock::now();
    const auto files = iodata::scan_natural_sorted(dir, "test1_lhadrons_", "_mesons_multi_src");
    auto t1 = Clock::now();
    size_t total_bytes = 0;
    for (const auto& f : files) total_bytes += static_cast<size_t>(std::filesystem::file_size(f));
    std::printf("[scan ] files=%zu  bytes=%.1f MB  time=%.4f s\n",
                files.size(), total_bytes / 1e6, secs(t0, t1));
    if (files.empty()) return 1;
    const size_t num_cfgs = files.size();

    // ---------- 1. 纯同步读取: 单线程 vs T 线程 ----------
    {
        auto a = Clock::now();
        size_t bytes = 0;
        std::string buf;
        for (const auto& f : files) { buf = iodata::read_file_to_string(f); bytes += buf.size(); }
        auto b = Clock::now();
        const double s = secs(a, b);
        std::printf("[read ] 1 thread : %.3f s  (%.0f MB/s)\n", s, bytes / 1e6 / s);
    }
    {
        auto a = Clock::now();
        const size_t chunk = (num_cfgs + T - 1) / T;
        std::vector<std::future<size_t>> fut;
        for (size_t t = 0; t < T; ++t) {
            const size_t s0 = t * chunk, e0 = std::min(s0 + chunk, num_cfgs);
            if (s0 >= e0) continue;
            fut.push_back(std::async(std::launch::async, [&, s0, e0]() {
                size_t bytes = 0;
                std::string buf;
                for (size_t i = s0; i < e0; ++i) { buf = iodata::read_file_to_string(files[i]); bytes += buf.size(); }
                return bytes;
            }));
        }
        size_t bytes = 0;
        for (auto& f : fut) bytes += f.get();
        auto b = Clock::now();
        const double s = secs(a, b);
        std::printf("[read ] %zu thread: %.3f s  (%.0f MB/s)\n", T, s, bytes / 1e6 / s);
    }

    // ---------- 2. 现有逐信道流水线 (Python 实际调用方式) ----------
    double sum_channel_time = 0.0;
    std::vector<std::vector<double>> ref_means, ref_errors;
    for (const auto& [ch, metas] : CHANNELS) {
        std::vector<lqcd::MesonChannelMapping> chs;
        for (const auto& m : metas) chs.push_back({m.first, m.second});
        auto a = Clock::now();
        auto res = lqcd::meson::run_meson_pipeline(dir, chs, binsize, num_lines, T, false);
        auto b = Clock::now();
        const double s = secs(a, b);
        sum_channel_time += s;
        ref_means.push_back(res.means);
        ref_errors.push_back(res.errors);
        std::printf("[pipe ] channel=%-4s metas=%zu  %zu threads : %.3f s\n", ch.c_str(), metas.size(), T, s);
    }
    std::printf("[pipe ] SUM of 6 channel calls        : %.3f s\n", sum_channel_time);

    // ---------- 2b. 单源模式逐信道流水线 (读取同样的 multi_src 大文件) ----------
    double sum_single_time = 0.0;
    for (const auto& [ch, metas] : CHANNELS) {
        std::vector<lqcd::MesonChannelMapping> chs;
        for (const auto& m : metas) chs.push_back({m.first, m.second});
        auto a = Clock::now();
        (void)lqcd::meson::run_meson_pipeline(dir, chs, binsize, num_lines, T, true);
        auto b = Clock::now();
        const double s = secs(a, b);
        sum_single_time += s;
        std::printf("[sing ] channel=%-4s metas=%zu  %zu threads : %.3f s\n", ch.c_str(), metas.size(), T, s);
    }
    std::printf("[sing ] SUM of 6 channel calls        : %.3f s\n", sum_single_time);

    // 单线程对照 (仅 AV)
    {
        std::vector<lqcd::MesonChannelMapping> chs;
        for (const auto& m : CHANNELS[0].second) chs.push_back({m.first, m.second});
        auto a = Clock::now();
        (void)lqcd::meson::run_meson_pipeline(dir, chs, binsize, num_lines, 1, false);
        auto b = Clock::now();
        std::printf("[pipe ] channel=AV   1 thread          : %.3f s\n", secs(a, b));
    }

    // ---------- 3. 单遍方案 A: 读一次 + 每信道各自扫描 (18 个 meta) ----------
    {
        std::vector<FlatMeta> metas;
        for (size_t c = 0; c < CHANNELS.size(); ++c)
            for (const auto& m : CHANNELS[c].second) metas.push_back({c, m.first, m.second, dir_axis(m.second)});
        const size_t num_ch = CHANNELS.size();
        std::vector<std::vector<double>> mat(num_ch, std::vector<double>(num_lines * num_cfgs, 0.0));
        auto a = Clock::now();
        const size_t chunk = (num_cfgs + T - 1) / T;
        std::vector<std::future<void>> fut;
        for (size_t t = 0; t < T; ++t) {
            const size_t s0 = t * chunk, e0 = std::min(s0 + chunk, num_cfgs);
            if (s0 >= e0) continue;
            fut.push_back(std::async(std::launch::async, [&, s0, e0]() {
                for (size_t ci = s0; ci < e0; ++ci) {
                    const std::string content = iodata::read_file_to_string(files[ci]);
                    if (content.empty()) continue;
                    for (size_t c = 0; c < num_ch; ++c) {
                        for (const auto& m : CHANNELS[c].second) {
                            const std::string tag = "--- " + m.first + " to " + m.first + " --- spatial:" + m.second + " ---";
                            const auto blk = lqcd::meson::extract_single_file_averaged_block(content, m.second, tag, num_lines);
                            auto& col = mat[c];
                            for (size_t r = 0; r < num_lines; ++r) col[r * num_cfgs + ci] += blk[r];
                        }
                    }
                }
            }));
        }
        for (auto& f : fut) f.get();
        auto b = Clock::now();
        std::printf("[onepass-A] 全部 6 信道, 每文件读 1 次     : %.3f s\n", secs(a, b));
    }

    // ---------- 4. 单遍方案 B: 单次线性扫描解析全部信道 + 数值一致性校验 ----------
    {
        std::vector<FlatMeta> metas;
        for (size_t c = 0; c < CHANNELS.size(); ++c)
            for (const auto& m : CHANNELS[c].second) metas.push_back({c, m.first, m.second, dir_axis(m.second)});
        const size_t num_ch = CHANNELS.size();
        std::vector<std::vector<double>> mat(num_ch, std::vector<double>(num_lines * num_cfgs, 0.0));
        auto a = Clock::now();
        const size_t chunk = (num_cfgs + T - 1) / T;
        std::vector<std::future<void>> fut;
        for (size_t t = 0; t < T; ++t) {
            const size_t s0 = t * chunk, e0 = std::min(s0 + chunk, num_cfgs);
            if (s0 >= e0) continue;
            fut.push_back(std::async(std::launch::async, [&, s0, e0]() {
                for (size_t ci = s0; ci < e0; ++ci) {
                    const std::string content = iodata::read_file_to_string(files[ci]);
                    if (content.empty()) continue;
                    const auto per_meta = parse_all_metas_one_pass(content, metas, num_lines);
                    for (size_t m = 0; m < metas.size(); ++m) {
                        auto& col = mat[metas[m].channel];
                        const auto& blk = per_meta[m];
                        for (size_t r = 0; r < num_lines; ++r) col[r * num_cfgs + ci] += blk[r];
                    }
                }
            }));
        }
        for (auto& f : fut) f.get();
        auto b = Clock::now();
        const double s = secs(a, b);
        std::printf("[onepass-B] 全部 6 信道, 单次线性扫描      : %.3f s\n", s);

        // 数值一致性 (注意: 流水线内部对同信道各 meta 求和后除以 meta 数)
        double max_dmean = 0.0, max_derr = 0.0;
        for (size_t c = 0; c < num_ch; ++c) {
            const double inv = 1.0 / static_cast<double>(CHANNELS[c].second.size());
            for (auto& v : mat[c]) v *= inv;
            std::vector<double> means, errors;
            run_stats(mat[c], num_lines, num_cfgs, binsize, means, errors);
            for (size_t i = 0; i < num_lines; ++i) {
                max_dmean = std::max(max_dmean, std::fabs(means[i] - ref_means[c][i]));
                max_derr = std::max(max_derr, std::fabs(errors[i] - ref_errors[c][i]));
            }
        }
        std::printf("[check ] onepass-B vs 现有流水线: max|dmean|=%.3e max|derr|=%.3e\n", max_dmean, max_derr);
    }

    // ---------- 5. 解析开销随 meta 数线性增长验证 ----------
    {
        const std::string content = iodata::read_file_to_string(files[0]);
        for (auto metas_count : {1u, 3u, 6u, 18u}) {
            std::vector<FlatMeta> metas;
            size_t added = 0;
            for (size_t c = 0; c < CHANNELS.size() && added < metas_count; ++c)
                for (const auto& m : CHANNELS[c].second) {
                    if (added >= metas_count) break;
                    metas.push_back({c, m.first, m.second, dir_axis(m.second)});
                    ++added;
                }
            auto a = Clock::now();
            const int reps = 10;
            for (int r = 0; r < reps; ++r) (void)parse_all_metas_one_pass(content, metas, num_lines);
            auto b = Clock::now();
            std::printf("[parse] 单文件 x %2zu meta = %.2f ms/pass (%.1f MB/s) file=%.2f MB\n",
                        metas.size(), secs(a, b) * 1000.0 / reps,
                        content.size() / 1e6 / (secs(a, b) / reps), content.size() / 1e6);
        }
    }
    return 0;
}

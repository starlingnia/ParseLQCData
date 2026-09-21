#include <CondensateAnalysis/CondensateExtractor.h>
#include <IOdata/FastParser.h>
#include <IOdata/FileReader.h>
#include <IOdata/DirectoryScanner.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <future>
#include <thread>

namespace lqcd::condensate {

[[nodiscard]] bool extract_pbp_from_xml(std::string_view xml_content, double& out_val) noexcept {
    const size_t tag_pos = xml_content.find("<pbp>(");
    if (tag_pos == std::string_view::npos) {
        return false;
    }
    const size_t num_start = tag_pos + 6; // len("<pbp>(") == 6
    const size_t comma_pos = xml_content.find(',', num_start);
    if (comma_pos == std::string_view::npos) {
        return false;
    }
    const std::string_view num_str = xml_content.substr(num_start, comma_pos - num_start);
    return iodata::parse_double(num_str, out_val);
}

[[nodiscard]] std::pair<double, double> extract_meas_pair(const std::filesystem::path& meas_dir) {
    const auto psibar_dir = meas_dir / "PsibarPsi";
    if (!std::filesystem::exists(psibar_dir)) {
        return {0.0, 0.0};
    }

    double sum_l = 0.0;
    size_t cnt_l = 0;
    double sum_s = 0.0;
    size_t cnt_s = 0;

    for (const auto& entry : std::filesystem::directory_iterator(psibar_dir)) {
        if (!entry.is_regular_file()) continue;
        const std::string fname = entry.path().filename().string();
        if (fname.ends_with(".xml")) {
            const std::string content = iodata::read_file_to_string(entry.path());
            double val = 0.0;
            if (extract_pbp_from_xml(content, val)) {
                if (fname.find("_pbp_l_") != std::string::npos) {
                    sum_l += val;
                    ++cnt_l;
                } else if (fname.find("_pbp_s_") != std::string::npos) {
                    sum_s += val;
                    ++cnt_s;
                }
            }
        }
    }

    const double avg_l = (cnt_l > 0) ? (sum_l / static_cast<double>(cnt_l)) : 0.0;
    const double avg_s = (cnt_s > 0) ? (sum_s / static_cast<double>(cnt_s)) : 0.0;

    return {avg_l, avg_s};
}

[[nodiscard]] CondensateResult process_chiral_condensate(
    const std::filesystem::path& base_dir,
    const double m_light,
    const double m_strange,
    const double m_residual,
    const double zm_factor) {

    CondensateResult res;

    // 扫描全部 meas.* 目录
    std::vector<std::filesystem::path> meas_dirs;
    if (!std::filesystem::exists(base_dir)) {
        return res;
    }

    for (const auto& entry : std::filesystem::directory_iterator(base_dir)) {
        if (entry.is_directory()) {
            const std::string name = entry.path().filename().string();
            if (name.starts_with("meas.")) {
                meas_dirs.push_back(entry.path());
            }
        }
    }

    std::sort(meas_dirs.begin(), meas_dirs.end(), [](const auto& a, const auto& b) {
        return iodata::extract_numeric_id(a.filename().string()) < iodata::extract_numeric_id(b.filename().string());
    });

    const size_t num_cfgs = meas_dirs.size();
    if (num_cfgs == 0) {
        return res;
    }
    res.num_cfgs = num_cfgs;

    std::vector<double> vals_l(num_cfgs, 0.0);
    std::vector<double> vals_s(num_cfgs, 0.0);

    const size_t thread_count = std::max(1u, std::thread::hardware_concurrency());
    const size_t chunk_size = (num_cfgs + thread_count - 1) / thread_count;
    std::vector<std::future<void>> futures;

    for (size_t t = 0; t < thread_count; ++t) {
        const size_t start_i = t * chunk_size;
        const size_t end_i = std::min(start_i + chunk_size, num_cfgs);
        if (start_i >= end_i) continue;

        futures.push_back(std::async(std::launch::async, [&, start_i, end_i]() {
            for (size_t i = start_i; i < end_i; ++i) {
                const auto [l, s] = extract_meas_pair(meas_dirs[i]);
                vals_l[i] = l;
                vals_s[i] = s;
            }
        }));
    }

    for (auto& f : futures) {
        f.get();
    }

    // Jackknife 留一重采样
    const double sum_l = std::accumulate(vals_l.begin(), vals_l.end(), 0.0);
    const double sum_s = std::accumulate(vals_s.begin(), vals_s.end(), 0.0);
    const double inv_n_minus_one = 1.0 / static_cast<double>(num_cfgs - 1);

    const double ml_eff = m_light + m_residual;
    const double ms_eff = m_strange + m_residual;
    const double mass_ratio = (ms_eff > 1e-15) ? (ml_eff / ms_eff) : 0.0;
    const double inv_zm = (std::abs(zm_factor) > 1e-15) ? (1.0 / zm_factor) : 1.0;

    res.jackknife_samples.resize(num_cfgs, 0.0);
    double jk_sum = 0.0;
    double jk_l_sum = 0.0;
    double jk_s_sum = 0.0;
    std::vector<double> jk_l_samples(num_cfgs, 0.0);
    std::vector<double> jk_s_samples(num_cfgs, 0.0);

    for (size_t i = 0; i < num_cfgs; ++i) {
        const double jk_l = (sum_l - vals_l[i]) * inv_n_minus_one;
        const double jk_s = (sum_s - vals_s[i]) * inv_n_minus_one;
        const double jk_sub = (jk_l - mass_ratio * jk_s) * inv_zm;
        res.jackknife_samples[i] = jk_sub;
        jk_sum += jk_sub;

        jk_l_samples[i] = jk_l;
        jk_l_sum += jk_l;

        jk_s_samples[i] = jk_s;
        jk_s_sum += jk_s;
    }

    const double mean = jk_sum / static_cast<double>(num_cfgs);
    res.mean = mean;

    double sq_diff_sum = 0.0;
    for (size_t i = 0; i < num_cfgs; ++i) {
        const double diff = res.jackknife_samples[i] - mean;
        sq_diff_sum += diff * diff;
    }

    const double var_ddof0 = sq_diff_sum / static_cast<double>(num_cfgs);
    res.error = std::sqrt(var_ddof0) * std::sqrt(static_cast<double>(num_cfgs - 1));

    // Bare light condensate mean & Jackknife error
    const double mean_l = jk_l_sum / static_cast<double>(num_cfgs);
    res.pbp_l_mean = mean_l;
    double sq_diff_l = 0.0;
    for (size_t i = 0; i < num_cfgs; ++i) {
        const double diff_l = jk_l_samples[i] - mean_l;
        sq_diff_l += diff_l * diff_l;
    }
    res.pbp_l_error = std::sqrt(sq_diff_l / static_cast<double>(num_cfgs)) * std::sqrt(static_cast<double>(num_cfgs - 1));

    // Bare strange condensate mean & Jackknife error
    const double mean_s = jk_s_sum / static_cast<double>(num_cfgs);
    res.pbp_s_mean = mean_s;
    double sq_diff_s = 0.0;
    for (size_t i = 0; i < num_cfgs; ++i) {
        const double diff_s = jk_s_samples[i] - mean_s;
        sq_diff_s += diff_s * diff_s;
    }
    res.pbp_s_error = std::sqrt(sq_diff_s / static_cast<double>(num_cfgs)) * std::sqrt(static_cast<double>(num_cfgs - 1));

    return res;
}

} // namespace lqcd::condensate

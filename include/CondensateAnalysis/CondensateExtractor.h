#pragma once

#include <filesystem>
#include <string_view>
#include <vector>
#include <utility>

namespace lqcd::condensate {

struct CondensateResult {
    double mean{0.0};
    double error{0.0};
    double pbp_l_mean{0.0};
    double pbp_l_error{0.0};
    double pbp_s_mean{0.0};
    double pbp_s_error{0.0};
    std::vector<double> jackknife_samples;
    size_t num_cfgs{0};
};

struct SusceptibilityResult {
    double mean_unscaled{0.0};
    double error_unscaled{0.0};
    double factor_vol{0.0};
    double mean_vol_scaled{0.0};
    double error_vol_scaled{0.0};
    double factor_scaled{0.0};
    double mean_scaled{0.0};
    double error_scaled{0.0};
    size_t num_cfgs{0};
    int ns{48};
    int nt{16};
    double temp{157.0};
};

// GSL 规范：小函数职责单一，[[nodiscard]] 强制调用方检查结果 (独立功能算子)
// 1. 从 XML 文本内容中直接提取 <pbp>(real, imag)</pbp> 的实数部分
[[nodiscard]] bool extract_pbp_from_xml(std::string_view xml_content, double& out_val) noexcept;

// 2. 从单个 meas 文件夹中并发或极速提取 light 和 strange 夸克的随机源均值
[[nodiscard]] std::pair<double, double> extract_meas_pair(const std::filesystem::path& meas_dir);

// 3. 从单个 meas 文件夹中并发或极速提取纯轻夸克的单构型无偏方均估计 (Obar, O2bar)
[[nodiscard]] std::pair<double, double> extract_meas_unbiased_quadratic(const std::filesystem::path& meas_dir);

// 4. 扫描整个 condensate 目录下的全部 meas.* 文件夹，执行 Jackknife 并完成物理重整化消除发散
[[nodiscard]] CondensateResult process_chiral_condensate(
    const std::filesystem::path& base_dir,
    double m_light,
    double m_strange,
    double m_residual,
    double zm_factor);

// 5. 扫描整个 condensate 目录下的全部 meas.* 文件夹，多线程高速提取轻夸克单构型无偏方均估计并执行 Jackknife
[[nodiscard]] SusceptibilityResult process_chiral_susceptibility(
    const std::filesystem::path& base_dir,
    int ns = 48,
    int nt = 16,
    double temp_mev = 157.0);

} // namespace lqcd::condensate


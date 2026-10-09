#include <CondensateAnalysis/CondensateExtractor.h>

extern "C" {

/**
 * @brief 手征凝聚快速提取与物理重整化 C ABI 导出接口
 */
int run_chiral_condensate_c_api(
    const char* base_dir,
    double m_light,
    double m_strange,
    double m_residual,
    double zm_factor,
    double* out_mean,
    double* out_error,
    int* out_num_cfgs) {

    if (!base_dir || !out_mean || !out_error) {
        return -1;
    }

    const auto res = lqcd::condensate::process_chiral_condensate(
        base_dir, m_light, m_strange, m_residual, zm_factor
    );

    *out_mean = res.mean;
    *out_error = res.error;
    if (out_num_cfgs) {
        *out_num_cfgs = static_cast<int>(res.num_cfgs);
    }
    return 0;
}

/**
 * @brief 手征凝聚全量测量导出接口 (含扣除残余质量与裸光/奇夸克手征凝聚)
 */
int run_chiral_condensate_full_c_api(
    const char* base_dir,
    double m_light,
    double m_strange,
    double m_residual,
    double zm_factor,
    double* out_pbp_sub_mean,
    double* out_pbp_sub_error,
    double* out_pbp_l_mean,
    double* out_pbp_l_error,
    double* out_pbp_s_mean,
    double* out_pbp_s_error,
    int* out_num_cfgs) {

    if (!base_dir || !out_pbp_sub_mean || !out_pbp_sub_error) {
        return -1;
    }

    const auto res = lqcd::condensate::process_chiral_condensate(
        base_dir, m_light, m_strange, m_residual, zm_factor
    );

    *out_pbp_sub_mean = res.mean;
    *out_pbp_sub_error = res.error;
    if (out_pbp_l_mean) *out_pbp_l_mean = res.pbp_l_mean;
    if (out_pbp_l_error) *out_pbp_l_error = res.pbp_l_error;
    if (out_pbp_s_mean) *out_pbp_s_mean = res.pbp_s_mean;
    if (out_pbp_s_error) *out_pbp_s_error = res.pbp_s_error;
    if (out_num_cfgs) {
        *out_num_cfgs = static_cast<int>(res.num_cfgs);
    }
    return 0;
}

/**
 * @brief 手征磁化率快速提取与物理标度 C ABI 导出接口
 */
int run_chiral_susceptibility_c_api(
    const char* base_dir,
    int ns,
    int nt,
    double temp_mev,
    double* out_mean_unscaled,
    double* out_error_unscaled,
    double* out_mean_vol_scaled,
    double* out_error_vol_scaled,
    double* out_mean_scaled,
    double* out_error_scaled,
    int* out_num_cfgs) {

    if (!base_dir || !out_mean_unscaled || !out_error_unscaled) {
        return -1;
    }

    const auto res = lqcd::condensate::process_chiral_susceptibility(
        base_dir, ns, nt, temp_mev
    );

    *out_mean_unscaled = res.mean_unscaled;
    *out_error_unscaled = res.error_unscaled;

    if (out_mean_vol_scaled) *out_mean_vol_scaled = res.mean_vol_scaled;
    if (out_error_vol_scaled) *out_error_vol_scaled = res.error_vol_scaled;

    if (out_mean_scaled) *out_mean_scaled = res.mean_scaled;
    if (out_error_scaled) *out_error_scaled = res.error_scaled;

    if (out_num_cfgs) {
        *out_num_cfgs = static_cast<int>(res.num_cfgs);
    }
    return 0;
}

/**
 * @brief 手征磁化率全量测量与样本序列导出接口 (含真实 Jackknife 样本序列及构型级 obar/o2bar)
 */
int run_chiral_susceptibility_full_c_api(
    const char* base_dir,
    int ns,
    int nt,
    double temp_mev,
    double zm_factor,
    double* out_mean_unscaled,
    double* out_error_unscaled,
    double* out_mean_vol_scaled,
    double* out_error_vol_scaled,
    double* out_mean_scaled,
    double* out_error_scaled,
    int* out_num_cfgs,
    int max_cfgs,
    double* out_jk_samples_renorm,
    double* out_obar,
    double* out_o2bar) {

    if (!base_dir || !out_mean_unscaled || !out_error_unscaled) {
        return -1;
    }

    const auto res = lqcd::condensate::process_chiral_susceptibility(
        base_dir, ns, nt, temp_mev
    );

    *out_mean_unscaled = res.mean_unscaled;
    *out_error_unscaled = res.error_unscaled;

    if (out_mean_vol_scaled) *out_mean_vol_scaled = res.mean_vol_scaled;
    if (out_error_vol_scaled) *out_error_vol_scaled = res.error_vol_scaled;

    if (out_mean_scaled) *out_mean_scaled = res.mean_scaled;
    if (out_error_scaled) *out_error_scaled = res.error_scaled;

    if (out_num_cfgs) {
        *out_num_cfgs = static_cast<int>(res.num_cfgs);
    }

    const double inv_zm2 = (std::abs(zm_factor) > 1e-15) ? (1.0 / (zm_factor * zm_factor)) : 1.0;
    const double renorm_factor = (res.factor_scaled / 1e6) * inv_zm2;

    const size_t copy_n = std::min(static_cast<size_t>(max_cfgs), res.num_cfgs);
    if (out_jk_samples_renorm && res.jackknife_samples.size() >= copy_n) {
        for (size_t i = 0; i < copy_n; ++i) {
            out_jk_samples_renorm[i] = res.jackknife_samples[i] * renorm_factor;
        }
    }
    if (out_obar && res.obar_list.size() >= copy_n) {
        for (size_t i = 0; i < copy_n; ++i) {
            out_obar[i] = res.obar_list[i];
        }
    }
    if (out_o2bar && res.o2bar_list.size() >= copy_n) {
        for (size_t i = 0; i < copy_n; ++i) {
            out_o2bar[i] = res.o2bar_list[i];
        }
    }
    return 0;
}

} // extern "C"


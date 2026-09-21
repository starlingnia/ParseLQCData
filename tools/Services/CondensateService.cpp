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

} // extern "C"

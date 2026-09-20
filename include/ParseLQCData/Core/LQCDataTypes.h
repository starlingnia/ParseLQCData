#pragma once

#include <string>
#include <vector>
#include <span>

namespace lqcd {

// 单个物理信道与空间方向的映射配置
struct MesonChannelMapping {
    std::string type;
    std::string dir;
};

// 完整的 Meson 分析统计计算结果
struct MesonAnalysisResult {
    std::vector<double> means;              // 最终对称均值 (长度 48)
    std::vector<double> errors;             // Jackknife 统计误差 (长度 48)
    std::vector<double> folded_jk_matrix;   // 折叠后的 Jackknife 样本矩阵 (48 x n_bins)
    size_t n_bins{0};
    size_t n_raw_cfgs{0};
};

} // namespace lqcd

#pragma once

#include <span>
#include <vector>
#include <cmath>
#include <cstddef>

namespace lqcd::stats {

// GSL 规范：小函数职责单一，[[nodiscard]] 强制检查结果
// 1. Binning 分块平滑：将 (rows, cols) 按 block_size 分组合并为 (rows, cols / block_size)
[[nodiscard]] std::vector<double> compute_block_binning(
    std::span<const double> matrix,
    size_t rows,
    size_t cols,
    size_t block_size);

// 2. Jackknife 留一重采样：计算 J_i = (RowSum - C_i) / (cols - 1)
[[nodiscard]] std::vector<double> compute_jackknife(
    std::span<const double> binned_matrix,
    size_t rows,
    size_t cols);

// 3. 对称折叠 (Time-reversal symmetry folding): C_fold(t) = (C(t) + C(T - t)) / 2
[[nodiscard]] std::vector<double> compute_folding(
    std::span<const double> jk_matrix,
    size_t rows,
    size_t cols);

// 4. 计算行均值与 Jackknife 统计误差：std_ddof0 * sqrt(N - 1)
void compute_row_mean_and_jackknife_error(
    std::span<const double> folded_jk_matrix,
    size_t rows,
    size_t cols,
    std::span<double> out_means,
    std::span<double> out_errors);

} // namespace lqcd::stats

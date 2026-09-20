#include <Statistics/Resampling.h>

#include <numeric>
#include <cmath>

namespace lqcd::stats {

[[nodiscard]] std::vector<double> compute_block_binning(
    std::span<const double> matrix,
    size_t rows,
    size_t cols,
    size_t block_size) {

    if (block_size <= 1) {
        return std::vector<double>(matrix.begin(), matrix.end());
    }

    const size_t new_cols = cols - (cols % block_size);
    const size_t num_bins = new_cols / block_size;
    std::vector<double> binned(rows * num_bins, 0.0);

    const double inv_block = 1.0 / static_cast<double>(block_size);

    for (size_t r = 0; r < rows; ++r) {
        for (size_t b = 0; b < num_bins; ++b) {
            double sum = 0.0;
            const size_t start_col = b * block_size;
            for (size_t k = 0; k < block_size; ++k) {
                sum += matrix[r * cols + start_col + k];
            }
            binned[r * num_bins + b] = sum * inv_block;
        }
    }

    return binned;
}

[[nodiscard]] std::vector<double> compute_jackknife(
    std::span<const double> binned_matrix,
    size_t rows,
    size_t cols) {

    std::vector<double> jk_matrix(rows * cols, 0.0);
    if (cols <= 1) {
        return jk_matrix;
    }

    const double inv_n_minus_one = 1.0 / static_cast<double>(cols - 1);

    for (size_t r = 0; r < rows; ++r) {
        double row_sum = 0.0;
        const size_t row_offset = r * cols;
        for (size_t c = 0; c < cols; ++c) {
            row_sum += binned_matrix[row_offset + c];
        }

        for (size_t c = 0; c < cols; ++c) {
            jk_matrix[row_offset + c] = (row_sum - binned_matrix[row_offset + c]) * inv_n_minus_one;
        }
    }

    return jk_matrix;
}

[[nodiscard]] std::vector<double> compute_folding(
    std::span<const double> jk_matrix,
    size_t rows,
    size_t cols) {

    std::vector<double> folded(rows * cols, 0.0);

    for (size_t r = 0; r < rows; ++r) {
        const size_t sym_r = (r == 0 || r == rows / 2) ? r : (rows - r);
        const size_t out_offset = r * cols;
        const size_t in_offset_r = r * cols;
        const size_t in_offset_sym = sym_r * cols;

        if (r == 0 || r == rows / 2) {
            for (size_t c = 0; c < cols; ++c) {
                folded[out_offset + c] = jk_matrix[in_offset_r + c];
            }
        } else {
            for (size_t c = 0; c < cols; ++c) {
                folded[out_offset + c] = (jk_matrix[in_offset_r + c] + jk_matrix[in_offset_sym + c]) * 0.5;
            }
        }
    }

    return folded;
}

void compute_row_mean_and_jackknife_error(
    std::span<const double> folded_jk_matrix,
    size_t rows,
    size_t cols,
    std::span<double> out_means,
    std::span<double> out_errors) {

    if (cols == 0) return;

    const double inv_cols = 1.0 / static_cast<double>(cols);
    const double sqrt_n_minus_one = std::sqrt(static_cast<double>(cols - 1));

    for (size_t r = 0; r < rows; ++r) {
        const size_t offset = r * cols;
        double sum = 0.0;
        for (size_t c = 0; c < cols; ++c) {
            sum += folded_jk_matrix[offset + c];
        }
        const double mean = sum * inv_cols;
        out_means[r] = mean;

        double sum_sq_diff = 0.0;
        for (size_t c = 0; c < cols; ++c) {
            const double diff = folded_jk_matrix[offset + c] - mean;
            sum_sq_diff += diff * diff;
        }

        const double var_ddof0 = sum_sq_diff * inv_cols;
        const double std_ddof0 = std::sqrt(var_ddof0);
        out_errors[r] = std_ddof0 * sqrt_n_minus_one;
    }
}

} // namespace lqcd::stats

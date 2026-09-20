#pragma once

#include <string_view>
#include <vector>
#include <cstddef>

namespace lqcd::meson {

// GSL 规范：小函数职责单一，[[nodiscard]] 强制检查错误
// 解析单文件内容并计算多源移位平均后的 48 维向量
[[nodiscard]] std::vector<double> extract_single_file_averaged_block(
    std::string_view file_content,
    std::string_view direction,
    std::string_view start_tag,
    size_t num_lines = 48,
    size_t max_blocks = 16);

} // namespace lqcd::meson

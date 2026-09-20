#pragma once

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace iodata {

// GSL 规范：小函数职责单一，[[nodiscard]] 强制检查结果
[[nodiscard]] std::vector<std::filesystem::path> scan_directory(
    const std::filesystem::path& dir_path,
    std::string_view extension_filter = "");

// 按照文件名中提取的数字做自然排序（与 Python natsort 严格等价）
[[nodiscard]] std::vector<std::filesystem::path> scan_natural_sorted(
    const std::filesystem::path& dir_path,
    std::string_view prefix_pattern = "",
    std::string_view suffix_pattern = "");

// 从形如 "test1_lhadrons_10500_mesons_multi_src" 中提取整数 ID
[[nodiscard]] long long extract_numeric_id(std::string_view filename) noexcept;

} // namespace iodata

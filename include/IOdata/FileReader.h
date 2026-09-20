#pragma once

#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>
#include <span>

namespace iodata {

// GSL 规范：小函数职责单一，[[nodiscard]] 强制检查错误
[[nodiscard]] std::string read_file_to_string(const std::filesystem::path& file_path);

// 零拷贝逐行遍历文本（避免反复 getline 堆分配）
template <typename LineCallback>
void for_each_line(std::string_view content, LineCallback&& callback) {
    size_t start = 0;
    const size_t n = content.size();
    while (start < n) {
        size_t end = content.find('\n', start);
        if (end == std::string_view::npos) {
            std::string_view line = content.substr(start);
            if (!line.empty() && line.back() == '\r') {
                line.remove_suffix(1);
            }
            callback(line);
            break;
        }
        std::string_view line = content.substr(start, end - start);
        if (!line.empty() && line.back() == '\r') {
            line.remove_suffix(1);
        }
        callback(line);
        start = end + 1;
    }
}

} // namespace iodata

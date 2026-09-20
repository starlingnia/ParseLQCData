#pragma once

#include <string_view>
#include <vector>
#include <charconv>
#include <cstdlib>
#include <cstring>
#include <span>

namespace iodata {

// GSL 规范：小函数职责单一，[[nodiscard]] 强制检查结果
[[nodiscard]] inline bool parse_double(std::string_view sv, double& out_val) noexcept {
    // 快速去除首尾空格
    while (!sv.empty() && (sv.front() == ' ' || sv.front() == '\t' || sv.front() == '\r' || sv.front() == '\n')) {
        sv.remove_prefix(1);
    }
    while (!sv.empty() && (sv.back() == ' ' || sv.back() == '\t' || sv.back() == '\r' || sv.back() == '\n')) {
        sv.remove_suffix(1);
    }
    if (sv.empty()) {
        return false;
    }

#if defined(__GNUC__) && !defined(__clang__)
    auto [ptr, ec] = std::from_chars(sv.data(), sv.data() + sv.size(), out_val);
    if (ec == std::errc() && ptr == sv.data() + sv.size()) {
        return true;
    }
#endif

    // 通用兼容方案（处理科学计数法与 Clang libc++ 浮点 from_chars 限制）
    char buf[64];
    if (sv.size() < sizeof(buf)) {
        std::memcpy(buf, sv.data(), sv.size());
        buf[sv.size()] = '\0';
        char* end = nullptr;
        out_val = std::strtod(buf, &end);
        return end == buf + sv.size();
    }

    return false;
}

[[nodiscard]] inline bool parse_int(std::string_view sv, int& out_val) noexcept {
    while (!sv.empty() && (sv.front() == ' ' || sv.front() == '\t' || sv.front() == '\r' || sv.front() == '\n')) {
        sv.remove_prefix(1);
    }
    while (!sv.empty() && (sv.back() == ' ' || sv.back() == '\t' || sv.back() == '\r' || sv.back() == '\n')) {
        sv.remove_suffix(1);
    }
    if (sv.empty()) {
        return false;
    }

    auto [ptr, ec] = std::from_chars(sv.data(), sv.data() + sv.size(), out_val);
    return ec == std::errc() && ptr == sv.data() + sv.size();
}

// 零拷贝空白字符分割视图提取
template <typename Callback>
inline void for_each_token(std::string_view sv, Callback&& callback) {
    size_t start = 0;
    const size_t n = sv.size();
    while (start < n) {
        while (start < n && (sv[start] == ' ' || sv[start] == '\t' || sv[start] == '\r' || sv[start] == '\n')) {
            ++start;
        }
        if (start >= n) break;
        size_t end = start;
        while (end < n && sv[end] != ' ' && sv[end] != '\t' && sv[end] != '\r' && sv[end] != '\n') {
            ++end;
        }
        callback(sv.substr(start, end - start));
        start = end;
    }
}

} // namespace iodata

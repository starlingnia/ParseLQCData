#include <IOdata/DirectoryScanner.h>

#include <algorithm>
#include <cctype>

namespace iodata {

[[nodiscard]] long long extract_numeric_id(std::string_view filename) noexcept {
    // 寻找最长连续数字序列
    long long result = -1;
    size_t i = 0;
    while (i < filename.size()) {
        if (std::isdigit(static_cast<unsigned char>(filename[i]))) {
            long long current = 0;
            while (i < filename.size() && std::isdigit(static_cast<unsigned char>(filename[i]))) {
                current = current * 10 + (filename[i] - '0');
                ++i;
            }
            result = current;
            // 继续寻找，如果有多段数字通常取中间配置编号（例如 test1_lhadrons_1000_... 取后面的 1000）
        } else {
            ++i;
        }
    }
    return result;
}

[[nodiscard]] std::vector<std::filesystem::path> scan_directory(
    const std::filesystem::path& dir_path,
    std::string_view extension_filter) {

    std::vector<std::filesystem::path> result;
    if (!std::filesystem::exists(dir_path) || !std::filesystem::is_directory(dir_path)) {
        return result;
    }

    for (const auto& entry : std::filesystem::directory_iterator(dir_path)) {
        if (!entry.is_regular_file()) continue;
        if (!extension_filter.empty()) {
            if (entry.path().extension().string() != extension_filter) {
                continue;
            }
        }
        result.push_back(entry.path());
    }

    return result;
}

[[nodiscard]] std::vector<std::filesystem::path> scan_natural_sorted(
    const std::filesystem::path& dir_path,
    std::string_view prefix_pattern,
    std::string_view suffix_pattern) {

    std::vector<std::filesystem::path> files;
    if (!std::filesystem::exists(dir_path) || !std::filesystem::is_directory(dir_path)) {
        return files;
    }

    for (const auto& entry : std::filesystem::directory_iterator(dir_path)) {
        if (!entry.is_regular_file()) continue;
        const std::string name = entry.path().filename().string();
        if (!prefix_pattern.empty() && name.find(prefix_pattern) == std::string::npos) {
            continue;
        }
        if (!suffix_pattern.empty() && name.find(suffix_pattern) == std::string::npos) {
            continue;
        }
        files.push_back(entry.path());
    }

    // 按照提取出的数值升序排序
    std::sort(files.begin(), files.end(), [](const std::filesystem::path& a, const std::filesystem::path& b) {
        const long long id_a = extract_numeric_id(a.filename().string());
        const long long id_b = extract_numeric_id(b.filename().string());
        if (id_a != id_b) {
            return id_a < id_b;
        }
        return a.filename().string() < b.filename().string();
    });

    return files;
}

} // namespace iodata

#include <ParseLQCData/MesonAnalysis/MesonExtractor.h>
#include <IOdata/FastParser.h>
#include <IOdata/FileReader.h>

#include <array>
#include <vector>
#include <string_view>

namespace lqcd::meson {

namespace {

// 快速从行尾解析 "x/y/z/t" 坐标 (例如 "0/0/0/0" 或 "12/0/0/0")
[[nodiscard]] std::array<int, 4> parse_coords(std::string_view line) noexcept {
    std::array<int, 4> nums = {0, 0, 0, 0};
    while (!line.empty() && (line.back() == ' ' || line.back() == '\t' || line.back() == '\r' || line.back() == '\n')) {
        line.remove_suffix(1);
    }
    const size_t last_space = line.find_last_of(" \t");
    std::string_view coord_part = (last_space != std::string_view::npos) ? line.substr(last_space + 1) : line;

    size_t idx = 0;
    size_t start = 0;
    while (idx < 4 && start < coord_part.size()) {
        size_t slash = coord_part.find('/', start);
        std::string_view token = (slash != std::string_view::npos) ? coord_part.substr(start, slash - start) : coord_part.substr(start);
        int val = 0;
        if (iodata::parse_int(token, val)) {
            nums[idx++] = val;
        } else {
            break;
        }
        if (slash == std::string_view::npos) break;
        start = slash + 1;
    }
    return nums;
}

} // namespace

[[nodiscard]] std::vector<double> extract_single_file_averaged_block(
    std::string_view file_content,
    std::string_view direction,
    std::string_view start_tag,
    const size_t num_lines,
    const size_t max_blocks) {

    std::vector<double> sum_rolled(num_lines, 0.0);
    std::vector<double> current_block;
    current_block.reserve(num_lines);

    std::array<int, 4> nums = {0, 0, 0, 0};
    bool in_block = false;
    size_t line_count = 0;
    size_t block_count = 0;

    iodata::for_each_line(file_content, [&](std::string_view line) {
        if (block_count >= max_blocks) {
            return;
        }

        if (line.find(start_tag) != std::string_view::npos) {
            nums = parse_coords(line);
            in_block = true;
            current_block.clear();
            line_count = 0;
            return;
        }

        if (in_block) {
            // 解析列：<index> <real_val> <imag_val>
            size_t token_idx = 0;
            double real_val = 0.0;
            bool got_real = false;

            iodata::for_each_token(line, [&](std::string_view token) {
                if (token_idx == 1) {
                    if (iodata::parse_double(token, real_val)) {
                        got_real = true;
                    }
                }
                ++token_idx;
            });

            if (got_real) {
                current_block.push_back(real_val);
                ++line_count;
            }

            if (line_count == num_lines) {
                int shift_offset = 0;
                if (direction == "DIRX") {
                    shift_offset = nums[0];
                } else if (direction == "DIRY") {
                    shift_offset = nums[1];
                } else if (direction == "DIRZ") {
                    shift_offset = nums[2];
                }

                // 环形循环移位并累加
                const int n = static_cast<int>(num_lines);
                for (int i = 0; i < n; ++i) {
                    const int src_idx = ((i + (shift_offset % n)) % n + n) % n;
                    sum_rolled[i] += current_block[src_idx];
                }

                ++block_count;
                in_block = false;
                line_count = 0;
                current_block.clear();
            }
        }
    });

    if (block_count == 0) {
        return std::vector<double>(num_lines, 0.0);
    }

    const double inv_count = 1.0 / static_cast<double>(block_count);
    for (size_t i = 0; i < num_lines; ++i) {
        sum_rolled[i] *= inv_count;
    }

    return sum_rolled;
}

} // namespace lqcd::meson

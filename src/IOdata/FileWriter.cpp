#include <IOdata/FileWriter.h>

#include <fstream>
#include <iomanip>
#include <print>

namespace iodata {

[[nodiscard]] bool write_matrix_to_csv(
    const std::filesystem::path& file_path,
    std::span<const double> matrix,
    size_t rows,
    size_t cols,
    bool include_header,
    std::string_view header_str) {

    if (matrix.size() < rows * cols) {
        return false;
    }

    if (file_path.has_parent_path()) {
        std::filesystem::create_directories(file_path.parent_path());
    }

    std::ofstream out(file_path);
    if (!out.is_open()) {
        return false;
    }

    if (include_header && !header_str.empty()) {
        out << header_str << "\n";
    }

    out << std::setprecision(16);
    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            out << matrix[r * cols + c];
            if (c + 1 < cols) {
                out << ',';
            }
        }
        out << '\n';
    }

    return true;
}

[[nodiscard]] bool write_pairs_to_csv(
    const std::filesystem::path& file_path,
    std::span<const double> col1,
    std::span<const double> col2,
    bool include_header,
    std::string_view header_str) {

    if (col1.size() != col2.size()) {
        return false;
    }

    if (file_path.has_parent_path()) {
        std::filesystem::create_directories(file_path.parent_path());
    }

    std::ofstream out(file_path);
    if (!out.is_open()) {
        return false;
    }

    if (include_header && !header_str.empty()) {
        out << header_str << "\n";
    }

    out << std::setprecision(17);
    for (size_t i = 0; i < col1.size(); ++i) {
        out << col1[i] << ',' << col2[i] << '\n';
    }

    return true;
}

} // namespace iodata

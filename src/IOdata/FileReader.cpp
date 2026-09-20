#include <IOdata/FileReader.h>

#include <fstream>
#include <sstream>

namespace iodata {

[[nodiscard]] std::string read_file_to_string(const std::filesystem::path& file_path) {
    std::ifstream file(file_path, std::ios::in | std::ios::binary);
    if (!file.is_open()) {
        return "";
    }
    const auto file_size = std::filesystem::file_size(file_path);
    std::string buffer;
    buffer.resize(file_size);
    file.read(buffer.data(), static_cast<std::streamsize>(file_size));
    return buffer;
}

} // namespace iodata

#include <ParseLQCData/MesonAnalysis/MesonPipeline.h>
#include <ParseLQCData/Core/LQCDataTypes.h>

#include <cstring>
#include <sstream>
#include <vector>

extern "C" {

// C 风格导出接口
int run_meson_pipeline_c_api(
    const char* input_dir,
    const char* channel_types_csv,
    const char* channel_dirs_csv,
    int binsize,
    int num_lines,
    int thread_count,
    double* out_means,
    double* out_errors,
    double** out_folded_jk,
    int* out_n_bins,
    int* out_n_raw_cfgs) {

    if (!input_dir || !channel_types_csv || !channel_dirs_csv || !out_means || !out_errors) {
        return -1;
    }

    // 解析逗号分隔的信道配置
    std::vector<std::string> types;
    std::vector<std::string> dirs;

    std::stringstream ss_types(channel_types_csv);
    std::string item;
    while (std::getline(ss_types, item, ',')) {
        if (!item.empty()) types.push_back(item);
    }

    std::stringstream ss_dirs(channel_dirs_csv);
    while (std::getline(ss_dirs, item, ',')) {
        if (!item.empty()) dirs.push_back(item);
    }

    if (types.size() != dirs.size() || types.empty()) {
        return -2;
    }

    std::vector<lqcd::MesonChannelMapping> channels;
    channels.reserve(types.size());
    for (size_t i = 0; i < types.size(); ++i) {
        channels.push_back({types[i], dirs[i]});
    }

    const auto result = lqcd::meson::run_meson_pipeline(
        input_dir,
        channels,
        static_cast<size_t>(binsize > 0 ? binsize : 4),
        static_cast<size_t>(num_lines > 0 ? num_lines : 48),
        static_cast<size_t>(thread_count > 0 ? thread_count : 0)
    );

    if (result.means.size() != static_cast<size_t>(num_lines)) {
        return -3;
    }

    std::memcpy(out_means, result.means.data(), num_lines * sizeof(double));
    std::memcpy(out_errors, result.errors.data(), num_lines * sizeof(double));

    if (out_n_bins) {
        *out_n_bins = static_cast<int>(result.n_bins);
    }
    if (out_n_raw_cfgs) {
        *out_n_raw_cfgs = static_cast<int>(result.n_raw_cfgs);
    }

    if (out_folded_jk && result.n_bins > 0) {
        const size_t matrix_size = static_cast<size_t>(num_lines) * result.n_bins;
        double* buffer = new double[matrix_size];
        std::memcpy(buffer, result.folded_jk_matrix.data(), matrix_size * sizeof(double));
        *out_folded_jk = buffer;
    }

    return 0;
}

void free_lqcd_buffer(void* ptr) {
    if (ptr) {
        delete[] static_cast<double*>(ptr);
    }
}

} // extern "C"

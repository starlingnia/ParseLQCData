#pragma once

#include <core/LQCDataTypes.h>
#include <span>
#include <filesystem>

namespace lqcd::meson {

[[nodiscard]] MesonAnalysisResult run_meson_pipeline(
    const std::filesystem::path& input_dir,
    std::span<const MesonChannelMapping> channels,
    size_t binsize = 4,
    size_t num_lines = 48,
    size_t thread_count = 0,
    bool is_single_source = false);

} // namespace lqcd::meson

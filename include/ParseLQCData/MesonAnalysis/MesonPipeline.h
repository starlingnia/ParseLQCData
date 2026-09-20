#pragma once

#include <ParseLQCData/Core/LQCDataTypes.h>
#include <filesystem>
#include <span>
#include <string_view>

namespace lqcd::meson {

// GSL 规范：[[nodiscard]] 强制检查结果
[[nodiscard]] MesonAnalysisResult run_meson_pipeline(
    const std::filesystem::path& input_dir,
    std::span<const MesonChannelMapping> channels,
    size_t binsize = 4,
    size_t num_lines = 48,
    size_t thread_count = 0);

} // namespace lqcd::meson

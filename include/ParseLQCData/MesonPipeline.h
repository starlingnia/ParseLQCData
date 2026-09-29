#pragma once

#include <core/LQCDataTypes.h>
#include <span>
#include <filesystem>
#include <map>
#include <string>

namespace lqcd::meson {

[[nodiscard]] MesonAnalysisResult run_meson_pipeline(
    const std::filesystem::path& input_dir,
    std::span<const MesonChannelMapping> channels,
    size_t binsize = 4,
    size_t num_lines = 48,
    size_t thread_count = 0,
    bool is_single_source = false);

[[nodiscard]] std::map<std::string, MesonAnalysisResult> run_meson_pipeline_batch(
    const std::filesystem::path& input_dir,
    const std::map<std::string, std::vector<MesonChannelMapping>>& channel_defs,
    size_t binsize = 4,
    size_t num_lines = 48,
    size_t thread_count = 0,
    bool is_single_source = false);

} // namespace lqcd::meson

#pragma once

#include <core/LQCDataTypes.h>
#include <span>
#include <filesystem>
#include <map>
#include <string>

namespace lqcd::meson {

// num_lines 为 0 时，从 input_dir 名称推导空间方向长度 Ns。
[[nodiscard]] MesonAnalysisResult run_meson_pipeline(
    const std::filesystem::path& input_dir,
    std::span<const MesonChannelMapping> channels,
    size_t binsize = 4,
    size_t num_lines = 0,
    size_t thread_count = 0,
    bool is_single_source = false);

[[nodiscard]] std::map<std::string, MesonAnalysisResult> run_meson_pipeline_batch(
    const std::filesystem::path& input_dir,
    const std::map<std::string, std::vector<MesonChannelMapping>>& channel_defs,
    size_t binsize = 4,
    size_t num_lines = 0,
    size_t thread_count = 0,
    bool is_single_source = false);

} // namespace lqcd::meson

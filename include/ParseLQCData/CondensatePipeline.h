#pragma once

#include <CondensateAnalysis/CondensateExtractor.h>
#include <filesystem>

namespace lqcd::condensate {

[[nodiscard]] CondensateResult run_condensate_pipeline(
    const std::filesystem::path& input_dir,
    double m_light,
    double m_strange,
    double m_residual,
    double zm_factor,
    size_t thread_count = 0);

} // namespace lqcd::condensate

#include <ParseLQCData/CondensatePipeline.h>
#include <CondensateAnalysis/CondensateExtractor.h>

#include <filesystem>

namespace lqcd::condensate {

[[nodiscard]] CondensateResult run_condensate_pipeline(
    const std::filesystem::path& input_dir,
    const double m_light,
    const double m_strange,
    const double m_residual,
    const double zm_factor,
    [[maybe_unused]] size_t thread_count) {

    std::filesystem::path base_path = input_dir;
    if (std::filesystem::exists(base_path / "test_condensate")) {
        base_path /= "test_condensate";
    }

    return process_chiral_condensate(base_path, m_light, m_strange, m_residual, zm_factor);
}

} // namespace lqcd::condensate

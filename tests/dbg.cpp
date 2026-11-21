#include <ParseLQCData/CondensatePipeline.h>
#include <CondensateAnalysis/CondensateExtractor.h>
#include <core/PhysicsSetup.h>
#include <iostream>

int main(int argc, char *argv[]) {
    // 1. 加载 docs/physics_setup.json 里的手征凝聚配置
    const auto configs = lqcd::load_condensate_configs_from_docs();
    std::cout << "Loaded configs count: " << configs.size() << "\n";

    // 2. 直接调用 Condensate 管道或算子 (在下面或 src/CondensateAnalysis/*.cpp 打红点断点)
    const std::filesystem::path test_dir = "data/readin/L32T12beta4.17/test_condensate";
    const auto res = lqcd::condensate::run_condensate_pipeline(
        test_dir, 0.001001, 0.0384, 0.000339722, 0.966247, 1
    );

    std::cout << "num_cfgs = " << res.num_cfgs
              << ", mean = " << res.mean
              << " +/- " << res.error << "\n";
    return 0;
}

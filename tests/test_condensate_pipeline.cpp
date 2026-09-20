#include <ParseLQCData/CondensatePipeline.h>
#include <core/PhysicsSetup.h>

#include <cassert>
#include <iostream>
#include <cmath>

int main() {
    std::cout << "=================================================\n";
    std::cout << "[TEST] C++ CondensatePipeline Direct Unit Test\n";
    std::cout << "=================================================\n";

    // 1. 测试 docs/physics_setup.json 动态加载 CONDENSATE_CONFIGS
    std::cout << ">>> [1/2] 测试 docs/physics_setup.json 加载手征凝聚配置...\n";
    const auto configs = lqcd::load_condensate_configs_from_docs();
    assert(!configs.empty());
    std::cout << "  -> 成功加载 " << configs.size() << " 组手征凝聚配置 (PASSED)\n";

    // 2. 测试直接运行核心管道 run_condensate_pipeline
    std::cout << ">>> [2/2] 测试直接运行 core/CondensatePipeline...\n";
    const std::filesystem::path test_dir = "data/readin/L32T12beta4.17/test_condensate";
    if (std::filesystem::exists(test_dir)) {
        const auto res = lqcd::condensate::run_condensate_pipeline(
            test_dir, 0.001001, 0.0384, 0.000339722, 0.966247, 4
        );
        std::cout << "  -> 处理构型数: " << res.num_cfgs << "\n";
        std::cout << "  -> 手征凝聚均值: " << res.mean << " +/- " << res.error << "\n";
        assert(res.num_cfgs > 0);
        assert(res.mean > 0.0);
        assert(res.error > 0.0);
        std::cout << "  -> CondensatePipeline 核心管道测试通过 (PASSED)\n";
    } else {
        std::cout << "  -> 未找到测试数据目录，跳过实际数据运算\n";
    }

    std::cout << "=================================================\n";
    std::cout << "[OK] Condensate 核心管道测试全部通过！\n";
    std::cout << "=================================================\n";
    return 0;
}

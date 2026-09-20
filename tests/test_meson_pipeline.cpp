#include <core/MesonPipeline.h>
#include <core/PhysicsSetup.h>
#include <Statistics/Resampling.h>
#include <IOdata/FastParser.h>

#include <cassert>
#include <iostream>
#include <vector>
#include <cmath>

int main() {
    std::cout << "=================================================\n";
    std::cout << "[TEST] C++ MesonPipeline & Components Unit Tests\n";
    std::cout << "=================================================\n";

    // 1. 测试 docs/physics_setup.json 动态加载
    std::cout << ">>> [1/3] 测试 docs/physics_setup.json 动态加载...\n";
    const auto channel_defs = lqcd::load_channel_configs_from_docs();
    assert(!channel_defs.empty());
    assert(channel_defs.contains("AV"));
    assert(channel_defs.contains("S"));
    assert(channel_defs.contains("Tt"));
    assert(channel_defs.contains("PS"));
    assert(channel_defs.contains("Xt"));
    assert(channel_defs.contains("Vec"));
    assert(channel_defs.at("AV").size() == 6);
    assert(channel_defs.at("S").size() == 3);
    std::cout << "  -> 成功加载 " << channel_defs.size() << " 个信道配置，AV 包含 " 
              << channel_defs.at("AV").size() << " 个空间分量 (PASSED)\n";

    // 2. 测试 FastParser 零拷贝快速解析
    std::cout << ">>> [2/3] 测试 FastParser 数值解析精度...\n";
    double val = 0.0;
    assert(iodata::parse_double("  -1.23456789e-05 \n", val));
    assert(std::abs(val - (-1.23456789e-05)) < 1e-12);
    int int_val = 0;
    assert(iodata::parse_int(" 10500 \r\n", int_val));
    assert(int_val == 10500);
    std::cout << "  -> FastParser double/int 解析测试通过 (PASSED)\n";

    // 3. 测试 Resampling Jackknife 留一重采样
    std::cout << ">>> [3/3] 测试 Jackknife 统计重采样...\n";
    // 构造简单样本矩阵: 1 行，4 列，值为 1.0, 2.0, 3.0, 4.0
    std::vector<double> sample_matrix = {1.0, 2.0, 3.0, 4.0};
    auto binned = lqcd::stats::compute_block_binning(sample_matrix, 1, 4, 1);
    auto jk = lqcd::stats::compute_jackknife(binned, 1, 4);
    assert(jk.size() == 4);
    // jk 样本分别为去掉一个后的均值: (2+3+4)/3 = 3, (1+3+4)/3 = 2.666..., etc.
    assert(std::abs(jk[0] - 3.0) < 1e-10);
    std::cout << "  -> Jackknife 重采样算法精度验证通过 (PASSED)\n";

    std::cout << "=================================================\n";
    std::cout << "[OK] 所有 C++ 单元测试全部通过！\n";
    std::cout << "=================================================\n";
    return 0;
}

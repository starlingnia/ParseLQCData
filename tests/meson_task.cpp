#include <ParseLQCData/Registry.h>
#include <ParseLQCData/MesonAnalysis/MesonPipeline.h>
#include <IOdata/FileWriter.h>

#include <chrono>
#include <filesystem>
#include <print>
#include <vector>

int run_meson_av_task(std::span<const std::string_view> args) {
    const std::string beta_str = args.empty() ? "17" : std::string(args[0]);
    const size_t binsize = args.size() > 1 ? static_cast<size_t>(std::stoul(std::string(args[1]))) : 4;

    const std::filesystem::path input_dir = "/Users/junxiongnie/code/ana/dat/readin/48x16b4." + beta_str + "/Output";
    const std::filesystem::path output_file = "output/b4." + beta_str + "/save_AV_cpp.csv";

    std::println("=================================================");
    std::println("启动 C++ 原生 Meson AV 分析任务 (beta: {}, binsize: {})", beta_str, binsize);
    std::println("输入目录: {}", input_dir.string());
    std::println("输出文件: {}", output_file.string());

    const std::vector<lqcd::MesonChannelMapping> av_channels = {
        {"AVector2", "DIRX"},
        {"AVector3", "DIRX"},
        {"AVector3", "DIRY"},
        {"AVector1", "DIRY"},
        {"AVector1", "DIRZ"},
        {"AVector2", "DIRZ"},
    };

    const auto start_time = std::chrono::steady_clock::now();
    const auto result = lqcd::meson::run_meson_pipeline(input_dir, av_channels, binsize, 48);
    const auto end_time = std::chrono::steady_clock::now();

    const double elapsed_s = std::chrono::duration<double>(end_time - start_time).count();
    std::println("C++ 端到端并发分析执行完毕！耗时: {:.3f} 秒 (处理配置数: {}, 生成 Bin 数: {})",
                 elapsed_s, result.n_raw_cfgs, result.n_bins);

    if (result.means.empty()) {
        std::println(stderr, "错误: 计算结果为空！");
        return 1;
    }

    if (iodata::write_pairs_to_csv(output_file, result.means, result.errors)) {
        std::println("成功将结果落盘至: {}", output_file.string());
    } else {
        std::println(stderr, "警告: 写入结果文件失败！");
    }

    std::println("前 5 行样例输出 (mean, error):");
    for (size_t i = 0; i < std::min(result.means.size(), size_t(5)); ++i) {
        std::println("  [{:02d}] {:.16e}, {:.16e}", i, result.means[i], result.errors[i]);
    }
    std::println("=================================================");

    return 0;
}

REGISTER_TASK("meson_av", run_meson_av_task, "Run Meson AV channel analysis on Lattice QCD dataset");

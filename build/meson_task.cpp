#include <ParseLQCData/Registry.h>
#include <core/PhysicsSetup.h>
#include <ParseLQCData/MesonPipeline.h>
#include <IOdata/FileWriter.h>

#include <chrono>
#include <filesystem>
#include <iostream>
#include <print>
#include <vector>
#include <string>

namespace {

// 辅助函数：根据 beta 智能定位 readin 目录下的强子关联函数数据目录
std::filesystem::path find_meson_input_dir(const std::string& beta_str) {
    const std::filesystem::path base_candidates[] = {
        "data/readin",
        "../data/readin",
        "../../data/readin",
        "/Users/junxiongnie/code/ana/dat/readin"
    };

    std::string b_clean = beta_str;
    if (b_clean.starts_with("4.")) b_clean = b_clean.substr(2);
    else if (b_clean.starts_with("b4.")) b_clean = b_clean.substr(3);

    for (const auto& base : base_candidates) {
        if (!std::filesystem::exists(base)) continue;

        // 1. 精确匹配 48x16b4.{beta}
        auto legacy_p = base / ("48x16b4." + b_clean) / "Output";
        if (std::filesystem::exists(legacy_p)) return legacy_p;
        legacy_p = base / ("48x16b4." + b_clean);
        if (std::filesystem::exists(legacy_p)) return legacy_p;

        // 2. 遍历扫描匹配 L48T16beta4.{beta}...
        for (const auto& entry : std::filesystem::directory_iterator(base)) {
            if (!entry.is_directory()) continue;
            const std::string fname = entry.path().filename().string();
            if (fname.find("beta4." + b_clean) != std::string::npos ||
                fname.find("b4." + b_clean) != std::string::npos) {
                if (std::filesystem::exists(entry.path() / "Output")) {
                    return entry.path() / "Output";
                }
                return entry.path();
            }
        }
    }

    // 默认回退
    return std::filesystem::path("data/readin/48x16b4." + b_clean + "/Output");
}

int execute_meson_measurement(
    std::string beta_str,
    std::string channel_name,
    size_t binsize,
    bool is_single_source) {

    if (beta_str.starts_with("4.")) beta_str = beta_str.substr(2);

    const auto channel_defs = lqcd::load_channel_configs_from_docs();
    auto it = channel_defs.find(channel_name);
    if (it == channel_defs.end()) {
        std::println(stderr, "错误: 未在 docs/physics_setup.json 中找到信道 '{}'", channel_name);
        std::print(stderr, "支持的信道包括: ");
        for (const auto& [ch, _] : channel_defs) std::print(stderr, "{} ", ch);
        std::println(stderr, "");
        return 1;
    }

    const auto input_dir = find_meson_input_dir(beta_str);
    const std::string mode_str = is_single_source ? "single" : "multi";
    const std::string out_dir = is_single_source ? "output/pickdata-singlesrc/b4." + beta_str
                                                 : "output/pickdata/b4." + beta_str;
    const std::filesystem::path output_file = out_dir + "/save_" + channel_name + "_cpp.csv";

    std::println("=================================================");
    std::println("启动 C++ Meson 关联函数测量 (模式: {}, Beta: 4.{}, 信道: {}, Binsize: {})",
                 mode_str, beta_str, channel_name, binsize);
    std::println("输入数据目录: {}", input_dir.string());
    std::println("结果保存路径: {}", output_file.string());

    if (!std::filesystem::exists(input_dir)) {
        std::println(stderr, "错误: 输入目录不存在: {}", input_dir.string());
        return 1;
    }

    const auto t0 = std::chrono::steady_clock::now();
    const auto result = lqcd::meson::run_meson_pipeline(
        input_dir, it->second, binsize, 48, 0, is_single_source
    );
    const auto t1 = std::chrono::steady_clock::now();

    const double elapsed_s = std::chrono::duration<double>(t1 - t0).count();
    std::println("C++ 流水线测量完成！耗时: {:.3f} 秒 (构型数: {}, 生成 Bin 数: {})",
                 elapsed_s, result.n_raw_cfgs, result.n_bins);

    if (result.means.empty()) {
        std::println(stderr, "错误: 计算结果为空！");
        return 1;
    }

    std::filesystem::create_directories(output_file.parent_path());
    if (iodata::write_pairs_to_csv(output_file, result.means, result.errors)) {
        std::println("成功落盘测量结果: {}", output_file.string());
    } else {
        std::println(stderr, "警告: 写入结果文件失败！");
    }

    std::println("关联函数前 5 步长数据 (mean, error):");
    for (size_t i = 0; i < std::min(result.means.size(), size_t(5)); ++i) {
        std::println("  t={:02d}: {:.16e} +/- {:.16e}", i, result.means[i], result.errors[i]);
    }
    std::println("=================================================");

    return 0;
}

} // namespace

// 1. 多源测量任务
int run_meson_multi_task(std::span<const std::string_view> args) {
    const std::string beta_str = args.empty() ? "17" : std::string(args[0]);
    const std::string channel = args.size() > 1 ? std::string(args[1]) : "AV";
    const size_t binsize = args.size() > 2 ? static_cast<size_t>(std::stoul(std::string(args[2]))) : 4;
    return execute_meson_measurement(beta_str, channel, binsize, false);
}

// 2. 单源测量任务
int run_meson_single_task(std::span<const std::string_view> args) {
    const std::string beta_str = args.empty() ? "17" : std::string(args[0]);
    const std::string channel = args.size() > 1 ? std::string(args[1]) : "S";
    const size_t binsize = args.size() > 2 ? static_cast<size_t>(std::stoul(std::string(args[2]))) : 4;
    return execute_meson_measurement(beta_str, channel, binsize, true);
}

// 3. 批量全信道测量任务: [beta] [multi|single]
int run_meson_all_task(std::span<const std::string_view> args) {
    const std::string beta_str = args.empty() ? "17" : std::string(args[0]);
    const std::string mode = args.size() > 1 ? std::string(args[1]) : "multi";
    const bool is_single = (mode == "single");

    const auto channel_defs = lqcd::load_channel_configs_from_docs();
    std::println(">>> 批量测量 Beta: 4.{} 下的所有 {} 种信道 (模式: {})",
                 beta_str, channel_defs.size(), mode);

    for (const auto& [ch, _] : channel_defs) {
        const size_t binsize = (beta_str == "18" && ch == "S" && !is_single) ? 5 : 4;
        execute_meson_measurement(beta_str, ch, binsize, is_single);
    }
    return 0;
}

REGISTER_TASK("meson_multi", run_meson_multi_task, "Run multi-source meson correlation measurement: [beta] [channel] [binsize]");
REGISTER_TASK("meson_single", run_meson_single_task, "Run single-source meson correlation measurement: [beta] [channel] [binsize]");
REGISTER_TASK("meson_all", run_meson_all_task, "Batch measure all meson channels for a given beta: [beta] [multi|single]");

#ifdef STANDALONE_MESON
int main(int argc, char* argv[]) {
    std::vector<std::string_view> args;
    for (int i = 1; i < argc; ++i) args.emplace_back(argv[i]);
    return run_meson_multi_task(args);
}
#endif

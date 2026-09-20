#include <ParseLQCData/Registry.h>
#include <core/PhysicsSetup.h>
#include <ParseLQCData/MesonPipeline.h>
#include <IOdata/FileWriter.h>

#include <chrono>
#include <filesystem>
#include <print>
#include <vector>
#include <string>

namespace {

std::filesystem::path find_single_meson_dir(const std::string& beta_str) {
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
        auto legacy_p = base / ("48x16b4." + b_clean) / "Output";
        if (std::filesystem::exists(legacy_p)) return legacy_p;
        legacy_p = base / ("48x16b4." + b_clean);
        if (std::filesystem::exists(legacy_p)) return legacy_p;

        for (const auto& entry : std::filesystem::directory_iterator(base)) {
            if (!entry.is_directory()) continue;
            const std::string fname = entry.path().filename().string();
            if (fname.find("beta4." + b_clean) != std::string::npos ||
                fname.find("b4." + b_clean) != std::string::npos) {
                if (std::filesystem::exists(entry.path() / "Output")) return entry.path() / "Output";
                return entry.path();
            }
        }
    }
    return std::filesystem::path("data/readin/48x16b4." + b_clean + "/Output");
}

int measure_channel_single(std::string beta_str, std::string ch, size_t binsize) {
    if (beta_str.starts_with("4.")) beta_str = beta_str.substr(2);
    const auto channel_defs = lqcd::load_channel_configs_from_docs();
    auto it = channel_defs.find(ch);
    if (it == channel_defs.end()) return 1;

    const auto input_dir = find_single_meson_dir(beta_str);
    const std::filesystem::path output_file = "output/pickdata-singlesrc/b4." + beta_str + "/save_" + ch + "_cpp.csv";

    const auto t0 = std::chrono::steady_clock::now();
    const auto res = lqcd::meson::run_meson_pipeline(input_dir, it->second, binsize, 48, 0, true);
    const auto t1 = std::chrono::steady_clock::now();

    const double elapsed_s = std::chrono::duration<double>(t1 - t0).count();
    std::println("Meson Single 测量完成 (Beta: 4.{}, 信道: {}, 耗时: {:.3f}s, 构型数: {}, Bins: {})",
                 beta_str, ch, elapsed_s, res.n_raw_cfgs, res.n_bins);

    if (!res.means.empty()) {
        std::filesystem::create_directories(output_file.parent_path());
        (void)iodata::write_pairs_to_csv(output_file, res.means, res.errors);
    }
    return 0;
}

} // namespace

int run_meson_single_standalone_task(std::span<const std::string_view> args) {
    const std::string beta_str = args.empty() ? "17" : std::string(args[0]);
    const std::string channel = args.size() > 1 ? std::string(args[1]) : "S";
    const size_t binsize = args.size() > 2 ? static_cast<size_t>(std::stoul(std::string(args[2]))) : 4;
    return measure_channel_single(beta_str, channel, binsize);
}

int run_meson_single_all_task(std::span<const std::string_view> args) {
    const std::string beta_str = args.empty() ? "17" : std::string(args[0]);
    const auto channel_defs = lqcd::load_channel_configs_from_docs();
    for (const auto& [ch, _] : channel_defs) {
        const size_t binsize = (beta_str == "18" || beta_str == "20" || beta_str == "23" || beta_str == "30") ? 5 : 4;
        measure_channel_single(beta_str, ch, binsize);
    }
    return 0;
}

#ifndef COMPILING_ALL_TASKS
REGISTER_TASK("meson_single", run_meson_single_standalone_task, "Run single-source meson measurement: [beta] [channel] [binsize]");
REGISTER_TASK("meson_single_all", run_meson_single_all_task, "Batch measure all single-source channels: [beta]");
#endif

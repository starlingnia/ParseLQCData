#include <ParseLQCData/Registry.h>
#include <core/PhysicsSetup.h>
#include <ParseLQCData/CondensatePipeline.h>
#include <IOdata/FileWriter.h>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <print>
#include <vector>
#include <string>

namespace {

// 智能定位手征凝聚数据集目录与物理参数
struct CondensateDatasetInfo {
    std::filesystem::path dir_path;
    double ml{0.001001};
    double ms{0.0384};
    double mres{0.000339722};
    double zm{0.966247};
    double beta{4.17};
    std::string name;
};

// 快速从目录名如 L32T12_beta4.17ms0.040m0.0020 或 L48T16beta4.13ms0.043547m0.000805 中提取物理参数
inline bool parse_condensate_dirname(std::string_view name, CondensateDatasetInfo& info) {
    size_t beta_pos = name.find("beta");
    if (beta_pos == std::string_view::npos) return false;
    size_t ms_pos = name.find("ms", beta_pos);
    if (ms_pos == std::string_view::npos) return false;
    size_t m_pos = name.find('m', ms_pos + 2);
    if (m_pos == std::string_view::npos) return false;

    // 提取 beta
    std::string_view beta_sv = name.substr(beta_pos + 4, ms_pos - (beta_pos + 4));
    (void)iodata::parse_double(beta_sv, info.beta);

    // 提取 ms
    std::string_view ms_sv = name.substr(ms_pos + 2, m_pos - (ms_pos + 2));
    (void)iodata::parse_double(ms_sv, info.ms);

    // 提取 ml
    size_t ml_start = m_pos + 1;
    if (ml_start < name.size() && name[ml_start] == '0' && ml_start + 1 < name.size() && name[ml_start + 1] != '.') {
        ml_start++;
    }
    size_t ml_end = name.find_first_of("_ \t\r\n", ml_start);
    if (ml_end == std::string_view::npos) ml_end = name.size();
    std::string_view ml_sv = name.substr(ml_start, ml_end - ml_start);
    (void)iodata::parse_double(ml_sv, info.ml);

    return true;
}

CondensateDatasetInfo resolve_condensate_info(std::string_view target_hint) {
    CondensateDatasetInfo info;

    const std::filesystem::path base_candidates[] = {
        "data/readin",
        "../data/readin",
        "../../data/readin",
        "/Users/junxiongnie/code/ana/dat/readin"
    };

    const auto configs = lqcd::load_condensate_configs_from_docs();

    std::filesystem::path base_dir = "data/readin";
    for (const auto& b : base_candidates) {
        if (std::filesystem::exists(b)) {
            base_dir = b;
            break;
        }
    }

    // 默认测试目录
    info.dir_path = base_dir / "L32T12beta4.17" / "test_condensate";
    info.name = "L32T12beta4.17";

    if (!target_hint.empty()) {
        std::string t(target_hint);
        // 1. 若直接是路径
        if (std::filesystem::exists(t)) {
            info.dir_path = t;
            info.name = std::filesystem::path(t).filename().string();
        } else if (std::filesystem::exists(base_dir / t)) {
            info.dir_path = base_dir / t;
            info.name = t;
        } else {
            // 2. 遍历 base_dir 模糊匹配
            for (const auto& entry : std::filesystem::directory_iterator(base_dir)) {
                if (!entry.is_directory()) continue;
                const std::string fname = entry.path().filename().string();
                if (fname.find(t) != std::string::npos) {
                    info.dir_path = entry.path();
                    info.name = fname;
                    break;
                }
            }
        }
    }

    // 优先从目录名解析 L, T, beta, ms, ml
    parse_condensate_dirname(info.name, info);

    // 从 docs/physics_setup.json 的 CONDENSATE_CONFIGS 补全残余质量 mres 与重整化常数 Zm
    for (const auto& c : configs) {
        if (std::abs(c.beta - info.beta) < 1e-4) {
            info.mres = c.mres;
            info.zm = c.zm;
            break;
        }
    }

    // 若未在列表中，使用 Domain Wall Fermion 拟合公式计算残余质量
    if (info.mres <= 1e-12 && info.beta > 0.1) {
        info.mres = 2.547e28 * std::exp(-17.559 * info.beta);
        info.zm = 1.0;
    }

    return info;
}

int execute_condensate_measurement(const CondensateDatasetInfo& info) {
    std::println("=================================================");
    std::println("启动 C++ 手征凝聚 (Chiral Condensate) 测量流水线");
    std::println("数据集: {}", info.name);
    std::println("目标路径: {}", info.dir_path.string());
    std::println("物理参数: ml={:.6f}, ms={:.6f}, mres={:.8f}, Zm={:.6f}",
                 info.ml, info.ms, info.mres, info.zm);

    if (!std::filesystem::exists(info.dir_path)) {
        std::println(stderr, "错误: 数据集路径不存在: {}", info.dir_path.string());
        return 1;
    }

    const auto t0 = std::chrono::steady_clock::now();
    const auto res = lqcd::condensate::run_condensate_pipeline(
        info.dir_path, info.ml, info.ms, info.mres, info.zm
    );
    const auto t1 = std::chrono::steady_clock::now();

    const double elapsed_s = std::chrono::duration<double>(t1 - t0).count();
    std::println("C++ 手征凝聚流水线执行完毕！耗时: {:.3f} 秒 (构型数: {})",
                 elapsed_s, res.num_cfgs);

    if (res.num_cfgs == 0) {
        std::println(stderr, "错误: 未找到可解析的 meas.* 或 XML 构型数据！");
        return 1;
    }

    std::println("测量物理结果:");
    std::println("  - 扣除残余质量并重整化后的手征凝聚: <pbp_sub> = {:.16e} +/- {:.16e}",
                 res.mean, res.error);

    // 写入输出文件
    const std::filesystem::path out_file = "output/condensate/results_condensate_cpp.csv";
    std::filesystem::create_directories(out_file.parent_path());
    std::vector<double> means = {res.mean};
    std::vector<double> errs = {res.error};
    if (iodata::write_pairs_to_csv(out_file, means, errs)) {
        std::println("结果已落盘至: {}", out_file.string());
    }
    std::println("=================================================");

    return 0;
}

} // namespace

// 1. 单数据集手征凝聚测量任务
int run_condensate_task(std::span<const std::string_view> args) {
    const std::string_view target = args.empty() ? "" : args[0];
    auto info = resolve_condensate_info(target);

    // 支持 CLI 自定义重载物理参数: [target] [ml] [ms] [mres] [zm]
    if (args.size() > 1) { (void)iodata::parse_double(args[1], info.ml); }
    if (args.size() > 2) { (void)iodata::parse_double(args[2], info.ms); }
    if (args.size() > 3) { (void)iodata::parse_double(args[3], info.mres); }
    if (args.size() > 4) { (void)iodata::parse_double(args[4], info.zm); }

    return execute_condensate_measurement(info);
}

// 2. 批量扫描全量手征凝聚测量任务
int run_condensate_all_task(std::span<const std::string_view> args) {
    std::filesystem::path readin_dir = "data/readin";
    if (!args.empty()) {
        readin_dir = args[0];
    } else if (!std::filesystem::exists(readin_dir)) {
        if (std::filesystem::exists("../data/readin")) readin_dir = "../data/readin";
        else if (std::filesystem::exists("/Users/junxiongnie/code/ana/dat/readin")) readin_dir = "/Users/junxiongnie/code/ana/dat/readin";
    }

    std::println(">>> 批量扫描 {} 下所有手征凝聚数据集...", readin_dir.string());
    int processed = 0;

    for (const auto& entry : std::filesystem::directory_iterator(readin_dir)) {
        if (!entry.is_directory()) continue;
        const auto p = entry.path();
        const bool has_meas = std::filesystem::exists(p / "test_condensate") || [p]() {
            for (const auto& sub : std::filesystem::directory_iterator(p)) {
                if (sub.is_directory() && sub.path().filename().string().starts_with("meas.")) return true;
            }
            return false;
        }();

        if (has_meas) {
            auto info = resolve_condensate_info(p.filename().string());
            info.dir_path = p;
            info.name = p.filename().string();
            execute_condensate_measurement(info);
            ++processed;
        }
    }

    std::println("批量测量结束，共完成 {} 组数据集分析。", processed);
    return 0;
}

REGISTER_TASK("condensate", run_condensate_task, "Run chiral condensate measurement: [dataset_path_or_beta] [ml] [ms] [mres] [zm]");
REGISTER_TASK("condensate_all", run_condensate_all_task, "Batch run chiral condensate measurement on all datasets: [readin_dir]");

#ifdef STANDALONE_CONDENSATE
int main(int argc, char* argv[]) {
    std::vector<std::string_view> args;
    for (int i = 1; i < argc; ++i) args.emplace_back(argv[i]);
    return run_condensate_task(args);
}
#endif

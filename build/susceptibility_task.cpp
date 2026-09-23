#include <ParseLQCData/Registry.h>
#include <CondensateAnalysis/CondensateExtractor.h>
#include <IOdata/FastParser.h>
#include <IOdata/FileWriter.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <format>
#include <iostream>
#include <map>
#include <print>
#include <vector>
#include <string>

namespace {

// 标准 Beta 与温度对应表 (Nt=16基准, MeV)
const std::map<std::string, double> TEMP_MAP = {
    {"4.13", 138.0},
    {"4.15", 145.0},
    {"4.17", 153.0},
    {"4.18", 157.0},
    {"4.20", 164.5},
    {"4.23", 176.0},
    {"4.30", 202.0},
    {"4.405", 241.6},
};

struct SusceptibilityMeta {
    std::filesystem::path dir_path;
    std::string name;
    int ns{48};
    int nt{16};
    std::string beta{"4.17"};
    double temp{153.31};
};

// 检查并识别是否为真正包含 PsibarPsi 随机源 XML 文件的有效手征凝聚目录
bool is_valid_condensate_dir(const std::filesystem::path& dir, std::filesystem::path& out_meas_root) {
    if (!std::filesystem::exists(dir) || !std::filesystem::is_directory(dir)) {
        return false;
    }

    // 若包含 Output 且无 meas.*，直接判定为介子关联函数目录，必须跳过
    if (std::filesystem::exists(dir / "Output") && !std::filesystem::exists(dir / "test_condensate")) {
        bool has_meas = false;
        for (const auto& sub : std::filesystem::directory_iterator(dir)) {
            if (sub.is_directory() && sub.path().filename().string().starts_with("meas.")) {
                has_meas = true;
                break;
            }
        }
        if (!has_meas) {
            return false;
        }
    }

    std::vector<std::filesystem::path> cand_roots = {dir, dir / "test_condensate"};
    for (const auto& r : cand_roots) {
        if (!std::filesystem::exists(r) || !std::filesystem::is_directory(r)) continue;
        for (const auto& sub : std::filesystem::directory_iterator(r)) {
            if (sub.is_directory() && sub.path().filename().string().starts_with("meas.")) {
                const auto psibar = sub.path() / "PsibarPsi";
                if (std::filesystem::exists(psibar)) {
                    out_meas_root = r;
                    return true;
                }
            }
        }
    }
    return false;
}

// 快速解析时空几何尺寸与 Beta
SusceptibilityMeta parse_meta(const std::filesystem::path& p) {
    SusceptibilityMeta meta;
    meta.dir_path = p;
    meta.name = p.filename().string();
    std::string name = meta.name;

    // 解析 L 与 T
    if (name.starts_with("L")) {
        size_t t_pos = name.find('T');
        if (t_pos != std::string::npos) {
            double ns_val = 0;
            if (iodata::parse_double(name.substr(1, t_pos - 1), ns_val)) {
                meta.ns = static_cast<int>(ns_val);
            }
            size_t under_pos = name.find_first_of("_b", t_pos + 1);
            if (under_pos != std::string::npos) {
                double nt_val = 0;
                if (iodata::parse_double(name.substr(t_pos + 1, under_pos - (t_pos + 1)), nt_val)) {
                    meta.nt = static_cast<int>(nt_val);
                }
            }
        }
    }

    // 解析 Beta
    size_t beta_pos = name.find("beta");
    if (beta_pos != std::string::npos) {
        size_t start = beta_pos + 4;
        size_t end = name.find_first_of("_ms", start);
        if (end == std::string::npos) end = name.size();
        meta.beta = name.substr(start, end - start);
    }

    // 计算温度
    auto it = TEMP_MAP.find(meta.beta);
    if (it != TEMP_MAP.end()) {
        meta.temp = it->second * (16.0 / static_cast<double>(meta.nt));
    } else {
        double beta_val = 4.17;
        iodata::parse_double(meta.beta, beta_val);
        meta.temp = 153.31 * (16.0 / static_cast<double>(meta.nt));
    }

    return meta;
}

int run_susceptibility_task(std::span<const std::string_view> args) {
    std::filesystem::path readin_dir = "data/readin";
    std::filesystem::path output_dir = "output/condensate";

    if (!args.empty()) {
        readin_dir = args[0];
    }
    if (args.size() > 1) {
        output_dir = args[1];
    }

    if (!std::filesystem::exists(readin_dir)) {
        if (std::filesystem::exists("../data/readin")) readin_dir = "../data/readin";
        else if (std::filesystem::exists("/Users/junxiongnie/code/ana/dat/readin")) readin_dir = "/Users/junxiongnie/code/ana/dat/readin";
    }

    std::println("==========================================================================");
    std::println("  ParseLQCData C++26 高性能手征磁化率抽取工具 (Chiral Susceptibility)");
    std::println("==========================================================================");
    std::println("扫描目录: {}", readin_dir.string());
    std::println("输出目录: {}", output_dir.string());

    std::filesystem::create_directories(output_dir);

    // 扫描所有有效目录，自动跳过强子关联函数目录 (如 48x16b4.13)
    std::vector<std::pair<SusceptibilityMeta, std::filesystem::path>> tasks;

    for (const auto& entry : std::filesystem::directory_iterator(readin_dir)) {
        if (!entry.is_directory()) continue;
        const auto p = entry.path();
        std::filesystem::path meas_root;
        if (is_valid_condensate_dir(p, meas_root)) {
            tasks.push_back({parse_meta(p), meas_root});
        } else {
            if (std::filesystem::exists(p / "Output")) {
                std::println("  [跳过介子目录] {} (包含 Output/ 强子关联函数，无 PsibarPsi 随机源)", p.filename().string());
            }
        }
    }

    if (tasks.empty()) {
        std::println("未在 {} 找到任何包含 PsibarPsi 测量的有效手征凝聚目录！", readin_dir.string());
        return 1;
    }

    // 按温度升序排序
    std::sort(tasks.begin(), tasks.end(), [](const auto& a, const auto& b) {
        return a.first.temp < b.first.temp;
    });

    std::println("\n共找到 {} 组有效手征凝聚格点数据集，开始多线程并发提取...", tasks.size());

    struct Record {
        std::string beta;
        double temp{0.0};
        int ns{48};
        int nt{16};
        double mean_unscaled{0.0};
        double error_unscaled{0.0};
        double mean_vol_scaled{0.0};
        double error_vol_scaled{0.0};
        double mean_scaled{0.0};
        double error_scaled{0.0};
        size_t cfgs{0};
    };

    std::vector<Record> records;
    const auto t_start = std::chrono::steady_clock::now();

    for (const auto& [meta, meas_root] : tasks) {
        std::println(">>> 正在分析 [{}] (Ns={}, Nt={}, β={}, T={:.1f} MeV)...",
                     meta.name, meta.ns, meta.nt, meta.beta, meta.temp);
        const auto res = lqcd::condensate::process_chiral_susceptibility(
            meas_root, meta.ns, meta.nt, meta.temp
        );

        if (res.num_cfgs > 0) {
            std::println("    -> 提取完成: 构型数={}, χ_unscaled={:.4e}±{:.1e}, χ_vol={:.4f}±{:.4f}, χ_scaled={:.4e}±{:.1e} MeV²",
                         res.num_cfgs, res.mean_unscaled, res.error_unscaled,
                         res.mean_vol_scaled, res.error_vol_scaled,
                         res.mean_scaled, res.error_scaled);

            records.push_back({
                meta.beta, meta.temp, meta.ns, meta.nt,
                res.mean_unscaled, res.error_unscaled,
                res.mean_vol_scaled, res.error_vol_scaled,
                res.mean_scaled, res.error_scaled,
                res.num_cfgs
            });
        } else {
            std::println("    [WARN] 未能提取到有效构型数据");
        }
    }

    const auto t_end = std::chrono::steady_clock::now();
    const double elapsed_s = std::chrono::duration<double>(t_end - t_start).count();
    std::println("\n全量计算完成！总耗时: {:.3f} 秒。", elapsed_s);

    // 写出 CSV
    const auto csv_path = output_dir / "results_susceptibility.csv";
    std::ofstream ofs(csv_path);
    if (ofs.is_open()) {
        ofs << "Beta,Temp,Ns,Nt,Mean_unscaled,Error_unscaled,Mean_vol_scaled,Error_vol_scaled,Mean_scaled,Error_scaled,Num_cfgs\n";
        for (const auto& r : records) {
            ofs << std::format("{},{:.2f},{},{},{:.14e},{:.14e},{:.14e},{:.14e},{:.14e},{:.14e},{}\n",
                               r.beta, r.temp, r.ns, r.nt,
                               r.mean_unscaled, r.error_unscaled,
                               r.mean_vol_scaled, r.error_vol_scaled,
                               r.mean_scaled, r.error_scaled,
                               r.cfgs);
        }
        ofs.close();
        std::println("结果已成功写出至: {}", csv_path.string());
    }

    return 0;
}

REGISTER_TASK("susceptibility", run_susceptibility_task, "High-performance C++26 chiral susceptibility calculation: [readin_dir] [output_dir]");


} // namespace

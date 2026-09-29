// 诊断：现有流水线(1线程 vs 8线程)与单遍解析器的数值差异定位
#include <IOdata/DirectoryScanner.h>
#include <IOdata/FileReader.h>
#include <IOdata/FastParser.h>
#include <MesonAnalysis/MesonExtractor.h>
#include <ParseLQCData/MesonPipeline.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

using Meta = std::pair<std::string, std::string>;
static const std::vector<Meta> AV = {
    {"AVector2","DIRX"},{"AVector3","DIRX"},{"AVector3","DIRY"},
    {"AVector1","DIRY"},{"AVector1","DIRZ"},{"AVector2","DIRZ"}};

int main(int argc, char** argv) {
    const std::filesystem::path dir = argv[1];
    const size_t num_lines = 48;
    const auto files = iodata::scan_natural_sorted(dir, "test1_lhadrons_", "_mesons_multi_src");

    std::vector<lqcd::MesonChannelMapping> chs;
    for (const auto& m : AV) chs.push_back({m.first, m.second});

    const auto r1 = lqcd::meson::run_meson_pipeline(dir, chs, 4, num_lines, 1, false);
    const auto r8 = lqcd::meson::run_meson_pipeline(dir, chs, 4, num_lines, 8, false);
    double d = 0.0;
    for (size_t i = 0; i < num_lines; ++i) d = std::max(d, std::fabs(r1.means[i] - r8.means[i]));
    std::printf("[race] AV 1 thread vs 8 threads: max|dmean| = %.6e\n", d);
    std::printf("  t=0  1T=%.12e  8T=%.12e\n", r1.means[0], r8.means[0]);
    std::printf("  t=5  1T=%.12e  8T=%.12e\n", r1.means[5], r8.means[5]);

    // 逐文件逐 meta：现有 extractor vs 手工统计该 meta 找到的 block 数
    for (size_t fi = 0; fi < std::min<size_t>(files.size(), 3); ++fi) {
        const std::string content = iodata::read_file_to_string(files[fi]);
        std::printf("--- file %zu (%s) ---\n", fi, files[fi].filename().c_str());
        for (const auto& m : AV) {
            const std::string tag = "--- " + m.first + " to " + m.first + " --- spatial:" + m.second + " ---";
            // 统计该 tag 在文件中出现次数
            size_t occ = 0, pos = 0;
            while ((pos = content.find(tag, pos)) != std::string::npos) { ++occ; pos += tag.size(); }
            const auto blk = lqcd::meson::extract_single_file_averaged_block(content, m.second, tag, num_lines);
            double sum = 0.0;
            for (double v : blk) sum += v;
            std::printf("   %-9s %-5s tag_occurrences=%2zu blk_sum=%.6e blk[0]=%.6e\n",
                        m.first.c_str(), m.second.c_str(), occ, sum, blk[0]);
        }
    }
    return 0;
}

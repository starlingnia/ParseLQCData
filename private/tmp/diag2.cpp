// 诊断 2：逐 (文件, meta) 对比现有 extractor 与单遍解析器
#include <IOdata/DirectoryScanner.h>
#include <IOdata/FileReader.h>
#include <IOdata/FastParser.h>
#include <MesonAnalysis/MesonExtractor.h>

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

struct FlatMeta { std::string type, dir; int axis; };

static std::vector<std::vector<double>> one_pass(std::string_view content,
                                                const std::vector<FlatMeta>& metas,
                                                size_t num_lines, size_t* counts_out) {
    std::vector<std::vector<double>> sums(metas.size(), std::vector<double>(num_lines, 0.0));
    std::vector<size_t> counts(metas.size(), 0);
    std::vector<double> blk(num_lines, 0.0);
    long cur = -1; int cur_shift = 0; size_t cur_fill = 0;

    iodata::for_each_line(content, [&](std::string_view line) {
        if (!line.empty() && line[0] == '-') {
            cur = -1; cur_fill = 0;
            const size_t sp = line.find(' ', 4);
            if (sp == std::string_view::npos) return;
            const std::string_view type = line.substr(4, sp - 4);
            const size_t spos = line.find("spatial:");
            if (spos == std::string_view::npos) return;
            const std::string_view dir = line.substr(spos + 8, 4);
            for (size_t m = 0; m < metas.size(); ++m) {
                if (metas[m].type == type && metas[m].dir == dir) { cur = (long)m; break; }
            }
            if (cur < 0) return;
            const size_t last = line.rfind("--- ");
            if (last == std::string_view::npos) { cur = -1; return; }
            std::string_view coord = line.substr(last + 4);
            int nums[4] = {0,0,0,0};
            size_t idx = 0, start = 0;
            while (idx < 4 && start < coord.size()) {
                const size_t slash = coord.find('/', start);
                const std::string_view tok = (slash == std::string_view::npos) ? coord.substr(start)
                                                                              : coord.substr(start, slash - start);
                int v = 0;
                if (iodata::parse_int(tok, v)) nums[idx++] = v; else break;
                if (slash == std::string_view::npos) break;
                start = slash + 1;
            }
            cur_shift = nums[metas[(size_t)cur].axis];
            return;
        }
        if (cur < 0 || cur_fill >= num_lines) return;
        size_t ti = 0; double rv = 0.0; bool got = false;
        iodata::for_each_token(line, [&](std::string_view tok) {
            if (ti == 1) { if (iodata::parse_double(tok, rv)) got = true; }
            ++ti;
        });
        if (!got) return;
        blk[cur_fill++] = rv;
        if (cur_fill == num_lines) {
            auto& s = sums[(size_t)cur];
            const int n = (int)num_lines;
            const int sh = ((cur_shift % n) + n) % n;
            for (int i = 0; i < n; ++i) s[(size_t)i] += blk[(size_t)((i + sh) % n)];
            ++counts[(size_t)cur];
            cur = -1; cur_fill = 0;
        }
    });
    for (size_t m = 0; m < metas.size(); ++m) {
        if (counts_out) counts_out[m] = counts[m];
        if (counts[m] == 0) continue;
        const double inv = 1.0 / (double)counts[m];
        for (auto& v : sums[m]) v *= inv;
    }
    return sums;
}

int main(int argc, char** argv) {
    const std::filesystem::path dir = argv[1];
    const size_t num_files = (argc > 2) ? std::stoul(argv[2]) : 20;
    const size_t num_lines = 48;
    const auto files = iodata::scan_natural_sorted(dir, "test1_lhadrons_", "_mesons_multi_src");

    std::vector<FlatMeta> metas;
    for (const auto& m : AV) metas.push_back({m.first, m.second, m.second == "DIRX" ? 0 : (m.second == "DIRY" ? 1 : 2)});

    double worst = 0.0; size_t worst_file = 0, worst_meta = 0;
    for (size_t fi = 0; fi < std::min(files.size(), num_files); ++fi) {
        const std::string content = iodata::read_file_to_string(files[fi]);
        std::vector<size_t> counts(metas.size(), 0);
        const auto op = one_pass(content, metas, num_lines, counts.data());
        for (size_t mi = 0; mi < metas.size(); ++mi) {
            const std::string tag = "--- " + AV[mi].first + " to " + AV[mi].first + " --- spatial:" + AV[mi].second + " ---";
            const auto ex = lqcd::meson::extract_single_file_averaged_block(content, AV[mi].second, tag, num_lines);
            double d = 0.0;
            for (size_t i = 0; i < num_lines; ++i) d = std::max(d, std::fabs(ex[i] - op[mi][i]));
            if (d > worst) { worst = d; worst_file = fi; worst_meta = mi; }
            if (d > 1e-12) {
                std::printf("file=%zu meta=%s/%s  blocks(one_pass)=%zu  max|diff|=%.6e\n",
                            fi, AV[mi].first.c_str(), AV[mi].second.c_str(), counts[mi], d);
                std::printf("    extractor[0..3]=");
                for (int i = 0; i < 4; ++i) std::printf(" %.12e", ex[i]);
                std::printf("\n    onepass  [0..3]=");
                for (int i = 0; i < 4; ++i) std::printf(" %.12e", op[mi][i]);
                std::printf("\n");
            }
        }
    }
    std::printf("扫描 %zu 个文件完成, 最大偏差 = %.6e (file=%zu meta=%s/%s)\n",
                std::min(files.size(), num_files), worst, worst_file,
                AV[worst_meta].first.c_str(), AV[worst_meta].second.c_str());
    return 0;
}

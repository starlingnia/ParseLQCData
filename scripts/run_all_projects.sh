#!/usr/bin/env bash
# ==============================================================================
# scripts/run_all_projects.sh
# ------------------------------------------------------------------------------
# 项目级综合调度编排脚本 (Master Project Runner):
# 支持独立执行或批量执行各个不同项目的分析流水线:
#   1. LCP 项目 (常物理线: 手征凝聚抽取与重整化、手征磁化率计算与标度、LCP 标准导出、绘图)
#   2. Meson 测量项目 (介子强子关联函数抽取、有效质量求解、平台拟合、能谱绘图)
#   3. CCOR 项目 (对称性破缺质量差分析与对照)
#
# 用法:
#   bash scripts/run_all_projects.sh --all            # 执行所有项目
#   bash scripts/run_all_projects.sh --lcp            # 仅执行 LCP 项目
#   bash scripts/run_all_projects.sh --meson          # 仅执行介子测量项目
#   bash scripts/run_all_projects.sh --ccor           # 仅执行 CCOR 项目
#   bash scripts/run_all_projects.sh --lcp --plot     # 执行 LCP 项目并绘图
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

RUN_MESON=false
RUN_LCP=false
RUN_CCOR=false
PLOT_FLAG=""

if [[ $# -eq 0 ]]; then
  RUN_LCP=true
  RUN_MESON=true
fi

while [[ $# -gt 0 ]]; do
  case $1 in
    --all)
      RUN_MESON=true
      RUN_LCP=true
      RUN_CCOR=true
      shift
      ;;
    --meson)
      RUN_MESON=true
      shift
      ;;
    --lcp)
      RUN_LCP=true
      shift
      ;;
    --ccor)
      RUN_CCOR=true
      shift
      ;;
    --plot)
      PLOT_FLAG="--plot"
      shift
      ;;
    -h|--help)
      echo "用法: bash scripts/run_all_projects.sh [选项]"
      echo "选项:"
      echo "  --all             执行所有项目流水线 (LCP, Meson, CCOR)"
      echo "  --lcp             执行常物理线 (LCP) 项目流水线"
      echo "  --meson           执行介子关联函数测量项目流水线"
      echo "  --ccor            执行 CCOR 对称性破缺项目流水线"
      echo "  --plot            各项目执行完成后自动生成图表"
      echo "  -h, --help        显示帮助信息"
      exit 0
      ;;
  esac
done

echo "=========================================================================="
echo ">>> [MASTER RUNNER] 格点 QCD 数据分析项目总调度器启动"
echo "  - 执行 LCP 项目:   ${RUN_LCP}"
echo "  - 执行 介子项目:   ${RUN_MESON}"
echo "  - 执行 CCOR 项目:  ${RUN_CCOR}"
echo "=========================================================================="

T_START=$(date +%s)

# 1. 调度 LCP 项目
if [ "$RUN_LCP" = true ]; then
  echo ""
  echo ">>>>>>>>>>>>>>>>>>>>>> [1] 执行常物理线 (LCP) 项目 <<<<<<<<<<<<<<<<<<<<<<"
  bash scripts/run_lcp_project.sh $PLOT_FLAG
fi

# 2. 调度介子测量项目
if [ "$RUN_MESON" = true ]; then
  echo ""
  echo ">>>>>>>>>>>>>>>>>>>>>> [2] 执行介子强子关联函数测量项目 <<<<<<<<<<<<<<<<<<<<<<"
  bash scripts/run_meson_project.sh $PLOT_FLAG
fi

# 3. 调度 CCOR 项目
if [ "$RUN_CCOR" = true ]; then
  echo ""
  echo ">>>>>>>>>>>>>>>>>>>>>> [3] 执行 CCOR 对称性破缺项目 <<<<<<<<<<<<<<<<<<<<<<"
  bash scripts/run_ccor_project.sh
fi

T_END=$(date +%s)
TOTAL_ELAPSED=$((T_END - T_START))

echo ""
echo "=========================================================================="
echo ">>> [MASTER RUNNER] 所选项目流水线全部顺利完成！总耗时: ${TOTAL_ELAPSED}s"
echo "=========================================================================="

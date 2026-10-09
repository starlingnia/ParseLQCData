#!/usr/bin/env bash
# ==============================================================================
# scripts/run_lcp_project.sh
# ------------------------------------------------------------------------------
# 常物理线项目 (LCP Project) 完整流水线 Shell 总控编排脚本:
# 依次调度执行 4 个原子化独立小任务，将手征凝聚与磁化率分析串联为一个统一项目任务:
#   - Task 1: 01_extract_condensate.py (C++ 核心库光/奇夸克手征凝聚抽取与 mres 减除)
#   - Task 2: 02_compute_susceptibility.py (C++ 无偏二次估计器手征磁化率计算与标度)
#   - Task 3: 03_export_lcp.py (Zm^2 重整化并导出 LCP 规范标准文本产物至 output/LCP/)
#   - Task 4: 04_plot_lcp.py (生成手征凝聚与磁化率热力学相变曲线)
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

# 参数默认值
DO_PLOT=false
EXTRA_ARGS=()
READIN_DIR=""
ENSEMBLE=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --plot)
      DO_PLOT=true
      shift
      ;;
    --readin-dir)
      READIN_DIR="$2"
      EXTRA_ARGS+=("--readin-dir" "$2")
      shift 2
      ;;
    --ensemble)
      ENSEMBLE="$2"
      EXTRA_ARGS+=("--ensemble" "$2")
      shift 2
      ;;
    -h|--help)
      echo "用法: bash scripts/run_lcp_project.sh [选项]"
      echo "选项:"
      echo "  --readin-dir <dir>  指定原始构型根目录 (默认 data/readin)"
      echo "  --ensemble <name>   指定单独分析某个格点系综"
      echo "  --plot              执行完成后自动绘制热力学图谱"
      echo "  -h, --help          显示帮助信息"
      exit 0
      ;;
    *)
      EXTRA_ARGS+=("$1")
      shift
      ;;
  esac
done

echo "=========================================================================="
echo ">>> [PROJECT: LCP] 启动常物理线 (LCP) 项目流水线总控编排"
echo "  - 包含任务: 手征凝聚抽取、手征磁化率计算、LCP 格式导出、热力学绘图"
echo "  - 绘图开启: ${DO_PLOT}"
echo "=========================================================================="

START_TIME=$(date +%s)

# Step 1: 手征凝聚抽取与残余质量减除
echo ""
echo ">>> [LCP STEP 1/4] 调度原子任务: 手征凝聚抽取与残余质量减除 (C++ CondensateExtractor)..."
uv run python scripts/tasks/lcp/01_extract_condensate.py "${EXTRA_ARGS[@]}"

# Step 2: 手征磁化率无偏二次估计与体积/温度标度
echo ""
echo ">>> [LCP STEP 2/4] 调度原子任务: 手征磁化率无偏二次交叉估计与标度..."
SUSC_ARGS=()
if [ -n "$READIN_DIR" ]; then
  SUSC_ARGS+=("--readin-dir" "$READIN_DIR")
fi
uv run python scripts/tasks/lcp/02_compute_susceptibility.py "${SUSC_ARGS[@]}"

# Step 3: LCP 规范化导出与 Zm^2 重整化
echo ""
echo ">>> [LCP STEP 3/4] 调度原子任务: LCP 标准产物导出与 Zm^2 质量重整化..."
uv run python scripts/tasks/lcp/03_export_lcp.py

# Step 4: 综合热力学相变曲线绘图
if [ "$DO_PLOT" = true ]; then
  echo ""
  echo ">>> [LCP STEP 4/4] 调度原子任务: 出版级热力学图谱绘制..."
  uv run python scripts/tasks/lcp/04_plot_lcp.py
else
  echo ""
  echo ">>> [LCP STEP 4/4] 跳过绘图步骤 (若需绘图请追加 --plot)"
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "=========================================================================="
echo ">>> [PROJECT: LCP] 常物理线项目全流水线执行完毕！总耗时: ${ELAPSED}s"
echo "=========================================================================="

#!/usr/bin/env bash
# ==============================================================================
# scripts/run_meson_project.sh
# ------------------------------------------------------------------------------
# 介子测量项目 (Meson Project) 完整流水线 Shell 总控编排脚本:
# 依次调度执行 4 个原子化独立小任务，将强子关联函数端到端分析串联为一个统一项目任务:
#   - Task 1: 01_extract_correlators.py (C++ 核心库空间关联函数抽取与折叠)
#   - Task 2: 02_compute_effective_mass.py (内存级/流式有效质量方程求解)
#   - Task 3: 03_fit_plateau.py (Jackknife 平台 cosh 最小二乘拟合与汇总)
#   - Task 4: 04_plot_meson.py (出版级物理图谱绘制)
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

# 参数默认值
SOURCE="all"
BETA="all"
DO_PLOT=false
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    --source)
      SOURCE="$2"
      shift 2
      ;;
    --beta)
      BETA="$2"
      shift 2
      ;;
    --plot)
      DO_PLOT=true
      shift
      ;;
    -h|--help)
      echo "用法: bash scripts/run_meson_project.sh [选项]"
      echo "选项:"
      echo "  --source [all|multi|single]  数据源模式 (默认 all)"
      echo "  --beta [all|17|...]          指定 Beta 编号 (默认 all)"
      echo "  --plot                       执行完成后自动绘制图谱"
      echo "  -h, --help                   显示帮助信息"
      exit 0
      ;;
    *)
      EXTRA_ARGS+=("$1")
      shift
      ;;
  esac
done

echo "=========================================================================="
echo ">>> [PROJECT: MESON] 启动介子测量项目流水线总控编排"
echo "  - 数据源模式: ${SOURCE}"
echo "  - Beta 设定:  ${BETA}"
echo "  - 绘图开启:   ${DO_PLOT}"
echo "=========================================================================="

START_TIME=$(date +%s)

# Step 1: 关联函数抽取
echo ""
echo ">>> [MESON STEP 1/4] 调度原子任务: 关联函数抽取与折叠 (C++ MesonExtractor)..."
uv run python scripts/tasks/meson/01_extract_correlators.py --source "$SOURCE" --beta "$BETA" "${EXTRA_ARGS[@]}"

# Step 2: 有效质量计算
echo ""
echo ">>> [MESON STEP 2/4] 调度原子任务: 有效质量超越方程求解..."
uv run python scripts/tasks/meson/02_compute_effective_mass.py --source "$SOURCE" --beta "$BETA"

# Step 3: 平台拟合
echo ""
echo ">>> [MESON STEP 3/4] 调度原子任务: 基态平台 cosh 拟合与 Jackknife 统计..."
uv run python scripts/tasks/meson/03_fit_plateau.py --source "$SOURCE" --beta "$BETA"

# Step 4: 绘图与可视化
if [ "$DO_PLOT" = true ]; then
  echo ""
  echo ">>> [MESON STEP 4/4] 调度原子任务: 出版级物理图谱绘制..."
  uv run python scripts/tasks/meson/04_plot_meson.py
else
  echo ""
  echo ">>> [MESON STEP 4/4] 跳过绘图步骤 (若需绘图请追加 --plot)"
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "=========================================================================="
echo ">>> [PROJECT: MESON] 介子测量项目全流水线执行完毕！总耗时: ${ELAPSED}s"
echo "=========================================================================="

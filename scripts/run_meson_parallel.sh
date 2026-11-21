#!/usr/bin/env bash
# ==============================================================================
# scripts/run_meson_parallel.sh
# ------------------------------------------------------------------------------
# 介子测量项目 (Meson Project) 系综级多进程并发流水线 Shell 脚本:
# 针对不同 Beta 系综 (13, 15, 17, 18, 20, 23, 30) 启动后台独立并发进程：
#   每个系综在独立子进程中按序执行:
#     1) 01_extract_correlators.py (关联函数抽取折叠)
#     2) 02_compute_effective_mass.py (有效质量求解)
#     3) 03_fit_plateau.py (平台拟合)
#   所有系综并发完成后，由主进程执行:
#     4) 03_fit_plateau.py --aggregate-only (统一合并生成总表)
#     5) 04_plot_meson.py (可选绘图)
# ==============================================================================

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

PYTHON_BIN="./.venv/bin/python"
if [ ! -f "$PYTHON_BIN" ]; then
  PYTHON_BIN="python3"
fi

SOURCE="all"
BETAS_INPUT="all"
DO_PLOT=false
MAX_JOBS=7

while [[ $# -gt 0 ]]; do
  case $1 in
    --source)
      SOURCE="$2"
      shift 2
      ;;
    --beta)
      BETAS_INPUT="$2"
      shift 2
      ;;
    --jobs|-j)
      MAX_JOBS="$2"
      shift 2
      ;;
    --plot)
      DO_PLOT=true
      shift
      ;;
    -h|--help)
      echo "用法: bash scripts/run_meson_parallel.sh [选项]"
      echo "选项:"
      echo "  --source [all|multi|single]  数据源模式 (默认: all)"
      echo "  --beta [all|13,15,17...]     指定 Beta 列表，逗号分隔或 'all'"
      echo "  --jobs, -j <num>             最大并发任务数 (默认: 7)"
      echo "  --plot                       执行完成后自动绘制图谱"
      echo "  -h, --help                   显示帮助信息"
      exit 0
      ;;
    *)
      echo "未知参数: $1"
      exit 1
      ;;
  esac
done

# 解析要跑的 Beta 列表
if [ "$BETAS_INPUT" = "all" ]; then
  BETAS=("13" "15" "17" "18" "20" "23" "30")
else
  IFS=',' read -r -a BETAS <<< "$BETAS_INPUT"
fi

echo "=========================================================================="
echo ">>> [PARALLEL MESON] 启动介子测量项目多系综并发调度"
echo "  - 数据源模式:   ${SOURCE}"
echo "  - 覆盖系综Beta: ${BETAS[*]}"
echo "  - 最大并发数:   ${MAX_JOBS}"
echo "  - 绘图开启:     ${DO_PLOT}"
echo "=========================================================================="

T_START=$(date +%s)

# 创建临时日志目录
LOG_DIR="${PROJECT_ROOT}/output/logs/meson"
mkdir -p "$LOG_DIR"

run_single_beta() {
  local b="$1"
  local log_file="${LOG_DIR}/beta_${b}.log"
  echo "  [LAUNCH] 启动 Beta 4.${b} 处理子进程 (日志: output/logs/meson/beta_${b}.log)..."

  {
    echo "=== [Beta 4.${b}] STEP 1: 抽取关联函数 ==="
    "$PYTHON_BIN" scripts/tasks/meson/01_extract_correlators.py --source "$SOURCE" --beta "$b"

    echo "=== [Beta 4.${b}] STEP 2: 求解有效质量 ==="
    "$PYTHON_BIN" scripts/tasks/meson/02_compute_effective_mass.py --source "$SOURCE" --beta "$b"

    echo "=== [Beta 4.${b}] STEP 3: 拟合基态平台 ==="
    "$PYTHON_BIN" scripts/tasks/meson/03_fit_plateau.py --source "$SOURCE" --beta "$b"
  } > "$log_file" 2>&1

  local status=$?
  if [ $status -eq 0 ]; then
    echo "  [DONE] Beta 4.${b} 完成！"
  else
    echo "  [FAIL] Beta 4.${b} 失败 (查看 ${log_file})"
    return $status
  fi
}

# 并发调度各 Beta
pids=()
for b in "${BETAS[@]}"; do
  # 统一格式，如 4.17 -> 17
  clean_b="${b//4./}"
  clean_b="${clean_b//b/}"

  run_single_beta "$clean_b" &
  pids+=($!)

  # 控制最大并发数
  if [ "${#pids[@]}" -ge "$MAX_JOBS" ]; then
    wait -n
  fi
done

# 等待全部 Beta 任务完成
has_error=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    has_error=1
  fi
done

if [ $has_error -ne 0 ]; then
  echo "[ERROR] 部分 Beta 系综任务执行失败，请检查 ${LOG_DIR}/ 下的日志！"
  exit 1
fi

echo ""
echo ">>> [PARALLEL MESON] 各系综计算完毕，开始汇总全量拟合大表..."
"$PYTHON_BIN" scripts/tasks/meson/03_fit_plateau.py --source "$SOURCE" --aggregate-only

if [ "$DO_PLOT" = true ]; then
  echo ""
  echo ">>> [PARALLEL MESON] 正在绘制介子物理图谱..."
  "$PYTHON_BIN" scripts/tasks/meson/04_plot_meson.py
fi

T_END=$(date +%s)
echo "=========================================================================="
echo ">>> [PARALLEL MESON] 介子全系综并发流水线顺利完成！总耗时: $((T_END - T_START))s"
echo "=========================================================================="

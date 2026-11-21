#!/usr/bin/env bash
# ==============================================================================
# scripts/run_all_parallel.sh
# ------------------------------------------------------------------------------
# 项目级端到端多进程并发总调度器 (Master Parallel Project Runner):
# 并发调度不同分析流水线，充分压榨多核算力:
#   1. 常物理线项目 (LCP): 凝聚抽取 -> 磁化率无偏二次估计 -> LCP导出 -> Padé拟合
#   2. 介子测量项目 (Meson): 各 Beta 系综后台独立并行 (抽取 -> 有效质量 -> 平台拟合)
#   3. CCOR 项目 (可选): 对称性破缺质量差复刻
#
# 用法:
#   bash scripts/run_all_parallel.sh                 # 默认同时并发 LCP 与 Meson
#   bash scripts/run_all_parallel.sh --plot          # 并发执行并在完成后出图
#   bash scripts/run_all_parallel.sh --all           # 并发执行全部 (LCP + Meson + CCOR)
#   bash scripts/run_all_parallel.sh --lcp           # 仅执行 LCP
#   bash scripts/run_all_parallel.sh --meson         # 仅执行 Meson (系综级并发)
# ==============================================================================

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

RUN_LCP=false
RUN_MESON=false
RUN_CCOR=false
PLOT_FLAG=""
EXTRA_ARGS=()

if [[ $# -eq 0 ]]; then
  RUN_LCP=true
  RUN_MESON=true
fi

while [[ $# -gt 0 ]]; do
  case $1 in
    --all)
      RUN_LCP=true
      RUN_MESON=true
      RUN_CCOR=true
      shift
      ;;
    --lcp)
      RUN_LCP=true
      shift
      ;;
    --meson)
      RUN_MESON=true
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
      echo "用法: bash scripts/run_all_parallel.sh [选项]"
      echo "选项:"
      echo "  --all             并发执行所有项目 (LCP + Meson + CCOR)"
      echo "  --lcp             仅执行 LCP 常物理线项目"
      echo "  --meson           仅执行 Meson 介子测量项目 (多系综并发)"
      echo "  --ccor            执行 CCOR 项目"
      echo "  --plot            执行完成后自动绘制图谱"
      echo "  -h, --help        显示帮助信息"
      exit 0
      ;;
    *)
      EXTRA_ARGS+=("$1")
      shift
      ;;
  esac
done

# 如果指定了 --ccor 但没指定其他，则只跑 CCOR；若都未指定则默认跑 LCP + Meson
if [ "$RUN_LCP" = false ] && [ "$RUN_MESON" = false ] && [ "$RUN_CCOR" = false ]; then
  RUN_LCP=true
  RUN_MESON=true
fi

LOG_DIR="${PROJECT_ROOT}/output/logs"
mkdir -p "$LOG_DIR"

echo "=========================================================================="
echo ">>> [MASTER PARALLEL RUNNER] 启动全项目流水线并发总调度"
echo "  - 常物理线 (LCP) 项目:    ${RUN_LCP}"
echo "  - 介子测量 (Meson) 项目:   ${RUN_MESON} (各Beta系综后台并发)"
echo "  - 对称性破缺 (CCOR) 项目:  ${RUN_CCOR}"
echo "  - 绘图标志:                ${PLOT_FLAG:-关闭}"
echo "=========================================================================="

T_START=$(date +%s)
pids=()
job_names=()

# 1. 后台并发启动 LCP 项目流水线
if [ "$RUN_LCP" = true ]; then
  echo "  [START] 后台并发启动: LCP 常物理线流水线 (日志: output/logs/project_lcp.log)..."
  (
    bash scripts/run_lcp_project.sh $PLOT_FLAG "${EXTRA_ARGS[@]}"
  ) > "${LOG_DIR}/project_lcp.log" 2>&1 &
  pids+=($!)
  job_names+=("LCP")
fi

# 2. 后台并发启动 Meson 项目流水线 (多系综并行)
if [ "$RUN_MESON" = true ]; then
  echo "  [START] 后台并发启动: Meson 介子测量流水线 (日志: output/logs/project_meson.log)..."
  (
    bash scripts/run_meson_parallel.sh $PLOT_FLAG "${EXTRA_ARGS[@]}"
  ) > "${LOG_DIR}/project_meson.log" 2>&1 &
  pids+=($!)
  job_names+=("Meson")
fi

# 3. 后台并发启动 CCOR 项目流水线 (若启用)
if [ "$RUN_CCOR" = true ]; then
  echo "  [START] 后台并发启动: CCOR 对称性破缺流水线 (日志: output/logs/project_ccor.log)..."
  (
    bash scripts/run_ccor_project.sh
  ) > "${LOG_DIR}/project_ccor.log" 2>&1 &
  pids+=($!)
  job_names+=("CCOR")
fi

echo ""
echo ">>> 所有子项目已进入后台并发执行，正在监控进程状态..."

failed_jobs=()
for i in "${!pids[@]}"; do
  pid="${pids[$i]}"
  name="${job_names[$i]}"
  if wait "$pid"; then
    echo "  [FINISHED] 项目 [${name}] 执行成功！"
  else
    echo "  [FAILED]   项目 [${name}] 执行失败！请检查日志: output/logs/project_${name,,}.log"
    failed_jobs+=("$name")
  fi
done

T_END=$(date +%s)
TOTAL_TIME=$((T_END - T_START))

echo ""
echo "=========================================================================="
if [ "${#failed_jobs[@]}" -eq 0 ]; then
  echo ">>> [MASTER PARALLEL RUNNER] 所有并发项目全部成功执行完毕！总耗时: ${TOTAL_TIME}s"
else
  echo ">>> [MASTER PARALLEL RUNNER] 部分项目失败: ${failed_jobs[*]}，总耗时: ${TOTAL_TIME}s"
  exit 1
fi
echo "=========================================================================="

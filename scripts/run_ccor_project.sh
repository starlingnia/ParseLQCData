#!/usr/bin/env bash
# ==============================================================================
# scripts/run_ccor_project.sh
# ------------------------------------------------------------------------------
# CCOR 对称性破缺质量差项目 (CCOR Project) Shell 总控编排脚本:
# 调度执行介子 block 关联函数抽取、逐样本拟合与 4 组质量差演化分析
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "=========================================================================="
echo ">>> [PROJECT: CCOR] 启动 CCOR 对称性破缺流水线总控编排"
echo "=========================================================================="

START_TIME=$(date +%s)

# Step 1: 关联函数抽取
echo ""
echo ">>> [CCOR STEP 1/3] 调度原子任务: block 关联函数抽取与折叠..."
uv run python scripts/tasks/ccor/01_extract_ccor.py "$@"

# Step 2: 逐样本平台拟合
echo ""
echo ">>> [CCOR STEP 2/3] 调度原子任务: 逐 Jackknife 样本 cosh 平台拟合..."
uv run python scripts/tasks/ccor/02_fit_ccor.py "$@"

# Step 3: 对称性质量差计算与比对
echo ""
echo ">>> [CCOR STEP 3/3] 调度原子任务: 4 组手征/自旋手征质量差 ΔM 计算与对照..."
uv run python scripts/tasks/ccor/03_compute_delta_m.py "$@"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "=========================================================================="
echo ">>> [PROJECT: CCOR] CCOR 项目全流水线执行完毕！总耗时: ${ELAPSED}s"
echo "=========================================================================="

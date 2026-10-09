#!/usr/bin/env bash
# scripts/lcp/run_lcp_all.sh
# 一键运行常物理线 (LCP) 全流程分析并出图

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
cd "$PROJECT_ROOT"

echo "=========================================================================="
echo "🚀 启动 ParseLQCData 常物理线 (LCP) 全流程分析流水线"
echo "=========================================================================="

uv run python scripts/lcp/run_lcp_all.py "$@"

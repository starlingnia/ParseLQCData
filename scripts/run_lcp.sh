#!/usr/bin/env bash
# ==============================================================================
# scripts/run_lcp.sh
# 快捷入口：执行常物理线 (LCP) 项目完整分析流水线
# ==============================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

exec bash "$SCRIPT_DIR/run_lcp_project.sh" "$@"

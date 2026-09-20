#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "[MESON-MULTI] Running multi-source meson data pipeline..."
uv run python scripts/reproduce_meson_multi.py
echo "[MESON-MULTI] Multi-source pipeline completed."

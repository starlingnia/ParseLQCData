#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "[MESON] Running multi-source pipeline..."
bash scripts/run_meson_multi.sh

echo "[MESON] Running single-source pipeline..."
bash scripts/run_meson_single.sh

echo "[MESON] Generating meson plots..."
uv run python scripts/plot_meson.py

echo "[MESON] Meson pipeline completed."

#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "[MESON-SINGLE] Running single-source meson data pipeline..."
uv run python scripts/reproduce_meson_single.py
echo "[MESON-SINGLE] Single-source pipeline completed."

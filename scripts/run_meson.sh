#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "[MESON] Running end-to-end meson pipeline (multi + single source)..."
uv run python scripts/run_meson.py --source all --plot

echo "[MESON] Meson pipeline completed successfully."

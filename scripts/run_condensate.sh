#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "[CONDENSATE] Running chiral condensate pipeline..."
uv run python scripts/reproduce_condensate.py

echo "[CONDENSATE] Generating chiral condensate plots..."
uv run python scripts/plot_condensate.py

echo "[CONDENSATE] Chiral condensate pipeline completed."

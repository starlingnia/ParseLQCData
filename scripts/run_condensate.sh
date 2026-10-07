#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "[CONDENSATE] Running end-to-end chiral condensate pipeline..."
uv run python scripts/run_condensate.py --plot

echo "[CONDENSATE] Chiral condensate pipeline completed successfully."

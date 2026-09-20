#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "[TEST] Running numerical regression tests vs ana..."
uv run python tests/compare_with_ana.py

echo "[TEST] Running concurrency and stress benchmarks..."
uv run python tests/benchmark_stress_test.py

echo "[TEST] Running binsize autocorrelation and error saturation test..."
uv run python tests/test_binsize_autocorr.py

echo "[TEST] All tests completed."

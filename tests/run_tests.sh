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

echo "[TEST] Running b4.17 Nt-scan task tests (correlators + meson masses)..."
uv run python tests/test_meson_b417_scan.py

echo "[TEST] Running ccor flow reproduction tests (vs ana/dat/ccor)..."
uv run python tests/test_ccor_flow.py

echo "[TEST] Running chiral susceptibility extraction and Zm scaling tests..."
uv run python tests/test_susceptibility.py

echo "[TEST] Running chiral susceptibility Padé rational approximation and pole tests..."
uv run python tests/test_pade_fitter.py

echo "[TEST] All tests completed."

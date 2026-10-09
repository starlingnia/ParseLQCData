#!/usr/bin/env bash
# scripts/run_meson_b417_nt_scan.sh
# ------------------------------------------------------------------------------
# 运行 beta=4.17 有限温度扫描任务 (关联函数 + 介子质量)
#   Ns^3 x Nt: 32^3x12 / 32^3x14 / 32^3x16 / 36^3x18 / 40^3x16 / 48^3x18
#   ms = 0.040, ml = 0.0020 / 0.0035 / 0.0070 / 0.0120
# 结果: output/meson_scan/b4.17/{multisrc,singlesrc}/...
# ------------------------------------------------------------------------------
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
cd "$PROJECT_ROOT"

PYTHON_BIN="${PYTHON_BIN:-$PROJECT_ROOT/.venv/bin/python}"
if [ ! -x "$PYTHON_BIN" ]; then
  PYTHON_BIN="python3"
fi

echo "[MESON-B417-SCAN] 开始: 关联函数 + 介子质量 (Ns^3 x Nt 扫描)"
"$PYTHON_BIN" scripts/b417_scan/reproduce_meson_b417_nt_scan.py "$@"
echo "[MESON-B417-SCAN] 完成, 结果位于 output/meson_scan/b4.17/"

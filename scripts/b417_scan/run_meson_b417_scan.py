#!/usr/bin/env python3
"""
scripts/run_meson_b417_scan.py
--------------------------------------------------------------------------------
Beta=4.17 有限温度介子扫描统一入口 (包含关联函数抽取、有效质量求解、平台拟合与出版级画图)
--------------------------------------------------------------------------------
"""

import argparse
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def main():
    parser = argparse.ArgumentParser(description="Beta=4.17 有限温度扫描统一分析入口")
    parser.add_argument("--source", choices=["multi", "single", "both"], default="both", help="数据源 (默认两者均跑)")
    parser.add_argument("--plot", action="store_true", help="分析完成后自动绘制 4 组质量的温度演化图")
    parser.add_argument("--workers", type=int, default=4, help="并发 worker 进程数")
    parser.add_argument("--force", action="store_true", help="强制覆盖已存在的结果")
    args, unknown = parser.parse_known_args()

    python_bin = sys.executable

    # 1. 运行核心扫描分析
    scan_script = PROJECT_ROOT / "scripts" / "reproduce_meson_b417_nt_scan.py"
    cmd = [python_bin, str(scan_script), "--workers", str(args.workers)]
    if args.source != "both":
        cmd.extend(["--source", args.source])
    if args.force:
        cmd.append("--force")
    cmd.extend(unknown)

    print(f"[INFO] 正在执行 b4.17 全量扫描: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # 2. 若指定 --plot，调用 measure_meson_b417.py 生成 4 幅演化图
    if args.plot:
        plot_script = PROJECT_ROOT / "scripts" / "measure_meson_b417.py"
        print(f"[INFO] 正在生成出版级介子质量随温度演化图...")
        subprocess.run([python_bin, str(plot_script)], check=True)

    print("[OK] b4.17 有限温度扫描与画图流程执行完毕！")


if __name__ == "__main__":
    main()

"""
端到端精度验证与性能评测自动化测试脚本
------------------------------------------
- 针对 ~/code/ana/dat/readin/48x16b4.17/Output 真实物理数据进行端到端并发计算
- 将计算产物保存到本项目独立目录 output/b4.17/save_AV.csv (绝不污染原项目)
- 与 ~/code/ana/dat/func/pickdata/b4.17/save_AV.csv 进行 100% 逐行浮点机器精度比对
"""

import os
import sys
import time
from pathlib import Path
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from src.parselqcdata.meson_pipeline import process_channel

def main():
    print("=================================================================")
    print("启动 ParseLQCData 高性能端到端验证测试 (Channel: AV, beta: 17)")
    print("=================================================================")

    beta_str = "17"
    channel = "AV"
    binsize = 4

    ground_truth_path = Path(f"/Users/junxiongnie/code/ana/dat/func/pickdata/b4.{beta_str}/save_{channel}.csv")
    if not ground_truth_path.exists():
        print(f"❌ 未找到原工程基准数据: {ground_truth_path}")
        sys.exit(1)

    print(f"1. 加载基准数据: {ground_truth_path}")
    gt_data = np.loadtxt(ground_truth_path, delimiter=",")
    gt_means = gt_data[:, 0]
    gt_stds = gt_data[:, 1]
    print(f"   基准数据行数: {len(gt_means)}")

    output_dir = PROJECT_ROOT / "output" / f"b4.{beta_str}"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"save_{channel}.csv"

    print(f"2. 驱动 C++ 多线程并发数据分析流水线...")
    t0 = time.perf_counter()
    means, errors = process_channel(channel, beta_str, binsize=binsize)
    t1 = time.perf_counter()
    elapsed = t1 - t0
    print(f"   ✨ C++ 并发端到端分析总耗时: {elapsed:.3f} 秒")

    print(f"3. 存储运算结果至本项目: {output_file}")
    save_data = np.column_stack([means, errors])
    np.savetxt(output_file, save_data, delimiter=",", fmt="%.17e")

    print("4. 进行全量数据精度严格比对...")
    mean_diff = np.abs(means - gt_means)
    std_diff = np.abs(errors - gt_stds)

    max_mean_diff = np.max(mean_diff)
    max_std_diff = np.max(std_diff)

    print(f"   均值最大绝对偏差 (Max Mean Diff):   {max_mean_diff:.3e}")
    print(f"   误差最大绝对偏差 (Max Std Diff):    {max_std_diff:.3e}")

    # 机器精度比对 (容差在 1e-12 以内)
    is_means_match = np.allclose(means, gt_means, rtol=1e-11, atol=1e-13)
    is_stds_match = np.allclose(errors, gt_stds, rtol=1e-11, atol=1e-13)

    if is_means_match and is_stds_match:
        print("\n🎉 测试成功！新架构计算结果与原项目基准数据 100% 精确吻合！")
        print("   前 5 行比对样例 (新值 vs 基准值):")
        for i in range(5):
            print(f"   Row {i:02d} | Calc Mean: {means[i]:.16e} | GT Mean: {gt_means[i]:.16e}")
            print(f"          | Calc Std:  {errors[i]:.16e} | GT Std:  {gt_stds[i]:.16e}")
        return 0
    else:
        print("\n❌ 精度比对不匹配，请检查计算差异！")
        return 1

if __name__ == "__main__":
    sys.exit(main())

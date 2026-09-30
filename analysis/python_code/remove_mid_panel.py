#!/usr/bin/env python3
"""θ–E 图: 每张 3 面板 → 2 面板(删除中间"每bin平均"面板 ax2); 右面板降采样 60000→50000."""
import json
from pathlib import Path

P = Path('/home/liuwch/R10_v2_clean/analysis/gamma_beam_analysis_R10.ipynb')
with open(P, encoding='utf-8') as f:
    nb = json.load(f)

cell = None
for c in nb['cells']:
    if 'def plot_theta_E_series' in ''.join(c.get('source', [])):
        cell = c
        break
assert cell is not None
src = ''.join(cell['source'])

# 中间面板(ax2)代码块: E_dep 版
edep_ax2 = (
    "        H, _, _ = np.histogram2d(th_s, ed_s, bins=(bins_theta, bins_edep))\n"
    "        Hw, _, _ = np.histogram2d(th_s, ed_s, bins=(bins_theta, bins_edep), weights=ei_s)\n"
    "        with np.errstate(divide='ignore', invalid='ignore'):\n"
    "            mean_other = np.where(H > 0, Hw / np.maximum(H, 1), np.nan)\n"
    "        im = ax2.pcolormesh(bins_theta, bins_edep, mean_other.T, cmap='viridis',\n"
    "                            norm=LogNorm(vmin=1000, vmax=8000))\n"
    "        ax2.set_xlabel(r'$\\theta$ [deg]'); ax2.set_ylabel(r'$E_{dep}$ [keV]')\n"
    "        ax2.set_title(f'{src_label}: mean $E_{{init}}$ per bin ({sample})')\n"
    "        plt.colorbar(im, ax=ax2, label='mean $E_{init}$ [keV]')\n"
    "\n"
)
# 中间面板(ax2)代码块: E_init 版
einit_ax2 = (
    "        H, _, _ = np.histogram2d(th_s, ei_s, bins=(bins_theta, bins_einit))\n"
    "        Hw, _, _ = np.histogram2d(th_s, ei_s, bins=(bins_theta, bins_einit), weights=ed_s)\n"
    "        with np.errstate(divide='ignore', invalid='ignore'):\n"
    "            mean_other = np.where(H > 0, Hw / np.maximum(H, 1), np.nan)\n"
    "        im = ax2.pcolormesh(bins_theta, bins_einit, mean_other.T, cmap='plasma', norm=LogNorm(vmin=1))\n"
    "        ax2.set_xlabel(r'$\\theta$ [deg]'); ax2.set_ylabel(r'$E_{init}$ [keV]')\n"
    "        ax2.set_title(f'{src_label}: mean $E_{{dep}}$ per bin ({sample})')\n"
    "        plt.colorbar(im, ax=ax2, label='mean $E_{dep}$ [keV]')\n"
    "\n"
)

pairs = [
    # 1) 两个 section 的 subplots: 3 面板 → 2 面板
    ("fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(22, 6.5))",
     "fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(15, 6.5))"),
    # 2) 删除 E_dep 的 ax2 块
    (edep_ax2, ""),
    # 3) 删除 E_init 的 ax2 块
    (einit_ax2, ""),
    # 4) 右面板降采样 N: 60000 → 50000
    ("size=min(60000, len(th_s))", "size=min(50000, len(th_s))"),
    # 5) 理论曲线循环只画 ax1, ax3
    ("for ax in [ax1, ax2, ax3]:", "for ax in [ax1, ax3]:"),
    # 6) docstring / 头部注释
    ('"""生成某源的全部 6 张 θ–E 图 (3 样本 × 2 纵轴), 每张 3 面板."""',
     '"""生成某源的全部 6 张 θ–E 图 (3 样本 × 2 纵轴), 每张 2 面板(计数密度/逐记录着色)."""'),
    ("# 每张图 3 面板: 计数密度 / 每 bin 平均 / 逐记录着色",
     "# 每张图 2 面板: 计数密度 / 逐记录着色"),
]

for old, new in pairs:
    n = src.count(old)
    print(f"  {'OK' if n>=1 else '!!'} 出现 {n} 次: {old[:55]!r}")
    src = src.replace(old, new)

cell['source'] = src.splitlines(keepends=True)
with open(P, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
    f.write('\n')
print("已写回")

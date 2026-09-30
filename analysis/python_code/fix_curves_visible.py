#!/usr/bin/env python3
"""修复 _compton_curve_edep:
1) 端点公式错误: 原 2·E0²·m/(m+2E0) (多乘了 m_ec²=511, 使 linspace 步长~2000keV,
   曲线只剩 1 个有效点 → 完全看不见)。正确: E_dep_max = 2·E0²/(m+2E0) [keV]。
2) 曲线加粗+黑色描边, 保证在亮/暗背景都清晰可见; 标注加黑底。"""
import json
from pathlib import Path

P = Path('/home/liuwch/R10_v2_clean/analysis/gamma_beam_analysis_R10.ipynb')
with open(P, encoding='utf-8') as f:
    nb = json.load(f)

cell = None
for c in nb['cells']:
    if 'def _compton_curve_edep' in ''.join(c.get('source', [])):
        cell = c
        break
assert cell is not None
src = ''.join(cell['source'])

pairs = [
    # 1) 修正端点公式 (去掉多余的 ×M_E_C2)
    ("e_dep_max = e0_keV**2 * 2 * M_E_C2 / (M_E_C2 + 2 * e0_keV)  # 2E0²/(m+2E0)",
     "e_dep_max = 2 * e0_keV**2 / (M_E_C2 + 2 * e0_keV)  # E_dep_max=2E0²/(m+2E0) [keV] (2026-08-07 修正: 去掉多余 m 因子)"),
    # 2) 曲线: 黑色描边 + 白色主线, 加粗
    ("ax.plot(th, e_dep_th[valid], 'w--', lw=1, alpha=0.6)",
     "ax.plot(th, e_dep_th[valid], 'k--', lw=2.4, alpha=0.95)   # 深色描边,亮/暗背景均可见\n        ax.plot(th, e_dep_th[valid], 'w--', lw=1.4, alpha=1.0)  # 白色主线"),
    # 3) 标注: 加粗 + 黑底
    ("ax.annotate(f'{e0_keV/1000:.0f} MeV', (th[idx], e_dep_th[valid][idx]),\n                        color='white', fontsize=8, alpha=0.7)",
     "ax.annotate(f'{e0_keV/1000:.0f} MeV', (th[idx], e_dep_th[valid][idx]),\n                        color='white', fontsize=9, fontweight='bold',\n                        bbox=dict(boxstyle='round,pad=0.15', fc='black', alpha=0.55))"),
]

for old, new in pairs:
    if old in src:
        src = src.replace(old, new)
        print(f"  OK: {old[:50]!r}")
    else:
        print(f"  !! NOT FOUND: {old[:60]!r}")

cell['source'] = src.splitlines(keepends=True)
with open(P, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
    f.write('\n')
print("已写回")

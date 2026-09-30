#!/usr/bin/env python3
"""给 θ–E E_dep 纵轴图加回理论曲线 E0=1/3/8 MeV (当前为 1/8 MeV).
同时更新说明 markdown。用标准 json 读写。"""
import json
from pathlib import Path

P = Path('/home/liuwch/R10_v2_clean/analysis/gamma_beam_analysis_R10.ipynb')
with open(P, encoding='utf-8') as f:
    nb = json.load(f)


def repl_in_cell(cell, pairs):
    src = ''.join(cell.get('source', []))
    for old, new in pairs:
        if old in src:
            src = src.replace(old, new)
            print(f"  OK: {old[:44]!r}")
        else:
            print(f"  !! NOT FOUND: {old[:50]!r}")
    cell['source'] = src.splitlines(keepends=True)


# ---- 1) θ–E 系列 code cell ----
code = [c for c in nb['cells'] if 'def _compton_curve_edep' in ''.join(c.get('source', []))]
assert len(code) == 1, f"code cells: {len(code)}"
print("== θ–E code cell ==")
repl_in_cell(code[0], [
    ("只画下边界 E0=1 MeV 与上边界 E0=8 MeV",
     "E0=1/3/8 MeV (1=下界, 8=上界, 3=带内参考)"),
    ("def _compton_curve_edep(ax, e0_list=(1000, 8000)):",
     "def _compton_curve_edep(ax, e0_list=(1000, 3000, 8000)):"),
    ("只画谱的上下边界: E0=1 MeV (Compton 带下边界) 和 E0=8 MeV (上边界).",
     "画 E0=1/3/8 MeV 三条: 1 MeV 下界、8 MeV 上界、3 MeV 带内参考."),
    ("# 只叠加谱边界理论曲线: E0=1 MeV (下界) 与 E0=8 MeV (上界)",
     "# 叠加理论曲线: E0=1/3/8 MeV"),
    ("_compton_curve_edep(ax, e0_list=(1000, 8000))",
     "_compton_curve_edep(ax, e0_list=(1000, 3000, 8000))"),
])

# ---- 2) 说明 markdown cell ----
md = [c for c in nb['cells'] if c.get('cell_type') == 'markdown'
      and '理论曲线（2026-08-07 精简）' in ''.join(c.get('source', []))]
assert len(md) == 1, f"md cells: {len(md)}"
print("== 说明 markdown cell ==")
repl_in_cell(md[0], [
    ("**仅 E_dep 纵轴图**叠加白色虚线，且只画**谱的上下边界**：$E_0=1$ MeV（Compton 带下界）与 $E_0=8$ MeV（上界），中间 2/3/5 MeV 三条已删除。",
     "**仅 E_dep 纵轴图**叠加白色虚线：$E_0=1/3/8$ MeV 三条（1 MeV = Compton 带下界，8 MeV = 上界，3 MeV = 带内参考）。"),
    (r"把上面公式对 $E_0 \in \{1, 8\}$ MeV 逐点求值，得到 2 条白色虚线。",
     r"把上面公式对 $E_0 \in \{1, 3, 8\}$ MeV 逐点求值，得到 3 条白色虚线。"),
    ("中间 2/3/5 MeV 曲线不提供新信息（只是带内细分），删除后图面更干净。",
     "3 MeV 作为带内参考，可直观看到数据按入射能量的分带。"),
    ("数据集中在 1/8 MeV 理论曲线围成的\"Compton 带\"上",
     "数据集中在 1/3/8 MeV 理论曲线围成的\"Compton 带\"上"),
])

with open(P, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
    f.write('\n')
print("已写回")

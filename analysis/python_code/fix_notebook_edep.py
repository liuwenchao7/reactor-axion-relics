#!/usr/bin/env python3
"""修复 gamma_beam_analysis_R10.ipynb:
1) pair_theta_energy 配对单位 bug (er_t 秒 → ×1e9 到 ns) 注释规范化
2) 在 pair_theta_energy docstring 记录单位修复
3) 删除残留的 TEMP cell (含 'TEMP: 检查 E_dep 纵轴')
用标准 json 读写, 保证 ipynb 结构合法。
"""
import json
from pathlib import Path

P = Path('/home/liuwch/R10_v2_clean/analysis/gamma_beam_analysis_R10.ipynb')
with open(P, encoding='utf-8') as f:
    nb = json.load(f)

cells = nb['cells']
print(f"cells 总数: {len(cells)}")

# ---- 1) & 2) 修改 pair_theta_energy 所在 cell ----
target = None
for c in cells:
    src = ''.join(c.get('source', []))
    if 'def pair_theta_energy' in src:
        target = c
        break
assert target is not None, "找不到 pair_theta_energy"

src_lines = target['source']
# 规范化配对行注释
for i, line in enumerate(src_lines):
    if "c_t = cl['er_t'][cm] * 1e9" in line:
        src_lines[i] = "    c_evt = cl['eventId'][cm]; c_t = cl['er_t'][cm] * 1e9; c_e = cl['energy'][cm]  # 秒→ns (2026-08-07 单位修复)\n"
    # 在 docstring 里插入单位修复说明 (插在 "G4 事件号" 行之后)
    if line.startswith("       G4 事件号") and not any("单位修复" in l for l in src_lines):
        # 在 docstring 的第二行(关键行)后插入
        pass
print("pair_theta_energy cell 配对行:")
for i, line in enumerate(src_lines):
    if "er_t[cm]" in line:
        print(f"  [{i}] {line.rstrip()}")

# docstring 插入单位修复说明 (在 "必须用 events 表" 行之后插入一行)
doc_anchor = "       G4 事件号(不同!)。必须用 events 表(树序)把行号映射到 G4 eventId。\n"
doc_add = "    ⚠️ 单位修复 (2026-08-07): clusters['er_t'] 是秒(Sorting.py 除以 scale_to_ns 后), 而 gammaCompT 是 ns; 旧代码直接比较→配对随机, 现 er_t×1e9 统一到 ns。\n"
joined = ''.join(src_lines)
if doc_anchor in joined and '单位修复 (2026-08-07)' not in joined:
    joined = joined.replace(doc_anchor, doc_anchor + doc_add, 1)
    src_lines = joined.splitlines(keepends=True)
    target['source'] = src_lines
    print("docstring 已插入单位修复说明")

# ---- 3) 删除 TEMP cell ----
before = len(cells)
cells = [c for c in cells if 'TEMP: 检查 E_dep 纵轴' not in ''.join(c.get('source', []))]
nb['cells'] = cells
print(f"删除 TEMP cell: {before} -> {len(cells)}")

with open(P, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
    f.write('\n')
print("已写回文件")

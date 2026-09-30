#!/usr/bin/env python3
"""θ–E 说明 markdown: 中/右图 → 右图; 并加注 2 面板改动."""
import json
from pathlib import Path

P = Path('/home/liuwch/R10_v2_clean/analysis/gamma_beam_analysis_R10.ipynb')
with open(P, encoding='utf-8') as f:
    nb = json.load(f)

cell = None
for c in nb['cells']:
    if c.get('cell_type') == 'markdown' and '样本定义' in ''.join(c.get('source', [])):
        cell = c
        break
assert cell is not None
src = ''.join(cell['source'])

n1 = src.count('中/右图')
n2 = src.count('中/ 右图')
print(f"中/右图: {n1} 次, 中/ 右图: {n2} 次")
src = src.replace('中/ 右图', '右图').replace('中/右图', '右图')

anchor = "**物理意义**：\n"
if anchor in src and '2 面板' not in src:
    src = src.replace(anchor, anchor + "- 2026-08-09 每张图改为 **2 面板**（计数密度 / 逐记录着色），已删除原中间\"每 bin 平均\"面板。\n", 1)

cell['source'] = src.splitlines(keepends=True)
with open(P, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
    f.write('\n')
print("已写回")

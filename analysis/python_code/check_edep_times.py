#!/usr/bin/env python3
"""检查 gammaCompT 与 cluster er_t 的实际单位/数值."""
import numpy as np
import awkward as ak
import uproot
import h5py
from pathlib import Path

RESULT = Path('/home/liuwch/R10_v2_clean/result')
root_path = RESULT / 'gamma_source2_1e6_clean.root'
h5_path = RESULT / 'gamma_source2_1e6_clean.h5'

with h5py.File(h5_path, 'r') as f:
    events = f['events'][:]
    clusters = f['clusters'][:]
    primaries = f['primaries'][:]

t = uproot.open(str(root_path))['mcTree']
tm = t['gammaCompT'].array()          # ns? 
td = t['td'].array()                  # hit 时间
print("== ROOT 分支单位探针 ==")
print(f"gammaCompT[0] = {np.asarray(ak.flatten(tm))[:3]}")
print(f"td[0] (前几个 hit 时间) = {np.asarray(ak.flatten(td))[:3]}")
print(f"td[1] = {np.asarray(ak.flatten(td))[np.asarray(ak.flatten(td))>0][:3] if np.any(np.asarray(ak.flatten(td))>0) else 'all<=0'}")

print("\n== h5 clusters er_t 探针 ==")
et = clusters['er_t']
print(f"er_t 非零: {(et>0).sum()} / {len(et)}")
print(f"er_t 前 5 个非零值: {et[et>0][:5]}")
print(f"er_t 统计: min={et.min():.3e} max={et.max():.3e} mean={np.nanmean(et):.3e}")
print(f"clusters['t']? 字段: {list(clusters.dtype.names)}")

# gammaCompT 统计
tm_flat = np.asarray(ak.flatten(tm))
print(f"\ngammaCompT 统计: min={tm_flat.min():.3e} max={tm_flat.max():.3e} mean={tm_flat.mean():.3e} (前5: {tm_flat[:5]})")

# 对照: events 表有没有时间字段
print(f"events 字段: {list(events.dtype.names)}")

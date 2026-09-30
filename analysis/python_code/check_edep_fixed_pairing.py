#!/usr/bin/env python3
"""验证配对单位 bug: er_t(秒) vs gammaCompT(ns)。
对比 ①现状配对(单位不匹配) ②修正配对(统一到 ns) 得到的簇能量与 γ 能量损失的关系。"""
import numpy as np
import awkward as ak
import uproot
import h5py
from pathlib import Path
from collections import defaultdict

RESULT = Path('/home/liuwch/R10_v2_clean/result')
root_path = RESULT / 'gamma_source2_1e6_clean.root'
h5_path = RESULT / 'gamma_source2_1e6_clean.h5'

with h5py.File(h5_path, 'r') as f:
    events = f['events'][:]
    clusters = f['clusters'][:]
    primaries = f['primaries'][:]

t = uproot.open(str(root_path))['mcTree']
nGC = np.asarray(t['nGammaComp'].array())
vol = t['gammaCompVol'].array()
prex = t['gammaCompPrePx'].array(); prey = t['gammaCompPrePy'].array(); prez = t['gammaCompPrePz'].array()
px = t['gammaCompPx'].array(); py = t['gammaCompPy'].array(); pz = t['gammaCompPz'].array()
tm = t['gammaCompT'].array()          # ns

in_tpc = vol == 'XenonDetector'
cnt = np.asarray(ak.sum(in_tpc, axis=1))
np_ = np.sqrt(prex**2 + prey**2 + prez**2)
np2 = np.sqrt(px**2 + py**2 + pz**2)
dot = prex*px + prey*py + prez*pz
den = np_ * np2
th = np.degrees(np.arccos(np.clip(np.asarray(ak.flatten(np.where(den > 0, dot/den, np.nan))), -1, 1)))
cnts_all = np.asarray(nGC)
fevt = np.repeat(np.arange(len(cnts_all)), cnts_all)
fin = np.asarray(ak.flatten(in_tpc))
E_pre = np.asarray(ak.flatten(np_))
E_post = np.asarray(ak.flatten(np2))
tm_all = np.asarray(ak.flatten(tm))   # ns

ok = fin & np.isfinite(th)
E_loss = (E_pre - E_post)[ok]
ev_all = fevt[ok]
tmf_ns = tm_all[ok]

# 簇
row_to_evtid = events['eventId']
vol_bytes = b'XenonDetector'
cm = (clusters['energy'] > 0) & (clusters['er_volume'] == vol_bytes)
c_evt = clusters['eventId'][cm]
c_t_s = clusters['er_t'][cm]           # 秒
c_e = clusters['energy'][cm]
evt_cl = defaultdict(list)
for i in range(len(c_evt)):
    evt_cl[int(c_evt[i])].append(i)
evt_cl = {k: np.asarray(v) for k, v in evt_cl.items()}

def pair(c_t_units_ns):
    e_loc = np.full(len(ev_all), np.nan)
    for i in range(len(ev_all)):
        eid = int(row_to_evtid[ev_all[i]])
        idxs = evt_cl.get(eid)
        if idxs is None or len(idxs) == 0:
            continue
        # 统一到 ns 比较
        dt = np.abs(c_t_units_ns[idxs] - tmf_ns[i])
        e_loc[i] = c_e[idxs[np.argmin(dt)]]
    return e_loc

# ① 现状(notebook 的 pair_theta_energy): er_t 秒 直接与 gammaCompT ns 比 (错误)
e_loc_broken = pair(c_t_s)
# ② 修正: er_t 秒 → ns (×1e9) 再比
e_loc_fixed = pair(c_t_s * 1e9)

for name, e_loc in [('现状配对(秒vs ns, bug)', e_loc_broken), ('修正配对(统一ns)', e_loc_fixed)]:
    v = np.isfinite(e_loc)
    dl = E_loss[v]; de = e_loc[v]; d = de - dl
    print(f"\n=== {name} ===")
    print(f"  配对率={100*v.mean():.1f}%  N={v.sum():,}")
    print(f"  γ损失 median={np.median(dl):.1f} | 簇 median={np.median(de):.1f} | corr={np.corrcoef(de,dl)[0,1]:.3f}")
    print(f"  簇≈γ损失(<5keV): {100*(np.abs(d)<5).mean():.1f}% | 簇>γ损失+5: {100*(d>5).mean():.1f}% | 簇<γ损失-5: {100*(d<-5).mean():.1f}%")

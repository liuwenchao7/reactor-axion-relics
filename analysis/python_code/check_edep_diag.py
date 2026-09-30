#!/usr/bin/env python3
"""对齐诊断: 对 multifirst 记录, 打印前 10 条 (γ损失, 配对簇) 与相关系数, 与示例脚本对照."""
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
tm = t['gammaCompT'].array()

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
tm_all = np.asarray(ak.flatten(tm))

ok = fin & np.isfinite(th)
E_loss = (E_pre - E_post)[ok]
ev_all = fevt[ok]
tmf = tm_all[ok]

cnt_rec = cnt[ev_all]
multi = cnt_rec >= 2
order = np.lexsort((tmf, ev_all))
s_ev = ev_all[order]
fm = np.concatenate([[True], s_ev[1:] != s_ev[:-1]])
first_mask = np.zeros(len(ev_all), bool)
first_mask[order[fm]] = cnt_rec[order[fm]] >= 2

# ---- 配对簇 (与 pair_cluster_energy 相同) ----
row_to_evtid = events['eventId']
vol_bytes = b'XenonDetector'
cm = (clusters['energy'] > 0) & (clusters['er_volume'] == vol_bytes)
c_evt = clusters['eventId'][cm]; c_t = clusters['er_t'][cm]; c_e = clusters['energy'][cm]
evt_cl = defaultdict(list)
for i in range(len(c_evt)):
    evt_cl[int(c_evt[i])].append(i)
evt_cl = {k: np.asarray(v) for k, v in evt_cl.items()}

e_loc = np.full(len(ev_all), np.nan)
for i in range(len(ev_all)):
    eid = int(row_to_evtid[ev_all[i]])
    idxs = evt_cl.get(eid)
    if idxs is None or len(idxs) == 0:
        continue
    best = idxs[np.argmin(np.abs(c_t[idxs] - tmf[i]))]
    e_loc[i] = c_e[best]

fm_idx = np.where(first_mask)[0]
print(f"multifirst N={len(fm_idx):,}, 配对率={100*np.isfinite(e_loc[fm_idx]).mean():.1f}%")
dl = E_loss[fm_idx]; de = e_loc[fm_idx]
print(f"γ损失 median={np.median(dl):.1f}, cluster median={np.median(de):.1f}, "
      f"corr={np.corrcoef(de, dl)[0,1]:.3f}, 配对簇≈γ损失(<5keV)={100*(np.abs(de-dl)<5).mean():.1f}%")
print("\n前 10 条 multifirst (idx, G4event, γ损失, 簇, 差):")
for i in fm_idx[:10]:
    eid = int(row_to_evtid[ev_all[i]])
    print(f"  idx={i:6d} evt={eid:8d} | γloss={E_loss[i]:8.1f} cluster={e_loc[i]:8.1f} diff={e_loc[i]-E_loss[i]:+8.1f}")

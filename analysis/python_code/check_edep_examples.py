#!/usr/bin/env python3
"""查看若干 multifirst 记录所在事件的全部簇, 理解 簇能量 vs γ 能量损失 的差异来源."""
import numpy as np
import awkward as ak
import uproot
import h5py
from pathlib import Path
from collections import defaultdict

RESULT = Path('/home/liuwch/R10_v2_clean/result')

# ---- 复用 inTPC / pair 逻辑 (Source 2) ----
h5_path = RESULT / 'gamma_source2_1e6_clean.h5'
root_path = RESULT / 'gamma_source2_1e6_clean.root'

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
X = t['gammaCompX'].array(); Y = t['gammaCompY'].array(); Z = t['gammaCompZ'].array()

in_tpc = vol == 'XenonDetector'
cnt = np.asarray(ak.sum(in_tpc, axis=1))
np_ = np.sqrt(prex**2 + prey**2 + prez**2)
np2 = np.sqrt(px**2 + py**2 + pz**2)
dot = prex*px + prey*py + prez*pz
den = np_ * np2
th = np.degrees(np.arccos(np.clip(np.asarray(ak.flatten(np.where(den > 0, dot/den, np.nan))), -1, 1)))
fevt = np.repeat(np.arange(len(np.asarray(nGC))), np.asarray(nGC))
fin = np.asarray(ak.flatten(in_tpc))
E_loss_all = (np.asarray(ak.flatten(np_)) - np.asarray(ak.flatten(np2)))[fin]
th_in = th[fin]
ev_in = fevt[fin]
tm_in = np.asarray(ak.flatten(tm))[fin]
x_in = np.asarray(ak.flatten(X))[fin]; y_in = np.asarray(ak.flatten(Y))[fin]; z_in = np.asarray(ak.flatten(Z))[fin]

ok = np.isfinite(th_in)
E_loss = E_loss_all[ok]; th_f = th_in[ok]; ev = ev_in[ok]; tmf = tm_in[ok]
xf = x_in[ok]; yf = y_in[ok]; zf = z_in[ok]

# multifirst 记录
cnt_ev = cnt[ev]
multi = cnt_ev >= 2
order = np.lexsort((tmf, ev))
s_ev = ev[order]
fm = np.concatenate([[True], s_ev[1:] != s_ev[:-1]])
first_mask = np.zeros(len(ev), bool)
first_mask[order[fm]] = cnt_ev[order[fm]] >= 2

# 簇数据
row_to_evtid = events['eventId']
vol_bytes = b'XenonDetector'
cm = (clusters['energy'] > 0) & (clusters['er_volume'] == vol_bytes)
c_evt = clusters['eventId'][cm]; c_t = clusters['er_t'][cm]; c_e = clusters['energy'][cm]
c_x = clusters['er_x'][cm]; c_y = clusters['er_y'][cm]; c_z = clusters['er_z'][cm]
evt_cl = defaultdict(list)
for i in range(len(c_evt)):
    evt_cl[int(c_evt[i])].append(i)

# 找前 5 个 multifirst 记录
idx = np.where(first_mask)[0][:5]
for i in idx:
    eid = int(row_to_evtid[ev[i]])
    print(f"\n=== multifirst 记录 #idx={i}  G4event={eid}  θ={th_f[i]:.1f}°  "
          f"γ损失={E_loss[i]:.0f} keV  位置=({xf[i]:.0f},{yf[i]:.0f},{zf[i]:.0f})mm  t={tmf[i]:.2f}ns ===")
    print(f"  TPC 内 γ 记录 (该事件): ", end='')
    for j in np.where(ev == ev[i])[0]:
        print(f"[θ={th_f[j]:.0f}°, γloss={E_loss[j]:.0f}keV, t={tmf[j]:.2f}ns]", end=' ')
    print()
    cl_idx = evt_cl.get(eid, [])
    print(f"  该事件 XenonDetector 簇 ({len(cl_idx)}):")
    for j in cl_idx:
        match = f"  <-- 与第{i}条配对" if j == int(np.argmin(np.abs(c_t[cl_idx] - tmf[i]))) else ""
        print(f"    E={c_e[j]:7.1f} keV  t={c_t[j]:8.2f}ns  pos=({c_x[j]:.0f},{c_y[j]:.0f},{c_z[j]:.0f})mm{match}")
    # primary 能量
    print(f"  primaryEnergy: {primaries['primaryEnergy'][np.where(primaries['eventId']==eid)[0][0]]:.0f} keV")

#!/usr/bin/env python3
"""验证: 输入谱峰位能量(1.4/1.8/2.7/6.1 MeV)的 γ 在 TPC 内散射角 θ 是否铺满 0-180°."""
import numpy as np
import awkward as ak
import uproot
import h5py

t = uproot.open('/home/liuwch/R10_v2_clean/result/gamma_source2_1e6_clean.root')['mcTree']
nGC = np.asarray(t['nGammaComp'].array())
vol = t['gammaCompVol'].array()
prex = t['gammaCompPrePx'].array(); prey = t['gammaCompPrePy'].array(); prez = t['gammaCompPrePz'].array()
px = t['gammaCompPx'].array(); py = t['gammaCompPy'].array(); pz = t['gammaCompPz'].array()
np_ = np.sqrt(prex**2 + prey**2 + prez**2)
np2 = np.sqrt(px**2 + py**2 + pz**2)
dot = prex*px + prey*py + prez*pz
th = np.degrees(np.arccos(np.clip(
    np.asarray(ak.flatten(np.where(np_*np2 > 0, dot/(np_*np2), np.nan))), -1, 1)))
fevt = np.repeat(np.arange(len(nGC)), nGC)
fin = np.asarray(ak.flatten(vol == 'XenonDetector'))
ok = fin & np.isfinite(th)
th = th[ok]; ev = fevt[ok]

with h5py.File('/home/liuwch/R10_v2_clean/result/gamma_source2_1e6_clean.h5', 'r') as f:
    row_to_eid = f['events']['eventId'][:]
    pe = f['primaries']['primaryEnergy'][:]
    pe_evt = f['primaries']['eventId'][:]
evt_to_p = {int(e): i for i, e in enumerate(pe_evt)}
e_init = np.array([pe[evt_to_p[int(row_to_eid[ev[i]])]] for i in range(len(ev))])

print("按 E_init 切片统计 θ 覆盖 (source2):")
for lo, hi in [(1000, 1600), (1300, 1400), (1700, 1850), (2600, 2850), (6000, 6200)]:
    m = (e_init >= lo) & (e_init < hi)
    print(f"  E_init∈[{lo/1000:.1f},{hi/1000:.1f}] MeV: N={m.sum():6,}  "
          f"θ∈[{th[m].min():5.1f}°, {th[m].max():5.1f}°]  (跨度 {th[m].max()-th[m].min():5.1f}°)")

print("\n6.1 MeV 切片的 θ 直方图 (0-180°, 30°/bin):",
      np.histogram(th[(e_init >= 6000) & (e_init < 6200)], bins=np.linspace(0, 180, 7))[0])
print("1.4 MeV 切片的 θ 直方图 (0-180°, 30°/bin):",
      np.histogram(th[(e_init >= 1300) & (e_init < 1400)], bins=np.linspace(0, 180, 7))[0])

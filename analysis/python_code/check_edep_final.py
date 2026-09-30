#!/usr/bin/env python3
"""
结论脚本: E_dep 轴 = 局部能量沉积(cluster energy), 不是 γ 能量损失。
且当前 pair_theta_energy 的"最近时间"配对有单位 bug (er_t=秒 vs gammaCompT=ns)。
对比 现状配对 vs 修正配对(统一 ns), 对 all/multi/multifirst 全部样本给出。
"""
import numpy as np
import awkward as ak
import uproot
import h5py
from pathlib import Path
from collections import defaultdict

RESULT = Path('/home/liuwch/R10_v2_clean/result')


def analyze(src_label, root_path, h5_path):
    with h5py.File(h5_path, 'r') as f:
        events = f['events'][:]
        clusters = f['clusters'][:]

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
    fevt = np.repeat(np.arange(len(np.asarray(nGC))), np.asarray(nGC))
    fin = np.asarray(ak.flatten(in_tpc))
    E_pre = np.asarray(ak.flatten(np_))
    E_post = np.asarray(ak.flatten(np2))
    tm_all = np.asarray(ak.flatten(tm))

    ok = fin & np.isfinite(th)
    E_loss = (E_pre - E_post)[ok]
    ev_all = fevt[ok]
    tmf_ns = tm_all[ok]

    cnt_rec = cnt[ev_all]
    multi = cnt_rec >= 2
    order = np.lexsort((tmf_ns, ev_all))
    s_ev = ev_all[order]
    fm = np.concatenate([[True], s_ev[1:] != s_ev[:-1]])
    first_mask = np.zeros(len(ev_all), bool)
    first_mask[order[fm]] = cnt_rec[order[fm]] >= 2

    # 簇
    row_to_evtid = events['eventId']
    vol_bytes = b'XenonDetector'
    cm = (clusters['energy'] > 0) & (clusters['er_volume'] == vol_bytes)
    c_evt = clusters['eventId'][cm]
    c_t_s = clusters['er_t'][cm]          # 秒
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
            e_loc[i] = c_e[idxs[np.argmin(np.abs(c_t_units_ns[idxs] - tmf_ns[i]))]]
        return e_loc

    e_broken = pair(c_t_s)          # 现状: 秒 vs ns
    e_fixed = pair(c_t_s * 1e9)     # 修正: 统一 ns

    print(f"\n################ {src_label} ################")
    for name, m in [('all', np.ones(len(ev_all), bool)),
                    ('multi', multi),
                    ('multifirst', first_mask)]:
        v = m & np.isfinite(e_broken)
        print(f"\n--- {src_label} [{name}]  N={v.sum():,} ---")
        print(f"  {'':10s} {'γ损失med':>9s} {'簇med':>9s} {'corr':>6s} {'≈(±5keV)':>8s} {'>γ+5':>6s} {'<γ-5':>6s}")
        for tag, e_loc in [('现状配对(单位bug)', e_broken), ('修正配对(统一ns)', e_fixed)]:
            dl = E_loss[v]; de = e_loc[v]; d = de - dl
            print(f"  {tag:20s} {np.median(dl):9.1f} {np.median(de):9.1f} {np.corrcoef(de,dl)[0,1]:6.3f} "
                  f"{100*(np.abs(d)<5).mean():8.1f} {100*(d>5).mean():6.1f} {100*(d<-5).mean():6.1f}")


analyze('Source 1', RESULT/'gamma_source1_1e6_clean.root', RESULT/'gamma_source1_1e6_clean.h5')
analyze('Source 2', RESULT/'gamma_source2_1e6_clean.root', RESULT/'gamma_source2_1e6_clean.h5')

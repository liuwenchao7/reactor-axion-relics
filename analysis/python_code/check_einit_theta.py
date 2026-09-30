#!/usr/bin/env python3
"""
检查 θ–E 图中 E_init 纵轴与 θ 横轴数据是否有 bug。
θ    = ∠(p_in, p_out)     (compute_inTPC_theta, 单位 deg)
E_init = primaryEnergy     (pair_theta_energy 按 eventId 映射, 单位 keV)
检查点:
 1) θ: NaN/越界/分布合理性
 2) E_init: 未配对率、离群点(120GeV)、范围
 3) 能量守恒: |p_in| <= E_init (E_init 是初级能量, 任何散射点的 γ 能量都不能超过它)
    —— 若出现 |p_in| > E_init, 说明 E_init 配错了事件!
 4) Compton 运动学: E_dep <= E_dep_max(E_init) = 2E_init²/(m+2E_init)
 5) 对齐: len(theta)==len(e_init), 抽查若干记录 (θ, E_init, |p_in|, |p_out|, E_dep)
"""
import numpy as np
import awkward as ak
import uproot
import h5py
from pathlib import Path
from collections import defaultdict

RESULT = Path('/home/liuwch/R10_v2_clean/result')
M_E_C2 = 511.0


def load_h5(h5_path):
    with h5py.File(h5_path, 'r') as f:
        return {k: f[k][:] for k in ['events', 'primaries', 'clusters']}


def check(src_label, root_path, h5_path):
    data = load_h5(h5_path)
    t = uproot.open(str(root_path))['mcTree']
    nGC = np.asarray(t['nGammaComp'].array())
    vol = t['gammaCompVol'].array()
    prex = t['gammaCompPrePx'].array(); prey = t['gammaCompPrePy'].array(); prez = t['gammaCompPrePz'].array()
    px = t['gammaCompPx'].array(); py = t['gammaCompPy'].array(); pz = t['gammaCompPz'].array()
    tm = t['gammaCompT'].array()

    in_tpc = vol == 'XenonDetector'
    np_ = np.sqrt(prex**2 + prey**2 + prez**2)   # |p_in| keV/c
    np2 = np.sqrt(px**2 + py**2 + pz**2)         # |p_out| keV/c
    dot = prex*px + prey*py + prez*pz
    den = np_ * np2
    th = np.degrees(np.arccos(np.clip(np.asarray(ak.flatten(np.where(den > 0, dot/den, np.nan))), -1, 1)))
    cnts_all = np.asarray(nGC)
    fevt = np.repeat(np.arange(len(cnts_all)), cnts_all)
    fin = np.asarray(ak.flatten(in_tpc))
    E_pre = np.asarray(ak.flatten(np_))
    E_post = np.asarray(ak.flatten(np2))
    ok = fin & np.isfinite(th)

    th_f = th[ok]; ev_f = fevt[ok]; E_pre_f = E_pre[ok]; E_post_f = E_post[ok]
    tm_f = np.asarray(ak.flatten(tm))[ok]   # 展平后的 gammaCompT (ns)
    E_loss = E_pre_f - E_post_f

    # ---- 复算 pair_theta_energy (与 notebook 一致, 含 er_t×1e9 修复) ----
    cl = data['clusters']; p = data['primaries']
    row_to_evtid = data['events']['eventId']
    vol_bytes = b'XenonDetector'
    cm = (cl['energy'] > 0) & (cl['er_volume'] == vol_bytes)
    c_evt = cl['eventId'][cm]; c_t = cl['er_t'][cm] * 1e9; c_e = cl['energy'][cm]
    evt_cl = defaultdict(list)
    for i in range(len(c_evt)):
        evt_cl[int(c_evt[i])].append(i)
    evt_cl = {k: np.asarray(v) for k, v in evt_cl.items()}
    evt_to_p = {int(p['eventId'][i]): i for i in range(len(p))}

    e_dep = np.full(len(th_f), np.nan)
    e_init = np.full(len(th_f), np.nan)
    for i in range(len(th_f)):
        eid = int(row_to_evtid[ev_f[i]])
        idxs = evt_cl.get(eid)
        if idxs is not None and len(idxs) > 0:
            e_dep[i] = c_e[idxs[np.argmin(np.abs(c_t[idxs] - tm_f[i]))]]
        pi = evt_to_p.get(eid)
        if pi is not None:
            e_init[i] = p['primaryEnergy'][pi]

    print(f"\n################ {src_label} ################")
    # 1) θ
    print(f"[θ] N={len(th_f):,}  NaN={np.isnan(th_f).sum()}  "
          f"min={th_f.min():.2f}  max={th_f.max():.2f}  越界(>180或<0)={((th_f<0)|(th_f>180)).sum()}"
          f"  median={np.median(th_f):.1f}°")

    # 2) E_init
    n_nan = np.isnan(e_init).sum()
    ei = e_init[np.isfinite(e_init)]
    n_hi8 = (ei > 8000).sum()      # 超出图纵轴上限 8000 keV
    n_hiG = (ei > 1e5).sum()       # 离谱大值(>100 MeV)
    print(f"[E_init] 未配对NaN={n_nan:,} ({100*n_nan/len(th_f):.2f}%)  范围=[{ei.min():.0f},{ei.max():.0f}] keV"
          f"  median={np.median(ei):.1f}  超出8000keV={n_hi8:,}   >100MeV={n_hiG:,}")

    # 3) 能量守恒 |p_in| <= E_init
    m3 = np.isfinite(e_init)
    viol = (E_pre_f[m3] > e_init[m3] + 1.0).sum()   # 容差 1 keV
    ratio = E_pre_f[m3] / e_init[m3]
    print(f"[守恒] |p_in| <= E_init: 违反 {viol:,} / {m3.sum():,} ({100*viol/m3.sum():.4f}%)")
    print(f"        |p_in|/E_init: p5={np.percentile(ratio,5):.3f} p50={np.median(ratio):.3f} p95={np.percentile(ratio,95):.3f}")

    # 4) Compton 运动学: E_dep <= 2E_init²/(m+2E_init)
    m4 = np.isfinite(e_init) & np.isfinite(e_dep)
    edep_max = 2 * e_init[m4]**2 / (M_E_C2 + 2 * e_init[m4])
    bad = (e_dep[m4] > edep_max + 5).sum()
    print(f"[运动学] E_dep <= 2E_init²/(m+2E_init): 违反 {bad:,} / {m4.sum():,} ({100*bad/m4.sum():.4f}%)")

    # 5) E_init 谱形 vs 输入谱
    hist, _ = np.histogram(ei, bins=np.linspace(1000, 8000, 71))
    peak_lo = hist[:6].sum()/hist.sum()*100   # 1-1.6 MeV 占比
    print(f"[谱形] E_init 在 1-1.6 MeV 占比={peak_lo:.1f}%  (输入谱该段约86.8%)")

    # 6) 抽查前 5 条记录
    print("[抽查] idx  θ(°)  E_init(keV) |p_in|(keV) |p_out|(keV) E_loss(keV) E_dep(keV)")
    for i in range(5):
        print(f"  {i:5d} {th_f[i]:6.1f} {e_init[i]:10.1f} {E_pre_f[i]:10.1f} {E_post_f[i]:10.1f} "
              f"{E_loss[i]:10.1f} {e_dep[i]:10.1f}")


check('Source 1', RESULT/'gamma_source1_1e6_clean.root', RESULT/'gamma_source1_1e6_clean.h5')
check('Source 2', RESULT/'gamma_source2_1e6_clean.root', RESULT/'gamma_source2_1e6_clean.h5')

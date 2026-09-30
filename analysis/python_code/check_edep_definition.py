#!/usr/bin/env python3
"""
检查 θ–E 图 E_dep 纵轴的数据定义:
  - γ 能量损失 (kinematic): |p_in| - |p_out|, 由 gammaCompPrePx/Py/Pz 与 gammaCompPx/Py/Pz 计算
    (光子 E=pc, 动量单位 keV/c → 差值即 keV)
  - 局部能量沉积 (local deposit): pair_theta_energy 配对的 clusters['energy']
    (Sorting.py DBSCAN 簇的总沉积能量, eps=3mm/~10ns)
对 all / multi / multifirst 三个样本逐一比较。
"""
import numpy as np
import awkward as ak
import uproot
import h5py
from pathlib import Path
from collections import defaultdict

RESULT = Path('/home/liuwch/R10_v2_clean/result')


def load_h5(h5_path):
    with h5py.File(h5_path, 'r') as f:
        return {k: f[k][:] for k in ['events', 'primaries', 'clusters']}


def inTPC(root_path, vol_filter='XenonDetector'):
    """与 notebook compute_inTPC_theta 完全一致的过滤/排序; 额外给出 γ 能量损失."""
    t = uproot.open(str(root_path))['mcTree']
    nGC = np.asarray(t['nGammaComp'].array())
    vol = t['gammaCompVol'].array()
    prex = t['gammaCompPrePx'].array(); prey = t['gammaCompPrePy'].array(); prez = t['gammaCompPrePz'].array()
    px = t['gammaCompPx'].array(); py = t['gammaCompPy'].array(); pz = t['gammaCompPz'].array()
    tm = t['gammaCompT'].array()

    in_tpc = vol == vol_filter
    cnt = np.asarray(ak.sum(in_tpc, axis=1))

    np_ = np.sqrt(prex**2 + prey**2 + prez**2)   # |p_in| keV/c ≈ 入射 γ 能量
    np2 = np.sqrt(px**2 + py**2 + pz**2)         # |p_out| keV/c ≈ 出射 γ 能量
    dot = prex * px + prey * py + prez * pz
    den = np_ * np2
    th = np.degrees(np.arccos(np.clip(
        np.asarray(ak.flatten(np.where(den > 0, dot / den, np.nan))), -1, 1)))

    cnts_all = np.asarray(nGC)
    fevt = np.repeat(np.arange(len(cnts_all)), cnts_all)
    fin = np.asarray(ak.flatten(in_tpc))
    E_pre = np.asarray(ak.flatten(np_))
    E_post = np.asarray(ak.flatten(np2))
    tm_all = np.asarray(ak.flatten(tm))

    ok = fin & np.isfinite(th)
    return dict(th_all=th[ok], ev_all=fevt[ok], tm_all=tm_all[ok],
                cnt=cnt, E_loss=(E_pre - E_post)[ok])


def pair_cluster_energy(res, data, vol_filter='XenonDetector'):
    """与 notebook pair_theta_energy 相同: 每条记录配对同事件时间最近 TPC 簇能量."""
    cl = data['clusters']
    row_to_evtid = data['events']['eventId']
    vol_bytes = vol_filter.encode()
    cm = (cl['energy'] > 0) & (cl['er_volume'] == vol_bytes)
    c_evt = cl['eventId'][cm]; c_t = cl['er_t'][cm]; c_e = cl['energy'][cm]
    evt_cl = defaultdict(list)
    for i in range(len(c_evt)):
        evt_cl[int(c_evt[i])].append(i)
    evt_cl = {k: np.asarray(v) for k, v in evt_cl.items()}

    th = res['th_all']; row = res['ev_all']; tt = res['tm_all']
    e_dep = np.full(len(th), np.nan)
    for i in range(len(th)):
        eid = int(row_to_evtid[row[i]])
        idxs = evt_cl.get(eid)
        if idxs is None or len(idxs) == 0:
            continue
        best = idxs[np.argmin(np.abs(c_t[idxs] - tt[i]))]
        e_dep[i] = c_e[best]
    return e_dep


def check(src_label, h5, root):
    data = load_h5(h5)
    res = inTPC(root)
    e_loc = pair_cluster_energy(res, data)
    E_loss = res['E_loss']
    ev_all = res['ev_all']
    cnt_rec = res['cnt'][ev_all]
    multi = cnt_rec >= 2

    # multifirst: 多重事件每事件时间最早一条
    order = np.lexsort((res['tm_all'], ev_all))
    s_ev = ev_all[order]
    fm_sorted = np.concatenate([[True], s_ev[1:] != s_ev[:-1]])
    first_mask = np.zeros(len(ev_all), bool)
    first_mask[order[fm_sorted]] = cnt_rec[order[fm_sorted]] >= 2

    valid = np.isfinite(e_loc)
    print(f"\n===== {src_label} =====")
    for name, m in [('all', np.ones(len(E_loss), bool)),
                    ('multi', multi),
                    ('multifirst', first_mask)]:
        vm = m & valid
        dl = E_loss[vm]; de = e_loc[vm]; d = de - dl
        print(f"--- {src_label} [{name}]  N={vm.sum():,} ---")
        print(f"  γ 能量损失 |p_in|-|p_out|:  median={np.median(dl):8.1f} keV  mean={dl.mean():8.1f}")
        print(f"  簇局部沉积 cluster energy:  median={np.median(de):8.1f} keV  mean={de.mean():8.1f}")
        print(f"  差 (簇 - γ损失):            median={np.median(d):+7.1f}  mean={d.mean():+7.1f}")
        print(f"  簇≈γ损失 (|差|<5 keV):      {100*(np.abs(d)<5).mean():5.1f}%")
        print(f"  簇>γ损失+5 keV (多余沉积):  {100*(d>5).mean():5.1f}%")
        print(f"  簇<γ损失-5 keV:             {100*(d<-5).mean():5.1f}%")
        print(f"  corr(簇, γ损失) = {np.corrcoef(de, dl)[0,1]:.3f}")
        extra = d[d > 5]
        if len(extra):
            print(f"  多余沉积量: median={np.median(extra):.1f} keV, "
                  f"占簇能量中位比={100*np.median(extra/de[d>5]):.0f}%")


if __name__ == '__main__':
    check('Source 1', RESULT / 'gamma_source1_1e6_clean.h5', RESULT / 'gamma_source1_1e6_clean.root')
    check('Source 2', RESULT / 'gamma_source2_1e6_clean.h5', RESULT / 'gamma_source2_1e6_clean.root')

import sys
import os
import json
import h5py
import numpy as np

# ==================== 全局参数（所有源共用） ====================
ER_MAX_0 = 200
ER_MAX = int(ER_MAX_0)
ER_MIN_0 = 0
ER_MIN = int(ER_MIN_0)
bins_ER_0 = 1000
bins_ER = int(bins_ER_0)

position_lxe = -37
water_shift_z = 0
solid_shift_x = -150
dis_gate = 0.25
dis_cathode = 15.3575
dis_botscreen = 0.2625
dis_fv_top_to_lxe = 0.35
fv_r = 80
fv_z1 = (position_lxe + water_shift_z) * 10
fv_z2 = (position_lxe + water_shift_z - dis_gate - dis_cathode) * 10

veto1 = 200
veto2 = 1000
veto_z1 = (position_lxe + water_shift_z) * 10
veto_z2 = (position_lxe + water_shift_z - dis_gate - dis_cathode - dis_botscreen) * 10

cut_nr2 = 0.1
cut_er1 = 0.05
cut_nr_er1 = 0.05

# 辅助函数（与原始代码一致）
def er_energy(events):
    return (events['max_er'] < ER_MAX) & (events['max_er'] >= ER_MIN)

def cut_fv(events, s):
    mask = np.sqrt((events[f'max_{s}_x'] - solid_shift_x)**2 + events[f'max_{s}_y']**2) < fv_r
    mask &= events[f'max_{s}_z'] < fv_z1
    mask &= events[f'max_{s}_z'] > fv_z2
    return mask

def muon_veto_er(events):
    return events['max_er_t'] < 350e-6

def rate_error_er(events, er_max):
    mask = (events['max_er'] < er_max) & (events['max_er'] >= ER_MIN)
    return len(events['max_er'][mask])

# # 归一化常数
# days = 24 * 3600
# kg = (fv_z1 - fv_z2)/10 * (fv_r/10) * (fv_r/10) * np.pi * 2.862 / 1000
# print(f'Common FV mass: {kg:.3f} kg')

# 读取归一化文件
with open('../../normalization.json', 'r') as f:
    genNormalize = json.loads(f.read())

# 归一化常数
days = 24 * 3600
kg = genNormalize['XenonDetector']['mass']
print(f'LXe mass: {kg} kg')
print('-----------------------')

# ------------------ 定义各源处理函数（返回 counts 和 error 数组） ------------------
def process_CRN(folder):
    bins_er = np.linspace(ER_MIN, ER_MAX, bins_ER+1)
    files = [os.path.join(folder, 'events_CRN.h5')]
    with h5py.File(files[0], 'r') as ipt:
        e = ipt['events'][:]
    weights = np.full(len(e), genNormalize['CRN']['activity']*days / genNormalize['CRN']['factor']) / np.diff(bins_er)[np.clip(np.digitize(e['max_er'], bins_er)-1, 0, bins_ER-1)]/kg
    counts = np.histogram(e['max_er'], weights=weights, bins=bins_er)[0]
    bins_diff = np.diff(bins_er)
    act = genNormalize['CRN']['factor'] / (genNormalize['CRN']['activity'] * days)
    c = 1 / (act * kg)
    error = np.array([(a * c / b)**0.5 for a, b in zip(counts, bins_diff)])
    return counts, error

def process_gamma_bottom(folder):
    bins_er = np.linspace(ER_MIN, ER_MAX, bins_ER+1)
    files = [os.path.join(folder, 'flux_gamma_OFF_BOTTOM_Z_OFF.h5')]
    with h5py.File(files[0], 'r') as ipt:
        e = ipt['events'][:]
    norm = genNormalize['gamma_E_BOTTOM_Z']['activity']*days / (genNormalize['gamma_E_BOTTOM_Z']['direct_factor'])
    weights = np.full(len(e), norm) / np.diff(bins_er)[np.clip(np.digitize(e['max_er'], bins_er)-1, 0, bins_ER-1)]/kg
    counts = np.histogram(e['max_er'], weights=weights, bins=bins_er)[0]
    bins_diff = np.diff(bins_er)
    act = (genNormalize['gamma_E_BOTTOM_Z']['direct_factor']) / (genNormalize['gamma_E_BOTTOM_Z']['activity']*days)
    c = 1 / (act * kg)
    error = np.array([(a * c / b)**0.5 for a, b in zip(counts, bins_diff)])
    return counts, error

def process_gamma_front(folder):
    bins_er = np.linspace(ER_MIN, ER_MAX, bins_ER+1)
    files = [os.path.join(folder, 'flux_gamma_OFF_FRONT_X_OFF.h5')]
    with h5py.File(files[0], 'r') as ipt:
        e = ipt['events'][:]
    norm = genNormalize['gamma_E_FRONT_X']['activity']*days / (genNormalize['gamma_E_FRONT_X']['direct_factor'])
    weights = np.full(len(e), norm) / np.diff(bins_er)[np.clip(np.digitize(e['max_er'], bins_er)-1, 0, bins_ER-1)]/kg
    counts = np.histogram(e['max_er'], weights=weights, bins=bins_er)[0]
    bins_diff = np.diff(bins_er)
    act = (genNormalize['gamma_E_FRONT_X']['direct_factor']) / (genNormalize['gamma_E_FRONT_X']['activity']*days)
    c = 1 / (act * kg)
    error = np.array([(a * c / b)**0.5 for a, b in zip(counts, bins_diff)])
    return counts, error

def process_gamma_side(folder):
    bins_er = np.linspace(ER_MIN, ER_MAX, bins_ER+1)
    files = [os.path.join(folder, 'flux_gamma_OFF_SIDE_Y_OFF.h5')]
    with h5py.File(files[0], 'r') as ipt:
        e = ipt['events'][:]
    norm = genNormalize['gamma_E_SIDE_Y']['activity']*days / (genNormalize['gamma_E_SIDE_Y']['direct_factor'])
    weights = np.full(len(e), norm) / np.diff(bins_er)[np.clip(np.digitize(e['max_er'], bins_er)-1, 0, bins_ER-1)]/kg
    counts = np.histogram(e['max_er'], weights=weights, bins=bins_er)[0]
    bins_diff = np.diff(bins_er)
    act = (genNormalize['gamma_E_SIDE_Y']['direct_factor']) / (genNormalize['gamma_E_SIDE_Y']['activity']*days)
    c = 1 / (act * kg)
    error = np.array([(a * c / b)**0.5 for a, b in zip(counts, bins_diff)])
    return counts, error

def process_material_main(folder):
    # 包含 steel, Teflon, PMTbase, PMTwindow, Flange
    bins_er = np.linspace(ER_MIN, ER_MAX, bins_ER+1)
    total_counts = np.zeros(bins_ER)
    total_error_sq = np.zeros(bins_ER)
    
    # 定义材料组
    material_groups = [
        (['U238','Th232','Co60','K40','Cs137'], ['steel']),
        (['U238','Th232','Co60','K40','Cs137'], ['Teflon']),
        (['U238','Th232','Co60','K40','Cs137'], ['copper']),
        (['U238','U235','Th232','Cs137','K40','Co60'], ['PMTbase']),
        (['U238','U235','Th232','Cs137','K40','Co60'], ['PMTwindow']),
        (['U238','Th232','Co60','K40','Cs137','U235'], ['Flange'])
    ]
    for isotopes, components in material_groups:
        for iso in isotopes:
            for comp in components:
                fpath = os.path.join(folder, f'events_{iso}_{comp}.h5')
                if not os.path.exists(fpath):
                    print(f"Warning: {fpath} not found, skip.")
                    continue
                with h5py.File(fpath, 'r') as ipt:
                    e = ipt['events'][:]
                key = iso+'_'+comp
                weights = np.full(len(e), genNormalize[key]['activity']*days / genNormalize[key]['factor']) / np.diff(bins_er)[np.clip(np.digitize(e['max_er'], bins_er)-1, 0, bins_ER-1)]/kg
                counts = np.histogram(e['max_er'], weights=weights, bins=bins_er)[0]
                total_counts += counts
                bins_diff = np.diff(bins_er)
                act = genNormalize[key]['factor'] / (genNormalize[key]['activity'] * days)
                c = 1 / (act * kg)
                err = np.array([(a * c / b)**0.5 for a, b in zip(counts, bins_diff)])
                total_error_sq += err**2
    total_error = np.sqrt(total_error_sq)
    return total_counts, total_error

def process_muon(folder):
    bins_er = np.linspace(ER_MIN, ER_MAX, bins_ER+1)
    files = [os.path.join(folder, 'events_muon.h5')]
    with h5py.File(files[0], 'r') as ipt:
        e = ipt['events'][:]
    # 注意：muon 代码中使用了 er_energy 和 muon_veto_er，但最终去掉了 veto，只保留原始全谱
    # mask = er_energy(e)   # 仅能量范围 cut
    weights = np.full(len(e), genNormalize['muon']['activity']*days / genNormalize['muon']['factor']) / np.diff(bins_er)[np.clip(np.digitize(e['max_er'], bins_er)-1, 0, bins_ER-1)]/kg
    counts = np.histogram(e['max_er'], weights=weights, bins=bins_er)[0]
    bins_diff = np.diff(bins_er)
    act = genNormalize['muon']['factor'] / (genNormalize['muon']['activity'] * days)
    c = 1 / (act * kg)
    error = np.array([(a * c / b)**0.5 for a, b in zip(counts, bins_diff)])
    return counts, error

def process_material_PE(folder):
    bins_er = np.linspace(ER_MIN, ER_MAX, bins_ER+1)
    total_counts = np.zeros(bins_ER)
    total_error_sq = np.zeros(bins_ER)
    isotopes = ['U238','Th232','K40','Cs137','Co60']
    for iso in isotopes:
        fpath = os.path.join(folder, f'events_{iso}_PE.h5')
        if not os.path.exists(fpath):
            print(f"Warning: {fpath} not found, skip.")
            continue
        with h5py.File(fpath, 'r') as ipt:
            e = ipt['events'][:]
        key = iso+'_PE'
        weights = np.full(len(e), genNormalize[key]['activity']*days / genNormalize[key]['factor']) / np.diff(bins_er)[np.clip(np.digitize(e['max_er'], bins_er)-1, 0, bins_ER-1)]/kg
        counts = np.histogram(e['max_er'], weights=weights, bins=bins_er)[0]
        total_counts += counts
        bins_diff = np.diff(bins_er)
        act = genNormalize[key]['factor'] / (genNormalize[key]['activity'] * days)
        c = 1 / (act * kg)
        err = np.array([(a * c / b)**0.5 for a, b in zip(counts, bins_diff)])
        total_error_sq += err**2
    total_error = np.sqrt(total_error_sq)
    return total_counts, total_error

def process_material_PS(folder):
    bins_er = np.linspace(ER_MIN, ER_MAX, bins_ER+1)
    total_counts = np.zeros(bins_ER)
    total_error_sq = np.zeros(bins_ER)
    isotopes = ['U238','Th232','K40']
    for iso in isotopes:
        fpath = os.path.join(folder, f'events_{iso}_PS.h5')
        if not os.path.exists(fpath):
            print(f"Warning: {fpath} not found, skip.")
            continue
        with h5py.File(fpath, 'r') as ipt:
            e = ipt['events'][:]
        key = iso+'_PS'
        weights = np.full(len(e), genNormalize[key]['activity']*days / genNormalize[key]['factor']) / np.diff(bins_er)[np.clip(np.digitize(e['max_er'], bins_er)-1, 0, bins_ER-1)]/kg
        counts = np.histogram(e['max_er'], weights=weights, bins=bins_er)[0]
        total_counts += counts
        bins_diff = np.diff(bins_er)
        act = genNormalize[key]['factor'] / (genNormalize[key]['activity'] * days)
        c = 1 / (act * kg)
        err = np.array([(a * c / b)**0.5 for a, b in zip(counts, bins_diff)])
        total_error_sq += err**2
    total_error = np.sqrt(total_error_sq)
    return total_counts, total_error

# ------------------ 设置所有源的数据路径（请根据实际情况修改） ------------------
paths = {
    # 'CRN': '../../../../result/CRN_500M_6m',
    'gamma_bottom': '../../../../result/gamma_OFF_direct/gammaOFF_Z_direct_30B',
    'gamma_front': '../../../../result/gamma_OFF_direct/gammaOFF_X_direct_30B',
    'gamma_side': '../../../../result/gamma_OFF_direct/gammaOFF_Y_direct_30B',
    # 'material_main': '/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/RelicsSim_v2_solid_origin_material_R_E/result/Material_10M_steel_window_1M_others',
    # 'muon': '../../../../result/muon_900M_6m',
    # 'PE': '../../../../result/Material/PE_150M',
    # 'PS': '../../../../result/Material/PS_50M_K40_500M'
}

# ------------------ 依次处理并累加 ------------------
total_counts = np.zeros(bins_ER)
total_error_sq = np.zeros(bins_ER)

# print("Processing CRN...")
# c, e = process_CRN(paths['CRN'])
# total_counts += c
# total_error_sq += e**2

print("Processing Gamma Bottom...")
c, e = process_gamma_bottom(paths['gamma_bottom'])
total_counts += c
total_error_sq += e**2

print("Processing Gamma Front...")
c, e = process_gamma_front(paths['gamma_front'])
total_counts += c
total_error_sq += e**2

print("Processing Gamma Side...")
c, e = process_gamma_side(paths['gamma_side'])
total_counts += c
total_error_sq += e**2

# print("Processing Main Materials (Steel/Teflon/PMT/Flange)...")
# c, e = process_material_main(paths['material_main'])
# total_counts += c
# total_error_sq += e**2

# print("Processing Muon...")
# c, e = process_muon(paths['muon'])
# total_counts += c
# total_error_sq += e**2

# print("Processing PE Material...")
# c, e = process_material_PE(paths['PE'])
# total_counts += c
# total_error_sq += e**2

# print("Processing PS Material...")
# c, e = process_material_PS(paths['PS'])
# total_counts += c
# total_error_sq += e**2

total_error = np.sqrt(total_error_sq)

# ------------------ 输出最终总谱 ------------------
out_dir = f"./ER_{ER_MIN}_{ER_MAX}keV_{bins_ER}bins"
os.makedirs(out_dir, exist_ok=True)

def convert_to_list(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj

bins_er = np.linspace(ER_MIN, ER_MAX, bins_ER+1)
with open(os.path.join(out_dir, "Gamma_er.json"), "w") as f:
    json.dump(convert_to_list(total_counts), f)
with open(os.path.join(out_dir, "Gamma_er_error.json"), "w") as f:
    json.dump(convert_to_list(total_error), f)
with open(os.path.join(out_dir, "Gamma_er_bin.json"), "w") as f:
    json.dump(convert_to_list(bins_er), f)

print(f"Gamma ER spectrum saved to {out_dir}")
print(f"Gamma counts sum: {np.sum(total_counts)/bins_ER} DRU")
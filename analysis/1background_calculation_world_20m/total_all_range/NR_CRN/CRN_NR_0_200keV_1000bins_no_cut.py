import sys
import os
import json
import h5py
import numpy as np

# ==================== 参数设置 ====================
ER_MAX = 20
ER_MIN = 0
bins_ER = 20
NR_MIN = 0
NR_MAX = 200
bins_NR = 1000

print(f'NR energy range: {NR_MIN} - {NR_MAX} keV, bins: {bins_NR}')

# 几何参数（仅用于kg计算）
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

# # 归一化常数
# days = 24 * 3600
# kg = (fv_z1 - fv_z2)/10 * (fv_r/10) * (fv_r/10) * np.pi * 2.862 / 1000
# print(f'LXe mass: {kg:.3f} kg')
# print('-----------------------')

# 读取归一化文件
with open('../../normalization.json', 'r') as f:
    genNormalize = json.loads(f.read())

# 归一化常数
days = 24 * 3600
kg = genNormalize['XenonDetector']['mass']
print(f'LXe mass: {kg} kg')
print('-----------------------')

# ------------------ 处理单个源（无任何cut，直接使用所有max_nr） ------------------
def process_source(folder, key, bins_nr):
    # 读取事件数据
    with h5py.File(os.path.join(folder, f'events_{key}.h5'), 'r') as ipt:
        e = ipt['events'][:]

    # 所有事件的 NR 能量（无任何 cut）
    energies = e['max_nr']

    # 计算加权直方图
    bin_edges = bins_nr
    bin_widths = np.diff(bin_edges)
    # 每个事件所属 bin 索引
    dig = np.clip(np.digitize(energies, bin_edges) - 1, 0, len(bin_widths)-1)
    # 权重因子 = activity * days / factor / kg
    weight_factor = genNormalize[key]['activity'] * days / genNormalize[key]['factor'] / kg
    # 每个事件的权重（已包含 bin 宽度归一化）
    weights = weight_factor / bin_widths[dig]

    # 直方图 counts (单位: 计数/kg/天/bin)
    counts, _ = np.histogram(energies, bins=bin_edges, weights=weights)

    # 误差计算（与参考代码完全一致）
    # c = 1 / (kg * factor / (activity * days))
    c = 1 / (kg * genNormalize[key]['factor'] / (genNormalize[key]['activity'] * days))
    error = np.sqrt(counts * c / bin_widths)

    return counts, error

# ------------------ 主程序 ------------------
bins_nr = np.linspace(NR_MIN, NR_MAX, bins_NR+1)

# 处理 muon
# folder_muon = '../../../../result/muon_900M_6m'
# muon_counts, muon_error = process_source(folder_muon, 'muon', bins_nr)

# 处理 CRN
folder_crn = '../../../../result/CRN_500M_6m'
crn_counts, crn_error = process_source(folder_crn, 'CRN', bins_nr)

# 合并
total_counts = crn_counts
total_error = np.sqrt( crn_error**2)

# 输出目录
out_dir = f"./NR_CRN_{NR_MIN}_{NR_MAX}keV_{bins_NR}bins"
os.makedirs(out_dir, exist_ok=True)

def convert_to_list(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj

with open(os.path.join(out_dir, "CRN_nr.json"), "w") as f:
    json.dump(convert_to_list(total_counts), f)
with open(os.path.join(out_dir, "CRN_nr_error.json"), "w") as f:
    json.dump(convert_to_list(total_error), f)
with open(os.path.join(out_dir, "CRN_nr_bin.json"), "w") as f:
    json.dump(convert_to_list(bins_nr), f)

print(f"CRN NR spectrum saved to {out_dir}")
print(f"Sum of CRN_counts: {np.sum(total_counts)/bins_NR} DRU")
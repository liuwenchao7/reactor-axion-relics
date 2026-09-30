import os
import json
import h5py
import numpy as np

print('Run')

# ======================================
# 【核心可配置参数区】
ER_MIN = 0.0        # 能量最小值 (keV)
ER_MAX = 200.0       # 能量最大值 (keV)
N_BINS_X = 800       # X轴分bin数量
# Z位置配置
Z_MIN = -160.0      # Z坐标最小值 (mm)
Z_MAX = 0.0         # Z坐标最大值 (mm)
N_BINS_Y = 60       # Y轴分bin数量
# 输出文件夹
OUTPUT_FOLDER_NAME = f"ER_{int(ER_MIN)}_{int(ER_MAX)}keV_{N_BINS_X}bins_material"
# ======================================

# ===================== 探测器几何参数 =====================
position_lxe = -2.5
water_shift_z = 0
solid_shift_x = -25

dis_gate = 0.25
dis_cathode = 15.3575
dis_botscreen = 0.2625

fv_r = 80
fv_z1 = (position_lxe + water_shift_z) * 10
fv_z2 = (position_lxe + water_shift_z - dis_gate - dis_cathode) * 10

# 切割阈值
cut_er1 = 0.05
cut_nr_er1 = 0.05

# ===================== 自动生成分bin边界 =====================
xedges = np.linspace(ER_MIN, ER_MAX, N_BINS_X + 1)
liquid_level = (position_lxe + water_shift_z) * 10
yedges_real = np.linspace(Z_MIN + liquid_level, Z_MAX + liquid_level, N_BINS_Y + 1)
yedges_target = np.linspace(Z_MIN, Z_MAX, N_BINS_Y + 1)

# ===================== 归一化参数 =====================
days = 24 * 3600
kg = (fv_z1 - fv_z2)/10 * (fv_r/10) * (fv_r/10) * np.pi * 2.862 / 1000
print('LXe mass:', kg, 'kg')
print('-----------------------')

# ===================== 切割函数 =====================
def er_energy(events):
    mask = events['max_er'] < ER_MAX
    mask &= events['max_er'] >= ER_MIN
    return mask

def cut_er_max_alt(events):
    mask = events['alt_er'] > cut_er1 * events['max_er']
    return ~mask

def cut_er(events):
    return cut_er_max_alt(events)

def cut_fv(events, s):
    mask = np.sqrt((events[f'max_{s}_x'] - solid_shift_x)**2 + events[f'max_{s}_y']**2) < fv_r
    mask &= events[f'max_{s}_z'] < fv_z1
    mask &= events[f'max_{s}_z'] > fv_z2
    return mask

# ===================== 核心约束：仅保留指定切割条件 =====================
def get_er_mask(events):
    e = events
    # 你指定的唯一cut条件
    mask_er = cut_fv(e, 'er') & cut_er(e) & er_energy(e)
    return mask_er

# ===================== 材料本底计算函数 =====================
def calc_material_background(folder, genNormalize, bins_er):
    total_er_1d = np.zeros(len(bins_er)-1)
    total_er_2d = np.zeros((N_BINS_Y, N_BINS_X))
    ER_rate_total = 0.0
    ER_error_total = 0.0

    # 定义所有材料和核素组合
    material_groups = [
        (['U238','Th232','Co60','K40','Cs137'], ['steel']),
        (['U238','Th232','Co60','K40','Cs137'], ['Teflon']),
        (['U238','U235','Th232','Cs137','K40','Co60'], ['PMTbase']),
        (['U238','U235','Th232','Cs137','K40','Co60'], ['PMTwindow']),
        (['U238','Th232','Co60','K40','Cs137','U235'], ['Flange']),
    ]

    for isotopes, components in material_groups:
        for iso in isotopes:
            for comp in components:
                file_path = os.path.join(folder, f'events_{iso}_{comp}.h5')
                if not os.path.exists(file_path):
                    continue
                
                with h5py.File(file_path, 'r', libver='latest', swmr=True) as ipt:
                    events = ipt['events'][:]
                
                # 获取筛选mask
                mask = get_er_mask(events)
                er_energy_vals = events['max_er'][mask]
                er_z_vals = events['max_er_z'][mask]

                # 权重计算（与原代码完全一致）
                norm_key = f"{iso}_{comp}"
                factor = genNormalize[norm_key]['activity'] * days / genNormalize[norm_key]['factor']
                weight_per = factor / kg
                weights = np.full(len(er_energy_vals), weight_per)

                # 一维能谱
                h1d, _ = np.histogram(er_energy_vals, bins=bins_er, weights=weights)
                total_er_1d += h1d

                # 二维Z-能量谱
                h2d, _, _ = np.histogram2d(er_energy_vals, er_z_vals, bins=[bins_er, yedges_real], weights=weights)
                total_er_2d += h2d.T

                # 总率统计
                ER_count = len(er_energy_vals)
                rate = ER_count * factor / kg
                err = np.sqrt(ER_count) * factor / kg
                ER_rate_total += rate
                ER_error_total += err**2

    ER_error_total = np.sqrt(ER_error_total)
    return total_er_1d, total_er_2d, ER_rate_total, ER_error_total

# ===================== 主流程 =====================
print('finish the function, then read the data')

# 读取归一化文件
with open('../normalization.json', 'r') as f:
    genNormalize = json.load(f)

# 数据路径（原代码完全保留）
folder = '/public/home/yyy/RELICS5/RELICS10_v2_0_basic_design_260304/RelicsSim_v2_solid_origin_material_R_E/result/Material_10M_steel_window_1M_others'

# 执行计算
er_1d, er_2d, total_er_rate, total_er_err = calc_material_background(folder, genNormalize, xedges)

# ===================== 保存输出文件 =====================
os.makedirs(OUTPUT_FOLDER_NAME, exist_ok=True)
os.makedirs(f"{OUTPUT_FOLDER_NAME}/temp", exist_ok=True)

# 保存 ER_total.json (标准格式)
er_total_dict = {
    "type": "ER_total_material",
    "unit": "kg^-1 day^-1",
    "x_label": "Electronic Recoil Energy [keVer]",
    "y_label": "z position [mm]",
    "xedges": xedges.tolist(),
    "yedges": yedges_target.tolist(),
    "rate_map": er_2d.tolist()
}
with open(os.path.join(OUTPUT_FOLDER_NAME, 'ER_total.json'), 'w') as f:
    json.dump(er_total_dict, f, indent=2)

# 保存一维谱、误差、bin
convert_to_list = lambda obj: obj.tolist() if isinstance(obj, np.ndarray) else obj
with open(f"{OUTPUT_FOLDER_NAME}/material_er.json", 'w') as f:
    json.dump(convert_to_list(er_1d), f)
with open(f"{OUTPUT_FOLDER_NAME}/material_er_bin.json", 'w') as f:
    json.dump(convert_to_list(xedges), f)

# ===================== 结果输出 =====================
er_rate_per_keV = total_er_rate / (ER_MAX - ER_MIN)
er_err_per_keV = total_er_err / (ER_MAX - ER_MIN)
er_rate_mdru = er_rate_per_keV * 1000
er_err_mdru = er_err_per_keV * 1000

print('='*50)
print(f'Material ER rate: {er_rate_per_keV:.6f} /kg/days/keV')
print(f'Material ER error: {er_err_per_keV:.6f}')
print(f'Material ER rate(mDRU): {er_rate_mdru:.6f}')
print(f'Material ER error(mDRU): {er_err_mdru:.6f}')
print(f'bins_er=np.linspace({ER_MIN},{ER_MAX},{N_BINS_X+1})')
print(f'FV mass: {kg:.2f} kg')
print('='*50)
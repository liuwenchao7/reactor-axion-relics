import sys
from scipy import optimize as opti
from scipy.stats import poisson
from scipy.interpolate import interp1d

import os
import copy
import math
from multiprocessing import Pool, cpu_count

from tqdm import tqdm
import h5py
import numpy as np
import pandas as pd

from scipy.integrate import quad
from statsmodels.stats.proportion import proportion_confint
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
# import radioactivedecay as rd

plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 100
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['font.size'] = 15

# plt.style.use('/public/home/yyy/relics_picture_style_huayu/relics_er.mplstyle')

import numpy as np
print('Run')

# ==================== 可修改参数 ====================
# ER 参数（若需要，保留原样）
ER_MAX_0 = 20
ER_MAX = int(ER_MAX_0)
bins_ER_0 = 20
bins_ER = int(bins_ER_0)

# NR 参数：新增下限和上限（单位 keV）
NR_MIN_0 = 4700          # NR 能量下限 (keV)
NR_MIN = int(NR_MIN_0)
NR_MAX_0 = 30000         # NR 能量上限 (keV)
NR_MAX = int(NR_MAX_0)
bins_NR_0 = 253         # bin 数量
bins_NR = int(bins_NR_0)

print(f'NR energy range: {NR_MIN} - {NR_MAX} keV, bins: {bins_NR}')
# ====================================================

position_lxe = -37 # 单位cm，包含了相对于固体屏蔽体的位移和相对于探测器外层结构的位移
water_shift_z = 0 # 单位cm，包含相对于水屏蔽体的位移
solid_shift_x = -150 # 单位mm

dis_gate = 0.25 # _to_lxe_top, cm
dis_cathode = 15.3575 # _to_gate, cm
dis_botscreen = 0.2625# _to_xedetector_bottom, cm
dis_fv_top_to_lxe = 0.35 #单位cm，因此距离gate为1mm
fv_r = 80              #主探测器中的半径
fv_z1= (position_lxe + water_shift_z) * 10              #主探测器的最上端，z_max
fv_z2= (position_lxe + water_shift_z - dis_gate - dis_cathode) *10                                           #主探测器的最小端，z_min

veto1=200          #veto中总共最大er能量
veto2=1000          #veto中总共最大nr能量

#veto有位置的，并不是所有LXe都是veto
veto_z1 = (position_lxe + water_shift_z) * 10
veto_z2 = (position_lxe + water_shift_z - dis_gate - dis_cathode - dis_botscreen) * 10

cut_nr2=0.1            #NR single中   最大的alt的下限（只需要第二个信号幅度过大就可以cut）
cut_er1=0.05          #ER single中   最大的ER与第二大的ER之间的关系
cut_nr_er1=0.05         #NR与ER    伴随NR信号同时产生的ER信号的上限

#veto_cut     针对veto层中，设定ER和NR上限
def cut_veto_er(events):
    mask = events['veto_max_er'] > veto1
    mask &= events['veto_max_er_z'] < veto_z1 
    mask &= events['veto_max_er_z'] > veto_z2
    return ~mask

def cut_veto_nr(events):
    mask = events['veto_max_nr'] > veto2
    mask &= events['veto_max_nr_z'] < veto_z1 
    mask &= events['veto_max_nr_z'] > veto_z2
    return ~mask

def cut_veto(events):
    mask = cut_veto_er(events) & cut_veto_nr(events)
    return mask

#NR_ER       主探测器中同时有NR与ER信号，ER信号有上限
def cut_nr_er(events):
    mask = events['max_er'] > cut_nr_er1
    return ~mask

#NR_SS       针对主探测器中的NR信号   设定第二大NR信号不能超过一定的幅值    
def cut_nr_alt_max(events):
    mask = events['alt_nr'] > cut_nr2
    return ~mask
def cut_nr(events):
    mask = cut_nr_alt_max(events)
    return mask

#ER_SS      针对主探测器中的ER信号   设定最大ER与第二大ER信号的关系   最大ER的范围
def cut_er_max_alt(events):
    mask = events['alt_er'] > cut_er1* events['max_er']
    return ~mask
def cut_er(events):
    mask = cut_er_max_alt(events) 
    return mask

def cut_fv(events,s):
    mask = np.sqrt((events[f'max_{s}_x']-solid_shift_x) ** 2 + events[f'max_{s}_y'] ** 2) < fv_r
    mask &= events[f'max_{s}_z'] < fv_z1
    mask &= events[f'max_{s}_z'] > fv_z2
    return mask 

def cut_nr_t(events):
    mask = events['max_nr_t'] < 1e8
    return mask

def cut_er_t(events):
    mask = events['max_er_t'] < 1e8
    return mask

def muon_veto_nr(events):
    mask = events['max_nr_t'] < 350*10**(-6)
    return mask

def muon_veto_er(events):
    mask = events['max_er_t'] < 350*10**(-6)
    return mask

def muon_rate_error(events,nr_max,nr_min,er_max):
    NR_number=0
    ER_number=0
    e=events
    mask=e['max_nr'] < nr_max
    mask&=e['max_nr'] > nr_min
    NR_number_veto=len(e['max_nr'][muon_veto_nr(e) & mask])
    NR_number = NR_number_veto
    return NR_number, ER_number

def rate_error(events, nr_max, nr_min, er_max):
    e = events
    mask = (e['max_nr'] < nr_max) & (e['max_nr'] > nr_min)
    NR_number = len(e['max_nr'][mask])
    ER_number = 0
    return NR_number, ER_number

#归一化，需要先跑好归一化python文件
days=24*3600
kg= (fv_z1 - fv_z2)/10* fv_r/10 *  fv_r/10 * np.pi *2.862 /1000
print('LXe mass:', kg, 'kg')
print('-----------------------')
    
def muon(folder,genNormalize,bins_nr,bins_er,nr_max,nr_min,er_max):
    kg= (fv_z1 - fv_z2)/10* fv_r/10 *  fv_r/10 * np.pi *2.862 /1000
    days=24*3600
    muon_veto=0.01   # 保留变量，但不再用于修正
    bins_nr=bins_nr
    bins_er=bins_er
    muon_data_nr=dict()
    muon_data_er=dict()
    muon_er_error=dict()
    NR_rate=dict()
    NR_error=dict()
    ER_rate=dict()
    ER_error=dict()

    files = [os.path.join(folder, f'events_muon.h5')]
    with h5py.File(files[0], 'r', libver='latest', swmr=True) as ipt:
        events = ipt['events'][:]
    e=events
    
    # 计算全谱（无 veto 压制）
    weights = np.full(len(e), genNormalize['muon']['activity']*days / genNormalize['muon']['factor']) / np.diff(bins_nr)[np.clip(np.digitize(e['max_nr'], bins_nr) - 1, 0, len(bins_nr) - 2)]/kg
    h = np.histogram(e['max_nr'], weights=weights, bins=bins_nr) 
    muon_nr_all = h[0]
    
    # 去掉 veto 修正，直接使用原始全谱
    muon_data_nr['muon'] = muon_nr_all
    
    # 总计数率（无 veto）
    NR_number_all, _ = rate_error(e, nr_max, nr_min, er_max)
    scale = (genNormalize['muon']['activity'] * days / genNormalize['muon']['factor']) / kg
    NR_rate['muon'] = NR_number_all * scale
    NR_error['muon'] = np.sqrt(NR_number_all) * scale
    
    # 占位返回（为了保持函数返回元组长度不变）
    muon_nr_error = np.zeros(len(bins_nr)-1)
    muon_er_error = {}
    
    return events, muon_data_nr, muon_data_er, muon_nr_error, muon_er_error, NR_rate, NR_error, ER_rate, ER_error

def cosmic_ray_neutron(folder,genNormalize,bins_nr,bins_er,nr_max,nr_min,er_max):
    kg= (fv_z1 - fv_z2)/10* fv_r/10 *  fv_r/10 * np.pi *2.862 /1000
    days=24*3600
    bins_nr=bins_nr
    bins_er=bins_er
    cosmic_data_nr=dict()
    cosmic_nr_error=dict()
    cosmic_data_er=dict()
    cosmic_er_error=dict()
    NR_rate=dict()
    NR_error=dict()
    ER_rate=dict()
    ER_error=dict()

    files = [os.path.join(folder, f'events_CRN.h5')]
    with h5py.File(files[0], 'r', libver='latest', swmr=True) as ipt:
        events = ipt['events'][:]
    e=events
    
    weights = np.full(len(e), genNormalize['CRN']['activity']*days / genNormalize['CRN']['factor']) / np.diff(bins_nr)[np.clip(np.digitize(e['max_nr'], bins_nr) - 1, 0, len(bins_nr) - 2)]/kg
    h = np.histogram(e['max_nr'], weights=weights, bins=bins_nr) 
    cosmic_nr_all = h[0]
    
    cosmic_data_nr['CRN'] = cosmic_nr_all
    
    NR_number_all, _ = rate_error(e, nr_max, nr_min, er_max)
    scale = (genNormalize['CRN']['activity'] * days / genNormalize['CRN']['factor']) / kg
    NR_rate['CRN'] = NR_number_all * scale
    NR_error['CRN'] = np.sqrt(NR_number_all) * scale
    
    cosmic_nr_error = np.zeros(len(bins_nr)-1)
    cosmic_er_error = {}
    
    return events, cosmic_data_nr, cosmic_data_er, cosmic_nr_error, cosmic_er_error, NR_rate, NR_error, ER_rate, ER_error

print('finish the function,then read the data')

# 归一化文件
import json
f = open('../../normalization.json', 'r')
content = f.read()
genNormalize = json.loads(content)

# 公共参数
bins_nr = np.linspace(NR_MIN, NR_MAX, bins_NR+1)
nr_max = 1.36
nr_min = 0.63
er_max = float(ER_MAX)
bins_er = np.linspace(0, er_max, bins_ER+1)

# 处理 muon
folder_muon = f'../../../../result/muon_900M_6m'
cosmic_events_muon, muon_data_nr, muon_data_er, muon_nr_error, muon_er_error, NR_rate_muon, NR_error_muon, ER_rate_muon, ER_error_muon = muon(
    folder_muon, genNormalize, bins_nr, bins_er, nr_max, nr_min, er_max)

# 处理 CRN
folder_crn = f'../../../../result/CRN_500M_6m'
cosmic_events_crn, crn_data_nr, crn_data_er, crn_nr_error, crn_er_error, NR_rate_crn, NR_error_crn, ER_rate_crn, ER_error_crn = cosmic_ray_neutron(
    folder_crn, genNormalize, bins_nr, bins_er, nr_max, nr_min, er_max)

# 合并能谱
total_counts = muon_data_nr['muon'] + crn_data_nr['CRN']
total_error = np.sqrt(muon_nr_error**2 + crn_nr_error**2)

# 合并总计数率
total_rate = NR_rate_muon['muon'] + NR_rate_crn['CRN']
total_rate_error = np.sqrt(NR_error_muon['muon']**2 + NR_error_crn['CRN']**2)

# 输出目录
out_dir = f"./NR_{NR_MIN}_{NR_MAX}keV_{bins_NR}bins"
os.makedirs(out_dir, exist_ok=True)
temp_dir = os.path.join(out_dir, "temp")
os.makedirs(temp_dir, exist_ok=True)

# 保存总谱
def convert_to_list(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: convert_to_list(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_list(item) for item in obj]
    else:
        return obj

json_str = json.dumps(convert_to_list(total_counts))
with open(os.path.join(out_dir, "total_nr.json"), "w") as f:
    f.write(json_str)

json_str = json.dumps(convert_to_list(total_error))
with open(os.path.join(out_dir, "total_nr_error.json"), "w") as f:
    f.write(json_str)

json_str = json.dumps(convert_to_list(bins_nr))
with open(os.path.join(out_dir, "total_nr_bin.json"), "w") as f:
    f.write(json_str)

# 可选：保存各组分（保留原代码输出）
json_str = json.dumps(convert_to_list(muon_data_nr))
with open(os.path.join(temp_dir, "muon_data.json"), "w") as f:
    f.write(json_str)
json_str = json.dumps(convert_to_list(NR_rate_muon))
with open(os.path.join(temp_dir, "muon_NR_rate.json"), "w") as f:
    f.write(json_str)
json_str = json.dumps(convert_to_list(NR_error_muon))
with open(os.path.join(temp_dir, "muon_NR_error.json"), "w") as f:
    f.write(json_str)

json_str = json.dumps(convert_to_list(crn_data_nr))
with open(os.path.join(temp_dir, "CRN_data.json"), "w") as f:
    f.write(json_str)
json_str = json.dumps(convert_to_list(NR_rate_crn))
with open(os.path.join(temp_dir, "CRN_NR_rate.json"), "w") as f:
    f.write(json_str)
json_str = json.dumps(convert_to_list(NR_error_crn))
with open(os.path.join(temp_dir, "CRN_NR_error.json"), "w") as f:
    f.write(json_str)

# 输出信息
print(f"Total NR rate: {total_rate:.3f} /kg/days (error: {total_rate_error:.3f})")
print(f"Output directory: {out_dir}")
print(f"bins_nr = np.linspace({NR_MIN}, {NR_MAX}, {bins_NR}+1)")
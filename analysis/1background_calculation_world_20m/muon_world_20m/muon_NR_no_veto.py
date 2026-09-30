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

plt.style.use('/public/home/yyy/relics_picture_style_huayu/relics_er.mplstyle')

import numpy as np
print('Run')

# 手动输入ER最大值
ER_MAX_0 = 20
print(type(ER_MAX_0))

ER_MAX = int(ER_MAX_0)
print(type(ER_MAX))

# 手动输入bins_ER数
bins_ER_0 = 20
print(type(ER_MAX_0))

bins_ER = int(bins_ER_0)
print(type(bins_ER))


# 手动输入bins_NR数
bins_NR_0 = 100

bins_NR = int(bins_NR_0)
print(type(bins_NR))

position_lxe = -37 # 单位cm，包含了相对于固体屏蔽体的位移和相对于探测器外层结构的位移

water_shift_z = 0 # 单位cm，包含相对于水屏蔽体的位移

solid_shift_x = -150 # 单位mm

dis_gate = 0.25 # _to_lxe_top, cm
dis_cathode = 15.3575 # _to_gate, cm
dis_botscreen = 0.2625# _to_xedetector_bottom, cm
dis_fv_top_to_lxe = 0.35 #单位cm，因此距离gate为1mm
fv_r = 80              #主探测器中的半径
fv_z1= (position_lxe + water_shift_z) * 10              #主探测器的最上端，z_max
fv_z2= (position_lxe + water_shift_z - dis_gate - dis_cathode) *10                                     #主探测器的最小端，z_min  （只有s2信号无法得到z方向上的位置重建）

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


#cut_fv         选取一个圆柱体，设定相应的r，z_max，z_min
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

# ===== 已移除 muon_veto_nr 和 muon_veto_er =====

def rate_error(events, nr_max, nr_min, er_max):
    """统计所有通过基本挑选的事件（不考虑 muon veto）"""
    e = events
    # NR 计数
    mask = (e['max_nr'] < nr_max) & (e['max_nr'] > nr_min)
    NR_number = len(e['max_nr'][cut_fv(e, 'nr') & cut_nr_er(e) & cut_nr(e) & cut_nr_t(e) & mask])
    # ER 计数（若需要可取消注释）
    # mask_er = e['max_er'] < er_max
    # ER_number = len(e['max_er'][cut_fv(e, 'er') & cut_er(e) & cut_er_t(e) & mask_er])
    ER_number = 0  # 暂不处理 ER
    return NR_number, ER_number

# 归一化参数（需提前运行归一化文件得到）
days = 24 * 3600
kg = (fv_z1 - fv_z2) / 10 * fv_r / 10 * fv_r / 10 * np.pi * 2.862 / 1000
print('LXe mass:', kg, 'kg')
print('-----------------------')
    
def muon(folder, genNormalize, bins_nr, bins_er, nr_max, nr_min, er_max):
    kg = (fv_z1 - fv_z2) / 10 * fv_r / 10 * fv_r / 10 * np.pi * 2.862 / 1000
    days = 24 * 3600
    bins_nr = bins_nr
    bins_er = bins_er

    muon_data_nr = dict()
    muon_data_er = dict()
    muon_nr_error = dict()
    muon_er_error = dict()
    NR_rate = dict()
    NR_error = dict()
    ER_rate = dict()
    ER_error = dict()

    files = [os.path.join(folder, f'events_muon.h5')]
    with h5py.File(files[0], 'r', libver='latest', swmr=True) as ipt:
        events = ipt['events'][:]
        # primaries = ipt['primaries'][:]
    e = events

    # 对 NR 能谱进行加权直方图（不应用 muon veto）
    weights = np.full(len(e), genNormalize['muon']['activity'] * days / genNormalize['muon']['factor']) / np.diff(bins_nr)[np.clip(np.digitize(e['max_nr'], bins_nr) - 1, 0, len(bins_nr) - 2)] / kg
    h = np.histogram(e['max_nr'][cut_fv(e, 'nr') & cut_nr_er(e) & cut_nr(e) & cut_nr_t(e)],
                     weights=weights[cut_fv(e, 'nr') & cut_nr_er(e) & cut_nr(e) & cut_nr_t(e)],
                     bins=bins_nr)
    muon_nr_all = h[0]

    # ER 部分（注释，如需可参照 NR 方式处理）
    # weights_er = ...
    # muon_er_all = ...

    # ===== 去除 muon veto 相关计算，直接使用全部事件 =====
    muon_data_nr['muon'] = muon_nr_all

    # 误差计算：仅基于泊松统计（sqrt(N)），乘以归一化因子再除以 bin 宽度
    bins_diff = np.diff(bins_nr)
    muon_activity = genNormalize['muon']['factor'] / (genNormalize['muon']['activity'] * days)
    # 注意：weights 中已经包含了归一化，但这里我们直接计算误差
    # 正确做法：每个 bin 的事件数（未加权）为 N_raw，但直方图已经加权，我们需要原始计数来计算误差
    # 更严谨：重新获取未加权直方图；但为简化，用加权后的值近似（因为权重相同）
    # 我们使用原始计数（未加权）来计算统计误差，然后乘以归一化因子/宽度
    # 重新计算未加权直方图
    h_raw = np.histogram(e['max_nr'][cut_fv(e, 'nr') & cut_nr_er(e) & cut_nr(e) & cut_nr_t(e)], bins=bins_nr)
    raw_counts = h_raw[0]
    norm_factor = genNormalize['muon']['activity'] * days / genNormalize['muon']['factor'] / kg
    muon_nr_error = np.sqrt(raw_counts) * norm_factor / bins_diff   # 除以 bin 宽度得到谱密度误差

    # 总计数率（不区分 veto）
    NR_number_all, _ = rate_error(e, nr_max, nr_min, er_max)
    NR_rate['muon'] = NR_number_all * (genNormalize['muon']['activity'] * days / genNormalize['muon']['factor']) / kg
    NR_error['muon'] = np.sqrt(NR_number_all) * (genNormalize['muon']['activity'] * days / genNormalize['muon']['factor']) / kg

    # ER 部分暂不处理，置为零
    ER_rate['muon'] = 0
    ER_error['muon'] = 0

    return events, muon_data_nr, muon_data_er, muon_nr_error, muon_er_error, NR_rate, NR_error, ER_rate, ER_error

print('finish the function,then read the data')
#归一化，需要先跑好归一化python文件

days = 24 * 3600
kg = (fv_z1 - fv_z2) / 10 * fv_r / 10 * fv_r / 10 * np.pi * 2.862 / 1000
#调用了活度，其中有0的需要筛选

import json
f = open('../normalization.json', 'r')
content = f.read()
genNormalize = json.loads(content)

#the position of cosmic data
folder = f'../../../result/muon_900M_6m'
bins_nr = np.logspace(-4, 7, bins_NR + 1)
nr_max = 1.36
nr_min = 0.63
er_max = float(ER_MAX)
bins_er = np.linspace(0, er_max, bins_ER + 1)
muon_events, muon_data_nr, muon_data_er, muon_nr_error, muon_er_error, NR_rate, NR_error, ER_rate, ER_error = muon(folder, genNormalize, bins_nr, bins_er, nr_max, nr_min, er_max)

print('finish reading data,then output the message')

#存数据用
# 将 Numpy 数组转换为列表的递归函数
def convert_to_list(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: convert_to_list(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_list(item) for item in obj]
    else:
        return obj
    
# 将字典转换为 JSON 字符串，并将其中的 Numpy 数组转换为列表
# 这里是直接存所有材料的数据，新建一个临时的temp_yyy，存用pandax的活度模拟的结果，如果要用xenon100的，单独存在temp中，记得在这里改写
json_str = json.dumps(convert_to_list(muon_data_nr))
with open(f"./NR_{bins_NR}bins/temp/muon_data.json", "w") as f:
    f.write(json_str) 
    
json_str = json.dumps(convert_to_list(NR_rate))
with open(f"./NR_{bins_NR}bins/temp/muon_NR_rate.json", "w") as f:
    f.write(json_str)
    
json_str = json.dumps(convert_to_list(NR_error))
with open(f"./NR_{bins_NR}bins/temp/muon_NR_error.json", "w") as f:
    f.write(json_str)

#进行数据的处理，得到汇总后的结果
#material data,figure
Material_Data = sum(muon_data_nr.values())
error = np.zeros(len(Material_Data))

bins = bins_nr
bins_diff = np.diff(bins)
kg = (fv_z1 - fv_z2) / 10 * fv_r / 10 * fv_r / 10 * np.pi * 2.862 / 1000 # 单位kg

for i in range(len(Material_Data)):
    isotopes = ['cosmic']
    components = ['muon']
    for j in isotopes:
        for k in components:
            c = 1 / (kg * genNormalize[k]['factor'] / (genNormalize[k]['activity'] * days))
            error[i] += (muon_data_nr[k][i] * c / bins_diff[i])
            #print(error)
   
error = np.sqrt(error)
Material_Data = sum(muon_data_nr.values())

#ER data，存储ER本底数据，xenon存在background的material文件夹下，pandax存在material_pandax文件夹下
json_str = json.dumps(convert_to_list(Material_Data))
with open(f"./NR_{bins_NR}bins/muon_nr.json", "w") as f:
    f.write(json_str)

json_str = json.dumps(convert_to_list(error))
with open(f"./NR_{bins_NR}bins/muon_nr_error.json", "w") as f:
    f.write(json_str)
    
json_str = json.dumps(convert_to_list(bins_nr))
with open(f"./NR_{bins_NR}bins/muon_nr_bin.json", "w") as f:
    f.write(json_str)

print('muon NR rate:', NR_rate, '/kg/days')
print('muon NR error:', NR_error)
print('bins_nr=np.logspace(-4, 7, 101)')
print('bins_er=np.linspace(0,20,21)')
print('max ER', er_max)
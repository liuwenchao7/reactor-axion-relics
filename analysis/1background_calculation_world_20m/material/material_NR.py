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


import numpy as np

# 手动输入ER最大值
ER_MAX_0 = 20
ER_MAX = int(ER_MAX_0)

# 手动输入ER最小值
ER_MIN_0 = 0
ER_MIN = int(ER_MIN_0)

# 手动输入bins_ER数
bins_ER_0 = 20
bins_ER = int(bins_ER_0)

# 手动输入bins_NR数
bins_NR_0 = 100
bins_NR = int(bins_NR_0)
print(type(bins_NR))

position_lxe = -2.5 # 单位cm，包含了相对于固体屏蔽体的位移和相对于探测器外层结构的位移

water_shift_z = 0 # 单位cm，包含相对于水屏蔽体的位移

solid_shift_x = -25 # 单位mm

dis_gate = 0.25 # _to_lxe_top, cm
dis_cathode = 15.3575 # _to_gate, cm
dis_botscreen = 0.2625# _to_xedetector_bottom, cm
dis_fv_top_to_lxe = 0.35 #单位cm，因此距离gate为1mm
fv_r = 70              #主探测器中的半径
fv_z1= (position_lxe + water_shift_z) * 10              #主探测器的最上端，z_max
fv_z2= (position_lxe + water_shift_z - dis_gate - dis_cathode) *10         

veto1=200          #veto中总共最大er能量
veto2=1000          #veto中总共最大nr能量

#veto有位置的，并不是所有LXe都是veto
veto_z1 = (position_lxe + water_shift_z) * 10
veto_z2 = (position_lxe + water_shift_z - dis_gate - dis_cathode - dis_botscreen) * 10

cut_nr2=0.1            #NR single中   最大的alt的下限（只需要第二个信号幅度过大就可以cut）


cut_er1=0.05          #ER single中   最大的ER与第二大的ER之间的关系

cut_nr_er1=0.05         #NR与ER    伴随NR信号同时产生的ER信号的上限

#energy cut   必须大于1keV
def er_energy(events):
    mask=events['max_er']<ER_MAX
    mask&=events['max_er']>= ER_MIN
    return mask

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

def rate_error(events,nr_max,nr_min,er_max):
    NR_number=0
    ER_number=0
    e=events
    # #NR,nr_max=1,nr_min=0.1
    mask=e['max_nr'][cut_fv(e,'nr') & cut_nr_er(e) & cut_nr(e)] < nr_max
    mask&=e['max_nr'][cut_fv(e,'nr') & cut_nr_er(e) & cut_nr(e)] > nr_min
    NR_number=len(e['max_nr'][cut_fv(e,'nr') & cut_nr_er(e) & cut_nr(e)][mask])
    
    # #ER,er_max=100
    # mask=e['max_er'][cut_fv(e,'er') & cut_er(e)] <er_max
    # mask&=e['max_er'][cut_fv(e,'er') & cut_er(e)] >= ER_MIN
    # ER_number=len(e['max_er'][cut_fv(e,'er') & cut_er(e)][mask])
    return NR_number,ER_number

#用于新模拟数据分析

def material(folder,genNormalize,bins_nr,nr_max,nr_min,er_max):
    kg= (fv_z1 - fv_z2)/10* fv_r/10 *  fv_r/10 * np.pi *2.862 /1000 # 单位kg
    days=24*3600
    bins = bins_nr
    material_data=dict()
    material_data_error=dict()
    events = dict()
    primaries = dict()
    NR_rate=dict()
    NR_error=dict()
    ER_rate=dict()
    ER_error=dict()

    isotopes=['U238','Th232','Co60','K40','Cs137']
    # components=['steel', 'copper']
    components=['steel']
    for i in isotopes:
        for j in components:
            files = [os.path.join(folder, f'events_{i}_{j}.h5')]
            with h5py.File(files[0], 'r', libver='latest', swmr=True) as ipt:
                events[i+j] = ipt['events'][:]
            e=events[i+j]

            weights = np.full(len(e), genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor']) / np.diff(bins)[np.clip(np.digitize(e['max_nr'], bins) - 1, 0, len(bins) - 2)]/kg
            h = np.histogram(e['max_nr'][cut_fv(e,'nr') & cut_nr_er(e) & cut_nr(e)], weights=weights[cut_fv(e,'nr') & cut_nr_er(e) & cut_nr(e)],bins=bins)
            material_data[i+j]=h[0]

            #to material error
            bins = bins_nr
            bins_diff=np.diff(bins)
            material_activity=genNormalize[i+'_'+j]['factor']/(genNormalize[i+'_'+j]['activity']*days)
            
            c=1/(material_activity*kg)
            material_data_error[i+j]= [ (a * c / b)**(1/2)  for a, b in zip(material_data[i+j], bins_diff)]
            
            NR_number,ER_number=rate_error(e,nr_max,nr_min,er_max)
            NR_rate[i+j]=NR_number*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            NR_error[i+j]=(NR_number)**(1/2)*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            # ER_rate[i+j]=ER_number*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            # ER_error[i+j]=(ER_number)**(1/2)*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg

    isotopes=['U238','Th232','Co60','K40','Cs137']
    components=['Teflon']
    for i in isotopes:
        for j in components:
            files = [os.path.join(folder, f'events_{i}_{j}.h5')]
            with h5py.File(files[0], 'r', libver='latest', swmr=True) as ipt:
                events[i+j] = ipt['events'][:]
            e=events[i+j]

            weights = np.full(len(e), genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor']) / np.diff(bins)[np.clip(np.digitize(e['max_nr'], bins) - 1, 0, len(bins) - 2)]/kg
            h = np.histogram(e['max_nr'][cut_fv(e,'nr') & cut_nr_er(e) & cut_nr(e)], weights=weights[cut_fv(e,'nr') & cut_nr_er(e) & cut_nr(e)],bins=bins)
            material_data[i+j]=h[0]

            #to material error
            bins = bins_nr
            bins_diff=np.diff(bins)
            material_activity=genNormalize[i+'_'+j]['factor']/(genNormalize[i+'_'+j]['activity']*days)
            
            c=1/(material_activity*kg)
            material_data_error[i+j]= [ (a * c / b)**(1/2)  for a, b in zip(material_data[i+j], bins_diff)]
            
            NR_number,ER_number=rate_error(e,nr_max,nr_min,er_max)
            NR_rate[i+j]=NR_number*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            NR_error[i+j]=(NR_number)**(1/2)*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            # ER_rate[i+j]=ER_number*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            # ER_error[i+j]=(ER_number)**(1/2)*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg

    isotopes = ['U238','U235','Th232','Cs137','K40','Co60']
    components = ['PMTbase']
    for i in isotopes:
        for j in components:
            files = [os.path.join(folder, f'events_{i}_{j}.h5')]
            with h5py.File(files[0], 'r', libver='latest', swmr=True) as ipt:
                events[i+j] = ipt['events'][:]
            e=events[i+j]

            weights = np.full(len(e), genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor']) / np.diff(bins)[np.clip(np.digitize(e['max_nr'], bins) - 1, 0, len(bins) - 2)]/kg
            h = np.histogram(e['max_nr'][cut_fv(e,'nr') & cut_nr_er(e) & cut_nr(e)], weights=weights[cut_fv(e,'nr') & cut_nr_er(e) & cut_nr(e)],bins=bins)
            material_data[i+j]=h[0]

            #to material error
            bins = bins_nr
            bins_diff=np.diff(bins)
            material_activity=genNormalize[i+'_'+j]['factor']/(genNormalize[i+'_'+j]['activity']*days)
            
            c=1/(material_activity*kg)
            material_data_error[i+j]= [ (a * c / b)**(1/2)  for a, b in zip(material_data[i+j], bins_diff)]
            
            NR_number,ER_number=rate_error(e,nr_max,nr_min,er_max)
            NR_rate[i+j]=NR_number*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            NR_error[i+j]=(NR_number)**(1/2)*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            # ER_rate[i+j]=ER_number*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            # ER_error[i+j]=(ER_number)**(1/2)*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg

    isotopes = ['U238','U235','Th232','Cs137','K40','Co60']
    components = ['PMTwindow']
    for i in isotopes:
        for j in components:
            files = [os.path.join(folder, f'events_{i}_{j}.h5')]
            with h5py.File(files[0], 'r', libver='latest', swmr=True) as ipt:
                events[i+j] = ipt['events'][:]
            e=events[i+j]

            weights = np.full(len(e), genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor']) / np.diff(bins)[np.clip(np.digitize(e['max_nr'], bins) - 1, 0, len(bins) - 2)]/kg
            h = np.histogram(e['max_nr'][cut_fv(e,'nr') & cut_nr_er(e) & cut_nr(e)], weights=weights[cut_fv(e,'nr') & cut_nr_er(e) & cut_nr(e)],bins=bins)
            material_data[i+j]=h[0]

            #to material error
            bins = bins_nr
            bins_diff=np.diff(bins)
            material_activity=genNormalize[i+'_'+j]['factor']/(genNormalize[i+'_'+j]['activity']*days)
            
            c=1/(material_activity*kg)
            material_data_error[i+j]= [ (a * c / b)**(1/2)  for a, b in zip(material_data[i+j], bins_diff)]
            
            NR_number,ER_number=rate_error(e,nr_max,nr_min,er_max)
            NR_rate[i+j]=NR_number*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            NR_error[i+j]=(NR_number)**(1/2)*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            # ER_rate[i+j]=ER_number*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            # ER_error[i+j]=(ER_number)**(1/2)*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
    
    isotopes=['U238','Th232','Co60','K40','Cs137','U235']
    components=['Flange']
    for i in isotopes:
        for j in components:
            files = [os.path.join(folder, f'events_{i}_{j}.h5')]
            with h5py.File(files[0], 'r', libver='latest', swmr=True) as ipt:
                events[i+j] = ipt['events'][:]
            e=events[i+j]

            weights = np.full(len(e), genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor']) / np.diff(bins)[np.clip(np.digitize(e['max_nr'], bins) - 1, 0, len(bins) - 2)]/kg
            h = np.histogram(e['max_nr'][cut_fv(e,'nr') & cut_nr_er(e) & cut_nr(e)], weights=weights[cut_fv(e,'nr') & cut_nr_er(e) & cut_nr(e)],bins=bins)
            material_data[i+j]=h[0]

            #to material error
            bins = bins_nr
            bins_diff=np.diff(bins)
            material_activity=genNormalize[i+'_'+j]['factor']/(genNormalize[i+'_'+j]['activity']*days)
            
            c=1/(material_activity*kg)
            material_data_error[i+j]= [ (a * c / b)**(1/2)  for a, b in zip(material_data[i+j], bins_diff)]
            
            NR_number,ER_number=rate_error(e,nr_max,nr_min,er_max)
            NR_rate[i+j]=NR_number*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            NR_error[i+j]=(NR_number)**(1/2)*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            # ER_rate[i+j]=ER_number*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg
            # ER_error[i+j]=(ER_number)**(1/2)*(genNormalize[i+'_'+j]['activity']*days / genNormalize[i+'_'+j]['factor'])/kg

    return events,material_data,material_data_error,NR_rate,NR_error,ER_rate,ER_error

print('finish the function,then read the data')
#归一化，需要先跑好归一化python文件

days=24*3600
kg= (fv_z1 - fv_z2)/10* fv_r/10 *  fv_r/10 * np.pi *2.862 /1000
#调用了活度，其中有0的需要筛选

import json
f = open('../normalization.json', 'r')
content = f.read()
genNormalize = json.loads(content)


folder = f'/public/home/yyy/RELICS5各个版本/RELICS5_v2/RelicsSim_v2_solid_origin_material_R_E/result/Material_10M_steel_window_1M_others'
bins_nr=np.logspace(-4, 7, bins_NR+1)
nr_max=1.36
nr_min=0.63
er_max=float(ER_MAX)
material_events,material_data,material_data_error,NR_rate,NR_error,ER_rate,ER_error=material(folder,genNormalize,bins_nr,nr_max,nr_min,er_max)

print('finish reading data,then output the message')

# 计算总的LXe ER率和误差
Total_NR_rate = sum(NR_rate.values())
Total_NR_error = np.sqrt(sum([err**2 for err in NR_error.values()]))

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
json_str = json.dumps(convert_to_list(material_data))
with open(f"./NR_{bins_NR}bins/temp/material_data.json", "w") as f:
    f.write(json_str) 
    
json_str = json.dumps(convert_to_list(NR_rate))
with open(f"./NR_{bins_NR}bins/temp/material_NR_rate.json", "w") as f:
    f.write(json_str)
    
json_str = json.dumps(convert_to_list(NR_error))
with open(f"./NR_{bins_NR}bins/temp/material_NR_error.json", "w") as f:
    f.write(json_str)

#进行数据的处理，得到汇总后的结果
#material data,figure
Material_Data=sum(material_data.values())
error=np.zeros(len(Material_Data))

bins = bins_nr
bins_diff=np.diff(bins)
kg= (fv_z1 - fv_z2)/10* fv_r/10 *  fv_r/10 * np.pi *2.862 /1000 # 单位kg

for i in range(len(Material_Data)):
    isotopes=['U238','Th232','Co60','K40','Cs137']
    # components=['steel', 'copper']
    components=['steel']
    for j in isotopes:
        for k in components:
            c=1/(kg*genNormalize[j+'_'+k]['factor']/(genNormalize[j+'_'+k]['activity']*days))
            error[i]+=(material_data[j+k][i] * c / bins_diff[i])
            #print(error)

    isotopes=['U238','Th232','Co60','K40','Cs137','U235']
    components=['Flange']
    for j in isotopes:
        for k in components:
            c=1/(kg*genNormalize[j+'_'+k]['factor']/(genNormalize[j+'_'+k]['activity']*days))
            error[i]+=(material_data[j+k][i] * c / bins_diff[i])
            #print(error)
            
    isotopes=['U238','Th232','Co60','K40','Cs137']
    components=['Teflon']
    for j in isotopes:
        for k in components:
            c=1/(kg*genNormalize[j+'_'+k]['factor']/(genNormalize[j+'_'+k]['activity']*days))
            error[i]+=(material_data[j+k][i] * c / bins_diff[i])
     
    isotopes=['U238','Th232','Co60','K40','Cs137','U235']
    components=['PMTwindow']
    for j in isotopes:
        for k in components:
            c=1/(kg*genNormalize[j+'_'+k]['factor']/(genNormalize[j+'_'+k]['activity']*days))
            error[i]+=(material_data[j+k][i] * c / bins_diff[i])
            
    isotopes=['U238','Th232','Co60','K40','Cs137','U235']
    components=['PMTbase']
    for j in isotopes:
        for k in components:
            c=1/(kg*genNormalize[j+'_'+k]['factor']/(genNormalize[j+'_'+k]['activity']*days))
            error[i]+=(material_data[j+k][i] * c / bins_diff[i])

        
error = np.sqrt(error)
Material_Data=sum(material_data.values())

#ER data，存储ER本底数据，xenon存在background的material文件夹下，pandax存在material_pandax文件夹下
json_str = json.dumps(convert_to_list(Material_Data))
with open(f"./NR_{bins_NR}bins/material_nr.json", "w") as f:
    f.write(json_str)

json_str = json.dumps(convert_to_list(error))
with open(f"./NR_{bins_NR}bins/material_nr_error.json", "w") as f:
    f.write(json_str)
    
json_str = json.dumps(convert_to_list(bins_nr))
with open(f"./NR_{bins_NR}bins/material_nr_bin.json", "w") as f:
    f.write(json_str)

#material_er and material_er_error
Material_NR=sum(NR_rate.values())
Material_Error=0
                
error_material_nr=0
for i in NR_error.keys():
    error_material_nr+=NR_error[i]**2
error_material_nr=error_material_nr**(1/2)

Material_Error=error_material_nr

print('CRN NR rate:',Material_NR,'/kg/days')
print('CRN NR error:',Material_Error)
print('bins_nr=np.logspace(-4, 7, 101)')
# # 输出总的LXe率
# print('\n--- Total LXe Results ---')
# print('Total LXe ER rate:', Total_ER_rate/(ER_MAX-ER_MIN), '/kg/days/keV')
# print('Total LXe ER error:', Total_ER_error/(ER_MAX-ER_MIN))
# print('Total LXe ER rate (mDRU):', Total_ER_rate/(ER_MAX-ER_MIN) * 1000)
# print('Total LXe ER error (mDRU):', Total_ER_error/(ER_MAX-ER_MIN) * 1000)
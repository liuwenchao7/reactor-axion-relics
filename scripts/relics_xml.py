#!/usr/bin/python3
import os
import argparse
import json
import copy
import shutil
from xml.dom.minidom import Document

import math
import numpy as np

from Reactor_gamma_relics.R10_v2_clean.scripts.arrangement import PMTArrangement, TankArrangement

parser = argparse.ArgumentParser(description='Automatically relics .xml geometry generator')

parser.add_argument('ipt',
                    action='store',
                    help='Input .json file, containing geometry parameters')
parser.add_argument('-o', dest='opt',
                    action='store', required=True,
                    help='Output .xml file name')
parser.add_argument('--run_number', dest='run_number',
                    action='store', default='0',
                    help='Run number in .xml file')
parser.add_argument('--gen', dest='generator',
                    choices=['muon', 'neutron', 'gamma', 'material', 'CRN'],
                    action='store', default='muon',
                    help='Generator used')
parser.add_argument('--cutlength', dest='cutlength', type=float,
                    action='store', default=0.1,  # in mm
                    help='Cut length used in simulation')
parser.add_argument('--enable_track', action='store_true',
                    help='Enable track of surface flux')
parser.add_argument('--force_sd', action='store_true',
                    help='Force all solid sensitive detector')
parser.add_argument('--check_overlap', action='store_true',
                    help='Whether check overlap of volumes')
parser.add_argument('--get_geo_mass', action='store_true',
                    help='Whether print geometry mass')
parser.add_argument('--save_txt', action='store_true',
                    help='Whether save .txt files of position')
parser.add_argument('--optical', action='store_true',
                    help='Whether do optical simulation')

args = parser.parse_args()

number = args.run_number
fipt = args.ipt
fopt = args.opt
gen = args.generator
cutlength = args.cutlength
enable_track = args.enable_track
force_sd = args.force_sd
check_overlap = args.check_overlap
get_geo_mass = args.get_geo_mass
save_txt = args.save_txt
optical = args.optical

with open(fipt) as f:
    params = json.load(f)

# For LXe position is 0cm
position_lxe = 0

# ========== 内部探测器参数（强制从JSON读取） ==========
dis_gate = params['dis_gate']
dis_cathode = params['dis_cathode']
dis_botscreen = params['dis_botscreen']
dis_botteflon = params['dis_botteflon']
thick_botteflon = params['thick_botteflon']
pmt_hole = params['pmt_hole']
radius_inn_lxe = params['radius_inn_lxe']
teflon_n = params['teflon_n']
dis_anode = params['dis_anode']
wire_diameter = params['wire_diameter']
wire_pitch = params['wire_pitch']
dis_topteflon = params['dis_topteflon']
thick_topteflon = params['thick_topteflon']
radius_outer_teflon = params['radius_outer_teflon']
PEEK_hole = params['PEEK_hole']
thick_botPEEK = params['thick_botPEEK']
radius_botPEEK = params['radius_botPEEK']
dis_botPEEK = params['dis_botPEEK']
thick_topPEEK = params['thick_topPEEK']
radius_topPEEK = params['radius_topPEEK']
dis_topPEEK = params['dis_topPEEK']
width_keel = params['width_keel']
length_keel = params['length_keel']
radius_inner_shapingrings = params['radius_inner_shapingrings']
radius_outer_shapingrings = params['radius_outer_shapingrings']
thick_shapingrings = params['thick_shapingrings']

# ========== 铅/聚乙烯屏蔽层参数（全部强制从JSON读取） ==========
lead_width_x     = params['lead_width_x']
lead_width_y     = params['lead_width_y']
lead_width_z     = params['lead_width_z']

lead_hole_width_x = params['lead_hole_width_x']
lead_hole_width_y = params['lead_hole_width_y']
lead_hole_width_z = params['lead_hole_width_z']

inner_pe_width_x = params['inner_pe_width_x']
inner_pe_width_y = params['inner_pe_width_y']
inner_pe_width_z = params['inner_pe_width_z']

pe_hole_width_x  = params['pe_hole_width_x']
pe_hole_width_y  = params['pe_hole_width_y']
pe_hole_width_z  = params['pe_hole_width_z']

air_width_x      = params['air_width_x']
air_width_y      = params['air_width_y']
air_width_z      = params['air_width_z']

# 偏移参数（强制从JSON读取）
lead_shift_x     = params['lead_shift_x']
lead_shift_y     = params['lead_shift_y']
lead_shift_z     = params['lead_shift_z']

lead_hole_shift_x = params['lead_hole_shift_x']
lead_hole_shift_y = params['lead_hole_shift_y']
lead_hole_shift_z = params['lead_hole_shift_z']

inner_pe_shift_x = params['inner_pe_shift_x']
inner_pe_shift_y = params['inner_pe_shift_y']
inner_pe_shift_z = params['inner_pe_shift_z']

pe_hole_shift_x  = params['pe_hole_shift_x']
pe_hole_shift_y  = params['pe_hole_shift_y']
pe_hole_shift_z  = params['pe_hole_shift_z']

air_shift_x      = params['air_shift_x']
air_shift_y      = params['air_shift_y']
air_shift_z      = params['air_shift_z']

# ========== 混凝土参数（全部强制从JSON读取） ==========
concrete_high = params['concrete_high']
concrete_width_x = params['concrete_width_x']
concrete_width_y = params['concrete_width_y']
concrete_shift_x = params['concrete_shift_x']
concrete_shift_y = params['concrete_shift_y']
concrete_shift_z = params['concrete_shift_z']

# ========== 缪子/CRN 参数 ==========
muon_x = params['muon_x']
muon_z = params['muon_z']
CRN_x = params['CRN_x']
CRN_z = params['CRN_z']

# ========== 派生参数计算 ==========
radius_outer_botteflon = radius_inn_lxe
height_inn_lxe = dis_gate + dis_cathode + dis_botteflon + 1/2 * thick_botteflon
radius_outer_topteflon = radius_inn_lxe
radius_inner_teflon = radius_inn_lxe
height_teflon_gaspart = dis_topteflon + dis_anode - 1/2 * thick_topteflon
height_teflon_lxepart = dis_gate + dis_cathode + dis_botteflon - 1/2 * thick_botteflon

# ========== 光学/物理模式设置 ==========
if optical:
    physics_name = 'RelicsOpticalPhysics'
    generator_name = 'ConfineGenerator'
    analysis_name = 'RelicsOpticalAnalysis'

    LXesd = False
    PMTsd = True
    opticalsd = True
    seed = int(number)

    materialList = {'LXe_SCINTILLATIONYIELD': '3000000.'}

    analysisList = {'EnableEnergyDeposition': 0,
                    'EnablePrimaryParticle': 1,
                    'SaveNullEvents': 1,
                    'OpticalSdName': 'bR8520',
                    'UserSeed': f'{seed}',
                    'GetGeoMass': int(get_geo_mass)}

    pmt_opsurface = {
        'base_type': 'dielectric_metal',
        'casing_finish': 'polished',
        'casing_type': 'dielectric_metal',
        'window_partner': 'parent',
        'window_finish': 'polished',
        'window_order': str(1),
        'photocathode_finish': 'polished'
    }
    surface = True
else:
    physics_name = 'PandaXPhysics'
    if gen == 'muon':
        generator_name = 'MuonGenerator'
    elif gen == 'neutron' or gen == 'gamma':
        generator_name = 'SimpleGPSGenerator'
    elif gen == 'material':
        generator_name = 'ConfineGenerator'
    elif gen == 'CRN':
        generator_name = 'NeutronGenerator'
    analysis_name = 'PandaXAnalysis'

    LXesd = True
    Plasticsd = False
    PMTsd = False
    opticalsd = False
    seed = 0

    materialList = dict()

    analysisList = {'EnableEnergyDeposition': 1,
                    'EnableSurfaceFlux': int(enable_track),
                    'EnablePrimaryParticle': 1,
                    'SaveNullEvents': 0,
                    'EnableDecayChainSplitting': 1,
                    'ChainSplittingLifeTime': '400*us',
                    'UserSeed': f'{seed}',
                    'GetGeoMass': int(get_geo_mass)}

    pmt_opsurface = dict()
    surface = False

# ========== XML 文档初始化 ==========
doc = Document()
djson = dict()

root = doc.createElement('BambooMC')
doc.appendChild(root)

run = doc.createElement('run')
run.setAttribute('number', number)
root.appendChild(run)
djson['run'] = int(number)

geometry = doc.createElement('geometry')
root.appendChild(geometry)
djson['geometry'] = dict()

material = doc.createElement('material')
geometry.appendChild(material)
djson['geometry']['material'] = dict()
material.setAttribute('name', 'Material')
djson['geometry']['material'].update({'name': 'Material'})

if len(list(materialList.keys())) > 0:
    djson['geometry']['material'].update({'parameters': dict()})
for k in materialList.keys():
    parameter = doc.createElement('parameter')
    material.appendChild(parameter)
    parameter.setAttribute('name', k)
    parameter.setAttribute('value', f'{materialList[k]}')
    djson['geometry']['material']['parameters'].update({k: materialList[k]})

detector = doc.createElement('detector')
geometry.appendChild(detector)
djson['geometry']['detectors'] = []

detector.setAttribute('type', 'World')
detector.setAttribute('name', 'World')
djson['geometry']['detectors'].append(dict())
djson['geometry']['detectors'][-1].update({'type': 'World'})
djson['geometry']['detectors'][-1].update({'name': 'World'})

Cuboid = 'Cuboid'
Cylinder = 'Cylinder'
Prism = 'Prism'
Polyhedra = 'Polyhedra'
HolesBoard = 'HolesBoard'
Array = 'Array'
PMTs = 'PMTsR8520'
LXe = 'LXe'
GXe = 'GXe'
Teflon = 'Teflon'
PEEK = 'PEEK'
PMTbasePEEK = 'PMTbasePEEK'
XenonVeto = 'XenonVeto'
XenonDetector = 'XenonDetector'
SS304LSteel = 'SS304LSteel'
Flange = 'Flange'
concrete = 'concrete'
topPMTs = 'data/topPMTs.txt'
botPMTs = 'data/botPMTs.txt'
XeInTopPMTTeflon = 'data/XeInTopPMTTeflon.txt'
XeInBotPMTTeflon = 'data/XeInBotPMTTeflon.txt'

# CAD 文件路径
CAD_OuterContainerRZ = 'CADdata/OuterContainerRZ.txt'
CAD_VacuumRZ = 'CADdata/VacuumRZ.txt'
CAD_InnerContainerRZ = 'CADdata/InnerContainerRZ.txt'
CAD_OuterGasXeRZ = 'CADdata/OuterGasXeRZ.txt'
CAD_XenonVetoRZ = 'CADdata/XenonVetoRZ.txt'
CAD_Lxe_OutContainer = 'CADdata/Lxe_OutContainer.txt'
CAD_Lxe_InnerContainer = 'CADdata/Lxe_InnerContainer.txt'
CAD_Top_PMT_position = 'CADdata/topPMTs.txt'
CAD_Bot_PMT_position = 'CADdata/botPMTs.txt'
CAD_Flange63 = 'CADdata/Flange63.txt'
CAD_Flange35 = 'CADdata/Flange35.txt'
CAD_shapingrings = 'CADdata/shapingrings.txt'

# 输出数据文件路径
shapingrings = 'data/shapingrings.txt'
OuterContainerRZ = 'data/OuterContainerRZ.txt'
VacuumRZ = 'data/VacuumRZ.txt'
InnerContainerRZ = 'data/InnerContainerRZ.txt'
OuterGasXeRZ = 'data/OuterGasXeRZ.txt'
XenonVetoRZ = 'data/XenonVetoRZ.txt'
Lxe_OutContainer = 'data/Lxe_OutContainer.txt'
Lxe_InnerContainer = 'data/Lxe_InnerContainer.txt'
Flange63 = 'data/Flange63.txt'
Flange35 = 'data/Flange35.txt'

shutil.copy(CAD_OuterContainerRZ, OuterContainerRZ)
shutil.copy(CAD_VacuumRZ, VacuumRZ)
shutil.copy(CAD_InnerContainerRZ, InnerContainerRZ)
shutil.copy(CAD_OuterGasXeRZ, OuterGasXeRZ)
shutil.copy(CAD_XenonVetoRZ, XenonVetoRZ)
shutil.copy(CAD_Lxe_OutContainer, Lxe_OutContainer)
shutil.copy(CAD_Lxe_InnerContainer, Lxe_InnerContainer)
shutil.copy(CAD_Flange63, Flange63)
shutil.copy(CAD_Flange35, Flange35)

# 法兰参数（固定）
flange63_height = 1.75 * 2
flange63_Outradius = 5.7
flange63_Innerradius = 3.35

flange35_height = 1.27 * 2
flange35_Outradius = 3.5
flange35_Innerradius = 2.075

# ========== PMT 位置处理函数 ==========
def replace_column_with_value(
    input_file: str,
    output_file: str,
    col_index: int,
    value: float,
    fmt: str = "{:.5f}",
):
    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        for line in fin:
            line_strip = line.strip()
            if not line_strip:
                fout.write(line)
                continue
            parts = line_strip.split()
            if len(parts) == 1:
                fout.write(line_strip + "\n")
                continue
            if len(parts) > col_index:
                parts[col_index] = fmt.format(value)
            fout.write(" ".join(parts) + "\n")

Top_PMT_positionz = position_lxe + dis_anode + dis_topteflon + 3.3/2 + thick_topteflon /2
Bot_PMT_positionz = position_lxe - (dis_gate + dis_cathode + dis_botteflon + 3.3/2 + thick_botteflon /2)

if save_txt:
    replace_column_with_value(CAD_Top_PMT_position, topPMTs, 3, Top_PMT_positionz)
    replace_column_with_value(CAD_Bot_PMT_position, botPMTs, 3, Bot_PMT_positionz)
    replace_column_with_value(CAD_Top_PMT_position, XeInTopPMTTeflon, 3, 0.00)
    replace_column_with_value(CAD_Bot_PMT_position, XeInBotPMTTeflon, 3, 0.00)
    print('Change the PMT position z and save txt files done!')

def subtract_from_column(
    input_file: str,
    output_file: str,
    col_index: int,
    position_lxe: float,
    height_inn_lxe: float,
    fmt: str = "{:.4f}",
):
    delta = position_lxe - 1/2 * height_teflon_lxepart
    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        for line in fin:
            line_strip = line.strip()
            if not line_strip:
                fout.write(line)
                continue
            parts = line_strip.split()
            if len(parts) == 1:
                fout.write(line_strip + "\n")
                continue
            if len(parts) > col_index:
                z = float(parts[col_index])
                z_new = z - delta
                parts[col_index] = fmt.format(z_new)
            fout.write(" ".join(parts) + "\n")
    print(f"处理完成，已生成文件: {output_file}")

if save_txt:
    subtract_from_column(CAD_shapingrings, shapingrings, 3, position_lxe, height_inn_lxe)

# ========== 世界尺寸 ==========
world_half_x = 1000
djson['geometry']['detectors'][-1].update({'parameters': dict()})
for x in ['x', 'y', 'z']:
    parameter = doc.createElement('parameter')
    parameter.setAttribute('name', f'half_{x}')
    parameter.setAttribute('value', f'{world_half_x}*cm')
    detector.appendChild(parameter)
    djson['geometry']['detectors'][-1]['parameters'].update({f'half_{x}': f'{world_half_x}*cm'})

def dict_merger(dict_a, dict_b):
    dict_c = copy.deepcopy(dict_a)
    dict_c.update(dict_b)
    return dict_c

# ========== 内部探测器位置计算 ==========
xenon_detector_shift_z = position_lxe - 1/2 * height_inn_lxe
radius_outer_lxe = np.sqrt((radius_outer_shapingrings + width_keel)**2 + length_keel**2) + 0.0001
siedteflon_lxe_offset = position_lxe - 1/2 * height_teflon_lxepart - xenon_detector_shift_z
bot_teflon_offset = position_lxe - (dis_gate + dis_cathode + dis_botteflon) - xenon_detector_shift_z
top_teflon_offset = position_lxe + dis_anode + dis_topteflon
siedteflon_gas_offset = position_lxe + 1/2 * height_teflon_gaspart

top_PEEK_offset = dis_topPEEK + position_lxe
bot_PEEK_offset = position_lxe - dis_botPEEK

keel_offsetxy = radius_outer_shapingrings + 1/2 * width_keel
keel_height_GXe = height_teflon_gaspart
keel_offsetz_GXe = position_lxe + 1/2 * keel_height_GXe
keel_height_LXe = height_teflon_lxepart
keel_offsetz_LXe = 0

def write_keel_positions_A_txt(output_file: str, keel_offsetxy: float, keel_offsetz: float):
    rows = [
        (0,  keel_offsetxy,  0.0,           keel_offsetz, "0.00",   "0.00"),
        (1,  0.0,            keel_offsetxy, keel_offsetz, "0.00", "90"),
        (2, -keel_offsetxy,  0.0,           keel_offsetz, "0.00",   "0.00"),
        (3,  0.0,           -keel_offsetxy, keel_offsetz, "0.00", "90"),
    ]
    with open(output_file, "w") as f:
        f.write("4\n")
        for idx, x, y, z, r1, r2 in rows:
            f.write(f"{idx} {x:.6f} {y:.6f} {z:.6f} {r1} {r2}\n")
    print("已生成:", output_file)

if save_txt:
    write_keel_positions_A_txt('./data/keel_GXe.txt', keel_offsetxy, keel_offsetz_GXe)
    write_keel_positions_A_txt('./data/keel_LXe.txt', keel_offsetxy, keel_offsetz_LXe)

# ========== 探测器列表（所有参数均来自JSON） ==========
detectorList = [
    # # 混凝土底板（尺寸和偏移全部由JSON控制）
    # {'type': Cuboid,
    #  'name': 'Concrete',
    #  'parent': 'World',
    #  'parameters': {
    #      'width_x': concrete_width_x,
    #      'width_y': concrete_width_y,
    #      'width_z': concrete_high,
    #      'material': concrete,
    #      'shift_x': concrete_shift_x,
    #      'shift_y': concrete_shift_y,
    #      'shift_z': concrete_shift_z}},

    # # 外层包裹
    # {'type': Cuboid,
    #  'name': 'AirOut',
    #  'parent': 'World',
    #  'parameters': {
    #      'width_x': 140,
    #      'width_y': 140,
    #      'width_z': 200,
    #      'material': 'G4_AIR',
    #      'shift_x': 0,
    #      'shift_y': 0,
    #      'shift_z': 0,
    #     #  'entrack' : '1',
    #     #  'sendname' : 'AirOut',
    #      }},

    # # 铅屏蔽层
    # {'type': Cuboid,
    #  'name': 'Lead',
    #  'parent': 'AirOut',
    #  'parameters': {
    #      'width_x': lead_width_x,
    #      'width_y': lead_width_y,
    #      'width_z': lead_width_z,
    #      'material': 'G4_Pb',
    #      'shift_x': lead_shift_x,
    #      'shift_y': lead_shift_y,
    #      'shift_z': lead_shift_z}},

    # # 铅孔洞
    # {'type': Cuboid,
    #  'name': 'LeadHole',
    #  'parent': 'Lead',
    #  'parameters': {
    #      'width_x': lead_hole_width_x,
    #      'width_y': lead_hole_width_y,
    #      'width_z': lead_hole_width_z,
    #      'material': 'G4_AIR',
    #      'shift_x': lead_hole_shift_x,
    #      'shift_y': lead_hole_shift_y,
    #      'shift_z': lead_hole_shift_z}},

    # # 聚乙烯内层
    # {'type': Cuboid,
    #  'name': 'InnerPE',
    #  'parent': 'Lead',
    #  'parameters': {
    #      'width_x': inner_pe_width_x,
    #      'width_y': inner_pe_width_y,
    #      'width_z': inner_pe_width_z,
    #      'material': 'G4_POLYETHYLENE',
    #      'shift_x': inner_pe_shift_x,
    #      'shift_y': inner_pe_shift_y,
    #      'shift_z': inner_pe_shift_z}},

    # # 聚乙烯孔洞
    # {'type': Cuboid,
    #  'name': 'PEHole',
    #  'parent': 'InnerPE',
    #  'parameters': {
    #      'width_x': pe_hole_width_x,
    #      'width_y': pe_hole_width_y,
    #      'width_z': pe_hole_width_z,
    #      'material': 'G4_AIR',
    #      'shift_x': pe_hole_shift_x,
    #      'shift_y': pe_hole_shift_y,
    #      'shift_z': pe_hole_shift_z}},

    # # 内部空腔 Air
    # {'type': Cuboid,
    #  'name': 'Air',
    #  'parent': 'InnerPE',
    #  'parameters': {
    #      'width_x': air_width_x,
    #      'width_y': air_width_y,
    #      'width_z': air_width_z,
    #      'material': 'G4_AIR',
    #      'shift_x': air_shift_x,
    #      'shift_y': air_shift_y,
    #      'shift_z': air_shift_z}},
   # # 混凝土底板 — 已禁用（无屏蔽模式）
    # {'type': Cuboid,
    #  'name': 'Concrete',
    #  'parent': 'World',
    #  'parameters': {
    #      'width_x': 500,
    #      'width_y': 500,
    #      'width_z': 5,
    #      'material': concrete,
    #      'shift_x': 0,
    #      'shift_y': 0,
    #      'shift_z': -137.5}},
    # # PE 外层 — 已禁用（无屏蔽模式）
    # {'type': Cuboid,
    #  'name': 'PE',
    #  'parent': 'World',
    #  'parameters': {
    #      'width_x': 230,
    #      'width_y': 200,
    #      'width_z': 270,
    #      'material': 'G4_POLYETHYLENE',
    #      'shift_x': 0,
    #      'shift_y': 0,
    #      'shift_z': 0}},
    # # Lead 铅层 — 已禁用（无屏蔽模式）
    # {'type': Cuboid,
    #  'name': 'Lead',
    #  'parent': 'PE',
    #  'parameters': {
    #      'width_x': 100,
    #      'width_y': 100,
    #      'width_z': 160,
    #      'material': 'G4_Pb',
    #      'shift_x': -15,
    #      'shift_y': 0,
    #      'shift_z': -45}},
    # # InnerPE 内层聚乙烯 — 已禁用（无屏蔽模式）
    # {'type': Cuboid,
    #  'name': 'InnerPE',
    #  'parent': 'Lead',
    #  'parameters': {
    #      'width_x': 90,
    #      'width_y': 90,
    #      'width_z': 150,
    #      'material': 'G4_POLYETHYLENE',
    #      'shift_x': 0,
    #      'shift_y': 0,
    #      'shift_z': 0}},
    # {'type': Cuboid,
    #  'name': 'Lead_FRONT',
    #  'parent': 'PE',
    #  'parameters': {
    #      'width_x': 5,
    #      'width_y': 200,
    #      'width_z': 265,
    #      'material': 'G4_Pb',
    #      'shift_x': 112.5,
    #      'shift_y': 0,
    #      'shift_z': 2.5,
    #     #  'entrack' : '1',
    #     #  'sendname' : 'AirOut',
    #      }},
    # {'type': Cuboid,
    #  'name': 'Lead_BACK',
    #  'parent': 'PE',
    #  'parameters': {
    #      'width_x': 5,
    #      'width_y': 200,
    #      'width_z': 265,
    #      'material': 'G4_Pb',
    #      'shift_x': -112.5,
    #      'shift_y': 0,
    #      'shift_z': 2.5,
    #     #  'entrack' : '1',
    #     #  'sendname' : 'AirOut',
    #      }},
    # {'type': Cuboid,
    #  'name': 'Lead_LEFT',
    #  'parent': 'PE',
    #  'parameters': {
    #      'width_x': 210,
    #      'width_y': 5,
    #      'width_z': 265,
    #      'material': 'G4_Pb',
    #      'shift_x': 0,
    #      'shift_y': 97.5,
    #      'shift_z': 2.5,
    #     #  'entrack' : '1',
    #     #  'sendname' : 'AirOut',
    #      }},
    # {'type': Cuboid,
    #  'name': 'Lead_RIGHT',
    #  'parent': 'PE',
    #  'parameters': {
    #      'width_x': 210,
    #      'width_y': 5,
    #      'width_z': 265,
    #      'material': 'G4_Pb',
    #      'shift_x': 0,
    #      'shift_y': -97.5,
    #      'shift_z': 2.5,
    #     #  'entrack' : '1',
    #     #  'sendname' : 'AirOut',
    #      }},
    # # PS 塑料闪烁体 — 已禁用（无屏蔽模式）
    # {'type': Cuboid,
    #  'name': 'PS',
    #  'parent': 'InnerPE',
    #  'parameters': {
    #      'width_x': 76,
    #      'width_y': 76,
    #      'width_z': 136,
    #      'material': 'PS',
    #      'shift_x': 0,
    #      'shift_y': 0,
    #      'shift_z': 0}},
    # Air 内腔 — 已禁用（v2 clean: 仅保留 LXe TPC，无前置材料）
    # {'type': Cuboid,
    #  'name': 'Air',
    #  'parent': 'World',
    #  'parameters': {
    #      'width_x': 60,
    #      'width_y': 60,
    #      'width_z': 120,
    #      'material': 'G4_AIR',
    #      'shift_x': 0,
    #      'shift_y': 0,
    #      'shift_z': 0,
    #      'entrack' : '1',
    #      'sendname' : 'Air',
    #      }},

    # ===== 容器/真空/法兰/气氙/外veto — 已禁用（v2 clean: gamma 直达 LXe TPC） =====
    # {'type': Polyhedra,
    #  'name': 'OuterContainer',
    #  'parent': 'Air',
    #  'parameters': {'source': OuterContainerRZ,
    #                 'material': SS304LSteel}},

    # {'type': Polyhedra,
    #  'name': 'Vacuum',
    #  'parent': 'OuterContainer',
    #  'parameters': {'source': VacuumRZ,
    #                 'material': 'Vacuum'}},

    # {'type': Array,
    #  'name': 'Flange63',
    #  'parent': 'Air',
    #  'parameters': {'objects': 'Prism',
    #                 'outer_radius': flange63_Outradius,
    #                 'inner_radius': flange63_Innerradius,
    #                 'sides': str(0),
    #                 'height': flange63_height,
    #                 'material': Flange,
    #                 'source': Flange63}},

    # {'type': Array,
    #  'name': 'Flange35',
    #  'parent': 'Air',
    #  'parameters': {'objects': 'Prism',
    #                 'outer_radius': flange35_Outradius,
    #                 'inner_radius': flange35_Innerradius,
    #                 'sides': str(0),
    #                 'height': flange35_height,
    #                 'material': Flange,
    #                 'source': Flange35}},

    # {'type': Polyhedra,
    #  'name': 'InnerContainer',
    #  'parent': 'Vacuum',
    #  'parameters': {'source': InnerContainerRZ,
    #                 'material': SS304LSteel}},

    # {'type': Polyhedra,
    #  'name': 'OuterGasXe',
    #  'parent': 'InnerContainer',
    #  'parameters': {'source': OuterGasXeRZ,
    #                 'material': GXe}},

    # {'type': Polyhedra,
    #  'name': XenonVeto,
    #  'parent': 'InnerContainer',
    #  'parameters': {'source': XenonVetoRZ,
    #                 'material': LXe},
    #  'sd': LXesd},

    # {'type': Polyhedra,
    #  'name': 'Lxe_OutContainer',
    #  'parent': XenonVeto,
    #  'parameters': {'source': Lxe_OutContainer,
    #                 'material': SS304LSteel}},

    # {'type': Polyhedra,
    #  'name': 'Lxe_InnerContainer',
    #  'parent': 'Lxe_OutContainer',
    #  'parameters': {'source': Lxe_InnerContainer,
    #                 'material': LXe}},

    # ===== 中央 LXe TPC（保留，直接挂 World） =====
    {'type': Prism,
     'name': XenonDetector,
     'parent': 'World',
     'parameters': {'outer_radius': radius_outer_lxe,
                    'sides': str(teflon_n),
                    'height': height_inn_lxe,
                    'material': LXe,
                    'shift_z': xenon_detector_shift_z},
     'sd': LXesd},

    {'type': Prism,
     'name': 'XenonVeto1',
     'parent': XenonDetector,
     'parameters': {'outer_radius': radius_outer_lxe,
                    'inner_radius': radius_inn_lxe,
                    'sides': str(teflon_n),
                    'height': height_teflon_lxepart,
                    'material': LXe,
                    'shift_z': siedteflon_lxe_offset},
     'sd': LXesd},

    {'type': HolesBoard,
     'name': 'BotPMTTeflon',
     'parent': XenonDetector,
     'parameters': {'radius': radius_outer_lxe,
                    'sides': str(teflon_n),
                    'height': thick_botteflon,
                    'width_x': pmt_hole,
                    'material': Teflon,
                    'shift_z': bot_teflon_offset,
                    'source': XeInBotPMTTeflon,
                    'surface': surface}},

    {'type': Prism,
     'name': 'XeSideTeflon1',
     'parent': 'XenonVeto1',
     'parameters': {'outer_radius': radius_outer_teflon,
                    'inner_radius': radius_inner_teflon,
                    'sides': str(teflon_n),
                    'height': height_teflon_lxepart,
                    'material': Teflon,
                    'surface': surface}},

    {'type': Array,
     'name': 'Shapingrings',
     'parent': 'XenonVeto1',
     'parameters': {'objects': 'Prism',
                    'outer_radius': radius_outer_shapingrings,
                    'inner_radius': radius_inner_shapingrings,
                    'sides': str(teflon_n),
                    'height': thick_shapingrings,
                    'material': 'Copper',
                    'surface': surface,
                    'source': shapingrings}},

    {'type': Array,
     'name': 'TopKee2',
     'parent': 'XenonVeto1',
     'parameters': {'objects': Cuboid,
                    'width_x': width_keel,
                    'width_y': length_keel,
                    'width_z': keel_height_LXe,
                    'material': Teflon,
                    'shift_z': keel_offsetz_LXe,
                    'source': 'data/keel_LXe.txt',
                    'surface': surface}},

    # 以下组件依赖已移除的 XenonVeto/OuterGasXe — v2 clean 中移除
    # {'type': Array,
    #  'name': 'bR8520',
    #  'parent': XenonVeto,
    #  'parameters': dict_merger({'objects': PMTs,
    #                             'physvol': 'R8520',
    #                             'source': botPMTs},
    #                             pmt_opsurface),
    #  'sd': PMTsd,
    #  'opticalsd': opticalsd},

    # {'type': HolesBoard,
    #  'name': 'BotBasePEEK',
    #  'parent': XenonVeto,
    #  'parameters': {'radius': radius_botPEEK,
    #                 'sides': str(teflon_n),
    #                 'height': thick_botPEEK,
    #                 'width_x': PEEK_hole,
    #                 'material': PMTbasePEEK,
    #                 'shift_z': bot_PEEK_offset,
    #                 'source': XeInBotPMTTeflon}},

    # {'type': HolesBoard,
    #  'name': 'TopPMTTeflon',
    #  'parent': 'OuterGasXe',
    #  'parameters': {'radius': radius_outer_lxe,
    #                 'sides': str(teflon_n),
    #                 'height': thick_topteflon,
    #                 'width_x': pmt_hole,
    #                 'material': Teflon,
    #                 'shift_z': top_teflon_offset,
    #                 'source': XeInTopPMTTeflon,
    #                 'surface': surface}},

    # {'type': HolesBoard,
    #  'name': 'TopBasePEEK',
    #  'parent': 'OuterGasXe',
    #  'parameters': {'radius': radius_topPEEK,
    #                 'sides': str(teflon_n),
    #                 'height': thick_topPEEK,
    #                 'width_x': PEEK_hole,
    #                 'material': PMTbasePEEK,
    #                 'shift_z': top_PEEK_offset,
    #                 'source': XeInTopPMTTeflon}},

    # {'type': Array,
    #  'name': 'TopKeel',
    #  'parent': 'OuterGasXe',
    #  'parameters': {'objects': Cuboid,
    #                 'width_x': width_keel,
    #                 'width_y': length_keel,
    #                 'width_z': keel_height_GXe,
    #                 'material': Teflon,
    #                 'shift_z': keel_offsetz_GXe,
    #                 'source': 'data/keel_GXe.txt',
    #                 'surface': surface}},

    # {'type': Array,
    #  'name': 'tR8520',
    #  'parent': 'OuterGasXe',
    #  'parameters': dict_merger({'objects': PMTs,
    #                             'physvol': 'R8520',
    #                             'source': topPMTs},
    #                             pmt_opsurface),
    #  'sd': PMTsd,
    #  'opticalsd': opticalsd},

    # {'type': Prism,
    #  'name': 'XeSideTeflon0',
    #  'parent': 'OuterGasXe',
    #  'parameters': {'outer_radius': radius_outer_teflon,
    #                 'inner_radius': radius_inner_teflon,
    #                 'sides': str(teflon_n),
    #                 'height': height_teflon_gaspart,
    #                 'material': Teflon,
    #                 'shift_z': siedteflon_gas_offset,
    #                 'surface': surface}},
]

# ========== 电极网格 ==========
electrodes = [
    {'name': 'Cathode',
     'parent': XenonDetector,
     'parameters': {'shift_z': position_lxe - dis_gate - dis_cathode - xenon_detector_shift_z}},
    {'name': 'BottomScreening',
     'parent': XenonDetector,
     'parameters': {'shift_z': position_lxe - dis_gate - dis_cathode - dis_botscreen - xenon_detector_shift_z}},
    {'name': 'Gate',
     'parent': XenonDetector,
     'parameters': {'shift_z': position_lxe - dis_gate - xenon_detector_shift_z}},
    # Anode 位于气氙区（OuterGasXe），v2 clean 中移除
    # {'name': 'Anode',
    #  'parent': 'OuterGasXe',
    #  'parameters': {'shift_z': position_lxe + dis_anode}}
]

for electrode in electrodes:
    electrode.update({'type': 'MeshGrid'})
    electrode['parameters'].update({'radius': radius_inn_lxe,
                                    'sides': str(teflon_n),
                                    'material': SS304LSteel,
                                    'type': 'dielectric_metal',
                                    'surface': surface})

detectorList += electrodes

# ========== 电极环 ==========
dis_electric_rings = 0.25
thick_electric_rings = 0.5

electrodes_rings = [
    {'name': 'Cathode_rings',
     'parent': 'XenonVeto1',
     'parameters': {'shift_z': position_lxe - dis_gate - dis_cathode - dis_electric_rings + 1/2 * height_teflon_lxepart}},
    {'name': 'BottomScreening_rings',
     'parent': 'XenonVeto1',
     'parameters': {'shift_z': position_lxe - dis_gate - dis_cathode - dis_botscreen + dis_electric_rings + 1/2 * height_teflon_lxepart}},
    {'name': 'Gate_rings',
     'parent': 'XenonVeto1',
     'parameters': {'shift_z': position_lxe - dis_gate - dis_electric_rings + 1/2 * height_teflon_lxepart}},
    # Anode_rings 位于气氙区（OuterGasXe），v2 clean 中移除
    # {'name': 'Anode_rings',
    #  'parent': 'OuterGasXe',
    #  'parameters': {'shift_z': position_lxe + dis_anode + dis_electric_rings}}
]

for electrodering in electrodes_rings:
    electrodering.update({'type': 'Prism'})
    electrodering['parameters'].update({'outer_radius': radius_outer_shapingrings,
                                        'inner_radius': radius_inner_shapingrings,
                                        'height': thick_electric_rings,
                                        'sides': 0,
                                        'material': SS304LSteel,
                                        'surface': surface})

detectorList += electrodes_rings

if force_sd:
    for d in detectorList:
        d.update({'sd': True})

volume = {'soliname': 'Solid',
          'logivol': 'Log',
          'physvol': ''}

for d in detectorList:
    detector = doc.createElement('detector')
    geometry.appendChild(detector)
    detector.setAttribute('type', d['type'])
    detector.setAttribute('name', d['name'])
    detector.setAttribute('parent', d['parent'])
    djson['geometry']['detectors'].append(dict())
    djson['geometry']['detectors'][-1].update({'type': d['type']})
    djson['geometry']['detectors'][-1].update({'name': d['name']})
    djson['geometry']['detectors'][-1].update({'parent': d['parent']})
    djson['geometry']['detectors'][-1].update({'parameters': dict()})
    for k in volume.keys():
        parameter = doc.createElement('parameter')
        detector.appendChild(parameter)
        parameter.setAttribute('name', k)
        name = d['parameters'].get('physvol', d['name'])
        parameter.setAttribute('value', f'{name}{volume[k]}')
        djson['geometry']['detectors'][-1]['parameters'].update({k: f'{name}{volume[k]}'})
    for x in d['parameters'].keys():
        if x == 'physvol':
            continue
        is_surface = d['parameters'].get('surface', False)
        if not is_surface and (x in ['model', 'finish', 'type', 'alpha']):
            continue
        parameter = doc.createElement('parameter')
        detector.appendChild(parameter)
        parameter.setAttribute('name', f'{x}')
        if isinstance(d['parameters'][x], str):
            parameter.setAttribute('value', f"{d['parameters'][x]}")
            djson['geometry']['detectors'][-1]['parameters'].update({f'{x}': f"{d['parameters'][x]}"})
        elif isinstance(d['parameters'][x], bool):
            parameter.setAttribute('value', f"{int(d['parameters'][x])}")
            djson['geometry']['detectors'][-1]['parameters'].update({f'{x}': f"{int(d['parameters'][x])}"})
        else:
            parameter.setAttribute('value', f"{d['parameters'][x]}*cm")
            djson['geometry']['detectors'][-1]['parameters'].update({f'{x}': f"{d['parameters'][x]}*cm"})
    if d.pop('sd', None):
        parameter = doc.createElement('parameter')
        detector.appendChild(parameter)
        parameter.setAttribute('name', 'sendname')
        name = d['name']
        parameter.setAttribute('value', f'{name}')
        djson['geometry']['detectors'][-1]['parameters'].update({'sendname': f'{name}'})
        if d.pop('opticalsd', False):
            parameter = doc.createElement('parameter')
            detector.appendChild(parameter)
            parameter.setAttribute('name', 'opticalsd')
            parameter.setAttribute('value', '1')
            djson['geometry']['detectors'][-1]['parameters'].update({'opticalsd': '1'})
        if enable_track:
            parameter = doc.createElement('parameter')
            detector.appendChild(parameter)
            parameter.setAttribute('name', 'entrack')
            parameter.setAttribute('value', '1')
            djson['geometry']['detectors'][-1]['parameters'].update({'entrack': '1'})
    if check_overlap:
        parameter = doc.createElement('parameter')
        detector.appendChild(parameter)
        parameter.setAttribute('name', 'check_overlap')
        parameter.setAttribute('value', '1')
        djson['geometry']['detectors'][-1]['parameters'].update({'check_overlap': '1'})

# ========== Physics, Generator, Analysis ==========
physics = doc.createElement('physics')
root.appendChild(physics)
djson['physics'] = dict()
physics.setAttribute('name', physics_name)
djson['physics'].update({'name': physics_name})

parameter = doc.createElement('parameter')
physics.appendChild(parameter)
djson['physics'].update({'parameters': dict()})
parameter.setAttribute('name', 'cutlength')
parameter.setAttribute('value', f'{cutlength}*mm')
djson['physics']['parameters'].update({'cutlength': f'{cutlength}*mm'})

generator = doc.createElement('generator')
root.appendChild(generator)
djson['generator'] = dict()
generator.setAttribute('name', generator_name)
djson['generator'].update({'name': generator_name})

djson['generator'].update({'parameters': dict()})

if gen == 'muon' and not optical:
    parameter = doc.createElement('parameter')
    generator.appendChild(parameter)
    parameter.setAttribute('name', 'shielding_x')
    parameter.setAttribute('value', f'{muon_x}*cm')
    djson['generator']['parameters'].update({'shielding_x': f'{muon_x}*cm'})

    parameter = doc.createElement('parameter')
    generator.appendChild(parameter)
    parameter.setAttribute('name', 'shielding_z')
    parameter.setAttribute('value', f'{muon_z}*cm')
    djson['generator']['parameters'].update({'shielding_z': f'{muon_z}*cm'})

    parameter = doc.createElement('parameter')
    generator.appendChild(parameter)
    parameter.setAttribute('name', 'E_low')
    parameter.setAttribute('value', '0.106*GeV')
    djson['generator']['parameters'].update({'E_low': '0.106*GeV'})

if gen == 'CRN' and not optical:
    parameter = doc.createElement('parameter')
    generator.appendChild(parameter)
    parameter.setAttribute('name', 'shielding_x')
    parameter.setAttribute('value', f'{CRN_x}*cm')
    djson['generator']['parameters'].update({'shielding_x': f'{CRN_x}*cm'})

    parameter = doc.createElement('parameter')
    generator.appendChild(parameter)
    parameter.setAttribute('name', 'shielding_z')
    parameter.setAttribute('value', f'{CRN_z}*cm')
    djson['generator']['parameters'].update({'shielding_z': f'{CRN_z}*cm'})

analysis = doc.createElement('analysis')
root.appendChild(analysis)
djson['analysis'] = dict()
analysis.setAttribute('name', analysis_name)
djson['analysis'].update({'name': analysis_name})

djson['analysis'].update({'parameters': dict()})
for k in analysisList.keys():
    parameter = doc.createElement('parameter')
    analysis.appendChild(parameter)
    parameter.setAttribute('name', k)
    parameter.setAttribute('value', f'{analysisList[k]}')
    djson['analysis']['parameters'].update({k: analysisList[k]})

fp = open(fopt, 'w')
doc.writexml(fp,
             indent='  ',
             addindent='  ',
             newl='\n',
             encoding="UTF-8")

with open(os.path.splitext(fopt)[0] + '.json', 'w') as fp:
    json.dump(djson, fp, indent='\t')
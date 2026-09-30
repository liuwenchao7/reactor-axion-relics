import os
import copy
import json
import argparse

import numpy as np
import pandas as pd
import shutil
from Reactor_gamma_relics.R10_v2_clean.Multi.energy.muonRate import muon_rate

parser = argparse.ArgumentParser(description='Generate normalization factor after simulation')

parser.add_argument('--InputFile', dest='InputFile',
                    action='store', required=True,
                    help='Input .json file name, '
                         'stored simulated events number')

parser.add_argument('--paramsFile', dest='paramsFile',
                    action='store', required=True,
                    help='params  .json file name')

parser.add_argument('--GeoFile', dest='GeoFile',
                    action='store', required=True,
                    help='Geometry setup .json file name')

parser.add_argument('--OutputFile', dest='OutputFile',
                    action='store', required=True,
                    help='Output .json file name')

args = parser.parse_args()

InputFilename = args.InputFile
paramsFilename = args.paramsFile
GeoFilename = args.GeoFile
OutputFilename = args.OutputFile

shutil.copy(paramsFilename, '../analysis/geo_params.json')

date_s = 24 * 3600  # seconds per day
pmtN = 56

Kr85 = 10  # natKr/Xe concentration in unit ppt
Rn222 = 40  # Rn222 concentration in unit uBq/kg

home = os.path.expanduser('~')

with open(paramsFilename) as f:
    params = json.load(f)

with open(GeoFilename) as f:
    geo = json.load(f)

with open(InputFilename) as f:
    norm = json.load(f)

normalization = dict()

shielding_x = params['muon_x'] / 100 # m
m_rate = muon_rate() * shielding_x ** 2 * 2 * np.pi
print(f'Total muon rate is {m_rate:.2f}Hz')
print(f'Total muon rate is {m_rate / (shielding_x ** 2):.2f}Hz/m^2')

shielding_x = params['CRN_x'] / 100 # m
CRN_rate = 248.59 / 25 * shielding_x * shielding_x
print(f'Total CRN rate is {CRN_rate:.2f}Hz')
print(f'Total CRN rate is {CRN_rate / (shielding_x ** 2):.2f}Hz/m^2')


detectors = copy.deepcopy(geo['geometry']['detectors'])

for i, det in enumerate(detectors):
    print(f'{i:d}'.rjust(2, ' '), det['name'])

# https://doi.org/10.1140/epjc/s10052-019-7160-2

power = 3 # GW

# in Hz
# normalization['neutron_ON'] = dict()
# normalization['neutron_ON']['rate'] = 8.13 * 1e4 * shielding_x ** 2 * power / 3600

# normalization['neutron_OFF'] = dict()
# normalization['neutron_OFF']['rate'] = 15.03 * 1e4 * shielding_x ** 2 / date_s

# normalization['gamma_ON'] = dict()
# normalization['gamma_ON']['rate'] = 26.689 * 1e4 * shielding_x ** 2 * power

# reactors = ['neutron_ON', 'neutron_OFF', 'gamma_ON']
# topsides = ['TOP', 'SIDE']
isotopes = ['U238', 'Th232', 'Co60', 'K40', 'Cs137', 'U235']
components = ['copper', 'steel', 'Teflon', 'PMTwindow', 'PMTbase', 'Flange']
# isotopes_lXe = ['Xe127', 'Kr85', 'Rn222']
isotopes_lXe = ['Kr85', 'Rn222']

reading = False
masses = []

with open('mass.log', 'r') as f:
    for line in f.readlines():
        if 'Printed Mass of Geometry' in line:
            reading = False
        if reading:
            masses.append([])
            line = line.split(':')
            for l in line:
                for i in l.split(' '):
                    if i != '' and i != '\n':
                        masses[-1].append(i)
        if 'Printing Mass of Geometry' in line:
            reading = True

masses = masses[:-1]
assert np.unique([len(m) for m in masses]).size == 1
masses = [m[:6] for m in masses]

for i, m in enumerate(masses):
    l = m[:2]
    if m[3] == 'mg':
        scale = 1e-6
    elif m[3] == 'g':
        scale = 1e-3
    elif m[3] == 'kg':
        scale = 1
    l.append(float(m[2]) * scale)
    if m[5] == 'mg/cm3':
        scale = 1e-3
    elif m[5] == 'g/cm3':
        scale = 1
    elif m[5] == 'kg/cm3':
        scale = 1e3
    l.append(float(m[4]) * scale)
    masses[i] = l

data = pd.DataFrame(masses, columns=['geometry', 'material', 'mass', 'density'])
materials = data.groupby('material').apply(lambda x:sum(x['mass'])).to_dict()
masses = dict([(m[0], m[1:]) for m in masses])

normalization['PMTwindow'] = dict()
normalization['PMTwindow']['density'] = masses['Window'][2]  # g/cm^3
normalization['PMTwindow']['mass'] = masses['Window'][1]
normalization['PMTwindow']['volume'] = 1e3 * normalization['PMTwindow']['mass'] / normalization['PMTwindow']['density']
normalization['PMTwindow']['num'] = pmtN

normalization['PMTcasing'] = dict()
normalization['PMTcasing']['density'] = masses['Casing'][2]  # g/cm^3
normalization['PMTcasing']['mass'] = masses['Casing'][1]
normalization['PMTcasing']['volume'] = 1e3 * normalization['PMTcasing']['mass'] / normalization['PMTcasing']['density']
normalization['PMTcasing']['num'] = 0  # Because already simulated in steel.mac

normalization['PMTbase'] = dict()
normalization['PMTbase']['density'] = masses['Base'][2]  # g/cm^3
normalization['PMTbase']['mass'] = masses['Base'][1]
normalization['PMTbase']['volume'] = 1e3 * normalization['PMTbase']['mass'] / normalization['PMTbase']['density']
normalization['PMTbase']['num'] = pmtN

normalization['PMT'] = dict()
normalization['PMT']['num'] = pmtN

normalization['copper'] = dict()
normalization['copper']['density'] = masses['Shapingrings'][2]  # g/cm^3
normalization['copper']['mass'] = materials['Copper']  # kg
normalization['copper']['volume'] = 1e3 * normalization['copper']['mass'] / normalization['copper']['density']

normalization['steel'] = dict()
normalization['steel']['density'] = masses['OuterContainer'][2]  # g/cm^3, 316Ti
normalization['outersteel'] = dict()
normalization['outersteel']['mass'] = masses['OuterContainer'][1]
normalization['outersteel']['volume'] = 1e3 * normalization['outersteel']['mass'] / normalization['steel']['density']
normalization['innersteel'] = dict()
normalization['innersteel']['mass'] = masses['InnerContainer'][1]
normalization['innersteel']['volume'] = 1e3 * normalization['innersteel']['mass'] / normalization['steel']['density']
# normalization['steel']['mass'] = materials['SS304LSteel'] + (pmtN - 1) * masses['Casing'][1]  # kg
normalization['steel']['total_mass'] = materials['SS304LSteel']
normalization['steel']['mass'] = masses['InnerContainer'][1] + masses['OuterContainer'][1] + masses['Flange35'][1] + masses['Flange63'][1] + masses['Anode'][1] + masses['Casing'][1] +masses['Lxe_OutContainer'][1] +masses['BottomScreening'][1] + masses['Cathode'][1] +masses['Gate'][1]   
normalization['steel']['volume'] = 1e3 * normalization['steel']['mass'] / normalization['steel']['density']

normalization['Flange'] = dict()
normalization['Flange']['density'] = masses['OuterContainer'][2]  # g/cm^3, 316Ti
normalization['Flange']['mass'] = masses['Flange35'][1] + masses['Flange63'][1]
normalization['Flange']['volume'] = 1e3 * normalization['Flange']['mass'] / normalization['Flange']['density']

normalization['Teflon'] = dict()
normalization['Teflon']['density'] = masses['XeSideTeflon1'][2]  # g/cm^3
normalization['Teflon']['mass'] = materials['Teflon']  # kg
normalization['Teflon']['volume'] = 1e3 * normalization['Teflon']['mass'] / normalization['Teflon']['density']

normalization['lXe'] = dict()
normalization['lXe']['density'] = masses['XenonDetector'][2]  # g/cm^3

normalization['XenonVeto'] = dict()
normalization['XenonVeto']['mass'] = masses['XenonVeto'][1] + masses['XenonVeto1'][1]
normalization['XenonVeto']['volume'] = 1e3 * normalization['XenonVeto']['mass'] / normalization['lXe']['density']
normalization['XenonDetector'] = dict()
normalization['XenonDetector']['mass'] = masses['XenonDetector'][1]
normalization['XenonDetector']['volume'] = 1e3 * normalization['XenonDetector']['mass'] / normalization['lXe']['density']

normalization['lXe']['volume'] = normalization['XenonVeto']['volume'] + normalization['XenonDetector']['volume']
normalization['lXe']['mass'] = (normalization['lXe']['density'] 
                                * normalization['lXe']['volume']
                                / 1e3)  # kg

for component in components + ['lXe', 'XenonVeto', 'XenonDetector']:
    print(f"{component} mass is {normalization[component]['mass']:.4f} kg")

# Factors might need to change corresponding to simulation, where we will simulate different numbers of primary particle. 

eventALL = norm['N']

string = 'muon'
normalization[string] = dict()
normalization[string]['factor'] = norm['_' + string] * eventALL
normalization[string]['activity'] = m_rate

string = 'CRN'
normalization[string] = dict()
normalization[string]['factor'] = norm['_' + string] * eventALL
normalization[string]['activity'] = CRN_rate

for isotope in isotopes:
    for component in components:
        string = f'{isotope}_{component}'
        normalization[string] = dict()
        normalization[string]['factor'] = norm['_' + string] * eventALL

string = f'Cs137_PMTwindow'
normalization[string] = dict()
normalization[string]['factor'] = norm['_' + string] * eventALL

component = 'lXe'
for isotope in isotopes_lXe:
    string = f'{isotope}_{component}'
    normalization[string] = dict()
    normalization[string]['factor'] = norm['_' + string] * eventALL

# http://arxiv.org/abs/1101.3866
# xe100PE = 1.6e3
# mass = normalization['PE']['mass']
# normalization['U238_PE']['concentration'] = 370 / xe100PE / 1e3
# # normalization['U235_PE']['concentration'] = 0 / xe100PE / 1e3
# # normalization['Ra226_PE']['concentration'] = 0 / xe100PE / 1e3
# normalization['Th232_PE']['concentration'] = 150 / xe100PE / 1e3
# # normalization['Th228_PE']['concentration'] = 0 / xe100PE / 1e3
# normalization['Co60_PE']['concentration'] = 0 / xe100PE / 1e3
# normalization['K40_PE']['concentration'] = 1100 / xe100PE / 1e3
# normalization['Pb210_PE']['concentration'] = 0 / xe100PE / 1e3
# normalization['Cs137_PE']['concentration'] = 0 / xe100PE / 1e3
# for isotope in isotopes:
#     string = f'{isotope}_PE'
#     normalization[string]['activity'] = normalization[string]['concentration'] * mass

# xe100lead = 27.2e3  # outer layer
# mass = normalization['lead']['mass']
# normalization['U238_lead']['concentration'] = 25000 / xe100lead / 1e3
# # normalization['U235_lead']['concentration'] = 0 / xe100lead / 1e3
# # normalization['Ra226_lead']['concentration'] = 0 / xe100lead / 1e3
# normalization['Th232_lead']['concentration'] = 19600 / xe100lead / 1e3
# # normalization['Th228_lead']['concentration'] = 0 / xe100lead / 1e3
# normalization['Co60_lead']['concentration'] = 3300 / xe100lead / 1e3
# normalization['K40_lead']['concentration'] = 380 / xe100lead / 1e3
# normalization['Pb210_lead']['concentration'] = 1.4e10 / xe100lead / 1e3
# normalization['Cs137_lead']['concentration'] = 0 / xe100lead / 1e3
# for isotope in isotopes:
#     string = f'{isotope}_lead'
#     normalization[string]['activity'] = normalization[string]['concentration'] * mass

xe100copper = 2.1e3
mass = normalization['copper']['mass']
normalization['U238_copper']['concentration'] = 170 / xe100copper / 1e3
# normalization['U235_copper']['concentration'] = 0 / xe100copper / 1e3
# normalization['Ra226_copper']['concentration'] = 0 / xe100copper / 1e3
normalization['Th232_copper']['concentration'] = 25 / xe100copper / 1e3
# normalization['Th228_copper']['concentration'] = 0 / xe100copper / 1e3
normalization['Co60_copper']['concentration'] = 82 / xe100copper / 1e3
normalization['K40_copper']['concentration'] = 7 / xe100copper / 1e3
# normalization['Pb210_copper']['concentration'] = 0 / xe100copper / 1e3
# normalization['Cs137_copper']['concentration'] = 0 / xe100copper / 1e3
for isotope in isotopes:
    string = f'{isotope}_copper'
    normalization[string]['activity'] = normalization[string]['concentration'] * mass

xe100steel = 73.61  # outer layer
mass = normalization['Flange']['mass']
normalization['U238_Flange']['concentration'] = 130 / xe100steel / 1e3
# normalization['U235_steel']['concentration'] = 0 / xe100steel / 1e3
# normalization['Ra226_steel']['concentration'] = 0 / xe100steel / 1e3
normalization['Th232_Flange']['concentration'] = 140 / xe100steel / 1e3
# normalization['Th228_steel']['concentration'] = 0 / xe100steel / 1e3
normalization['Co60_Flange']['concentration'] = 400 / xe100steel / 1e3
normalization['K40_Flange']['concentration'] = 660 / xe100steel / 1e3
# normalization['Pb210_Flange']['concentration'] = 0 / xe100steel / 1e3
# normalization['Cs137_Flange']['concentration'] = 0 / xe100steel / 1e3
for isotope in isotopes:
    string = f'{isotope}_Flange'
    normalization[string]['activity'] = normalization[string]['concentration'] * mass


xe100steel = 73.61  # outer layer
mass = normalization['steel']['mass']
normalization['U238_steel']['concentration'] = 130 / xe100steel / 1e3
# normalization['U235_steel']['concentration'] = 0 / xe100steel / 1e3
# normalization['Ra226_steel']['concentration'] = 0 / xe100steel / 1e3
normalization['Th232_steel']['concentration'] = 140 / xe100steel / 1e3
# normalization['Th228_steel']['concentration'] = 0 / xe100steel / 1e3
normalization['Co60_steel']['concentration'] = 400 / xe100steel / 1e3
normalization['K40_steel']['concentration'] = 660 / xe100steel / 1e3
# normalization['Pb210_steel']['concentration'] = 0 / xe100steel / 1e3
# normalization['Cs137_steel']['concentration'] = 0 / xe100steel / 1e3
for isotope in isotopes:
    string = f'{isotope}_steel'
    normalization[string]['activity'] = normalization[string]['concentration'] * mass


xe100Teflon = 11.86  # outer layer
mass = normalization['Teflon']['mass']
normalization['U238_Teflon']['concentration'] = 0.7 / xe100Teflon / 1e3
# normalization['U235_Teflon']['concentration'] = 0 / xe100Teflon / 1e3
# normalization['Ra226_Teflon']['concentration'] = 0 / xe100Teflon / 1e3
normalization['Th232_Teflon']['concentration'] = 1.2 / xe100Teflon / 1e3
# normalization['Th228_Teflon']['concentration'] = 0 / xe100Teflon / 1e3
normalization['Co60_Teflon']['concentration'] = 0.4 / xe100Teflon / 1e3
normalization['K40_Teflon']['concentration'] = 8.9 / xe100Teflon / 1e3
# normalization['Pb210_Teflon']['concentration'] = 0 / xe100Teflon / 1e3
# normalization['Cs137_Teflon']['concentration'] = 0 / xe100Teflon / 1e3
for isotope in isotopes:
    string = f'{isotope}_Teflon'
    normalization[string]['activity'] = normalization[string]['concentration'] * mass

xe100PMTwindow = 242
num = normalization['PMTwindow']['num']
# in Bq/pieces
normalization['U238_PMTwindow']['concentration'] = 36 / xe100PMTwindow / 1e3
# normalization['U235_PMTwindow']['concentration'] = 0 / xe100PMTwindow / 1e3
# normalization['Ra226_PMTwindow']['concentration'] = 0 / xe100PMTwindow / 1e3
normalization['Th232_PMTwindow']['concentration'] = 41 / xe100PMTwindow / 1e3
# normalization['Th228_PMTwindow']['concentration'] = 0 / xe100PMTwindow / 1e3
normalization['Co60_PMTwindow']['concentration'] = 150 / xe100PMTwindow / 1e3
normalization['K40_PMTwindow']['concentration'] = 2700 / xe100PMTwindow / 1e3
# normalization['Pb210_PMTwindow']['concentration'] = 0 / xe100PMTwindow / 1e3
normalization['Cs137_PMTwindow']['concentration'] = 190 / xe100PMTwindow / 1e3
for isotope in isotopes:
    string = f'{isotope}_PMTwindow'
    normalization[string]['activity'] = normalization[string]['concentration'] * num

# xe100PMTcasing = 242
# num = normalization['PMTcasing']['num']
# # in Bq/pieces
# normalization['U238_PMTcasing']['concentration'] = 39 / xe100PMTcasing / 1e3
# # normalization['U235_PMTcasing']['concentration'] = 0 / xe100PMTcasing / 1e3
# # normalization['Ra226_PMTcasing']['concentration'] = 0 / xe100PMTcasing / 1e3
# normalization['Th232_PMTcasing']['concentration'] = 17 / xe100PMTcasing / 1e3
# # normalization['Th228_PMTcasing']['concentration'] = 0 / xe100PMTcasing / 1e3
# normalization['Co60_PMTcasing']['concentration'] = 2.4 / xe100PMTcasing / 1e3
# normalization['K40_PMTcasing']['concentration'] = 39 / xe100PMTcasing / 1e3
# normalization['Pb210_PMTcasing']['concentration'] = 0 / xe100PMTcasing / 1e3
# normalization['Cs137_PMTcasing']['concentration'] = 0 / xe100PMTcasing / 1e3
# for isotope in isotopes:
#     string = f'{isotope}_PMTcasing'
#     normalization[string]['activity'] = normalization[string]['concentration'] * num

xe100PMTbase = 242
num = normalization['PMTbase']['num']
# in Bq/pieces
normalization['U238_PMTbase']['concentration'] = 39.0 / xe100PMTbase / 1e3
# normalization['U235_PMTbase']['concentration'] = 0 / xe100PMTbase / 1e3
# normalization['Ra226_PMTbase']['concentration'] = 0 / xe100PMTbase / 1e3
normalization['Th232_PMTbase']['concentration'] = 17.0 / xe100PMTbase / 1e3
# normalization['Th228_PMTbase']['concentration'] = 0 / xe100PMTbase / 1e3
normalization['Co60_PMTbase']['concentration'] = 2.4 / xe100PMTbase / 1e3
normalization['K40_PMTbase']['concentration'] = 39.0 / xe100PMTbase / 1e3
# normalization['Pb210_PMTbase']['concentration'] = 0 / xe100PMTbase / 1e3
# normalization['Cs137_PMTbase']['concentration'] = 0 / xe100PMTbase / 1e3
for isotope in isotopes:
    string = f'{isotope}_PMTbase'
    normalization[string]['activity'] = normalization[string]['concentration'] * num

# http://arxiv.org/abs/2007.08796
mass = normalization['lXe']['mass']
# normalization['Xe127_lXe']['concentration'] = 0
# normalization['Xe127_lXe']['activity'] = normalization['Xe127_lXe']['concentration'] * mass
# normalization['Kr85_lXe']['concentration'] = Kr85 / 1.66053906660e-27 / 131.293 * 83.798 * 1e-12 * 1.7e-11 * np.log(2) / (10.76 * 365 * date_s)
normalization['Kr85_lXe']['concentration'] = Kr85 * 1e3 / 131.293 * 6.022e23 * 1e-12 * 1.7e-11 * np.log(2) / (10.76 * 365 * date_s)
normalization['Kr85_lXe']['activity'] = normalization['Kr85_lXe']['concentration'] * mass
normalization['Rn222_lXe']['concentration'] = Rn222 / 1e6
normalization['Rn222_lXe']['activity'] = normalization['Rn222_lXe']['concentration'] * mass

# for reactor in reactors:
#     for topside in topsides:
#         if 'gamma' in reactor and topside == 'TOP':
#             continue
#         string = f'{reactor}_{topside}'
#         normalization[string]['activity'] = normalization[reactor]['rate']

with open(OutputFilename, 'w') as fp:
    json.dump(normalization, fp, indent='\t')

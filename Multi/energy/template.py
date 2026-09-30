import os
import json

import h5py
import numpy as np
import pandas as pd
from numpy.lib import recfunctions
from tqdm import tqdm

import relicsapt

from Reactor_gamma_relics.R10_v2_clean.Multi.energy.cut import *

class TemplateGenerate():
    # energy bin shall not include 0
    bins = dict()
    bins['nr'] = np.logspace(-4, 2, 31)
    bins['er'] = np.linspace(1e-4, 3.5e3, 31)
    bins['x'] = np.linspace(-150, 150, 51)
    bins['y'] = np.linspace(-150, 150, 51)
    bins['z'] = np.linspace(-210, 100, 51)
    # bins['r'] = np.linspace(0, 150, 51)
    bins['r'] = np.linspace(0, 150 ** 2, 51)

    # isotopes = ['U238', 'U235', 'Ra226', 'Th232', 'Th228', 'Co60', 'K40', 'Pb210', 'Cs137']
    isotopes = ['U238', 'Th232', 'Co60', 'K40', 'Pb210', 'Cs137']
    components = ['PE', 'lead', 'copper', 'steel', 'Teflon', 'PMTwindow', 'PMTbase']#, 'PMTcasing'
    # isotopes_lXe = ['Xe127', 'Kr85', 'Rn222']
    isotopes_lXe = ['Kr85', 'Rn222']
    reactors = ['ON', 'OFF']
    topsides = ['TOP', 'SIDE']
    long_decay = ['Pb210', 'Bi210', 'Hg206', 'Po210', 'Tl206', 'Pb206']
    tagged = ['Bi214', 'Po214']

    dtype = np.dtype([('cut_er_fiducial_volume', bool),
                      ('cut_nr_fiducial_volume', bool),
                      ('cut_veto_trigger_4pi', bool),
                      ('cut_nr_single_scatter', bool),
                      ('cut_er_single_scatter', bool),
                      ('cut_er_tagging', bool),
                      ('cut_nr_roi', bool),
                      ('cut_er_roi', bool),
                      ('weights', float)])

    cut_combine = {'nr': cut_combine_nr, 'er': cut_combine_er}

    def __init__(self, folder, norm, veto=1e-2) -> None:
        self.events = dict()
        self.primaries = dict()
        self.folder = folder
        self.veto = veto
        assert self.veto >= 0. and self.veto <= 1.

        with open(norm) as f:
            self.normalization = json.load(f)

        self.mass = self.normalization['XenonDetector']['mass']
        self.mass_infv = self.mass_cylinder(r_threshold, z_max, z_min)
        self.bad_primaries = np.array(self.long_decay + self.tagged)

    def mass_cylinder(self, r_max, z_max, z_min):
        density = 2.862
        mass = (r_max / 10) ** 2 * np.pi * (z_max - z_min) / 10 * density / 1e3  # in kg
        return mass

    def change_bytes(self, events):
        for field in events.dtype.names:
            if np.issubdtype(events.dtype[field], object):
                events[field] = events[field].astype(str)
            else:
                continue

    def read_single(self, folder, key):
        self.events[key] = dict()
        self.primaries[key] = dict()
        file = os.path.join(folder, f'events_{key}.h5')
        if key == 'muon':
            file = '/home/xudc/results/muon/events_muon.h5'
        with h5py.File(file, 'r', libver='latest', swmr=True) as ipt:
            events = ipt['events'][:]
            primaries = ipt['primaries'][:]
        # self.change_bytes(events)
        # self.change_bytes(primaries)
        # Sorting like this is a bit buggy...
        # events = np.sort(events, kind='stable', order=['runId', 'eventId'])
        # primaries = np.sort(primaries, kind='stable', order=['runId', 'eventId'])
        self.events[key]['events'] = events
        self.primaries[key]['primaries'] = primaries
        if len(self.events[key]['events']) == 0:
            print(f'{key} is empty!')

    def apply_weight(self, key):
        # count per day
        date_s = 24 * 3600
        self.normalization[key]['perday'] = self.normalization[key]['activity'] * date_s

        events = self.events[key]['events']
        factor = self.normalization[key]['factor']
        # apply muon veto factor
        if key == 'muon':
            factor /= self.veto

        weight = np.full(len(events), self.normalization[key]['perday'] / factor)
        self.events[key]['weights'] = weight / self.mass  # in /day/kg

    def apply_cut(self, key):
        events = self.events[key]['events']
        cuts = np.zeros(len(events), dtype=self.dtype)
        cuts['cut_nr_fiducial_volume'] = cut_fv(events, 'nr')
        cuts['cut_er_fiducial_volume'] = cut_fv(events, 'er')
        cuts['cut_veto_trigger_4pi'] = cut_veto(events)
        cuts['cut_nr_single_scatter'] = cut_nrms(events)
        cuts['cut_er_single_scatter'] = cut_erms(events)
        cuts['cut_er_tagging'] = cut_er_tag(events)
        cuts['cut_nr_roi'] = cut_roi_nr(events)
        cuts['cut_er_roi'] = cut_roi_er(events)
        # Some primary events are not considered
        if key == 'Rn222_lXe':
            bad_events = np.isin(self.events[key]['primaries']['primaryType'].values.astype(str), self.bad_primaries)
            cuts['weights'] = np.where(bad_events, 0, self.events[key]['weights'])
        else:
            cuts['weights'] = self.events[key]['weights']
        assert len(events) == len(cuts)
        self.events[key]['cuts'] = cuts
        max_nr_r2 = events['max_nr_x'] ** 2 + events['max_nr_y'] ** 2
        max_er_r2 = events['max_er_x'] ** 2 + events['max_er_y'] ** 2
        names = cuts.dtype.names + ('max_nr_r2', 'max_er_r2', 'type')
        data = [cuts[name] for name in self.dtype.names] + [max_nr_r2, max_er_r2] + [[key] * len(events)]
        events = pd.DataFrame(events)
        for name, datum in zip(names, data):
            events[name] = datum
        self.events[key]['events'] = events

    def read(self):
        # Muon
        muonfolder = f'{self.folder}/muon'
        key = 'muon'
        print(f'Reading {key}')
        self.read_single(muonfolder, key)
        assert len(self.events[key]['events']) > 0

        # Material
        materialfolder = f'{self.folder}/material'

        for isotope in self.isotopes:
            for component in self.components:
                key = f'{isotope}_{component}'
                print(f'Reading {key}')
                self.read_single(materialfolder, key)

        component = 'lXe'
        for isotope in self.isotopes_lXe:
            key = f'{isotope}_{component}'
            print(f'Reading {key}')
            self.read_single(materialfolder, key)

        # Neutron
        neutronfolder = f'{self.folder}/neutron'
        for reactor in self.reactors:
            for topside in self.topsides:
                key = f'{reactor}_{topside}'
                print(f'Reading {key}')
                self.read_single(neutronfolder, key)

    def get_weights(self):
        for key in self.events.keys():
            self.apply_weight(key)

    def get_cuts(self):
        for key in tqdm(self.events.keys()):
            self.apply_cut(key)

    def get_hist(self):
        for key in self.events.keys():
            for e in ['nr', 'er']:
                mask = self.cut_combine[e](self.events[key]['events'])
                event = self.events[key]['events']
                weights = self.events[key]['events']['weights'] * mask
                self.events[key][e] = np.histogram(event[f'max_{e}'], bins=self.bins[e], 
                                                weights=weights)[0]

                self.events[key][f'{e}_xy'] = np.histogram2d(event[f'max_{e}_x'], event[f'max_{e}_y'], 
                                                    bins=(self.bins['x'], self.bins['y']), 
                                                    weights=weights)[0]
                r2 = (event[f'max_{e}_x'] ** 2 + event[f'max_{e}_y'] ** 2)
                # r = r2 ** 0.5
                self.events[key][f'{e}_rz'] = np.histogram2d(r2, event[f'max_{e}_z'], 
                                                    bins=(self.bins['r'], self.bins['z']), 
                                                    weights=weights)[0]

    def get_primaries(self):
        for key in self.events.keys():
            events = self.events[key]['events']
            primaries = self.primaries[key]['primaries']
            print(f'Reading primary particles of {key}...')
            primaries_index = relicsapt.get_primaries_indices(events, primaries)
            self.events[key]['primaries'] = pd.DataFrame(primaries[primaries_index])

    def big_combine(self):
        # self.combine = np.hstack([self.events[key]['events'] for key in self.events.keys()])
        self.combine = pd.concat([self.events[key]['events'] for key in self.events.keys() if 'OFF' not in key])
        self.combine_primaries = pd.concat([self.events[key]['primaries'] for key in self.events.keys() if 'OFF' not in key])

        bad_events = np.isin(self.combine_primaries['primaryType'], self.bad_primaries)
        bad_events &= self.combine['type'].values.astype(str) == 'Rn222_lXe'
        assert np.all(self.combine['weights'][bad_events] == 0)

    def add_hist(self, exclude=[]):
        self.nr_hist = np.zeros_like(self.events['muon']['nr'])
        self.er_hist = np.zeros_like(self.events['muon']['er'])
        self.nrxy_hist = np.zeros_like(self.events['muon']['nr_xy'])
        self.nrrz_hist = np.zeros_like(self.events['muon']['nr_rz'])
        self.erxy_hist = np.zeros_like(self.events['muon']['er_xy'])
        self.errz_hist = np.zeros_like(self.events['muon']['er_rz'])

        for key in self.events.keys():
            # if key in exclude:
            if np.any([ex in key for ex in exclude]):
                continue
            if 'OFF' in key:
                continue
            self.nr_hist += self.events[key]['nr']
            self.er_hist += self.events[key]['er']
            self.nrxy_hist += self.events[key]['nr_xy']
            self.nrrz_hist += self.events[key]['nr_rz']
            self.erxy_hist += self.events[key]['er_xy']
            self.errz_hist += self.events[key]['er_rz']

    def save(self, file):
        # normalize to /keV/kg/y
        np.savez(file,
                NR=self.nr_hist * 365 * self.mass / self.mass_infv / np.diff(self.bins['nr']),
                NRbins=self.bins['nr'],
                ER=self.er_hist * 365 * self.mass / self.mass_infv / np.diff(self.bins['er']),
                ERbins=self.bins['er'])

        print(f'Saved to {file}')

# Only consider events within ChainSplittingLifeTime(default 400us)
# In XenonVeto:
# 1. Using max ER energy
# 2. Considering NR quenching factor
# In XenonDetector:
# 1. Using max ER
# 2. Tagging multiple scatter
# 3. Fiducial volume cut on radius

import numpy as np

quenching = 5  # quenching factor
r_threshold = 120  # radius cut in mm
z_max = 80  # z upper cut in mm
z_min = -160  # z lower cut in mm should be the height of cathode
veto_max_er_threshold = 100  # ER threshold in veto in keV
veto_max_nr_threshold = veto_max_er_threshold * quenching  # NR threshold in veto in keV
max_er_threshold = 10  # ER threshold in detector in keV
alt_nr_threshold = 1  # alternative NR threshold in detector in keV
alt_er_ratio = 5e-2  # alternative ER threshold ratio in detector
nr_lower_threshold = 1e-1  # NR ROI lower limit in keV
nr_upper_threshold = 3  # NR ROI upper limit in keV
er_lower_threshold = 1  # ER ROI lower limit in keV
er_upper_threshold = 200  # ER ROI upper limit in keV

def cut_fv(events, s):
    mask = np.sqrt(events[f'max_{s}_x'] ** 2 + events[f'max_{s}_y'] ** 2) < r_threshold
    mask &= events[f'max_{s}_z'] > z_min
    mask &= events[f'max_{s}_z'] < z_max  # TODO: also add for S2Only?
    return mask

def cut_veto_er(events):
    mask = events['veto_max_er'] > veto_max_er_threshold
    return ~mask

def cut_veto_nr(events):
    mask = events['veto_max_nr'] > veto_max_nr_threshold
    return ~mask

def cut_veto(events):
    mask = cut_veto_er(events) & cut_veto_nr(events)
    return mask

def cut_er_tag(events):
    mask = events['max_er'] > max_er_threshold
    return ~mask

def cut_nrms(events):
    mask = events['alt_nr'] > alt_nr_threshold
    return ~mask

def cut_roi_nr(events):
    mask = events['max_nr'] < nr_upper_threshold
    mask &= events['max_nr'] > nr_lower_threshold
    return mask

def cut_erms(events):
    mask = events['alt_er'] > alt_er_ratio * events['max_er']
    return ~mask

def cut_roi_er(events):
    mask = events['max_er'] < er_upper_threshold
    mask &= events['max_er'] > er_lower_threshold
    return mask

def cut_combine_nr(events, roi=False):
    mask = cut_veto(events) & cut_fv(events, 'nr') & cut_er_tag(events) & cut_nrms(events)
    if roi:
        mask &= cut_roi_nr(events)
    return mask

def cut_combine_er(events, roi=False):
    mask = cut_veto(events) & cut_fv(events, 'er') & cut_erms(events)
    if roi:
        mask &= cut_roi_er(events)
    return mask

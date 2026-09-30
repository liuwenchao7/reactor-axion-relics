import argparse

import h5py
import numpy as np
from tqdm import tqdm
from Reactor_gamma_relics.R10_v2_clean.Multi.energy.dtypes import event_dtype, primary_dtype, cluster_dtype

parser = argparse.ArgumentParser(description='Combine h5 files together')

parser.add_argument('--InputFiles', dest='InputFiles',
                    action='store', required=True, type=str,
                    help='Input .h5 file names parser')

parser.add_argument('--OutputFile', dest='OutputFile',
                    action='store', required=True,
                    help='Output .h5 file name')

parser.add_argument('-N', dest='N', type=int,
                    action='store', required=True,
                    help='Number of .h5 file name')

args = parser.parse_args()

InputFiles = args.InputFiles
OutputFile = args.OutputFile
N = args.N

events = []
primaries = []
# clusters = []

for n in tqdm(range(1, N + 1)):
    file = InputFiles % n
    with h5py.File(file, 'r', libver='latest', swmr=True) as ipt:
        events.append(ipt['events'][:])
        primaries.append(ipt['primaries'][:])
        # clusters.append(ipt['clusters'][:])

events = np.hstack(events).astype(event_dtype)
primaries = np.hstack(primaries).astype(primary_dtype)
# clusters = np.hstack(clusters).astype(cluster_dtype)

with h5py.File(OutputFile, 'w') as opt:
    opt.create_dataset('events', data=events, compression='gzip')
    opt.create_dataset('primaries', data=primaries, compression='gzip')
    # opt.create_dataset('clusters', data=primaries, compression='gzip')

print(f'{OutputFile} saved')
print(f'{len(events)} events in the file')

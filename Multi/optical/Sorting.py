###############################################################
## Fetching optical photon hits number for different vertex
################################################################

import argparse

from ROOT import TFile
from tqdm import tqdm
import numpy as np
import h5py

parser = argparse.ArgumentParser(description='Cluster energy deposition based on position')

parser.add_argument('--InputFile', dest='InputFile',
                    action='store', required=True,
                    help='Input Geant4 ROOT file')

parser.add_argument('--OutputFile', dest='OutputFile',
                    action='store', required=True,
                    help='Output .h5 file name')

parser.add_argument('--MaxPMT', dest='MaxPMT', type=int,
                    action='store', default=175,
                    help='Max PMT number')

parser.add_argument('--MaxInnerPMT', dest='MaxInnerPMT', type=int,
                    action='store', default=127,
                    help='Max inner PMT number')

args = parser.parse_args()

InputFile = args.InputFile
OutputFile = args.OutputFile
MaxPMT = args.MaxPMT
MaxInnerPMT = args.MaxInnerPMT

# Load input ROOT file

Infile = TFile(InputFile)
tree = Infile.Get('mcTree')

eventN = int(tree.GetEntries())

print(f'Loading: {InputFile}')
print(f'Total event number = {eventN}')

def set_nan_defaults(result):
    for field in result.dtype.names:
        if np.issubdtype(result.dtype[field], np.integer):
            result[field][:] = 0
        elif np.issubdtype(result.dtype[field], float):
            result[field][:] = np.nan

primary_dtype = np.dtype([('runId', np.uint32),
                          ('eventId', np.uint32),
                          ('primaryEnergy', np.float64),
                          ('primaryX', np.float64),
                          ('primaryY', np.float64),
                          ('primaryZ', np.float64),
                          ('nScintillation', np.uint32),
                          ('nCerenkov', np.uint32),
                          ('nHits', np.uint32),
                          ('nHitsInner', np.uint32),
                          ('nHitsOuter', np.uint32),
                          ('hits_per_channel', np.uint32, MaxPMT + 1)])

# Loop over & clustering

nPrimaries = []
for ii in range(eventN):
    tree.GetEntry(ii)
    nPrimaries.append(int(tree.nPrimaries))
    assert nPrimaries[-1] == 1

primary = np.zeros(sum(nPrimaries), dtype=primary_dtype)
set_nan_defaults(primary)

for ii in tqdm(range(eventN)):
    tree.GetEntry(ii)
    primary[ii]['runId'] = int(tree.runId)
    primary[ii]['eventId'] = int(tree.eventId)
    primary[ii]['nScintillation'] = int(tree.nScintillation)
    primary[ii]['nCerenkov'] = int(tree.nCerenkov)
    primary[ii]['nHits'] = int(tree.nHits)
    pmtNumber = np.array(list(tree.pmtNumber))
    primary[ii]['hits_per_channel'] = np.histogram(pmtNumber, np.arange(MaxPMT + 2))[0]
    primary[ii]['nHitsInner'] = (pmtNumber <= MaxInnerPMT).sum()
    primary[ii]['nHitsOuter'] = (pmtNumber > MaxInnerPMT).sum()
    primary[ii]['primaryEnergy'] = list(tree.primaryEnergy)[0]
    primary[ii]['primaryX'] = list(tree.primaryX)[0]
    primary[ii]['primaryY'] = list(tree.primaryY)[0]
    primary[ii]['primaryZ'] = list(tree.primaryZ)[0]

assert np.all(primary['nHits'] == primary['nHitsInner'] + primary['nHitsOuter'])
assert np.all(primary['hits_per_channel'].sum(axis=1) == primary['nHits'])

with h5py.File(OutputFile, 'w') as opt:
    opt.create_dataset('primaries', data=primary, compression='gzip')

print(f'Saving: {OutputFile}')

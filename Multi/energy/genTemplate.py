import argparse

from Reactor_gamma_relics.R10_v2_clean.Multi.energy.template import TemplateGenerate

parser = argparse.ArgumentParser(description='Generate normalization factor after simulation')

parser.add_argument('--InputFolder', dest='InputFolder',
                    action='store', required=True,
                    help='Input folder stores the simulated files')

parser.add_argument('--InputNorm', dest='InputNorm',
                    action='store', required=True,
                    help='Normalization factor .json file')

parser.add_argument('--OutputFile', dest='OutputFile',
                    action='store', required=True,
                    help='Output .npz file name')

parser.add_argument('--Veto', dest='Veto',
                    action='store', default=1e-2,
                    type=float, help='Muon Veto factor')

args = parser.parse_args()

InputFolder = args.InputFolder
InputNorm = args.InputNorm
OutputFile = args.OutputFile
Veto = args.Veto

tempgen = TemplateGenerate(InputFolder, InputNorm, veto=Veto)
tempgen.read()
tempgen.get_primaries()
tempgen.get_weights()
tempgen.get_cuts()
tempgen.get_hist()
tempgen.add_hist()
tempgen.save(OutputFile)

import numpy as np

from cosmologies import ORIGINAL_COSMO, TARGET_COSMO
from cfg import OUT_PATH, SIMTYPE, SIGNAL

# the powerlaw exponents (from halo model)
POWERLAWS = dict(y={'h': 1.59,
                    'Om': 0.92,
                    'Ob': 0.84,
                    'ns': 1.52,
                    's8': 3.86},
                 T={'h': -0.86,
                    'Om': 0.36,
                    'Ob': 0.09,
                    'ns': -0.84,
                    's8': 1.93})

def GetTruth(simtype=None, signal=None) :
    simtype = simtype if simtype is not None else SIMTYPE
    signal = signal if signal is not None else SIGNAL

    sim_ident = 'IllustrisTNG300-1' if simtype=='IllustrisTNG' else 'SIMBA100'
    y_true = float(open('%s/mean_%s_data/%s/%s_mean.dat'%(OUT_PATH, signal, simtype, sim_ident), 'r').readline())

    for name, power in POWERLAWS[signal].items() :
        factor = (TARGET_COSMO[name] / ORIGINAL_COSMO[simtype][name])**power
        y_true *= factor
        print('from %s: %.2f'%(name, factor))

    return y_true

if __name__ == '__main__' :
    from sys import argv
    simtype = argv[1]
    signal = argv[2]
    print(GetTruth(simtype, signal))

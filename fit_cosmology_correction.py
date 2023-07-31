"""Small script to get cosmology correction fitting functions.

Command line argument:
    [1] simulation type -- IllustrisTNG or SIMBA

"""

from sys import argv

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

from cfg import SIM_PATH, OUT_PATH, COSMO_ASTRO_SEED

SIMTYPE = argv[1]

BLOCK_LEN = 11
VAR_PARAMS = ['Om', 's8', 'SN1', 'AGN1', 'SN2', 'AGN2', ]

class FitFun :
    def __init__(self, pivot_x, pivot_y) :
        self.pivot_x = pivot_x
        self.pivot_y = pivot_y

    def __call__(self, x, pwr) :
        return pivot_y * (x/self.pivot_x)**pwr
    


fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10,10))

par_values = np.loadtxt(COSMO_ASTRO_SEED,
                        skiprows=1000, max_rows=len(VAR_PARAMS)*BLOCK_LEN, usecols=(1,2,3,4,5,6))

for param_idx, param in enumerate(['s8', 'Om', ]) :

    par_idx = VAR_PARAMS.index(param)
    these_par_values = par_values[par_idx*BLOCK_LEN : (par_idx+1)*BLOCK_LEN, par_idx]

    for signal_idx, signal in enumerate(['y', 'T', ]) :

        _, sim_signal = np.loadtxt('%s/mean_%s_data/%s/1P_mean.dat'%(OUT_PATH, signal, SIMTYPE),
                                   unpack=True)
        sim_signal = sim_signal[par_idx*BLOCK_LEN : (par_idx+1)*BLOCK_LEN]

        ax[param_idx][signal_idx].scatter(these_par_values, sim_signal)

        ax[param_idx][signal_idx].set_xlabel(param)
        ax[param_idx][signal_idx].set_ylabel(signal)

        idx_pivot = BLOCK_LEN // 2
        pivot_x = these_par_values[idx_pivot]
        pivot_y = sim_signal[idx_pivot]

        f = FitFun(pivot_x, pivot_y)
        popt, _ = curve_fit(f, these_par_values, sim_signal)
        x = np.linspace(np.min(these_par_values), np.max(these_par_values), num=100)
        ax[param_idx][signal_idx].plot(x, f(x, *popt))

for a in ax.flatten() :
#    a.set_xscale('log')
#    a.set_yscale('log')
    pass

fig.suptitle(SIMTYPE)
plt.show()

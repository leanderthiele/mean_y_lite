"""Tries to factor out the cosmology dependence in the LH set (no mass function correction).

Uses the 1P set to get fitting functions for the dependencies on s8 and Om.
"""

import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

from cfg import SIM_PATH, OUT_PATH, SIMTYPE, SIGNAL, COSMO_ASTRO_SEED
from get_truth import GetTruth
from cosmologies import TARGET_COSMO

PLOT = True

BLOCK_LEN = 11
VAR_PARAMS = ['Om', 's8', 'SN1', 'AGN1', 'SN2', 'AGN2', ]

idx_LH, x_LH = np.loadtxt('%s/mean_%s_data/%s/LH_mean.dat'%(OUT_PATH, SIGNAL, SIMTYPE), unpack=True)

x_mean_true = GetTruth()

x_CV = np.loadtxt('%s/mean_%s_data/%s/CV_mean.dat'%(OUT_PATH, SIGNAL, SIMTYPE), usecols=1)
x_mean_CV = np.mean(x_CV)

# construct the interpolators
interpolators = {}
par_values_1P = np.loadtxt(COSMO_ASTRO_SEED,
                           skiprows=1000, max_rows=len(VAR_PARAMS)*BLOCK_LEN, usecols=(1,2,3,4,5,6))
signal_1P = np.loadtxt('%s/mean_%s_data/%s/1P_mean.dat'%(OUT_PATH, SIGNAL, SIMTYPE), usecols=1)
for param in ['s8', 'Om', ] :
    par_idx = VAR_PARAMS.index(param)
    these_par_values = par_values_1P[par_idx*BLOCK_LEN : (par_idx+1)*BLOCK_LEN, par_idx]
    these_signal_values = signal_1P[par_idx*BLOCK_LEN : (par_idx+1)*BLOCK_LEN]
    interpolators[param] = interp1d(these_par_values, these_signal_values, bounds_error=True)

    if PLOT :
        x = np.linspace(np.min(these_par_values), np.max(these_par_values), num=100)
        plt.scatter(these_par_values, these_signal_values)
        plt.plot(x, interpolators[param](x))
        plt.show()

# now rescale the LH measurements
Om_LH, s8_LH = np.loadtxt(COSMO_ASTRO_SEED,
                          max_rows=1000, usecols=(1,2), unpack=True)

x_LH *= ( interpolators['s8'](TARGET_COSMO['s8']) / interpolators['s8'](s8_LH) ) \
       *( interpolators['Om'](TARGET_COSMO['Om']) / interpolators['Om'](Om_LH) ) \
       * x_mean_true / x_mean_CV

np.savetxt('%s/mean_%s_data/%s/LH_mean_rescaled_cosmology.dat'%(OUT_PATH, SIGNAL, SIMTYPE),
           np.stack((np.array(idx_LH), np.array(x_LH))).T,
           header='rescaled to ground truth, including cosmology (but no mass function) correction\n'\
                  'Cosmology rescalings computed using simple interpolators through 1P measurements\n'\
                  'idx, %s'%SIGNAL,
           fmt=['%.3d', '%.18e'])

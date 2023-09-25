"""Tries to factor out the cosmology dependence in the LH set as well as the mass function.

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

idx_CV, x_CV = np.loadtxt('%s/mean_%s_data/%s/CV_mean.dat'%(OUT_PATH, SIGNAL, SIMTYPE), unpack=True)

# NOTE we also need to do the rescaling on the CV set to get the mean right!
idx_CV_hm, x_CV_hm_measured, x_CV_hm_Tinker \
    = np.loadtxt('%s/mean_%s_data/%s/CV_custom_hmf_mean_sigma0.19.dat'%(OUT_PATH, SIGNAL, SIMTYPE),
                 unpack=True)
assert np.allclose(idx_CV, idx_CV_hm)
x_CV *= x_CV_hm_Tinker / x_CV_hm_measured

x_mean_CV = np.mean(x_CV)

# construct the interpolators
interpolators = {}
par_values_1P = np.loadtxt(COSMO_ASTRO_SEED,
                           skiprows=1000, max_rows=len(VAR_PARAMS)*BLOCK_LEN, usecols=(1,2,3,4,5,6))
idx_1P, x_1P = np.loadtxt('%s/mean_%s_data/%s/1P_mean.dat'%(OUT_PATH, SIGNAL, SIMTYPE),
                          unpack=True)
idx_1P_hm, x_1P_hm_measured, x_1P_hm_Tinker \
    = np.loadtxt('%s/mean_%s_data/%s/1P_custom_hmf_mean_sigma0.19.dat'%(OUT_PATH, SIGNAL, SIMTYPE),
                 unpack=True)
assert np.allclose(idx_1P, idx_1P_hm)
x_1P *= x_1P_hm_Tinker / x_1P_hm_measured
for param in ['s8', 'Om', ] :
    par_idx = VAR_PARAMS.index(param)
    these_par_values = par_values_1P[par_idx*BLOCK_LEN : (par_idx+1)*BLOCK_LEN, par_idx]
    these_signal_values = x_1P[par_idx*BLOCK_LEN : (par_idx+1)*BLOCK_LEN]
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

idx_LH_hm, x_LH_hm_measured, x_LH_hm_Tinker \
    = np.loadtxt('%s/mean_%s_data/%s/LH_custom_hmf_mean_sigma0.19.dat'%(OUT_PATH, SIGNAL, SIMTYPE),
                 unpack=True)
assert np.allclose(idx_LH, idx_LH_hm)
x_LH *= x_LH_hm_Tinker / x_LH_hm_measured

np.savetxt('%s/mean_%s_data/%s/LH_mean_rescaled_cosmology_hmf.dat'%(OUT_PATH, SIGNAL, SIMTYPE),
           np.stack((np.array(idx_LH), np.array(x_LH))).T,
           header='rescaled to ground truth, including cosmology and mass function correction\n'\
                  'Cosmology rescalings computed using simple interpolators through 1P measurements\n'\
                  'idx, %s'%SIGNAL,
           fmt=['%.3d', '%.18e'])

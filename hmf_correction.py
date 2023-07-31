"""Command line arguments:

    [1] chunk size -- how many simulations this process should handle
    [2] chunk index -- [ 0, #Total/chunk size - 1 ]

Can also omit, then all simulations are being done.
"""

from sys import argv

import numpy as np
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter1d

import h5py

from cosmologies import TARGET_COSMO
from cfg import SIM_PATH, OUT_PATH, SIMTYPE, SUITE, NSIMS, NSNAPS, SIGNAL, COSMO_ASTRO_SEED
from halomodel_meany import meany, c_str, c_char

if len(argv) > 1 :
    CHUNK_SIZE = int(argv[1])
    CHUNK_IDX = int(argv[2])
else :
    CHUNK_SIZE = NSIMS
    CHUNK_IDX = 0

M_MIN = 1e11
M_MAX = 1e16
M_NBINS = 1024

# width of the Gaussian convolution (in dex)
SIGMA_LOGM_DEX = 0.19

# whether we try to compensate for the change in cosmology by running
# the Tinker mass function with a fiducial cosmology
COSMO_CONST = False

# load the class parameter template file
CLASS_TEMPLATE = open('template_class.ini', 'r').read()

# -------------------------------------------------------------------------#
# -------------------------------------------------------------------------#

# bins for the log(mass) histogram
LOGM_BINS = np.linspace(np.log(M_MIN), np.log(M_MAX), num=M_NBINS+1)
LOGM_CENTERS = 0.5 * (LOGM_BINS[1:] + LOGM_BINS[:-1])
DLOGM = LOGM_BINS[1] - LOGM_BINS[0]
SIGMA_LOGM = SIGMA_LOGM_DEX * np.log(10)

sim_idx_arr = []
custom_hmf_mean_arr = []
Tinker_hmf_mean_arr = []

h = TARGET_COSMO['h']
Ob = TARGET_COSMO['Ob']

# loop over realizations
for sim_idx in range(CHUNK_IDX * CHUNK_SIZE, (CHUNK_IDX+1) * CHUNK_SIZE) :
    
    if sim_idx >= NSIMS :
        break

    print('sim_idx = ', sim_idx)
    
    # figure out the cosmology
    # NOTE that for some reason SIMBA doesn't have OmegaBaryon in the snap_*.hdf5 files
    with h5py.File('%s/%s/%s_%d/snap_000.hdf5'%(SIM_PATH, SIMTYPE, SUITE, sim_idx), 'r') as f :
        Om = f['Header'].attrs['Omega0'] # we only load this as a useful cross check
        BoxSize = f['Header'].attrs['BoxSize'] # kpc/h

        if SIMTYPE == 'IllustrisTNG' :
            assert abs(f['Header'].attrs['HubbleParam']-h) < 1e-5
            assert abs(f['Header'].attrs['OmegaBaryon']-Ob) < 1e-5

    # lines in parameter table
    start = 0 if SUITE == 'LH' else 1066 if SUITE == 'CV' else 1000 if SUITE == '1P' else None
    
    idx_list, Om_list, s8_list = np.loadtxt(COSMO_ASTRO_SEED,
                                            skiprows=start, max_rows=NSIMS, usecols=(0,1,2), unpack=True,
                                            converters={0: lambda _s: float(_s.decode().split('_')[-1])})
    assert int(idx_list[sim_idx]) == sim_idx
    assert abs(Om - Om_list[sim_idx]) < 1e-5
    
    # write the CLASS .ini file
    class_fname = '%s/mean_%s_data/class_ini_files/%s/%s_%d.ini'%(OUT_PATH, SIGNAL, SIMTYPE, SUITE, sim_idx)
    open(class_fname, 'w').write(CLASS_TEMPLATE.format(h=h, Ob=Ob, Omega_cdm=Om-Ob, s8=s8_list[sim_idx]))

    z_arr = []
    hmf_arr = []

    # now figure out the mass function
    for snap_idx in range(NSNAPS) :
        
        with h5py.File('%s/%s/%s_%d/fof_subhalo_tab_%.3d.hdf5'%(SIM_PATH,
                                                                SIMTYPE,
                                                                SUITE,
                                                                sim_idx,
                                                                snap_idx), 'r') as f :
            z = f['Header'].attrs['Redshift']
            M = f['Group/Group_M_Mean200'][...] # in 1e10 Msun/h

        # convert into Msun
        M *= 1e10 / h

        # filter out zero masses
        M = M[M>1]

        # get the histogram
        hist, _ = np.histogram(np.log(M), bins=LOGM_BINS, density=False)

        # normalize properly
        hist = hist.astype(np.float64) / ( DLOGM * (1e-3*BoxSize/h)**3 )

        z_arr.append(z)
        hmf_arr.append(hist)

    # convert into numpy
    z_arr = np.array(z_arr, dtype=np.float64)
    hmf_arr = np.array(hmf_arr, dtype=np.float64) # [z x M]

    # sort the redshifts, the GSL interpolator requires this!
    sorter = np.argsort(z_arr)
    z_arr = z_arr[sorter]
    hmf_arr = hmf_arr[sorter]

    # apply convolution with Gaussian to make HMF a bit smoother
    hmf_arr_smoothed = gaussian_filter1d(hmf_arr, sigma=SIGMA_LOGM/DLOGM, axis=-1, mode='nearest')

    # FIXME debugging
    if False :
        fig, ax = plt.subplots(ncols=2)
        ax[0].matshow(np.log(hmf_arr), aspect='auto')
        ax[1].matshow(np.log(hmf_arr_smoothed), aspect='auto')
        fig.savefig('hmf_correction_debugging_%d.pdf'%sim_idx, bbox_inches='tight')

    # now call the halo model compiled code
    custom_hmf_mean = meany(c_str(class_fname),
                            c_char(SIGNAL),
                            1, # we want custom mass function
                            len(z_arr), # number of redshift sample points
                            len(LOGM_CENTERS), # number of mass sample points
                            z_arr, LOGM_CENTERS,
                            # NOTE the GSL uses column-major convention for some reason
                            #      (only in the 2D interpolation routine!)
                            #      That's why we need to transpose the array
                            hmf_arr_smoothed.T.flatten())

    Tinker_hmf_mean = meany(c_str(class_fname if not COSMO_CONST else 'CAMELS_fid.ini'),
                            c_char(SIGNAL),
                            0, # we want Tinker mass function
                            0, 0, np.zeros(1), np.zeros(1), np.zeros(1))

    sim_idx_arr.append(sim_idx)
    custom_hmf_mean_arr.append(custom_hmf_mean)
    Tinker_hmf_mean_arr.append(Tinker_hmf_mean)

np.savetxt('%s/mean_%s_data/%s/%s_custom_hmf_mean_chunk%d.dat'%(OUT_PATH, SIGNAL, SIMTYPE, SUITE, CHUNK_IDX), 
           np.stack((np.array(sim_idx_arr),
                     np.array(custom_hmf_mean_arr),
                     np.array(Tinker_hmf_mean_arr))).T,
           header='Configs: M_MIN=%.2e, M_MAX=%.2e, M_NBINS=%d, SIGMA_LOGM_DEX=%f\n\n'\
                  'idx, %s (custom HMF), %s (Tinker HMF)'%(M_MIN, M_MAX, M_NBINS, SIGMA_LOGM_DEX, SIGNAL, SIGNAL),
           fmt=['%.3d', '%.18e', '%.18e'])

# either SIMBA or IllustrisTNG
SIMTYPE = 'IllustrisTNG'

# 1P, CV, LH
SUITE = 'LH'

# either 'y' or 'T'
SIGNAL = 'y'

# ======================================= #

# to switch between Rusty and Tiger
SIM_PATH = '/mnt/ceph/users/camels/Sims'
#OUT_PATH = '/mnt/ceph/users/lthiele'
OUT_PATH = 'data'

# difference between Tiger and Rusty layout
COSMO_ASTRO_SEED = '%s/%s/CosmoAstroSeed_params.txt'%(SIM_PATH, SIMTYPE) \
                   if 'mnt/ceph' in SIM_PATH \
                   else '%s/CosmoAstroSeed_params_%s.txt'%(SIM_PATH, SIMTYPE) \
                   if 'projects/QUIJOTE' in SIM_PATH \
                   else None

NSIMS = 27 if SUITE == 'CV' \
        else 66 if SUITE == '1P' \
        else 1000 if SUITE == 'LH' \
        else None

NSNAPS = 34

assert SIMTYPE in ['IllustrisTNG', 'SIMBA', ]
assert SUITE in ['CV', '1P', 'LH', ]
assert SIGNAL in ['y', 'T', ]

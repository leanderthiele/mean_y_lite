"""collect some cosmological parameters"""

# the cosmology the `truth' sim was run with
ORIGINAL_COSMO = dict(IllustrisTNG={'h': 0.6774,
                                    'Om': 0.3089,
                                    'Ob': 0.0486,
                                    'ns': 0.9667,
                                    's8': 0.8159},
                      SIMBA={'h': 0.68,
                             'Om': 0.3,
                             'Ob': 0.048,
                             'ns': 0.97,
                             's8': 0.82})

# the cosmology we want to rescale to
TARGET_COSMO = {'h': 0.6711,
                'Om': 0.3,
                'Ob': 0.049,
                'ns': 0.9624,
                's8': 0.8}


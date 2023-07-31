import numpy as np

import h5py

from cfg import SIM_PATH, OUT_PATH, SIMTYPE, SUITE, NSIMS, NSNAPS, SIGNAL

for sim_idx in range(NSIMS) :

    print('sim_idx = ', sim_idx)

    z_arr = []
    S_avg_arr = []

    for snap_idx in range(NSNAPS) :

        with h5py.File('%s/%s/%s_%d/snap_%.3d.hdf5'%(SIM_PATH,
                                                     SIMTYPE,
                                                     SUITE,
                                                     sim_idx,
                                                     snap_idx), 'r') as f :

            a = f['Header'].attrs['Time']
            z = f['Header'].attrs['Redshift']
            h = f['Header'].attrs['HubbleParam']

            # 1e10 Msun/h / (ckpc/h)**3
            rho = f['PartType0/Density'][...]

            # dimensionless
            x = f['PartType0/ElectronAbundance'][...]

            # 1e10 Msun/h
            m = f['PartType0/Masses'][...]

            # (km/s)**2
            e = f['PartType0/InternalEnergy'][...]

        gamma = 5/3
        XH = 0.76
        c0 = 299792.458 # speed of light, km/s
        mp = 938280 # proton mass, keV

        # electron pressure times particle volume
        # 1e10 Msun/h * (km/s)**2
        S_times_V = 4.0 * x * XH / (1.0 + 3.0*XH + 4.0*XH*x) * (gamma-1.0) * m * e
        if SIGNAL == 'T' :
            # *electron* temperature in keV
            T = (gamma-1.0) * e * 4.0 / (1.0 + 3.0*XH + 4.0*XH*x) * mp / c0**2
            S_times_V *= T

        # particle volumes
        # (ckpc/h)**3
        V = m / rho

        # mean pressure in some units
        # 'y': 1e10 Msun/h * (km/s)**2 / (ckpc/h)**3
        # 'T': 1e10 Msun/h * (km/s)**2 / (ckpc/h)**3 * keV
        S_avg = np.sum(S_times_V) / np.sum(V)

        # get rid of little h
        # 'y': 1e10 Msun * (km/s)**2 / ckpc**3
        # 'T': 1e10 Msun * (km/s)**2 / ckpc**3 * keV
        S_avg *= h**2

        z_arr.append(z)
        S_avg_arr.append(S_avg)

    # save to file
    np.savetxt('%s/mean_%s_data/%s/%s_%d.dat'%(OUT_PATH, SIGNAL, SIMTYPE, SUITE, sim_idx), np.stack((np.array(z_arr), np.array(S_avg_arr))).T,
               header='z, <%s> [1e10 Msun * (km/s)**2 / ckpc**3%s]'%('P' if SIGNAL=='y' else 'P*T', '' if SIGNAL=='y' else ' * keV'))

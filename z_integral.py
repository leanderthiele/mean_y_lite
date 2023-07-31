import numpy as np
from scipy.integrate import simps

from matplotlib import pyplot as plt

def zIntegral(fname, h, Om, plot_integrand=False) :

    z, P = np.loadtxt(fname, unpack=True)

    sorter = np.argsort(z)
    z = z[sorter]
    P = P[sorter]

    # (km/s) / kpc
    H = h * 0.1 * np.sqrt( Om*(1+z)**3 + (1-Om) )

    # 1e10 Msun * (km/s) / kpc^2
    integrand = (1+z)**2 / H * P

    if plot_integrand :
        plt.plot(z, integrand)
        plt.xlabel('redshift $z$')
        plt.ylabel('integrand')
        plt.show()

    integral = simps(integrand, x=z)

    # now get the units right

    # kpc^2
    sigma_T = 6.98684e-68

    # 1e10 Msun
    m_e = 4.58110e-71

    # km/s
    c_0 = 299792.458

    y = sigma_T / m_e / c_0 * integral

    return y

avg_signal.py
  measure the average signals in individual boxes

z_integral.py
  perform the quadrature over redshift, Eq. (6) in the paper

fit_cosmology_correction.py:
  uses the 1P set to fit power laws in the cosmology parameters.
  CAUTION: I recall that at some point the 1P convention was changed,
           but I don't remember if it was before or after this project.
  NOTE: In the final version we actually didn't use these fits but
        rather used halo-model fits.
        These don't depend on the feedback implementation,
        so one can universally use Eq. (B1) in the paper.
        At the moment, I can't recall the exact reasoning why we
        preferred this over fitting to the 1P set, but I have a plot
        showing that the halo model power laws work very well for the 1P
        set.

halomodel_meany.c:
  this is a wrapper around my hmpdf code to compute the halo model
  prediction for the monopoles.
  It is compiled into a dynamic library to be used from python.

halomodel_meany.py:
  just a small python wrapper around the above dynamic library.
  The function <meany> is the only relevant object there.
  It can compute the prediction for a custom halo mass function,
  to be passed via the trailing arguments.

hmf_correction.py:
  uses the provided FOF catalogs to compute the HMF correction factors

cfg.py:
  this has some globals which are imported by other scripts.
  When running code I would just change this file for different cases
  (not good practice, I know...)


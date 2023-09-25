avg_signal.py
  measure the average signals in individual boxes

z_integral.py
  perform the quadrature over redshift, Eq. (6) in the paper

fit_cosmology_correction.py:
  uses the 1P set to fit power laws in the cosmology parameters.
  CAUTION: I recall that at some point the 1P convention was changed,
           but I don't remember if it was before or after this project.

rescale_LH_cosmology.py:
  rescales the LH realizations to a reference cosmology, using
  interpolators through the 1P set

rescale_LH_cosmology_hmf.py:
  additional rescaling, accounting for the sample variance through the
  halo mass function.

halomodel_meany.c:
  this is a wrapper around my hmpdf code to compute the halo model
  prediction for the monopoles.
  It is compiled into a dynamic library to be used from python.
  NOTE: hmpdf mean_y branch, do not use python wrapper (broken),
        may need to update source to work with newer gcc

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

get_truth.py:
  gets the predictions from big boxes. These need to be rescaled to the
  fiducial CAMELS cosmology, which is being done with the halo model
  cosmology fits.
  We do not use the 1P fits here since they do not allow accounting for
  all differences in parameters.

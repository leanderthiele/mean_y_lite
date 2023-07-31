/* Compile as
 *
 * gcc --std=c99 --shared -fPIC -Wall -Wextra -O3 -I/home/lthiele/hmpdf/include -o libhalomodel_meany.so halomodel_meany.c -L/home/lthiele/hmpdf -lgsl -lgslcblas -lm -lhmpdf
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>

#include "hmpdf.h"

typedef struct
{//{{{
    int Nz;
    int NM;

    const double *z_grid;
    const double *logM_grid;
    const double *data; // Nz x NM

    gsl_interp2d *interpolator;
    gsl_interp_accel *zacc;
    gsl_interp_accel *Macc;
}//}}}
hmf_interpolator_params;

double
hmf_interpolator(double z, double M, void *params)
{// {{{
    hmf_interpolator_params *p = (hmf_interpolator_params *)params;

    double logM = log(M);

    // handle out-of-range issues
    if (z < p->z_grid[0])
        z = p->z_grid[0];
    else if (z > p->z_grid[p->Nz-1])
        z = p->z_grid[p->Nz-1];
    
    if (logM < p->logM_grid[0])
        logM = p->logM_grid[0];
    else if (logM > p->logM_grid[p->NM-1])
        return 0.0;

    // make sure we always return positive numbers (regardless of what interpolator says)
    return GSL_MAX(0.0, gsl_interp2d_eval(p->interpolator,
                                          p->z_grid, p->logM_grid, p->data,
                                          z, logM,
                                          p->zacc, p->Macc));
}//}}}

void
init_interpolator(hmf_interpolator_params *p)
{//{{{
    p->interpolator = gsl_interp2d_alloc(gsl_interp2d_bilinear, p->Nz, p->NM);

    p->zacc = gsl_interp_accel_alloc();
    p->Macc = gsl_interp_accel_alloc();

    gsl_interp2d_init(p->interpolator, p->z_grid, p->logM_grid, p->data, p->Nz, p->NM);
}//}}}

void
free_interpolator(hmf_interpolator_params *p)
{//{{{
    gsl_interp2d_free(p->interpolator);
    gsl_interp_accel_free(p->zacc);
    gsl_interp_accel_free(p->Macc);
}//}}}

double
meany(char *class_ini, char outtype, int custom_hmf, int Nz, int NM, const double *z_grid, const double *logM_grid, const double *hmf_data)
{//{{{
    if (outtype != 'y' && outtype != 'T')
        return -42;

    hmf_interpolator_params p;

    if (custom_hmf)
    {
        p.Nz = Nz; p.NM = NM; p.z_grid = z_grid; p.logM_grid = logM_grid; p.data = hmf_data;
        init_interpolator(&p);
    }

    hmpdf_obj *d = hmpdf_new();
    if (!d)
        return -1.0;

    if (hmpdf_init(d, class_ini, hmpdf_tsz,
                   hmpdf_rout_scale, 2.5,
                   hmpdf_custom_hmf, (custom_hmf) ? &hmf_interpolator : NULL,
                   hmpdf_custom_hmf_params, (custom_hmf) ? (void *)&p : NULL,
                   hmpdf_N_M, 100, // this is extremely well converged
                   hmpdf_N_z, 20, // this is converged
                   hmpdf_N_theta, 64))
        return -2.0;
    
    double out;
    if (((outtype == 'y') ? hmpdf_get_mean : hmpdf_get_mean_T)(d, &out))
        return -3.0;

    if (hmpdf_delete(d))
        return -4.0;

    if (custom_hmf)
        free_interpolator(&p);

    return out;
}//}}}

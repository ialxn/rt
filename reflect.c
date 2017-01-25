/*	reflect.c
 *
 * Copyright (C) 2014,2015,2016,2017 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#define _GNU_SOURCE		/* for sincos() */

#include <cblas.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "likely.h"
#include "math_utils.h"
#include "off.h"
#include "reflect.h"


void reflect_specular(ray_t * r, const double N[3],
		      const double P[3], const gsl_rng
		      __attribute__ ((__unused__)) * rng, void
		      __attribute__ ((__unused__)) * model_params)
/*
 * return specularly reflected ray 'r'. surface normal of reflecting surface
 * is 'N'. incoming ray has been determined before to intersect at 'P'. thus
 * origin or reflected ray will be 'P'.
 */
{
    const double t = cblas_ddot(3, N, 1, r->dir, 1);	/* 'N' dot 'r' */

    cblas_daxpy(3, -2.0 * t, N, 1, r->dir, 1);	/* 'r' - 2 * 'N' dot 'r' * 'N' */
    memcpy(r->orig, P, 3 * sizeof(double));	/* update origin */

    if (unlikely(r->n_refl == UCHAR_MAX))	/* too many reflections is unlikely */
	fprintf(stderr,
		"                INFO: maximum path length exceeded (wrap around of counter occurs)\n");
    ++r->n_refl;
}

void reflect_lambertian(ray_t * r, const double N[3], const double P[3],
			const gsl_rng * rng, void
			__attribute__ ((__unused__)) * model_params)
    /*
     * return diffusely reflected ray 'r'. surface normal of reflecting surface
     * is 'N'. incoming ray has been determined before to intersect at 'P'. thus
     * origin or reflected ray will be 'P'.
     */
{
    get_uniform_random_vector_hemisphere(r->dir, 1.0, N, rng);
    memcpy(r->orig, P, 3 * sizeof(double));

    if (unlikely(r->n_refl == UCHAR_MAX))	/* too many reflections is unlikely */
	fprintf(stderr,
		"                INFO: maximum path length exceeded (wrap around of counter occurs)\n");
    ++r->n_refl;
}

void reflect_microfacet_gaussian(ray_t * r, const double N[3],
				 const double P[3], const gsl_rng * rng,
				 void *model_params)
{
/*
 * surface consists of uniformly distributed random facets. surface normal of
 * these microfacets follows a gaussian distribution with zero mean (i.e. the
 * distribution is centered at 'N', the macroscopic surface normal, and has a
 * standard deviation of 'sigma'.
 *
 * algorithm:
 * - get random theta (gaussian distributed with zero mean and standard
 *   deviation sigma)	-inf < theta < +inf
 * - discard values |theta| >= pi/2 (surface normal of micro facet is tilted
 *   by maximally 90 degrees relative to the macroscopic surface normal)
 *   90 degrees, otherwise reflected ray passes through surface.)
 * - tilt 'N' by random_theta and rotate around original 'N' by (uniform)
 *   random phi ( 0 <= phi < pi ).
 * - reflect 'r' at new surface normal
 * - check that ray is reflected and not transmitted ('N' dot 'r' > 0)
 */
    double alpha, beta;
    double dummy[3];
    double original_ray_dir[3], original_ray_orig[3];
    double *sigma = (double *) model_params;

    /*
     * determine alpha / beta to transform 'N' into local
     * system where it becomes 0,0,1 (returned as 'dummy' and discarded).
     */
    g2l_off_rot(N, dummy, &alpha, &beta);

    memcpy(original_ray_dir, r->dir, 3 * sizeof(double));	/* save ray */
    memcpy(original_ray_orig, r->orig, 3 * sizeof(double));

    do {
	double theta = GSL_DBL_MAX;
	double phi;
	double cos_theta, sin_theta;
	double sin_phi, cos_phi;
	double random_N[3], new_N[3];

	memcpy(r->dir, original_ray_dir, 3 * sizeof(double));	/* restore ray */
	memcpy(r->orig, original_ray_orig, 3 * sizeof(double));

	do {			/* gaussian theta */
	    theta = gsl_ran_gaussian(rng, *sigma);
	} while (fabs(theta) > M_PI_2);

	phi = M_PI * gsl_rng_uniform(rng);

	sincos(theta, &sin_theta, &cos_theta);
	sincos(phi, &sin_phi, &cos_phi);

	random_N[0] = sin_theta * cos_phi;
	random_N[1] = sin_theta * sin_phi;
	random_N[2] = cos_theta;

	/*
	 * transform 'random_N' into global system
	 */
	l2g_off_rot(random_N, new_N, alpha, beta);

	reflect_specular(r, new_N, P, NULL, NULL);
	--r->n_refl;		/* reflect_specular increases count */

    } while (cblas_ddot(3, r->dir, 1, N, 1) < 0.0);	/* 'r' transmitted (not reflected) */

    if (unlikely(r->n_refl == UCHAR_MAX))	/* too many reflections is unlikely */
	fprintf(stderr,
		"                INFO: maximum path length exceeded (wrap around of counter occurs)\n");
    ++r->n_refl;
}

void reflect_ray(ray_t * r, const double N[3], const double P[3],
		 const gsl_rng * rng, const refl_model_t * model)
{
    int i = 0;
    const double threshold = gsl_rng_uniform(rng);

    /*
     * find out which model to apply.
     */
    while (threshold > model->defn[i]->threshold)
	++i;

    (model->defn[i]->f) (r, N, P, rng, model->defn[i]->par);
}

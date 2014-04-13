/*	reflect.c
 *
 * Copyright (C) 2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#define _GNU_SOURCE		/* for sincos() */


#include <math.h>
#include <string.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "off.h"
#include "reflect.h"


static void reflect_specular(ray_t * r, const double N[3],
			     const double P[3])
/*
 * return specularly reflected ray 'r'. surface normal of reflecting surface
 * is 'N'. incoming ray has been determined before to intersect at 'P'. thus
 * origin or reflected ray will be 'P'.
 */
{
    const double t = cblas_ddot(3, N, 1, r->dir, 1);	/* 'N' dot 'r' */

    cblas_daxpy(3, -2.0 * t, N, 1, r->dir, 1);	/* 'r' - 2 * 'N' dot 'r' * 'N' */
    memcpy(r->orig, P, 3 * sizeof(double));	/* update origin */
}

static void reflect_lambertian(ray_t * r, const double N[3],
			       const double P[3], const gsl_rng * rng)
 /*
  * return diffusely reflected ray 'r'. surface normal of reflecting surface
  * is 'N'. incoming ray has been determined before to intersect at 'P'. thus
  * origin or reflected ray will be 'P'.
  */
{
    double t, random_ray[3];
    double cos_theta, sin_theta;
    double phi, sin_phi, cos_phi;

    t = gsl_rng_uniform(rng);
    cos_theta = 1.0 - 2.0 * t;
    sin_theta = sin(acos(cos_theta));

    t = gsl_rng_uniform(rng);
    phi = 2.0 * M_PI * t;
    sincos(phi, &sin_phi, &cos_phi);

    random_ray[0] = sin_theta * cos_phi;
    random_ray[1] = sin_theta * sin_phi;
    random_ray[2] = cos_theta;

    if (cblas_ddot(3, N, 1, random_ray, 1) < 0.0) {
	/*
	 * 'N' and 'random_ray' point into opposite directions (are anti-
	 * parallel). thus ray is not reflected but transmitted. use inverted
	 * 'random_ray'.
	 */
	size_t i;

	for (i = 0; i < 3; i++)
	    r->dir[i] = -random_ray[i];
    } else
	memcpy(r->dir, random_ray, 3 * sizeof(double));

    memcpy(r->orig, P, 3 * sizeof(double));

}

static void reflect_microfacet_gaussian(ray_t * r, const double N[3],
					const double P[3],
					const gsl_rng * rng,
					const double sigma)
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
 *   random phi ( 0 <= phi < 2pi ).
 * - reflect 'r' at new surface normal
 * - check that ray is reflected and not transmitted ('N' dot 'r' > 0)
 */
    double dot_product = -1.0;
    double alpha, beta;
    const double O[] = { 0.0, 0.0, 0.0 };
    double dummy[3];
    double original_ray_dir[3];

    /*
     * determine alpha / beta to transform 'N' into local
     * system where it becomes 0,0,1 (returned as 'dummy' and discarded).
     */
    g2l_off(O, N, dummy, &alpha, &beta);

    memcpy(original_ray_dir, r->dir, 3 * sizeof(double));	/* save */

    do {
	double theta = GSL_DBL_MAX;
	double phi;
	double cos_theta, sin_theta;
	double sin_phi, cos_phi;
	double random_N[3], new_N[3];

	memcpy(r->dir, original_ray_dir, 3 * sizeof(double));	/* restore */

	do {			/* gaussian theta */
	    theta = gsl_ran_gaussian(rng, sigma);
	} while (theta > M_PI_2);

	phi = 2.0 * M_PI * gsl_rng_uniform(rng);

	sincos(theta, &sin_theta, &cos_theta);
	sincos(phi, &sin_phi, &cos_phi);

	random_N[0] = sin_theta * cos_phi;
	random_N[1] = sin_theta * sin_phi;
	random_N[2] = cos_theta;

	/*
	 * transform 'random_N' into global system
	 */
	l2g_off(O, random_N, new_N, alpha, beta);

	reflect_specular(r, new_N, P);
	dot_product = cblas_ddot(3, r->dir, 1, N, 1);

    } while (dot_product < 0.0);	/* 'r' transmitted (not reflected) */

}

void reflect(ray_t * r, const double N[3], const double P[3],
	     const char model, const gsl_rng * rng, void *model_params)
{

    switch (model) {

    case SPECULAR:
	reflect_specular(r, N, P);
	break;

    case LAMBERTIAN:
	reflect_lambertian(r, N, P, rng);
	break;

    case MICROFACET_GAUSSIAN:
	{
	    double *sigma = (double *) model_params;

	    reflect_microfacet_gaussian(r, N, P, rng, *sigma);
	    break;
	}

    }

}

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
#include <gsl/gsl_rng.h>

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
    } else {
	size_t i;

	for (i = 0; i < 3; i++)
	    r->dir[i] = random_ray[i];
    }

    memcpy(r->orig, P, 3 * sizeof(double));

}


void reflect(ray_t * r, const double N[3], const double P[3],
	     const char model, const gsl_rng * rng)
{

    switch (model) {

    case SPECULAR:
	reflect_specular(r, N, P);
	break;

    case LAMBERTIAN:
	reflect_lambertian(r, N, P, rng);
	break;

    }

}

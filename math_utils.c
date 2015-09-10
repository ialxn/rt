/*	math_utils.c
 *
 * Copyright (C) 2011,2012,2013,2014,2015 Ivo Alxneit
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

#include "io_utils.h"
#include "math_utils.h"


inline void a_plus_cb(double result[3], const double a[3], const double c,
		      const double b[3])
{
    result[0] = a[0] + c * b[0];
    result[1] = a[1] + c * b[1];
    result[2] = a[2] + c * b[2];
}

inline void a_times_const(double result[3], const double a[3],
			  const double c)
{
    result[0] = a[0] * c;
    result[1] = a[1] * c;
    result[2] = a[2] * c;
}

inline void diff(double result[3], const double a[3], const double b[3])
{
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
}

double d_sqr(const double a[3], const double b[3])
/*
 * return squared distance between point 'a' and 'b'
 */
{
    size_t i;
    double d;

    for (i = 0, d = 0.0; i < 3; i++) {
	const double t = a[i] - b[i];
	d += t * t;
    }

    return (d);

}

inline void cross_product(const double a[3], const double b[3],
			  double result[3])
{
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

double normalize(double a[3])
{
    double norm = cblas_dnrm2(3, a, 1);

    cblas_dscal(3, 1.0 / norm, a, 1);
    return norm;
}

int orthonormalize(double x[3], double y[3], double z[3])
/*
 * make sure 'x', 'y', 'z' form a orthonormal basis.
 * 1): 'z' is a given (maybe normalize).
 * 2): 'y' = 'z' cross 'x' and normalize
 * 3): verify that 'x' dot 'z' ('x' dot 'z' must be zero).
 *     if this is not the case 'y' = 'z' cross 'x'.
 *
 * return NO_ERR / ERR if 'x' was not / was changed
 */
{
    int status = NO_ERR;

    normalize(z);

    cross_product(z, x, y);
    normalize(y);

    if (cblas_ddot(3, x, 1, z, 1)) {	/* 'x' is not perpendicular to 'z' */
	cross_product(y, z, x);
	status = ERR;
    }

    return status;
}

void g2l(const double *mat, const double *origin, const double *g,
	 double *l)
/*
 * expresses vector 'vec' (global) in local coordinates
 *     l(x, y, z) = M (g(x, y, z) - o(x, y, z))
 */
{
    int i;
    double t[3];

    diff(t, g, origin);

    for (i = 0; i < 3; i++)
	l[i] = cblas_ddot(3, t, 1, &mat[3 * i], 1);
}

void l2g(const double *mat, const double *origin, const double *l,
	 double *g)
/*
 * expresses vector 'vec' (local) in global coordinates
 *     g(x, y, z) = MT l(x, y, z) + o(x, y, z)
 */
{
    int i;

    for (i = 0; i < 3; i++)
	g[i] = cblas_ddot(3, &mat[i], 3, l, 1) + origin[i];

}

void g2l_rot(const double *mat, const double *g, double *l)
/*
 * global -> local: performs rotation part only
 *     l(x, y, z) = M (g(x, y, z))
 */
{
    int i;

    for (i = 0; i < 3; i++)
	l[i] = cblas_ddot(3, g, 1, &mat[3 * i], 1);
}

void l2g_rot(const double *mat, const double *l, double *g)
/*
 * local -> global: performs rotation part only
 */
{
    int i;

    for (i = 0; i < 3; i++)
	g[i] = cblas_ddot(3, &mat[i], 3, l, 1);
}


void get_uniform_random_vector(double *result, const double l,
			       const gsl_rng * r)
{
/*
 * return random vector of length 'l'
 */
    double sin_theta, cos_theta;
    double phi, sin_phi, cos_phi;
    double R_sin_theta;

    cos_theta = 1.0 - 2.0 * gsl_rng_uniform(r);
    sin_theta = sin(acos(cos_theta));
    phi = 2.0 * M_PI * gsl_rng_uniform(r);
    sincos(phi, &sin_phi, &cos_phi);
    R_sin_theta = l * sin_theta;

    result[0] = R_sin_theta * cos_phi;
    result[1] = R_sin_theta * sin_phi;
    result[2] = l * cos_theta;
}

void get_uniform_random_vector_hemisphere(double *result,
					  const double radius,
					  const double *normal,
					  const gsl_rng * r)
{
    double random_ray[3];

    get_uniform_random_vector(random_ray, radius, r);

    if (cblas_ddot(3, normal, 1, random_ray, 1) < 0.0) {
	/*
	 * 'normal' and 'random_ray' point into opposite directions (are anti-
	 * parallel). use inverted 'random_ray'.
	 */
	size_t i;

	for (i = 0; i < 3; i++)
	    result[i] = -random_ray[i];
    } else
	memcpy(result, random_ray, 3 * sizeof(double));
}

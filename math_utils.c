/*	math_utils.c
 *
 * Copyright (C) 2011,2012,2013 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <string.h>
#include <gsl/gsl_cblas.h>

#include "io_utils.h"
#include "math_utils.h"


void a_plus_cb(double result[3], const double a[3], const double c,
	       const double b[3])
/*
 * 'result' = 'a' + 'c' * 'b'
 */
{
    size_t i;

    for (i = 0; i < 3; i++)
	result[i] = a[i] + c * b[i];
}

void diff(double result[3], const double a[3], const double b[3])
/*
 * calculate difference of vectors 'a' and 'b'
 * 'result' = 'a' - 'b'
 */
{
    size_t i;

    for (i = 0; i < 3; i++)
	result[i] = a[i] - b[i];
}

double d_sqr(const double a[3], const double b[3])
/*
 * return squared distance between point 'a' and 'b'
 */
{
    size_t i;
    double d = 0.0;

    for (i = 0; i < 3; i++) {
	const double t = a[i] - b[i];
	d += t * t;
    }

    return (d);

}

void cross_product(const double a[3], const double b[3], double result[3])
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

void reflect(ray_t * r, const double N[3], const double P[3])
{
    const double t = cblas_ddot(3, N, 1, r->direction, 1);	/* 'N' dot 'r' */

    cblas_daxpy(3, -2.0 * t, N, 1, r->direction, 1);	/* 'r' - 2 * 'N' dot 'r' * 'N' */
    memcpy(r->origin, P, 3 * sizeof(double));	/* update origin */
}

void g2l(const double *mat, const double *origin, const double *g,
	 double *l)
/*
 * expresses vector 'vec' (global) in local coordinates
 *     l(x, y, z) = MT (g(x, y, z) - o(x, y, z))
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
 *     g(x, y, z) = M l(x, y, z) + o(x, y, z)
 */
{
    int i;

    for (i = 0; i < 3; i++)
	g[i] = cblas_ddot(3, &mat[i], 3, l, 1) + origin[i];

}

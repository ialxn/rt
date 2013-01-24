/*	vector_math.h
 *
 * Copyright (C) 2011,2012,2013 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __VECTOR_MATH_H__
#define __VECTOR_MATH_H__

#include "ray.h"

extern void v_a_plus_cb(double result[3], const double a[3],
			const double c, const double b[3]);
extern void v_diff(double result[3], const double a[3], const double b[3]);
extern void cross_product(const double a[3], const double b[3],
			  double result[3]);
extern double normalize(double a[3]);
extern int orthonormalize(double x[3], double y[3], double z[3]);
extern void reflect(ray_t * r, const double N[3], const double P[3]);

extern void g2l(const double *mat, const double *origin, const double *g,
		double *l);
extern void l2g(const double *mat, const double *origin, const double *l,
		double *g);

#endif				/* __VECTOR_MATH_H__ */

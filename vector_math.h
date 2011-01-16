/*	vector_math.h
 *
 * Copyright (C) 2011 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __VECTOR_MATH_H__
#define __VECTOR_MATH_H__

#include "ray.h"

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

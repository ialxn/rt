/*	math_utils.h
 *
 * Copyright (C) 2011 - 2018 Ivo Alxneit, Paul Scherrer Institute
 *
 * This file is part of rt
 *
 * rt is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rt is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rt. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef __MATH_UTILS_H__
#define __MATH_UTILS_H__

#include <gsl/gsl_rng.h>

#include "ray.h"


#define SWAP(x, y) do { __typeof__((x)) temp = (x); (x) = (y); (y) = temp; } while (0)


extern void a_plus_cb(double result[3], const double a[3],
		      const double c, const double b[3]);
extern void a_times_const(double result[3], const double a[3],
			  const double c);
extern void diff(double result[3], const double a[3], const double b[3]);
extern double d_sqr(const double a[3], const double b[3]);
extern void cross_product(const double a[3], const double b[3],
			  double result[3]);
extern double normalize(double a[3]);
extern int orthonormalize(double x[3], double y[3], double z[3]);
extern void g2l(const double *mat, const double *origin, const double *g,
		double *l);
extern void l2g(const double *mat, const double *origin, const double *l,
		double *g);
extern void g2l_rot(const double *mat, const double *g, double *l);
extern void l2g_rot(const double *mat, const double *l, double *g);

extern void get_uniform_random_vector(double *result, const double radius,
				      const gsl_rng * r);
extern void get_uniform_random_vector_hemisphere(double *result,
						 const double radius,
						 const double *normal,
						 const gsl_rng * r);

#endif				/* __MATH_UTILS_H__ */

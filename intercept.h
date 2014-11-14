/*	intercept.h
 *
 * Copyright (C) 2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __INTERCEPT_H__
#define __INTERCEPT_H__

#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include <libconfig.h>

#include "ray.h"
#include "math_utils.h"
#include "reflect.h"
#include "targets.h"


#define LAST_WAS_HIT (1<<0)	/* target was hit by last ray */
#define ABSORBED     (1<<1)	/* ray was absorbed on target */
#define ICPT_ON_CONVEX_SIDE	(1<<2)	/* ray intercepted by convex side */

/*
 * surface normals of non-planar targets (pointing away from convex side
 */
extern void cyl_surf_normal(double *const icpt, const double *C,
			    const double *a, const double r,
			    double *const normal);
void ell_surf_normal(const double *point, const double *axes,
		     double *const normal);


/*
 * functions to return intercepts of ray with target
 */
extern double *intercept_cylinder(const ray_t * ray, const double *c,
				  const double *a, const double r,
				  const double l, int *hits_outside);
extern double *intercept_ellipsoid(const ray_t * ray, const double *M,
				   const double center[3],
				   const double axes[3],
				   const double z_min, const double z_max);
extern double *intercept_plane(const ray_t * ray,
			       const double *plane_normal,
			       const double *plane_point, int *hits_front);
extern double *intercept_sphere(const ray_t * ray, const double *center,
				const double radius);

#endif				/* __INTERCEPT_H__ */

/*	reflect.h
 *
 * Copyright (C) 2014,2015,2016 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __REFLECT_H__
#define __REFLECT_H__

#include <gsl/gsl_rng.h>

#include "ray.h"


extern void reflect_specular(ray_t * r, const double N[3],
			     const double P[3], const gsl_rng * rng,
			     void *model_params);
extern void reflect_lambertian(ray_t * r, const double N[3],
			       const double P[3], const gsl_rng * rng,
			       void *model_params);
extern void reflect_microfacet_gaussian(ray_t * r, const double N[3],
					const double P[3],
					const gsl_rng * rng,
					void *model_params);

#endif				/* __REFLECT_H__ */

/*	reflect.h
 *
 * Copyright (C) 2014 - 2018 Ivo Alxneit, Paul Scherrer Institute
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

#ifndef __REFLECT_H__
#define __REFLECT_H__

#include <gsl/gsl_rng.h>

#include "ray.h"


typedef void (*refl_func_t) (ray_t * r, const double N[3],
			     const double P[3],
			     const gsl_rng * rng, void *model_params);

typedef struct model_def_t {
    refl_func_t f;		/* reflect function */
    void *par;			/* (possible) parameters for refl_func */
    double threshold;		/* relative weights of refl_funcs */
} model_def_t;

typedef struct refl_model_t {
    int n_models;		/* number of specific reflectivity models used */
    model_def_t **defn;		/* array with n_models specific definitions */
} refl_model_t;


extern void reflect_ray(ray_t * r, const double N[3], const double P[3],
			const gsl_rng * rng, const refl_model_t * model);

extern void reflect_microfacet_gaussian(ray_t * r, const double N[3],
					const double P[3],
					const gsl_rng * rng,
					void *model_params);
extern void reflect_lambertian(ray_t * r, const double N[3],
			       const double P[3], const gsl_rng * rng,
			       void *model_params);
extern void reflect_specular(ray_t * r, const double N[3],
			     const double P[3], const gsl_rng * rng,
			     void *model_params);

#endif				/* __REFLECT_H__ */

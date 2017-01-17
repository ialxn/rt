/*	reflect.h
 *
 * Copyright (C) 2014,2015,2016,2017 Ivo Alxneit
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

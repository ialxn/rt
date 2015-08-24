/*	reflect.h
 *
 * Copyright (C) 2014,2015 Ivo Alxneit
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

/*
 * reflectivity models
 * use bits 1-15 of an int
 */
#define MODEL_NONE		0
#define SPECULAR		(1<<0)
#define LAMBERTIAN		(1<<1)
#define MICROFACET_GAUSSIAN	(1<<2)


extern void reflect(ray_t * r, const double N[3], const double P[3],
		    const int model, const gsl_rng * rng,
		    void *model_params);

#endif				/* __REFLECT_H__ */

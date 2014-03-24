/*	reflect.h
 *
 * Copyright (C) 2014 Ivo Alxneit
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
 */
#define SPECULAR	1
#define LAMBERTIAN	2


extern void reflect(ray_t * r, const double N[3], const double P[3],
		    const char model, const gsl_rng * rng);

#endif				/* __REFLECT_H__ */

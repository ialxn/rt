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

#include "ray.h"

extern void reflect(ray_t * r, const double N[3], const double P[3]);

#endif				/* __REFLECT_H__ */

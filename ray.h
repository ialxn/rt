/*	ray.h
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __RAY_H__
#define __RAY_H__

typedef struct vec_t {
    double x;
    double y;
    double z;
} vec_t;

typedef struct ray_t {
    vec_t origin;
    vec_t direction;
    double power;
} ray_t;

#endif				/* __RAY_H__ */

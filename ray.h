/*	ray.h
 *
 * Copyright (C) 2010, 2011 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __RAY_H__
#define __RAY_H__

typedef struct ray_t {
    double origin[3];
    double direction[3];
    double power;
} ray_t;

#endif				/* __RAY_H__ */

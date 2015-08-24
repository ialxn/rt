/*	ray.h
 *
 * Copyright (C) 2010,2011,2012,2013,2014,2015 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __RAY_H__
#define __RAY_H__


typedef struct ray_t {
    double orig[3];
    double dir[3];
    double lambda;		/* wavelength */
    unsigned char n_refl;	/* no more than 255 reflections */
} ray_t;

#endif				/* __RAY_H__ */

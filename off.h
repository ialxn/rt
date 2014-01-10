/*	off.h
 *
 * Copyright (C) 2010,2011,2012,2013,2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __OFF_H__
#define __OFF_H__

#include <libconfig.h>

extern void output_geometry(config_t * cfg);

extern void g2l_off(const double *P, const double *N, double *L,
		    double *alpha, double *beta);
extern void l2g_off(const double *P, const double *L, double *G,
		    const double alpha, const double beta);


#endif				/* __OFF_H__ */

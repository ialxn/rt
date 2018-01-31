/*	off.h
 *
 * Copyright (C) 2010 - 2018 Ivo Alxneit, Paul Scherrer Institute
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

#ifndef __OFF_H__
#define __OFF_H__

#include <libconfig.h>
#include <stdio.h>


extern void output_geometry(config_t * cfg);
extern FILE *open_off(const char *name);
extern void write_ray(FILE * f, const double *start, const double *stop,
		      const int r, const int g, const int b);
extern void g2l_off(const double *P, const double *N, double *L,
		    double *alpha, double *beta);
extern void l2g_off(const double *P, const double *L, double *G,
		    const double alpha, const double beta);
extern void g2l_off_rot(const double *N, double *L, double *alpha,
			double *beta);
extern void l2g_off_rot(const double *L, double *G, const double alpha,
			const double beta);



#endif				/* __OFF_H__ */

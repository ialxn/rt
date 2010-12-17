/*	off.h
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __OFF_H__
#define __OFF_H__

extern void off_sphere(const char *name, double *O, const double radius,
		       const double r, const double g, const double b);
extern void off_axes(const double size);
extern void off_plane(const char *name, const double *P, const double *N,
		      const double x, const double y, const double rf,
		      const double gf, const double bf, const double rb,
		      const double gb, const double bb, const double dz);

#endif				/* __OFF_H__ */

/*	off.h
 *
 * Copyright (C) 2010, 2011 Ivo Alxneit
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
extern void off_cone(const char *name, double *O, double *dir,
		     const double l, const double r, const double g,
		     const double b);
extern void off_axes(const char *name, const double *origin,
		     const double *X, const double *y, const double *Z);
extern void off_plane(const char *name, const double *P, const double *N,
		      const double x, const double y, const double rf,
		      const double gf, const double bf, const double rb,
		      const double gb, const double bb, const double dz);
extern void off_rectangle(const char *name, const double *P,
			  const double *X, const double *Y,
			  const double rf, const double gf,
			  const double bf, const double rb,
			  const double gb, const double bb,
			  const double dz);
extern void off_triangle(const char *name, const double *P1,
			 const double *P2, const double *P3,
			 const double *N, const double rf, const double gf,
			 const double bf, const double rb, const double gb,
			 const double bb, const double dz);
extern void off_ellipsoid(const char *name, const double *origin,
			  const double *Z, const double *axes,
			  const double z_min, const double z_max,
			  const double ri, const double gi,
			  const double bi, const double ro,
			  const double go, const double bo,
			  const double dz);

extern void g2l_off(const double *P, const double *N, double *L,
		    double *alpha, double *beta);
extern void l2g_off(const double *P, const double *L, double *G,
		    const double alpha, const double beta);


#endif				/* __OFF_H__ */

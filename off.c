/*	off.c
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
#define _GNU_SOURCE		/* for sincos() */

#include <math.h>
#include <string.h>

#include <gsl/gsl_math.h>

#include "io_utils.h"
#include "math_utils.h"
#include "off.h"

#define INSIDE 0		/* reflecting surface of non-planar targets */
#define OUTSIDE 1

#define DZ 0.005		/* offset between inside/outside or front/back surface */

#define RAY_DIAMETER 0.05	/* diameter of rays logged to OFF file */

#define AXES_LENGTH 5.0
#define AXES_WIDTH 0.05

#define S_SCREEN 10.0		/* size of screens */

#define N_TRANS 11		/* rotationally symmertric targets */
#define N_ROT 24

#define R_REFL 1.0		/* r g b values of reflecting surface (targets) */
#define G_REFL 1.0
#define B_REFL 1.0

#define R_ABS 0.0		/* r g b values of absorbing surface (targets) */
#define G_ABS 0.0
#define B_ABS 0.0

#define R_SCREEN 1.0		/* r g b values of screens (targets) */
#define G_SCREEN 0.0
#define B_SCREEN 0.0

#define R_TRNSP 0.390		/* r g b values of transparent targets */
#define G_TRNSP 0.785
#define B_TRNSP 0.785

#define R_SRC 1.0		/* r g b value of (yellow) sources */
#define G_SRC 1.0
#define B_SRC 0.0

#define SPOT_LENGTH 1.2		/* length of cone for spot source */
#define SPOT_RADIUS 0.3		/* radius of cone for spot source */

#define R_POINT_SRC 0.2		/* radius of point source */


/*
 * macros and functions to write vertices
 */
#define WRITE_VERTEX(_F, _X, _Y, _Z, _P, _ALPHA, _BETA) do { \
        double _L[3], _G[3]; \
        _L[0] = _X; \
        _L[1] = _Y; \
        _L[2] = _Z; \
        l2g_off(_P, _L, _G, _ALPHA, _BETA); \
        fprintf(_F, "%12.6g\t%12.6g\t%12.6g\n", _G[0], _G[1], _G[2]); \
} while(0);

static void write_ring_vertices(FILE * outf, const double l2,
				const double r, const double *origin,
				const double alpha, const double beta)
{
    int j;
    double l[3];
    const double delta_phi = 2.0 * M_PI / N_ROT;

    l[2] = l2;
    for (j = 0; j < N_ROT; j++) {
	const double phi = j * delta_phi;
	double g[3];

	l[0] = r * sin(phi);
	l[1] = r * cos(phi);

	l2g_off(origin, l, g, alpha, beta);

	fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", g[0], g[1], g[2]);
    }
}

static void write_ann_vertices(FILE * outf, const double R, const double r,
			       const double *origin, const double alpha,
			       const double beta, const double offset)
{
    const double delta_phi = 2.0 * M_PI / N_ROT;
    double P[3], g_P[3];
    int i;

    P[2] = offset;
    for (i = 0; i < N_ROT; i++) {
	const double arg = i * delta_phi;
	const double sa = sin(arg);
	const double ca = cos(arg);
	P[0] = R * sa;
	P[1] = R * ca;
	l2g_off(origin, P, g_P, alpha, beta);
	fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", g_P[0], g_P[1], g_P[2]);
	P[0] = r * sa;
	P[1] = r * ca;
	l2g_off(origin, P, g_P, alpha, beta);
	fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", g_P[0], g_P[1], g_P[2]);
    }
}

static void block_vertices(FILE * f, const double *P, const double *N,
			   const double x, const double y)
/*
 * writes the vertices of a block of length 'N' - 'P'
 * and diameter 'x' times 'y' to file 'f'
 */
{
    double L[3];
    double alpha, beta;
    const double x2 = x / 2.0;
    const double y2 = y / 2.0;

/*
 * determine alpha, beta for transformation from local to global system.
 * length is z component of 'L'.
 */
    g2l_off(P, N, L, &alpha, &beta);

    WRITE_VERTEX(f, x2, y2, 0.0, P, alpha, beta);
    WRITE_VERTEX(f, x2, -y2, 0.0, P, alpha, beta);
    WRITE_VERTEX(f, -x2, -y2, 0.0, P, alpha, beta);
    WRITE_VERTEX(f, -x2, y2, 0.0, P, alpha, beta);
    WRITE_VERTEX(f, x2, y2, L[2], P, alpha, beta);
    WRITE_VERTEX(f, x2, -y2, L[2], P, alpha, beta);
    WRITE_VERTEX(f, -x2, -y2, L[2], P, alpha, beta);
    WRITE_VERTEX(f, -x2, y2, L[2], P, alpha, beta);
}

static void write_ell_vertices(FILE * outf, const double s,
			       const double z_min, const double z_max,
			       const double *axes, const double *origin,
			       const double alpha, const double beta)
{
    const double delta_z = (z_max - z_min) / (N_TRANS - 1);
    int i;

    for (i = 0; i < N_TRANS; i++) {
	double l = z_min + i * delta_z;
	double arg = 1.0 - l * l / (axes[2] * axes[2]);
	double r;

	if (arg <= 0)
	    r = 0;
	else
	    r = s * axes[0] * sqrt(arg);

	write_ring_vertices(outf, l, r, origin, alpha, beta);
    }
}

static void write_par_vertices(FILE * outf, const double s,
			       const double z_min, const double z_max,
			       const double foc, const double *origin,
			       const double alpha, const double beta)
{
    const double delta_z = (z_max - z_min) / (N_TRANS - 1);
    int i;

    for (i = 0; i < N_TRANS; i++) {
	double l = z_min + i * delta_z;
	double r = s * sqrt(fabs(4.0 * foc * l));

	write_ring_vertices(outf, l, r, origin, alpha, beta);
    }
}


/*
 * macros and functions to write faces
 */
#define WRITE_ANN_FACE(_F, _OFFSET, _R, _G, _B) do { \
    int _I, _BASE; \
    for (_I = _OFFSET; _I < _OFFSET + N_ROT - 1; _I++) { \
	_BASE = 2 * _I; \
	fprintf(_F, "5 %d %d %d %d %d %12.6g\t%12.6g\t%12.6g 1.0\n", \
	              _BASE, _BASE + 1, _BASE + 3, _BASE + 2, _BASE, \
	              _R, _G, _B); \
    } \
    _BASE += 2; \
    fprintf(_F, "5 %d %d %d %d %d %12.6g\t%12.6g\t%12.6g 1.0\n", \
                   _BASE, _BASE + 1, 2 * _OFFSET + 1, 2 * _OFFSET, _BASE, \
                   _R, _G, _B); \
} while(0);

#define WRITE_C_WALL(_F, _OFFSET, _R, _G, _B) do { \
    int _I; \
    for (_I = _OFFSET; _I < N_ROT - 1 + _OFFSET; _I++) \
	fprintf(_F, "5 %d %d %d %d %d %12.6g %12.6g %12.6g 1.0\n", \
	        _I, _I + 1, _I + N_ROT + 1, _I + N_ROT, _I, _R, _G, _B); \
    fprintf(_F, "5 %d %d %d %d %d %12.6g %12.6g %12.6g 1.0\n", \
            _OFFSET + N_ROT - 1, _OFFSET, _OFFSET + N_ROT, \
            _OFFSET + 2 * N_ROT - 1, _OFFSET + N_ROT - 1, _R, _G, _B); \
} while(0);

#define WRITE_W_FACE(_F, _OFFSET, _R, _G, _B) do { \
    int _I; \
    fprintf(_F, "%d ", N_ROT + 1); \
    for (_I = _OFFSET; _I < N_ROT + _OFFSET; _I++) \
	fprintf(_F, "%d ", _I); \
    fprintf(_F, "%d %12.6g\t%12.6g\t%12.6g\t 1.0\n", _OFFSET, _R, _G, _B); \
} while(0);

static void block_faces(FILE * f, const int i, const double r,
			const double g, const double b)
/*
 * prints the faces of a block (in OFF format) to 'f'.
 * 'i' denotes offset (index) of first vertex of block.
 * NOTE: vertices have to be written before to OFF file
 *       in the correct order! (use 'block_vertices()'
 */
{
    fprintf(f, "5 %d %d %d %d %d\t%12.6g %12.6g %12.6g 1.0\n", i, i + 1,
	    i + 2, i + 3, i, r, g, b);
    fprintf(f, "5 %d %d %d %d %d\t%12.6g %12.6g %12.6g 1.0\n", i, i + 1,
	    i + 5, i + 4, i, r, g, b);
    fprintf(f, "5 %d %d %d %d %d\t%12.6g %12.6g %12.6g 1.0\n", i, i + 3,
	    i + 7, i + 4, i, r, g, b);
    fprintf(f, "5 %d %d %d %d %d\t%12.6g %12.6g %12.6g 1.0\n", i + 2,
	    i + 1, i + 5, i + 6, i + 2, r, g, b);
    fprintf(f, "5 %d %d %d %d %d\t%12.6g %12.6g %12.6g 1.0\n", i + 2,
	    i + 3, i + 7, i + 6, i + 2, r, g, b);
    fprintf(f, "5 %d %d %d %d %d\t%12.6g %12.6g %12.6g 1.0\n", i + 4,
	    i + 5, i + 6, i + 7, i + 4, r, g, b);
}

static void write_ell_faces(FILE * outf, const int offset,
			    const double R, const double G, const double B)
{
    int i;

    for (i = 0; i < (N_TRANS - 1); i++)
	WRITE_C_WALL(outf, offset + i * N_ROT, R, G, B);
}


/*
 * individual off_XXX() functions for axes, sources and targets
 */
static void off_axes(const char *name, const double *origin,
		     const double *X, const double *Y, const double *Z)
{
/*
 * draw x,y,z (r,g,b) axes of global coordinate system as colored bars
 */
    double P[3];
    int i;
    FILE *outf;

    if (!name)			/* axes of global system */
	outf = open_off("axes");
    else {			/* axes of local system */
	char f_name[256];

	snprintf(f_name, 256, "axes_%s", name);
	outf = open_off(f_name);
    }

    fprintf(outf, "OFF\n");
    fprintf(outf, "24 18 0\n");

    /* output vertices of axes */
    for (i = 0; i < 3; i++)
	P[i] = origin[i] + X[i];
    block_vertices(outf, origin, P, AXES_WIDTH, AXES_WIDTH);
    for (i = 0; i < 3; i++)
	P[i] = origin[i] + Y[i];
    block_vertices(outf, origin, P, AXES_WIDTH, AXES_WIDTH);
    for (i = 0; i < 3; i++)
	P[i] = origin[i] + Z[i];
    block_vertices(outf, origin, P, AXES_WIDTH, AXES_WIDTH);

    /* output faces of axes */
    block_faces(outf, 0, 1.0, 0.0, 0.0);
    block_faces(outf, 8, 0.0, 1.0, 0.0);
    block_faces(outf, 16, 0.0, 0.0, 1.0);

    fclose(outf);
}

static void off_sphere_src(const char *name, double *O,
			   const double radius)
{
    int i;
    double axes[3];
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "%d %d 0\n", N_TRANS * N_ROT, (N_TRANS - 1) * N_ROT);

    for (i = 0; i < 3; i++)
	axes[i] = radius;

    write_ell_vertices(outf, 1.0, -radius, radius, axes, O, 0.0, 0.0);
    write_ell_faces(outf, 0, R_SRC, G_SRC, B_SRC);

    fclose(outf);
}

static void off_spot_src(const char *name, double *origin, double *dir)
{
    double P[3];
    double alpha, beta;
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "%d %d 0\n", 2 * N_ROT, 1 * N_ROT + 1);

    g2l_off_rot(dir, P, &alpha, &beta);

    write_ring_vertices(outf, SPOT_LENGTH, SPOT_RADIUS, origin, alpha,
			beta);
    write_ring_vertices(outf, 0.0, 0.0, origin, alpha, beta);

    WRITE_C_WALL(outf, 0, R_SRC, G_SRC, B_SRC);
    WRITE_W_FACE(outf, 0, R_SRC, G_SRC, B_SRC);

    fclose(outf);
}


static void off_annulus(const char *name, const double *origin,
			const double *dir, const double R, const double r)
/*
 * writes 'OFF' file to 'name.off' for annulus defined by its center 'origin',
 * normal vector 'dir' and radii 'R' and 'r'.
 */
{
    double P[3];
    double alpha, beta;
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "%d %d 0\n", 2 * 2 * N_ROT, 2 * N_ROT);

    g2l_off_rot(dir, P, &alpha, &beta);

    write_ann_vertices(outf, R, r, origin, alpha, beta, 0.0);
    write_ann_vertices(outf, R, r, origin, alpha, beta, DZ);

    WRITE_ANN_FACE(outf, 0, R_ABS, G_ABS, B_ABS);
    WRITE_ANN_FACE(outf, N_ROT, R_REFL, G_REFL, B_REFL);

    fclose(outf);
}

static void off_cone(const char *name, const double *origin,
		     const double *dir, const double h, double Radius,
		     double radius, const int r_surface)
{
/*
 * writes 'OFF' file to 'name.off' for cone with radius of base disk 'Radius',
 * radius of top disk 'radius', heigth (between the two disks) 'h'.
 */
    double P[3];
    double alpha, beta;
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "%d %d 0\n", 4 * N_ROT, 2 * N_ROT);

    /*
     * determine alpha, beta for transformation from local to global system.
     * discard 'P'
     */
    g2l_off_rot(dir, P, &alpha, &beta);

    /*
     * write vertices: outside/bottom, outside/top, inside/bottom, inside/top
     */
    write_ring_vertices(outf, 0.0, Radius, origin, alpha, beta);
    write_ring_vertices(outf, h, radius, origin, alpha, beta);
    write_ring_vertices(outf, 0.0, Radius * (1 - DZ), origin, alpha, beta);
    write_ring_vertices(outf, h, radius * (1 - DZ), origin, alpha, beta);

    if (r_surface == OUTSIDE) {
	WRITE_C_WALL(outf, 0, R_REFL, G_REFL, B_REFL);	/* outside */
	WRITE_C_WALL(outf, 2 * N_ROT, R_ABS, R_ABS, R_ABS);	/* inside */
    } else {
	WRITE_C_WALL(outf, 0, R_ABS, R_ABS, R_ABS);	/* outside */
	WRITE_C_WALL(outf, 2 * N_ROT, R_REFL, G_REFL, B_REFL);	/* inside */
    }

    fclose(outf);
}

static void off_ellipsoid(const char *name, const double *origin,
			  const double *Z, const double *axes,
			  const double z_min, const double
			  z_max, const int r_surface)
{
    double alpha, beta;
    double l[3];
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "%d %d 0\n", 2 * N_TRANS * N_ROT,
	    2 * (N_TRANS - 1) * N_ROT);

    g2l_off_rot(Z, l, &alpha, &beta);

    write_ell_vertices(outf, 1.0, z_min, z_max, axes, origin, alpha, beta);
    write_ell_vertices(outf, 1.0 - DZ, z_min, z_max, axes, origin, alpha,
		       beta);

    if (r_surface == OUTSIDE) {
	write_ell_faces(outf, 0, R_REFL, G_REFL, B_REFL);	/* outside */
	write_ell_faces(outf, N_TRANS * N_ROT, R_ABS, R_ABS, R_ABS);	/* inside */
    } else {
	write_ell_faces(outf, 0, R_ABS, R_ABS, R_ABS);	/* outside */
	write_ell_faces(outf, N_TRANS * N_ROT, R_REFL, G_REFL, B_REFL);	/* inside */
    }

    fclose(outf);
}

static void off_paraboloid(const char *name, const double *origin,
			   const double *Z, const double foc,
			   const double z_min, const double z_max,
			   const int r_surface)
{
    double alpha, beta;
    double l[3];
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "%d %d 0\n", 2 * N_TRANS * N_ROT,
	    2 * (N_TRANS - 1) * N_ROT);

    g2l_off_rot(Z, l, &alpha, &beta);

    write_par_vertices(outf, 1.0, z_min, z_max, foc, origin, alpha, beta);
    write_par_vertices(outf, 1.0 - DZ, z_min, z_max, foc, origin, alpha,
		       beta);

    if (r_surface == OUTSIDE) {
	write_ell_faces(outf, 0, R_REFL, G_REFL, B_REFL);	/* outside */
	write_ell_faces(outf, N_TRANS * N_ROT, R_ABS, R_ABS, R_ABS);	/* inside */
    } else {
	write_ell_faces(outf, 0, R_ABS, R_ABS, R_ABS);	/* outside */
	write_ell_faces(outf, N_TRANS * N_ROT, R_REFL, G_REFL, B_REFL);	/* inside */
    }

    fclose(outf);
}

static void off_rectangle(const char *name, const double *P,
			  const double *X, const double *Y)
{
    double N[3];
    FILE *outf = open_off(name);

    cross_product(X, Y, N);
    normalize(N);

    fprintf(outf, "OFF\n");
    fprintf(outf, "8 2 0\n");	/* 8 vertices, 2 faces */

    /* back side */
    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", P[0], P[1], P[2]);
    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", P[0] + X[0], P[1] + X[1],
	    P[2] + X[2]);
    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", P[0] + X[0] + Y[0],
	    P[1] + X[1] + Y[1], P[2] + X[2] + Y[2]);
    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", P[0] + Y[0], P[1] + Y[1],
	    P[2] + Y[2]);

    /* front side: offset by 'DZ' times 'N' */
    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", P[0] + DZ * N[0],
	    P[1] + DZ * N[1], P[2] + DZ * N[2]);
    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", P[0] + DZ * N[0] + X[0],
	    P[1] + DZ * N[1] + X[1], P[2] + DZ * N[2] + X[2]);
    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n",
	    P[0] + DZ * N[0] + X[0] + Y[0], P[1] + DZ * N[1] + X[1] + Y[1],
	    P[2] + DZ * N[2] + X[2] + Y[2]);
    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", P[0] + DZ * N[0] + Y[0],
	    P[1] + DZ * N[1] + Y[1], P[2] + DZ * N[2] + Y[2]);

    /*
     * output faces (back, front)
     */
    fprintf(outf, "5 0 1 2 3 0 %12.6g\t%12.6g\t%12.6g 1.0\n", R_ABS, G_ABS,
	    B_ABS);
    fprintf(outf, "5 4 5 6 7 4 %12.6g\t%12.6g\t%12.6g 1.0\n", R_REFL,
	    G_REFL, B_REFL);

    fclose(outf);

}

static void off_screen(const char *name, const double *P, const double *N,
		       const double rf, const double gf, const double bf,
		       const double rb, const double gb, const double bb)
/*
 * writes 'OFF' file to 'name.off' for screens defined
 * by the point 'P' and its normal vector 'N'.
 */
{
    int i;
    double P2[3];
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "8 2 0\n");	/* 8 vertices, 2 faces */

    for (i = 0; i < 3; i++)	/* 'P2' is point on front side */
	P2[i] = P[i] + DZ * N[i];

    block_vertices(outf, P, P2, S_SCREEN, S_SCREEN);

    /*
     * print back face ('rb', 'gb', 'bb'), front face ('rf', 'gf', 'bf') 
     */
    fprintf(outf, "5 0 1 2 3 0 %12.6g\t%12.6g\t%12.6g 1.0\n", rb, gb, bb);
    fprintf(outf, "5 4 5 6 7 4 %12.6g\t%12.6g\t%12.6g 1.0\n", rf, gf, bf);

    fclose(outf);
}

static void off_triangle(const char *name, const double *P1,
			 const double *P2, const double *P3,
			 const double *N)
/*
 * writes 'OFF' file to 'name.off' for triangle defined
 * by the point 'P1' and the two sides 'P2' and 'P3'. its normal vector is 'N'.
 */
{
    const double dx = N[0] * DZ;
    const double dy = N[1] * DZ;
    const double dz = N[2] * DZ;
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "6 2 0\n");	/* 6 vertices, 2 faces */

    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", P1[0], P1[1], P1[2]);
    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", P2[0] + P1[0], P2[1] + P1[1],
	    P2[2] + P1[2]);
    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", P3[0] + P1[0], P3[1] + P1[1],
	    P3[2] + P1[2]);
    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", P1[0] + dx, P1[1] + dy,
	    P1[2] + dz);
    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", P2[0] + P1[0] + dx,
	    P2[1] + P1[1] + dy, P2[2] + P1[2] + dz);
    fprintf(outf, "%12.6g\t%12.6g\t%12.6g\n", P3[0] + P1[0] + dx,
	    P3[1] + P1[1] + dy, P3[2] + P1[2] + dz);

    /*
     * output faces (back, front)
     */
    fprintf(outf, "4 0 1 2 0 %12.6g\t%12.6g\t%12.6g 1.0\n", R_ABS, G_ABS,
	    B_ABS);
    fprintf(outf, "4 3 4 5 0 %12.6g\t%12.6g\t%12.6g 1.0\n", R_REFL, G_REFL,
	    B_REFL);

    fclose(outf);
}

static void off_window(const char *name, const double *origin,
		       const double *dir, const double R, const double d)
{
/*
 * writes 'OFF' file to 'name.off' for window. the window is defined by the
 * center 'origin' of one face, the direction vector 'dir' of the cylinder,
 * its radius 'r', and its thickness.
 */
    double P[3];
    double alpha, beta;
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "%d %d 0\n", 2 * N_ROT, N_ROT + 2);
    /* 2 faces with 1 patch each plus cylinder wall with 'N_ROT' patches */

    /*
     * determine alpha, beta for transformation from local to global system.
     * discard 'P'
     */
    g2l_off_rot(dir, P, &alpha, &beta);

    write_ring_vertices(outf, 0.0, R, origin, alpha, beta);
    write_ring_vertices(outf, d, R, origin, alpha, beta);

    WRITE_W_FACE(outf, 0, R_TRNSP, G_TRNSP, B_TRNSP);
    WRITE_W_FACE(outf, N_ROT, R_TRNSP, G_TRNSP, B_TRNSP);

    WRITE_C_WALL(outf, 0, R_ABS, G_ABS, B_ABS);

    fclose(outf);
}



static void output_sources(const config_t * cfg)
{
    int i;
    const config_setting_t *s = config_lookup(cfg, "sources");
    const int n_sources = config_setting_length(s);

    for (i = 0; i < n_sources; ++i) {	/* iterate through all sources */
	const char *type, *name;
	const config_setting_t *this_s =
	    config_setting_get_elem(s, (unsigned int) i);

	config_setting_lookup_string(this_s, "name", &name);
	config_setting_lookup_string(this_s, "type", &type);

	if (!strcmp(type, "arc")) {
	    double origin[3];
	    double dir[3];
	    double length;
	    double radius;
	    double P[3];
	    double alpha, beta;

	    FILE *outf = open_off(name);

	    read_vector(this_s, "origin", origin);
	    read_vector(this_s, "direction", dir);
	    config_setting_lookup_float(this_s, "radius", &radius);
	    config_setting_lookup_float(this_s, "length", &length);

	    fprintf(outf, "OFF\n");
	    fprintf(outf, "%d %d 0\n", 2 * N_ROT, N_ROT + 2);
	    /* 2 faces with 1 patch each plus cylinder wall with 'N_ROT' patches */

	    /*
	     * determine alpha, beta for transformation from local to global system.
	     * discard 'P'
	     */
	    g2l_off_rot(dir, P, &alpha, &beta);

	    write_ring_vertices(outf, 0.0, radius, origin, alpha, beta);
	    write_ring_vertices(outf, length, radius, origin, alpha, beta);

	    WRITE_W_FACE(outf, 0, R_SRC, G_SRC, B_SRC);
	    WRITE_W_FACE(outf, N_ROT, R_SRC, G_SRC, B_SRC);
	    WRITE_C_WALL(outf, 0, R_SRC, G_SRC, B_SRC);

	    fclose(outf);

	} /* end 'arc' */
	else if (!strcmp(type, "solid cone")) {
	    double origin[3];
	    double z[3];
	    double h;
	    double R, r;
	    double P[3];
	    double alpha, beta;
	    int ans;

	    FILE *outf = open_off(name);

	    read_vector(this_s, "origin", origin);
	    read_vector(this_s, "z", z);
	    config_setting_lookup_float(this_s, "R", &R);
	    config_setting_lookup_float(this_s, "r", &r);
	    config_setting_lookup_float(this_s, "h", &h);

	    fprintf(outf, "OFF\n");
	    fprintf(outf, "%d %d 0\n", 2 * N_ROT, N_ROT + 2);
	    /* 2 faces with 1 patch each plus cylinder wall with 'N_ROT' patches */

	    /*
	     * determine alpha, beta for transformation from local to global system.
	     * discard 'P'
	     */
	    g2l_off_rot(z, P, &alpha, &beta);

	    write_ring_vertices(outf, 0.0, R, origin, alpha, beta);
	    write_ring_vertices(outf, h, r, origin, alpha, beta);

	    config_setting_lookup_bool(this_s, "base_face_emits", &ans);
	    if (!ans) {
		WRITE_W_FACE(outf, 0, R_ABS, G_ABS, B_ABS);
	    } else {
		WRITE_W_FACE(outf, 0, R_SRC, G_SRC, B_SRC);
	    }

	    config_setting_lookup_bool(this_s, "top_face_emits", &ans);
	    if (!ans) {
		WRITE_W_FACE(outf, N_ROT, R_ABS, G_ABS, B_ABS);
	    } else {
		WRITE_W_FACE(outf, N_ROT, R_SRC, G_SRC, B_SRC);
	    }

	    WRITE_C_WALL(outf, 0, R_SRC, G_SRC, B_SRC);

	    fclose(outf);

	} /* end 'solid cone' */
	else if (!strcmp(type, "solid cylinder")) {
	    double origin[3];
	    double dir[3];
	    double length;
	    double radius;
	    double P[3];
	    double alpha, beta;
	    int ans;

	    FILE *outf = open_off(name);

	    read_vector(this_s, "origin", origin);
	    read_vector(this_s, "direction", dir);
	    config_setting_lookup_float(this_s, "radius", &radius);
	    config_setting_lookup_float(this_s, "length", &length);

	    fprintf(outf, "OFF\n");
	    fprintf(outf, "%d %d 0\n", 2 * N_ROT, N_ROT + 2);
	    /* 2 faces with 1 patch each plus cylinder wall with 'N_ROT' patches */

	    /*
	     * determine alpha, beta for transformation from local to global system.
	     * discard 'P'
	     */
	    g2l_off_rot(dir, P, &alpha, &beta);

	    write_ring_vertices(outf, 0.0, radius, origin, alpha, beta);
	    write_ring_vertices(outf, length, radius, origin, alpha, beta);

	    config_setting_lookup_bool(this_s, "base_face_emits", &ans);
	    if (!ans) {
		WRITE_W_FACE(outf, 0, R_ABS, G_ABS, B_ABS);
	    } else {
		WRITE_W_FACE(outf, 0, R_SRC, G_SRC, B_SRC);
	    }

	    config_setting_lookup_bool(this_s, "top_face_emits", &ans);
	    if (!ans) {
		WRITE_W_FACE(outf, N_ROT, R_ABS, G_ABS, B_ABS);
	    } else {
		WRITE_W_FACE(outf, N_ROT, R_SRC, G_SRC, B_SRC);
	    }

	    WRITE_C_WALL(outf, 0, R_SRC, G_SRC, B_SRC);

	    fclose(outf);

	} /* end 'solid cylinder' */
	else if (!strcmp(type, "sphere") || !strcmp(type, "solid sphere")) {
	    double O[3];
	    double radius;

	    read_vector(this_s, "origin", O);
	    config_setting_lookup_float(this_s, "radius", &radius);

	    off_sphere_src(name, O, radius);
	} /* end 'sphere' or 'solid_sphere' */
	else if (!strcmp(type, "spot source")) {
	    double O[3], dir[3];

	    read_vector(this_s, "origin", O);
	    read_vector(this_s, "direction", dir);

	    off_spot_src(name, O, dir);

	} /* end 'spot source' */
	else if (!strcmp(type, "uniform point source")) {
	    double O[3];

	    read_vector(this_s, "origin", O);

	    off_sphere_src(name, O, R_POINT_SRC);
	}			/* end 'uniform point source' */
    }				/* end all sources */
}

static void output_targets(const config_t * cfg)
{
    int i;
    const config_setting_t *t = config_lookup(cfg, "targets");
    const int n_targets = config_setting_length(t);

    for (i = 0; i < n_targets; ++i) {	/* iterate through all targets */
	const char *type, *name;
	const config_setting_t *this_t =
	    config_setting_get_elem(t, (unsigned int) i);

	config_setting_lookup_string(this_t, "name", &name);
	config_setting_lookup_string(this_t, "type", &type);

	if (!strcmp(type, "annulus")) {
	    double O[3], X[3], Y[3], Z[3];
	    double R, r;

	    read_vector(this_t, "P", O);
	    read_vector_normalize(this_t, "x", X);
	    read_vector_normalize(this_t, "N", Z);
	    config_setting_lookup_float(this_t, "R", &R);
	    config_setting_lookup_float(this_t, "r", &r);

	    orthonormalize(X, Y, Z);

	    off_axes(name, O, X, Y, Z);	/*local system */

	    off_annulus(name, O, Z, R, r);

	} else if (!strcmp(type, "cone")) {
	    double O[3], X[3], Y[3], Z[3];
	    double R, r;
	    double h;
	    const char *S;

	    read_vector(this_t, "origin", O);
	    read_vector_normalize(this_t, "x", X);
	    read_vector_normalize(this_t, "axis", Z);
	    orthonormalize(X, Y, Z);

	    off_axes(name, O, X, Y, Z);	/*local system */

	    config_setting_lookup_float(this_t, "R", &R);
	    config_setting_lookup_float(this_t, "r", &r);
	    config_setting_lookup_float(this_t, "h", &h);

	    config_setting_lookup_string(this_t, "reflecting_surface", &S);
	    if (!strcmp(S, "inside"))
		off_cone(name, O, Z, h, R, r, INSIDE);
	    else
		off_cone(name, O, Z, h, R, r, OUTSIDE);

	} else if (!strcmp(type, "cpc")) {
	    double O[3], X[3], Y[3], Z[3];
	    double R, r;
	    double phi_a, theta_t;
	    double h;
	    double foc2, term;

	    read_vector(this_t, "origin", O);
	    read_vector_normalize(this_t, "x", X);
	    read_vector_normalize(this_t, "axis", Z);
	    orthonormalize(X, Y, Z);

	    off_axes(name, O, X, Y, Z);	/*local system */

	    config_setting_lookup_float(this_t, "acceptance_angle",
					&phi_a);
	    phi_a *= M_PI / 180.0;
	    config_setting_lookup_float(this_t, "truncation_angle",
					&theta_t);
	    theta_t *= M_PI / 180.0;
	    config_setting_lookup_float(this_t, "exit_radius", &r);

	    foc2 = 2.0 * r * (1.0 + sin(phi_a));
	    term = foc2 / (1 - cos(theta_t + phi_a));

	    R = sin(theta_t) * term - r;
	    h = cos(theta_t) * term;

	    /*
	     * FIXME: use better surface such as paraboloid or ellipsoid
	     */
	    off_cone(name, O, Z, h, r, R, INSIDE);

	} else if (!strcmp(type, "cylinder")) {
	    double O[3], X[3], Y[3], Z[3];
	    double r, l;
	    const char *S;

	    read_vector(this_t, "C", O);
	    read_vector_normalize(this_t, "x", X);
	    read_vector_normalize(this_t, "a", Z);

	    orthonormalize(X, Y, Z);

	    off_axes(name, O, X, Y, Z);	/*local system */

	    config_setting_lookup_float(this_t, "r", &r);
	    config_setting_lookup_float(this_t, "l", &l);

	    config_setting_lookup_string(this_t, "reflecting_surface", &S);
	    if (!strcmp(S, "inside"))
		off_cone(name, O, Z, l, r, r, INSIDE);
	    else
		off_cone(name, O, Z, l, r, r, OUTSIDE);

	} else if (!strcmp(type, "disk")) {
	    double O[3], X[3], Y[3], Z[3];
	    double r;

	    read_vector(this_t, "P", O);
	    read_vector_normalize(this_t, "x", X);
	    read_vector_normalize(this_t, "N", Z);
	    config_setting_lookup_float(this_t, "r", &r);

	    orthonormalize(X, Y, Z);

	    off_axes(name, O, X, Y, Z);	/*local system */

	    off_annulus(name, O, Z, r, 0.0);

	} else if (!strcmp(type, "ellipsoid")) {
	    double O[3], X[3], Y[3], Z[3];
	    double axes[3];
	    double z_min, z_max;
	    const char *S;

	    read_vector(this_t, "center", O);
	    read_vector_normalize(this_t, "x", X);
	    read_vector_normalize(this_t, "z", Z);
	    read_vector(this_t, "axes", axes);

	    config_setting_lookup_float(this_t, "z_min", &z_min);
	    config_setting_lookup_float(this_t, "z_max", &z_max);

	    orthonormalize(X, Y, Z);

	    off_axes(name, O, X, Y, Z);	/*local system */

	    config_setting_lookup_string(this_t, "reflecting_surface", &S);
	    if (!strcmp(S, "inside"))
		off_ellipsoid(name, O, Z, axes, z_min, z_max, INSIDE);
	    else
		off_ellipsoid(name, O, Z, axes, z_min, z_max, OUTSIDE);

	} else if (!strcmp(type, "paraboloid")) {
	    double vertex[3], X[3], Y[3], Z[3];
	    double foc;
	    double z_min, z_max;
	    const char *S;

	    read_vector(this_t, "vertex", vertex);
	    read_vector_normalize(this_t, "x", X);
	    read_vector_normalize(this_t, "z", Z);
	    orthonormalize(X, Y, Z);

	    off_axes(name, vertex, X, Y, Z);	/*local system */

	    config_setting_lookup_float(this_t, "focal_length", &foc);

	    config_setting_lookup_float(this_t, "z_min", &z_min);
	    config_setting_lookup_float(this_t, "z_max", &z_max);

	    if (z_max < z_min)	/* safety */
		SWAP(z_max, z_min);

	    config_setting_lookup_string(this_t, "reflecting_surface", &S);
	    if (!strcmp(S, "inside"))
		off_paraboloid(name, vertex, Z, foc, z_min, z_max, INSIDE);
	    else
		off_paraboloid(name, vertex, Z, foc, z_min, z_max,
			       OUTSIDE);

	} else if (!strcmp(type, "one-sided plane screen")) {
	    double P[3], N[3];
	    double X[3], Y[3];

	    read_vector(this_t, "point", P);
	    read_vector_normalize(this_t, "normal", N);
	    read_vector_normalize(this_t, "x", X);

	    cross_product(N, X, Y);

	    off_axes(name, P, X, Y, N);	/* local system */

	    off_screen(name, P, N, R_SCREEN, G_SCREEN, B_SCREEN, R_ABS,
		       G_ABS, B_ABS);

	} else if (!strcmp(type, "two-sided plane screen")) {
	    double P[3], N[3];
	    double X[3], Y[3];

	    read_vector(this_t, "point", P);
	    read_vector_normalize(this_t, "normal", N);
	    read_vector_normalize(this_t, "x", X);

	    cross_product(N, X, Y);

	    off_axes(name, P, X, Y, N);	/*local system */

	    off_screen(name, P, N, R_SCREEN, G_SCREEN, B_SCREEN, R_SCREEN,
		       G_SCREEN, B_SCREEN);

	} else if (!strcmp(type, "rectangle")) {
	    int j;
	    double P[3], N[3];
	    double X[3], Y[3];

	    read_vector(this_t, "P1", P);
	    read_vector(this_t, "P2", X);
	    read_vector(this_t, "P3", Y);
	    for (j = 0; j < 3; j++) {
		X[j] -= P[j];
		Y[j] -= P[j];
	    }

	    off_rectangle(name, P, X, Y);

	    /* make 'P1' point to center of rectangle */
	    for (j = 0; j < 3; j++)
		P[j] += (X[j] + Y[j]) / 2.0;

	    normalize(X);
	    normalize(Y);
	    /* surface normal N = X cross Y */
	    cross_product(X, Y, N);

	    off_axes(name, P, X, Y, N);	/*local system */

	} else if (!strcmp(type, "sphere")) {
	    double O[3], X[3], Y[3], Z[3];
	    double axes[3];
	    double radius;
	    double z_min, z_max;
	    const char *S;

	    read_vector(this_t, "origin", O);
	    read_vector_normalize(this_t, "x", X);
	    read_vector_normalize(this_t, "z", Z);
	    orthonormalize(X, Y, Z);

	    off_axes(name, O, X, Y, Z);	/*local system */

	    config_setting_lookup_float(this_t, "radius", &radius);

	    axes[0] = radius;
	    axes[1] = radius;
	    axes[2] = radius;

	    config_setting_lookup_float(this_t, "z_min", &z_min);
	    config_setting_lookup_float(this_t, "z_max", &z_max);

	    /*
	     * draw sphere as ellipsoid with all axes=radius^2
	     */
	    config_setting_lookup_string(this_t, "reflecting_surface", &S);
	    if (!strcmp(S, "inside"))
		off_ellipsoid(name, O, Z, axes, z_min, z_max, INSIDE);
	    else
		off_ellipsoid(name, O, Z, axes, z_min, z_max, OUTSIDE);

	} else if (!strcmp(type, "triangle")) {
	    int j;
	    double P1[3];
	    double P2[3];	/* redefined 'P2'='P2'-'P1' */
	    double P3[3];	/* redefined 'P3'='P3'-'P1' */
	    double N[3];	/* 'P2' cross 'P3' */
	    double Y[3];	/* local system. 'X' parallel 'P2' */
	    config_setting_t *this;

	    read_vector(this_t, "P1", P1);

	    this = config_setting_get_member(this_t, "P2");
	    for (j = 0; j < 3; j++)
		P2[j] = config_setting_get_float_elem(this, j) - P1[j];

	    this = config_setting_get_member(this_t, "P3");
	    for (j = 0; j < 3; j++)
		P3[j] = config_setting_get_float_elem(this, j) - P1[j];

	    cross_product(P2, P3, N);
	    normalize(N);

	    off_triangle(name, P1, P2, P3, N);

	    normalize(P2);
	    cross_product(N, P2, Y);

	    off_axes(name, P1, P2, Y, N);	/*local system */

	} else if (!strcmp(type, "window")) {
	    double O[3], X[3], Y[3], Z[3];
	    double d, r;

	    read_vector(this_t, "C", O);
	    read_vector_normalize(this_t, "x", X);
	    read_vector_normalize(this_t, "a", Z);
	    config_setting_lookup_float(this_t, "d", &d);
	    config_setting_lookup_float(this_t, "r", &r);

	    orthonormalize(X, Y, Z);

	    off_axes(name, O, X, Y, Z);	/*local system */

	    off_window(name, O, Z, r, d);

	}
    }				/* end all targets */
}



void output_geometry(config_t * cfg)
{
    const double O[] = { 0.0, 0.0, 0.0 };
    const double X[] = { AXES_LENGTH, 0.0, 0.0 };
    const double Y[] = { 0.0, AXES_LENGTH, 0.0 };
    const double Z[] = { 0.0, 0.0, AXES_LENGTH };

    off_axes(NULL, O, X, Y, Z);

    output_sources(cfg);
    output_targets(cfg);
}

FILE *open_off(const char *name)
{
    char f_name[256];

    snprintf(f_name, 256, "%s.off", name);
    return fopen(f_name, "w");
}

void write_ray(FILE * f, const double *start, const double *stop,
	       const int r, const int g, const int b)
{
    fprintf(f, "{ = OFF\n");
    fprintf(f, "8 6 0\n");

    block_vertices(f, start, stop, RAY_DIAMETER, RAY_DIAMETER);
    block_faces(f, 0, r, g, b);
    fprintf(f, "}\n");
}

void g2l_off(const double *P, const double *N, double *L,
	     double *alpha, double *beta)
/*
 * convert vector N from global coordinate system to local one by
 *
 * 1) translation by -P
 * 2a) determine beta
 * 2b) rotation around z axis into x-z plane (beta)
 * 3a) determine alpha
 * 3b) rotation around y axis onto z axis (alpha)
 *
 * return alpha, beta, and resulting vector in L
 */
{
    double norm;

    diff(L, N, P);

    norm = sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);

    *beta = atan2(L[1], L[0]);
    *alpha = atan2(L[0] * cos(-*beta) - L[1] * sin(-*beta), L[2]);

    /*
     * after translation, rotation L lies on z axis
     * and therefore becomes {0, 0, norm}
     */
    L[0] = 0.0;
    L[1] = 0.0;
    L[2] = norm;

}

void l2g_off(const double *P, const double *L, double *G,
	     const double alpha, const double beta)
/*
 * convert vector L from local coordinate system to global one by
 *
 * 1) rotation around y axis by -alpha
 * 2) rotation around z axis by -beta
 * 3) translation by P
 *
 * return alpha, beta, and resulting vector in G
 */
{
    double sa, ca, sb, cb;
    double x;

    sincos(-alpha, &sa, &ca);
    sincos(-beta, &sb, &cb);
    x = L[0] * ca - L[2] * sa;

    G[0] = x * cb - L[1] * sb + P[0];
    G[1] = -x * sb - L[1] * cb + P[1];
    G[2] = L[0] * sa + L[2] * ca + P[2];
}

void g2l_off_rot(const double *N, double *L, double *alpha, double *beta)
/*
 * convert vector N from global coordinate system to local
 * rotation only by
 *
 * 1a) determine beta
 * 1b) rotation around z axis into x-z plane (beta)
 * 2a) determine alpha
 * 2b) rotation around y axis onto z axis (alpha)
 *
 * return alpha, beta, and resulting vector in L
 */
{
    const double norm = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);

    *beta = atan2(N[1], N[0]);
    *alpha = atan2(N[0] * cos(-*beta) - N[1] * sin(-*beta), N[2]);

    /*
     * after translation, rotation L lies on z axis
     * and therefore becomes {0, 0, norm}
     */
    L[0] = 0.0;
    L[1] = 0.0;
    L[2] = norm;
}

void l2g_off_rot(const double *L, double *G, const double alpha,
		 const double beta)
/*
 * convert vector L from local coordinate system to global
 * rotation only by
 *
 * 1) rotation around y axis by -alpha
 * 2) rotation around z axis by -beta
 *
 * return alpha, beta, and resulting vector in G
 */
{
    double sa, ca, sb, cb;
    double x;

    sincos(-alpha, &sa, &ca);
    sincos(-beta, &sb, &cb);
    x = L[0] * ca - L[2] * sa;

    G[0] = x * cb - L[1] * sb;
    G[1] = -x * sb - L[1] * cb;
    G[2] = L[0] * sa + L[2] * ca;
}

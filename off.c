/*	off.c
 *
 * Copyright (C) 2010,2011,2012,2013 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>

#include "io_utils.h"
#include "math_utils.h"
#include "off.h"

#define DZ 0.005
#define AXES_LENGTH 5.0


static FILE *open_off(const char *name)
{
    char f_name[256];

    snprintf(f_name, 256, "%s.off", name);
    return fopen(f_name, "w");
}

static void block_vertices(FILE * f, const double *P, const double *N,
			   const double x, const double y)
/*
 * writes the vertices of a block of length 'N' - 'P'
 * and diameter 'x' times 'y' to file 'f'
 */
{
    double L[3], G[3];
    double alpha, beta;
    double l;
    const double x2 = x / 2.0;
    const double y2 = y / 2.0;

/*
 * determine alpha, beta for transformation from local to global system.
 * copy length from z component of 'L'.
 */
    g2l_off(P, N, L, &alpha, &beta);
    l = L[2];

/*
 * write the vertices
 *	+, +, 0
 *	+, -, 0
 *	-, -, 0
 *	-, +, 0
 *	+, +, l
 *	+, -, l
 *	-, -, l
 *	-, +, l
 *
 */
    L[0] = x2;
    L[1] = y2;
    L[2] = 0.0;
    l2g_off(P, L, G, alpha, beta);
    fprintf(f, "%f\t%f\t%f\n", G[0], G[1], G[2]);

    L[1] = -y2;
    l2g_off(P, L, G, alpha, beta);
    fprintf(f, "%f\t%f\t%f\n", G[0], G[1], G[2]);

    L[0] = -x2;
    l2g_off(P, L, G, alpha, beta);
    fprintf(f, "%f\t%f\t%f\n", G[0], G[1], G[2]);

    L[1] = y2;
    l2g_off(P, L, G, alpha, beta);
    fprintf(f, "%f\t%f\t%f\n", G[0], G[1], G[2]);

    L[0] = x2;
    L[2] = l;
    l2g_off(P, L, G, alpha, beta);
    fprintf(f, "%f\t%f\t%f\n", G[0], G[1], G[2]);

    L[1] = -y2;
    l2g_off(P, L, G, alpha, beta);
    fprintf(f, "%f\t%f\t%f\n", G[0], G[1], G[2]);

    L[0] = -x2;
    l2g_off(P, L, G, alpha, beta);
    fprintf(f, "%f\t%f\t%f\n", G[0], G[1], G[2]);

    L[1] = y2;
    l2g_off(P, L, G, alpha, beta);
    fprintf(f, "%f\t%f\t%f\n", G[0], G[1], G[2]);

}

static void block_faces(FILE * f, const int i, const double r,
			const double g, const double b)
/*
 * prints the faces of a block (in OFF format) to 'f'.
 * all faces are colored 'r','g','b'.
 * 'i' denotes offset (index) of first vertex of block.
 * NOTE: vertices have to be written before to OFF file
 *       in the correct order! (use 'block_vertices()'
 */
{
    fprintf(f, "5 %d %d %d %d %d\t%g %g %g 1.0\n", i, i + 1, i + 2, i + 3,
	    i, r, g, b);
    fprintf(f, "5 %d %d %d %d %d\t%g %g %g 1.0\n", i, i + 1, i + 5, i + 4,
	    i, r, g, b);
    fprintf(f, "5 %d %d %d %d %d\t%g %g %g 1.0\n", i, i + 3, i + 7, i + 4,
	    i, r, g, b);
    fprintf(f, "5 %d %d %d %d %d\t%g %g %g 1.0\n", i + 2, i + 1, i + 5,
	    i + 6, i + 2, r, g, b);
    fprintf(f, "5 %d %d %d %d %d\t%g %g %g 1.0\n", i + 2, i + 3, i + 7,
	    i + 6, i + 2, r, g, b);
    fprintf(f, "5 %d %d %d %d %d\t%g %g %g 1.0\n", i + 4, i + 5, i + 6,
	    i + 7, i + 4, r, g, b);
}



static void off_axes(const char *name, const double *origin,
		     const double *X, const double *Y, const double *Z)
{
/*
 * draw axes of global coordinate system as colored bars
 * x-axis (X): red
 * y-axis (Y): green
 * z-axis (Z): blue
 */
    const double s = 0.05;
    double P[3];
    int i;
    FILE *outf;

    if (!name)
	outf = open_off("axes");
    else {
	char f_name[256];

	snprintf(f_name, 256, "axes_%s", name);
	outf = open_off(f_name);
    }

    fprintf(outf, "OFF\n");
    fprintf(outf, "24 18 0\n");

    /* output vertices of axes */
    for (i = 0; i < 3; i++)
	P[i] = origin[i] + X[i];
    block_vertices(outf, origin, P, s, s);
    for (i = 0; i < 3; i++)
	P[i] = origin[i] + Y[i];
    block_vertices(outf, origin, P, s, s);
    for (i = 0; i < 3; i++)
	P[i] = origin[i] + Z[i];
    block_vertices(outf, origin, P, s, s);

    /* output faces of axes */
    block_faces(outf, 0, 1.0, 0.0, 0.0);
    block_faces(outf, 8, 0.0, 1.0, 0.0);
    block_faces(outf, 16, 0.0, 0.0, 1.0);

    fclose(outf);
}



static void off_sphere(const char *name, double *O, const double radius,
		       const double r, const double g, const double b)
{
/*
 * print octaeder as source
 */
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "6 8 0\n");	/* 6 vertices, 8 faces */

    /*
     * list of vertices
     * 0 ... 5
     */
    fprintf(outf, "%f\t%f\t%f\n", O[0], O[1], O[2] + radius);
    fprintf(outf, "%f\t%f\t%f\n", O[0] + radius, O[1], O[2]);
    fprintf(outf, "%f\t%f\t%f\n", O[0], O[1] + radius, O[2]);
    fprintf(outf, "%f\t%f\t%f\n", O[0] - radius, O[1], O[2]);
    fprintf(outf, "%f\t%f\t%f\n", O[0], O[1] - radius, O[2]);
    fprintf(outf, "%f\t%f\t%f\n", O[0], O[1], O[2] - radius);

    /*
     * list of (triangular) faces, each defined by 3 vertices, as
     * a->b->c->a
     */
    fprintf(outf, "4 0 1 2 0 %f\t%f\t%f 1.0\n", r, g, b);
    fprintf(outf, "4 0 2 3 0 %f\t%f\t%f 1.0\n", r, g, b);
    fprintf(outf, "4 0 3 4 0 %f\t%f\t%f 1.0\n", r, g, b);
    fprintf(outf, "4 0 4 1 0 %f\t%f\t%f 1.0\n", r, g, b);
    fprintf(outf, "4 5 1 2 5 %f\t%f\t%f 1.0\n", r, g, b);
    fprintf(outf, "4 5 2 3 5 %f\t%f\t%f 1.0\n", r, g, b);
    fprintf(outf, "4 5 3 4 5 %f\t%f\t%f 1.0\n", r, g, b);
    fprintf(outf, "4 5 4 1 5 %f\t%f\t%f 1.0\n", r, g, b);

    fclose(outf);
}

static void off_cone(const char *name, double *origin, double *dir,
		     const double l, const double r, const double g,
		     const double b)
{
    double P[3], g_P[3];
    const double O[] = { 0.0, 0.0, 0.0 };
    double alpha, beta;
    int i;

    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "7 7 0\n");	/* 7 vertices, 7 faces */

    /*
     * determine alpha, beta for transformation from local to global system.
     * discard 'P'
     */
    g2l_off(O, dir, P, &alpha, &beta);

    fprintf(outf, "%f\t%f\t%f\n", origin[0], origin[1], origin[2]);
    /*
     * vertices at hexagonal base of cone
     */
    for (i = 0; i < 6; i++) {
	const double arg = i / 3.0 * M_PI;	/* i * 60 / 360 * 2 * M_PI */
	const double radius = l / 4.0;
	P[0] = radius * sin(arg);
	P[1] = radius * cos(arg);
	P[2] = l;
	l2g_off(origin, P, g_P, alpha, beta);
	fprintf(outf, "%f\t%f\t%f\n", g_P[0], g_P[1], g_P[2]);
    }

    /*
     * list of (triangular) faces, each defined by 3 vertices, as
     * a->b->c->a
     */
    fprintf(outf, "4 0 1 2 0 %f\t%f\t%f 1.0\n", r, g, b);
    fprintf(outf, "4 0 2 3 0 %f\t%f\t%f 1.0\n", r, g, b);
    fprintf(outf, "4 0 3 4 0 %f\t%f\t%f 1.0\n", r, g, b);
    fprintf(outf, "4 0 4 5 0 %f\t%f\t%f 1.0\n", r, g, b);
    fprintf(outf, "4 0 5 6 0 %f\t%f\t%f 1.0\n", r, g, b);
    fprintf(outf, "4 0 6 1 0 %f\t%f\t%f 1.0\n", r, g, b);

    /* hexagonal face */
    fprintf(outf, "7 1 2 3 4 5 6 1 %f\t%f\t%f 1.0\n", r, g, b);

    fclose(outf);
}


static void off_plane(const char *name, const double *P, const double *N,
		      const double x, const double y, const double rf,
		      const double gf, const double bf, const double rb,
		      const double gb, const double bb, const double dz)
/*
 * writes 'OFF' file to 'name.off' for two sided plane. the plane is defined
 * by the point 'P' and its normal vector 'N'.
 * two planes are actually draw:
 *	- front side ('P2', offset by 'dz'*'N') colored 'rf','gf','bf'
 *	- back side ('P')colored rb,gb,bb
 * the dimension of the planes drawn is 'x' by 'y'.
 * the structure of the vertices for two parallel
 * planes is identical to the one for a block. we
 * thus can abuse 'block_vertices()' for this task
 */
{
    int i;
    double P2[3];
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "8 2 0\n");	/* 8 vertices, 2 faces */

    for (i = 0; i < 3; i++)	/* 'P2' is point on front side */
	P2[i] = P[i] + dz * N[i];

    block_vertices(outf, P, P2, x, y);

    /*
     * print back face ('rb', 'gb', 'bb'), front face ('rf', 'gf', 'bf') 
     */
    fprintf(outf, "5 0 1 2 3 0 %f\t%f\t%f 1.0\n", rb, gb, bb);
    fprintf(outf, "5 4 5 6 7 4 %f\t%f\t%f 1.0\n", rf, gf, bf);

    fclose(outf);
}

static void off_rectangle(const char *name, const double *P,
			  const double *X, const double *Y,
			  const double rf, const double gf,
			  const double bf, const double rb,
			  const double gb, const double bb,
			  const double dz)
{
    double N[3];
    FILE *outf = open_off(name);

    cross_product(X, Y, N);
    normalize(N);

    fprintf(outf, "OFF\n");
    fprintf(outf, "8 2 0\n");	/* 6 vertices, 2 faces */

    /* front side */
    fprintf(outf, "%f\t%f\t%f\n", P[0], P[1], P[2]);
    fprintf(outf, "%f\t%f\t%f\n", P[0] + X[0], P[1] + X[1], P[2] + X[2]);
    fprintf(outf, "%f\t%f\t%f\n", P[0] + X[0] + Y[0], P[1] + X[1] + Y[1],
	    P[2] + X[2] + Y[2]);
    fprintf(outf, "%f\t%f\t%f\n", P[0] + Y[0], P[1] + Y[1], P[2] + Y[2]);

    /* backside side: offset by 'dz' times 'N' */
    fprintf(outf, "%f\t%f\t%f\n", P[0] + dz * N[0], P[1] + dz * N[1],
	    P[2] + dz * N[2]);
    fprintf(outf, "%f\t%f\t%f\n", P[0] + dz * N[0] + X[0],
	    P[1] + dz * N[1] + X[1], P[2] + dz * N[2] + X[2]);
    fprintf(outf, "%f\t%f\t%f\n", P[0] + dz * N[0] + X[0] + Y[0],
	    P[1] + dz * N[1] + X[1] + Y[1],
	    P[2] + dz * N[2] + X[2] + Y[2]);
    fprintf(outf, "%f\t%f\t%f\n", P[0] + dz * N[0] + Y[0],
	    P[1] + dz * N[1] + Y[1], P[2] + dz * N[2] + Y[2]);

    /*
     * print back face ('rb', 'gb', 'bb'), front face ('rf', 'gf', 'bf') 
     */
    fprintf(outf, "5 0 1 2 3 0 %f\t%f\t%f 1.0\n", rb, gb, bb);
    fprintf(outf, "5 4 5 6 7 4 %f\t%f\t%f 1.0\n", rf, gf, bf);
    fclose(outf);

}

static void off_triangle(const char *name, const double *P1,
			 const double *P2, const double *P3,
			 const double *N, const double rf, const double gf,
			 const double bf, const double rb, const double gb,
			 const double bb, const double d)
/*
 * writes 'OFF' file to 'name.off' for two sided triangle. the triangle is defined
 * by the point 'P1' and the two sides 'P2' and 'P3'. its normal vector is 'N'.
 * two triangles are actually draw:
 *	- front side ('Pi', offset by 'dz'*'N') colored 'rf','gf','bf'
 *	- back side ('Pi')colored rb,gb,bb
 */
{
    const double dx = N[0] * d;
    const double dy = N[1] * d;
    const double dz = N[2] * d;
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "6 2 0\n");	/* 6 vertices, 2 faces */

    fprintf(outf, "%f\t%f\t%f\n", P1[0], P1[1], P1[2]);
    fprintf(outf, "%f\t%f\t%f\n", P2[0] + P1[0], P2[1] + P1[1],
	    P2[2] + P1[2]);
    fprintf(outf, "%f\t%f\t%f\n", P3[0] + P1[0], P3[1] + P1[1],
	    P3[2] + P1[2]);
    fprintf(outf, "%f\t%f\t%f\n", P1[0] + dx, P1[1] + dy, P1[2] + dz);
    fprintf(outf, "%f\t%f\t%f\n", P2[0] + P1[0] + dx, P2[1] + P1[1] + dy,
	    P2[2] + P1[2] + dz);
    fprintf(outf, "%f\t%f\t%f\n", P3[0] + P1[0] + dx, P3[1] + P1[1] + dy,
	    P3[2] + P1[2] + dz);

    /*
     * print back face ('rb', 'gb', 'bb'), front face ('rf', 'gf', 'bf') 
     */
    fprintf(outf, "4 0 1 2 0 %f\t%f\t%f 1.0\n", rb, gb, bb);
    fprintf(outf, "4 3 4 5 0 %f\t%f\t%f 1.0\n", rf, gf, bf);

    fclose(outf);
}

static void off_ellipsoid(const char *name, const double *origin,
			  const double *Z, const double *axes,
			  const double z_min, const double z_max,
			  const double ri, const double gi,
			  const double bi, const double ro,
			  const double go, const double bo,
			  const double dz)
{
#define N_TRANS 5
#define N_ROT 12

    int i, j;
    double alpha, beta;
    double O[] = { 0.0, 0.0, 0.0 };
    const double delta_z = (z_max - z_min) / (N_TRANS - 1);
    const double delta_phi = 2.0 * M_PI / N_ROT;
    double l[3], g[3];
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "%d %d 0\n", 2 * N_TRANS * N_ROT,
	    2 * (N_TRANS - 1) * N_ROT);

    /*
     * determine alpha, beta for transformation from local to global system.
     * discard 'l'
     */
    g2l_off(O, Z, l, &alpha, &beta);

    /*
     * calculate and output vertices
     * outside surface ('ro', 'go', 'bo')
     */
    for (i = 0; i < N_TRANS; i++) {
	double r;

	l[2] = z_min + i * delta_z;
	r = 1.0 - l[2] * l[2] / (axes[2] * axes[2]);
	if (r < GSL_SQRT_DBL_EPSILON)
	    r = GSL_SQRT_DBL_EPSILON;

	for (j = 0; j < N_ROT; j++) {
	    const double phi = j * delta_phi;

	    l[0] = r * axes[0] * sin(phi);
	    l[1] = r * axes[1] * cos(phi);

	    l2g_off(origin, l, g, alpha, beta);

	    fprintf(outf, "%f\t%f\t%f\n", g[0], g[1], g[2]);
	}
    }

    /*
     * calculate and output vertices
     * inside surface scaled by 'dz' ('ri', 'gi', 'bi')
     */
    for (i = 0; i < N_TRANS; i++) {
	double r;

	l[2] = z_min + i * delta_z;
	r = dz * (1.0 - l[2] * l[2] / (axes[2] * axes[2]));
	if (r < GSL_SQRT_DBL_EPSILON)
	    r = GSL_SQRT_DBL_EPSILON;

	for (j = 0; j < N_ROT; j++) {
	    const double phi = j * delta_phi;

	    l[0] = r * axes[0] * sin(phi);
	    l[1] = r * axes[1] * cos(phi);

	    l2g_off(origin, l, g, alpha, beta);

	    fprintf(outf, "%f\t%f\t%f\n", g[0], g[1], g[2]);
	}
    }

    /*
     * outside surface
     * output faces (N_ROT=12, N_TRANS=5)
     *   0       1      13      12       0
     *   1       2      14      13       1
     *           .
     *  10      11      23      22      10
     *  11       0      12      23      11
     *  12      13      25      24      12
     *           .
     *           .
     *           .
     *  46      47      59      58      46
     *  47      36      48      59      47
     */
    for (i = 0; i < (N_TRANS - 1) * N_ROT; i++)
	if ((i + 1) % N_ROT)
	    fprintf(outf, "5 %d\t%d\t%d\t%d\t%d\t%f\t%f\t%f 1.0\n",
		    i, i + 1, i + 1 + N_ROT, i + N_ROT, i, ro, go, bo);
	else
	    fprintf(outf, "5 %d\t%d\t%d\t%d\t%d\t%f\t%f\t%f 1.0\n",
		    i, i + 1 - N_ROT, i + 1, i + N_ROT, i, ro, go, bo);

    /*
     * same for inside surface
     */
    for (i = 0; i < (N_TRANS - 1) * N_ROT; i++)
	if ((i + 1) % N_ROT)
	    fprintf(outf, "5 %d\t%d\t%d\t%d\t%d\t%f\t%f\t%f 1.0\n",
		    i + N_TRANS * N_ROT, i + 1 + N_TRANS * N_ROT,
		    i + 1 + N_ROT + N_TRANS * N_ROT,
		    i + N_ROT + N_TRANS * N_ROT, i + N_TRANS * N_ROT, ri,
		    gi, bi);
	else
	    fprintf(outf, "5 %d\t%d\t%d\t%d\t%d\t%f\t%f\t%f 1.0\n",
		    i + N_TRANS * N_ROT, i + N_TRANS * N_ROT + 1 - N_ROT,
		    i + N_TRANS * N_ROT + 1, i + N_TRANS * N_ROT + N_ROT,
		    i + N_TRANS * N_ROT, ri, gi, bi);

    fclose(outf);

}

static void off_disk(const char *name, const double *origin,
		     const double *dir, const double r, const double rf,
		     const double gf, const double bf, const double rb,
		     const double gb, const double bb, const double dz)
/*
 * writes 'OFF' file to 'name.off' for disk. the disk is defined by its
 * center 'origin', normal vector 'dir' and radius 'r'.
 * two planes are actually draw:
 *	- front side (offset by 'dz') colored 'rf','gf','bf'
 *	- back side colored rb,gb,bb
 */
{
    double P[3], g_P[3];
    const double O[] = { 0.0, 0.0, 0.0 };
    double alpha, beta;
    int i;

    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "24 2 0\n");	/* 2*12 vertices, 2 faces */

    /*
     * determine alpha, beta for transformation from local to global system.
     * discard 'P'
     */
    g2l_off(O, dir, P, &alpha, &beta);

    /*
     * vertices at center
     */
    for (i = 0; i < 12; i++) {
	const double arg = i / 6.0 * M_PI;	/* i * 30 / 360 * 2 * M_PI */
	P[0] = r * sin(arg);
	P[1] = r * cos(arg);
	P[2] = 0.0;
	l2g_off(origin, P, g_P, alpha, beta);
	fprintf(outf, "%f\t%f\t%f\n", g_P[0], g_P[1], g_P[2]);
    }

    /*
     * vertices at center + 'dz'
     */
    for (i = 0; i < 12; i++) {
	const double arg = i / 6.0 * M_PI;	/* i * 30 / 360 * 2 * M_PI */
	P[0] = r * sin(arg);
	P[1] = r * cos(arg);
	P[2] = dz;
	l2g_off(origin, P, g_P, alpha, beta);
	fprintf(outf, "%f\t%f\t%f\n", g_P[0], g_P[1], g_P[2]);
    }
    /*
     * print back face ('rb', 'gb', 'bb'), front face ('rf', 'gf', 'bf') 
     */
    fprintf(outf, "13 ");
    for (i = 0; i < 12; i++)
	fprintf(outf, "%d ", i);
    fprintf(outf, "0 %f\t%f\t%f 1.0\n", rb, gb, bb);

    fprintf(outf, "13 ");
    for (i = 12; i < 24; i++)
	fprintf(outf, "%d ", i);
    fprintf(outf, "12 %f\t%f\t%f 1.0\n", rf, gf, bf);

    fclose(outf);
}

static void off_annulus(const char *name, const double *origin,
			const double *dir, const double R, const double r,
			const double rf, const double gf, const double bf,
			const double rb, const double gb, const double bb,
			const double dz)
/*
 * writes 'OFF' file to 'name.off' for disk. the disk is defined by its
 * center 'origin', normal vector 'dir' and radius 'r'.
 * two planes are actually draw:
 *	- front side (offset by 'dz') colored 'rf','gf','bf'
 *	- back side colored rb,gb,bb
 */
{
    double P[3], g_P[3];
    const double O[] = { 0.0, 0.0, 0.0 };
    double alpha, beta;
    int i;

    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "48 24 0\n");	/* 2*12 vertices,
				   2 faces with 12 patches each */

    /*
     * determine alpha, beta for transformation from local to global system.
     * discard 'P'
     */
    g2l_off(O, dir, P, &alpha, &beta);

    /*
     * vertices at center
     */
    for (i = 0; i < 12; i++) {
	const double arg = i / 6.0 * M_PI;	/* i * 30 / 360 * 2 * M_PI */
	const double sa = sin(arg);
	const double ca = cos(arg);
	P[0] = R * sa;
	P[1] = R * ca;
	P[2] = 0.0;
	l2g_off(origin, P, g_P, alpha, beta);
	fprintf(outf, "%f\t%f\t%f\n", g_P[0], g_P[1], g_P[2]);
	P[0] = r * sa;
	P[1] = r * ca;
	P[2] = 0.0;
	l2g_off(origin, P, g_P, alpha, beta);
	fprintf(outf, "%f\t%f\t%f\n", g_P[0], g_P[1], g_P[2]);
    }

    /*
     * vertices at center + 'dz'
     */
    for (i = 0; i < 12; i++) {
	const double arg = i / 6.0 * M_PI;	/* i * 30 / 360 * 2 * M_PI */
	const double sa = sin(arg);
	const double ca = cos(arg);
	P[0] = R * sa;
	P[1] = R * ca;
	P[2] = dz;
	l2g_off(origin, P, g_P, alpha, beta);
	fprintf(outf, "%f\t%f\t%f\n", g_P[0], g_P[1], g_P[2]);
	P[0] = r * sa;
	P[1] = r * ca;
	P[2] = dz;
	l2g_off(origin, P, g_P, alpha, beta);
	fprintf(outf, "%f\t%f\t%f\n", g_P[0], g_P[1], g_P[2]);
    }
    /*
     * print back face ('rb', 'gb', 'bb'), front face ('rf', 'gf', 'bf') 
     */
    for (i = 0; i < 11; i++) {
	const int base = 2 * i;
	fprintf(outf, "5 %d %d %d %d %d %f\t%f\t%f 1.0\n", base, base + 1,
		base + 3, base + 2, base, rb, gb, bb);
    }
    fprintf(outf, "5 22 23 1 0 22 %f\t%f\t%f 1.0\n", rb, gb, bb);

    for (i = 12; i < 23; i++) {
	const int base = 2 * i;
	fprintf(outf, "5 %d %d %d %d %d %f\t%f\t%f 1.0\n", base, base + 1,
		base + 3, base + 2, base, rf, gf, bf);
    }
    fprintf(outf, "5 46 47 25 24 46 %f\t%f\t%f 1.0\n", rf, gf, bf);

    fclose(outf);
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

	if (!strcmp(type, "one-sided plane screen")) {
	    double P[3], N[3];
	    double X[3], Y[3];

	    read_vector(this_t, "point", P);
	    read_vector_normalize(this_t, "normal", N);
	    read_vector_normalize(this_t, "x", X);

	    /* Y = N cross X */
	    cross_product(N, X, Y);

	    off_axes(name, P, X, Y, N);	/* local system */

	    /*
	     * draw plane:
	     *   front (red, counter) at 'P'+'DZ'
	     *   back side (black) at 'P'
	     */
	    off_plane(name, P, N, 10.0, 10.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		      DZ);

	} else if (!strcmp(type, "two-sided plane screen")) {
	    double P[3], N[3];
	    double X[3], Y[3];

	    read_vector(this_t, "point", P);
	    read_vector_normalize(this_t, "normal", N);
	    read_vector_normalize(this_t, "x", X);

	    /* Y = N cross X */
	    cross_product(N, X, Y);

	    off_axes(name, P, X, Y, N);	/*local system */

	    /*
	     * draw plane:
	     *   front (red, counter) at 'P'+'DZ'
	     *   back side (red, counter) at 'P'
	     */
	    off_plane(name, P, N, 10.0, 10.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
		      DZ);

	} else if (!strcmp(type, "rectangle")) {
	    int j;
	    double P[3], N[3];
	    double X[3], Y[3];

	    /*
	     * read three corner points:
	     * P1:
	     * P2: -> P2-P1 defines x axis
	     * P3: -> P3-P1 defines y axis
	     */
	    read_vector(this_t, "P1", P);
	    read_vector(this_t, "P2", X);
	    read_vector(this_t, "P3", Y);
	    for (j = 0; j < 3; j++) {
		X[j] -= P[j];
		Y[j] -= P[j];
	    }

	    /*
	     * draw plane:
	     *   front (white, mirror) at 'P'+'DZ'
	     *   rear side (black, absorbs) at 'P'
	     */
	    off_rectangle(name, P, X, Y, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, DZ);

	    /* make 'P1' point to center of rectangle */
	    for (j = 0; j < 3; j++)
		P[j] += (X[j] + Y[j]) / 2.0;

	    normalize(X);
	    normalize(Y);
	    /* surface normal N = X cross Y */
	    cross_product(X, Y, N);

	    off_axes(name, P, X, Y, N);	/*local system */

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


	    /* N = P2 cross P3 */
	    cross_product(P2, P3, N);
	    normalize(N);

	    /*
	     * draw triangle:
	     *   front (white, mirror) at 'P1' + 'N' * 'DZ'
	     *   rear side (black, absorbs) at 'P1'
	     */
	    off_triangle(name, P1, P2, P3, N, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
			 DZ);

	    /* 'Y' = 'N' cross 'P2' */
	    normalize(P2);
	    cross_product(N, P2, Y);

	    off_axes(name, P1, P2, Y, N);	/*local system */

	} else if (!strcmp(type, "ellipsoid")) {
	    double O[3], X[3], Y[3], Z[3];
	    double axes[3];
	    double z_min, z_max;

	    read_vector(this_t, "center", O);
	    read_vector_normalize(this_t, "x", X);
	    read_vector_normalize(this_t, "z", Z);
	    read_vector(this_t, "axes", axes);

	    config_setting_lookup_float(this_t, "z_min", &z_min);
	    config_setting_lookup_float(this_t, "z_max", &z_max);

	    orthonormalize(X, Y, Z);

	    off_axes(name, O, X, Y, Z);	/*local system */

	    /*
	     * draw ellipsoid:
	     *   inside: (white, reflecting) scaled by 1-'DZ'
	     *   outside (black, absorbing)
	     */
	    off_ellipsoid(name, O, Z, axes, z_min, z_max, 1.0, 1.0, 1.0,
			  0.0, 0.0, 0.0, 1.0 - DZ);

	} else if (!strcmp(type, "disk")) {
	    double O[3], X[3], Y[3], Z[3];
	    double r;

	    read_vector(this_t, "P", O);
	    read_vector_normalize(this_t, "x", X);
	    read_vector_normalize(this_t, "N", Z);
	    config_setting_lookup_float(this_t, "r", &r);

	    orthonormalize(X, Y, Z);

	    off_axes(name, O, X, Y, Z);	/*local system */

	    /*
	     * draw disk:
	     *   front (white, mirror) at 'P'+'DZ'
	     *   rear side (black, absorbs) at 'P'
	     */
	    off_disk(name, O, Z, r, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, DZ);

	} else if (!strcmp(type, "annulus")) {
	    double O[3], X[3], Y[3], Z[3];
	    double R, r;

	    read_vector(this_t, "P", O);
	    read_vector_normalize(this_t, "x", X);
	    read_vector_normalize(this_t, "N", Z);
	    config_setting_lookup_float(this_t, "R", &R);
	    config_setting_lookup_float(this_t, "r", &r);

	    orthonormalize(X, Y, Z);

	    off_axes(name, O, X, Y, Z);	/*local system */

	    /*
	     * draw disk:
	     *   front (white, mirror) at 'P'+'DZ'
	     *   rear side (black, absorbs) at 'P'
	     */
	    off_annulus(name, O, Z, R, r, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
			DZ);

	}
    }				/* end all targets */
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

	if (!strcmp(type, "uniform point source")) {
	    double O[3];

	    read_vector(this_s, "origin", O);

	    /*
	     * draw yellow (rgb=1.0,1.0,0.0) octahedron
	     * with size 1.2
	     * at origin 'O'
	     */
	    off_sphere(name, O, 1.2, 1.0, 1.0, 0.0);
	} /* end 'uniform point source' */
	else if (!strcmp(type, "sphere")) {
	    double O[3];
	    double radius;

	    read_vector(this_s, "origin", O);
	    config_setting_lookup_float(this_s, "radius", &radius);

	    /*
	     * draw yellow (rgb=1.0,1.0,0.0) octahedron
	     * with size 'radius'
	     * at origin 'O'
	     */
	    off_sphere(name, O, radius, 1.0, 1.0, 0.0);
	} /* end 'sphere' */
	else if (!strcmp(type, "spot source")) {
	    double O[3], dir[3];

	    read_vector(this_s, "origin", O);
	    read_vector(this_s, "direction", dir);

	    /*
	     * draw yellow (rgb=1.0,1.0,0.0) cone
	     * with size 1.2
	     * at origin 'O'
	     * in direction 'dir'
	     */
	    off_cone(name, O, dir, 1.2, 1.0, 1.0, 0.0);
	}
	/* end 'spot source' */
    }				/* end all sources */
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
    int i;
    double norm;

    for (i = 0; i < 3; i++)	/* translate to origin */
	L[i] = N[i] - P[i];

    norm = cblas_dnrm2(3, L, 1);	/* length of L */

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
    int i;
    const double ca = cos(-alpha);
    const double cb = cos(-beta);
    const double sa = sin(-alpha);
    const double sb = sin(-beta);

    const double x = L[0] * ca - L[2] * sa;

    G[0] = x * cb - L[1] * sb;
    G[1] = -x * sb - L[1] * cb;
    G[2] = L[0] * sa + L[2] * ca;;

    for (i = 0; i < 3; i++)
	G[i] += P[i];

}

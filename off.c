/*	off.c
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <math.h>
#include <stdio.h>

#include <gsl/gsl_cblas.h>

#include "off.h"

static FILE *open_off(const char *name)
{
    char f_name[256];

    snprintf(f_name, 256, "%s.off", name);
    return fopen(f_name, "w");
}

static void g2l_off(const double *P, const double *N, double *L,
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

static void l2g_off(const double *P, const double *L, double *G,
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

extern void off_axes(const double size)
{
/*
 * draw axes of global coordinate system as colored bars
 * x-axis: red
 * y-axis: green
 * z-axis: blue
 */
    const double s = size / 25.0;
    FILE *outf = open_off("axes");

    fprintf(outf, "OFF\n");
    fprintf(outf, "24 18 0\n\n");

    /*
     * vertices of x axis
     */
    fprintf(outf, "0\t%f\t%f\n", s, s);
    fprintf(outf, "0\t%f\t%f\n", s, -s);
    fprintf(outf, "0\t%f\t%f\n", -s, -s);
    fprintf(outf, "0\t%f\t%f\n", -s, s);
    fprintf(outf, "%f\t%f\t%f\n", size, s, s);
    fprintf(outf, "%f\t%f\t%f\n", size, s, -s);
    fprintf(outf, "%f\t%f\t%f\n", size, -s, -s);
    fprintf(outf, "%f\t%f\t%f\n", size, -s, s);
    /*
     * vertices of z axis
     */
    fprintf(outf, "%f\t0\t%f\n", s, s);
    fprintf(outf, "%f\t0\t%f\n", s, -s);
    fprintf(outf, "%f\t0\t%f\n", -s, -s);
    fprintf(outf, "%f\t0\t%f\n", -s, s);
    fprintf(outf, "%f\t%f\t%f\n", s, size, s);
    fprintf(outf, "%f\t%f\t%f\n", s, size, -s);
    fprintf(outf, "%f\t%f\t%f\n", -s, size, -s);
    fprintf(outf, "%f\t%f\t%f\n", -s, size, s);
    /*
     * vertices of z axis
     */
    fprintf(outf, "%f\t%f\t0\n", s, s);
    fprintf(outf, "%f\t%f\t0\n", s, -s);
    fprintf(outf, "%f\t%f\t0\n", -s, -s);
    fprintf(outf, "%f\t%f\t0\n", -s, s);
    fprintf(outf, "%f\t%f\t%f\n", s, s, size);
    fprintf(outf, "%f\t%f\t%f\n", s, -s, size);
    fprintf(outf, "%f\t%f\t%f\n", -s, -s, size);
    fprintf(outf, "%f\t%f\t%f\n", -s, s, size);

    /*
       faces of x axis in red
     */
    fprintf(outf, "5 0 1 2 3 0 1.0 0.0 0.0\n");
    fprintf(outf, "5 0 1 5 4 0 1.0 0.0 0.0\n");
    fprintf(outf, "5 0 3 7 4 0 1.0 0.0 0.0\n");
    fprintf(outf, "5 2 1 5 6 2 1.0 0.0 0.0\n");
    fprintf(outf, "5 2 3 7 6 2 1.0 0.0 0.0\n");
    fprintf(outf, "5 4 5 6 7 4 1.0 0.0 0.0\n");
    /*
       faces of x axis in green
     */
    fprintf(outf, "5 8 9 10 11 8 0.0 1.0 0.0\n");
    fprintf(outf, "5 8 9 13 12 8 0.0 1.0 0.0\n");
    fprintf(outf, "5 8 11 15 12 8 0.0 1.0 0.0\n");
    fprintf(outf, "5 10 9 13 14 10 0.0 1.0 0.0\n");
    fprintf(outf, "5 10 11 15 14 10 0.0 1.0 0.0\n");
    fprintf(outf, "5 12 13 14 15 12 0.0 1.0 0.0\n");
    /*
       faces of x axis in blue
     */
    fprintf(outf, "5 16 17 18 19 16 0.0 0.0 1.0\n");
    fprintf(outf, "5 16 17 21 20 16 0.0 0.0 1.0\n");
    fprintf(outf, "5 16 19 23 20 16 0.0 0.0 1.0\n");
    fprintf(outf, "5 18 17 21 22 18 0.0 0.0 1.0\n");
    fprintf(outf, "5 18 19 23 22 18 0.0 0.0 1.0\n");
    fprintf(outf, "5 20 21 22 23 20 0.0 0.0 1.0\n");

    fclose(outf);
}


extern void off_sphere(const char *name, double *O, const double radius,
		       const double r, const double g, const double b)
{
/*
 * print octaeder as source
 */
    FILE *outf = open_off(name);

    fprintf(outf, "OFF\n");
    fprintf(outf, "6 8 0\n\n");	/* 6 vertices, 8 faces */

    /*
     * list of vertices
     * 0 ... 5
     */
    fprintf(outf, "%f\t%f\t%f\n", O[0], O[1], O[3] + radius);
    fprintf(outf, "%f\t%f\t%f\n", O[0] + radius, O[1], O[3]);
    fprintf(outf, "%f\t%f\t%f\n", O[0], O[1] + radius, O[3]);
    fprintf(outf, "%f\t%f\t%f\n", O[0] - radius, O[1], O[3]);
    fprintf(outf, "%f\t%f\t%f\n", O[0], O[1] - radius, O[3]);
    fprintf(outf, "%f\t%f\t%f\n", O[0], O[1], O[3] - radius);

    /*
     * list of (triangular) faces, each defined by 3 vertices, as
     * a->b->c->a
     */
    fprintf(outf, "4 0 1 2 0 %f\t%f\t%f\n", r, g, b);
    fprintf(outf, "4 0 2 3 0 %f\t%f\t%f\n", r, g, b);
    fprintf(outf, "4 0 3 4 0 %f\t%f\t%f\n", r, g, b);
    fprintf(outf, "4 0 4 1 0 %f\t%f\t%f\n", r, g, b);
    fprintf(outf, "4 5 1 2 5 %f\t%f\t%f\n", r, g, b);
    fprintf(outf, "4 5 2 3 5 %f\t%f\t%f\n", r, g, b);
    fprintf(outf, "4 5 3 4 5 %f\t%f\t%f\n", r, g, b);
    fprintf(outf, "4 5 4 1 5 %f\t%f\t%f\n", r, g, b);

    fclose(outf);
}

/*	off.c
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <stdio.h>

#include "off.h"

static FILE *open_off(const char *name)
{
    char f_name[256];

    snprintf(f_name, 256, "%s.off", name);
    return fopen(f_name, "w");
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

/*	get_flux.c
 *
 * Copyright (C) 2011 - 2018 Ivo Alxneit, Paul Scherrer Institute
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
#include <getopt.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_math.h>

#include "math_utils.h"
#include "io_utils.h"


static void help(void)
{
    fprintf(stdout, "\nget_flux Version: %s  %s\n", RELEASE, RELEASE_INFO);
    fprintf(stdout, "Usage: get_flux\n");
    fprintf(stdout,
	    "       --global, -g      report flux distribution in global\n");
    fprintf(stdout,
	    "                         coordinate system [local]\n");
    fprintf(stdout, "       --help, -h        Print this help message\n");
    fprintf(stdout, "       --Version, -V     Print version number\n");
    fprintf(stdout, "\n");
}


int main(int argc, char **argv)
{
    double lambda_max = GSL_DBL_MAX;
    double lambda_min = 0.0;
    float t[MAX_FLOAT_ITEMS];
    unsigned char tmp_8;
    size_t idx_l;
    size_t n_float_items_read, n_uchar_items_read;
    int coordinates = LOCAL;
    double M[9];
    double origin[3];

    while (1) {
	int c;
	int option_index = 0;
	static struct option long_options[] = {
	    {"global", no_argument, 0, 'g'},
	    {"minl", required_argument, 0, 'l'},
	    {"maxl", required_argument, 0, 'L'},
	    {"Version", no_argument, 0, 'V'},
	    {"help", no_argument, 0, 'h'},
	    {0, 0, 0, 0}
	};

	c = getopt_long(argc, argv, "gl:L:Vh", long_options,
			&option_index);

	if (c == -1)
	    break;

	switch (c) {

	case 'g':
	    coordinates = GLOBAL;
	    break;

	case 'l':
	    lambda_min = atof(optarg);
	    break;

	case 'L':
	    lambda_max = atof(optarg);
	    break;

	case 'V':
	    fprintf(stdout, "get_flux Version: %s  %s\n", RELEASE,
		    RELEASE_INFO);
	    exit(EXIT_SUCCESS);
	    break;

	case 'h':
	    help();
	    exit(EXIT_SUCCESS);
	    break;

	default:
	    exit(EXIT_FAILURE);
	}			/* end 'switch (c)' */
    }				/* end 'while(1) */

    if (optind < argc) {
	fprintf(stderr, "ERROR: non-option ARGV-elements:");

	while (optind < argc)
	    fprintf(stderr, "%s", argv[optind++]);
	putc('\n', stderr);

	exit(EXIT_FAILURE);
    }

    if (get_idx(stdin, &idx_l) == ERR)
	exit(EXIT_FAILURE);

    if (coordinates == LOCAL) {
	if (skip_N_comments(stdin, HEADER_LINES) == ERR)
	    return ERR;
    } else
	read_transformation(stdin, M, origin);

    /*
     * read data. (x,y,[z,]lambda)
     */
    n_float_items_read = fread(t, sizeof(float), idx_l, stdin);

    if (!n_float_items_read) {
	fprintf(stderr, "No data found\n");
	exit(EXIT_FAILURE);
    }
    n_uchar_items_read = fread(&tmp_8, sizeof(unsigned char), 1, stdin);

    while (n_float_items_read && n_uchar_items_read) {
	double g_xyz[3];
	double l_xyz[3];


	if ((n_float_items_read < idx_l) || !n_uchar_items_read) {	/* insufficient data read */
#if __STDC_VERSION__ >= 199901L
	    fprintf(stderr,
		    "Incomplete data set read (%zu floats instead of %zu and %zu chars instead of 1)\n",
		    n_float_items_read, idx_l, n_uchar_items_read);
#else
	    fprintf(stderr,
		    "Incomplete data set read (%lu floats instead of %lu and %lu chars instead of 1)\n",
		    (unsigned long) n_float_items_read,
		    (unsigned long) idx_l,
		    (unsigned long) n_uchar_items_read);
#endif
	    exit(EXIT_FAILURE);
	}

	if ((t[idx_l - 1] >= lambda_min) && (t[idx_l - 1] <= lambda_max)) {
	    /*
	     * lambda is within selected interval
	     */
	    if (idx_l == MAX_FLOAT_ITEMS) {	/* non-planar target */
		if (coordinates == GLOBAL) {
		    l_xyz[0] = t[0];
		    l_xyz[1] = t[1];
		    l_xyz[2] = t[2];

		    l2g(M, origin, l_xyz, g_xyz);
		    fprintf(stdout, "%12.6g\t%12.6g\t%12.6g\t%12.6g\t%u\n",
			    g_xyz[0], g_xyz[1], g_xyz[2], t[idx_l - 1],
			    tmp_8);

		} else
		    fprintf(stdout, "%12.6g\t%12.6g\t%12.6g\t%12.6g\t%u\n",
			    t[0], t[1], t[2], t[idx_l - 1], tmp_8);

	    } else {		/* planar target */
		if (coordinates == GLOBAL) {
		    l_xyz[0] = t[0];
		    l_xyz[1] = t[1];
		    l_xyz[2] = 0.0;

		    l2g(M, origin, l_xyz, g_xyz);
		    fprintf(stdout, "%12.6g\t%12.6g\t%12.6g\t%12.6g\t%u\n",
			    g_xyz[0], g_xyz[1], g_xyz[2], t[idx_l - 1],
			    tmp_8);

		} else
		    fprintf(stdout, "%12.6g\t%12.6g\t%12.6g\t%u\n", t[0],
			    t[1], t[idx_l - 1], tmp_8);
	    }
	}

	n_float_items_read = fread(t, sizeof(float), idx_l, stdin);
	n_uchar_items_read =
	    fread(&tmp_8, sizeof(unsigned char), 1, stdin);

    }

    exit(EXIT_SUCCESS);
}

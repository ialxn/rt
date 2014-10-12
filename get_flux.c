/*	spectrum.c
 *
 * Copyright (C) 2011,2012,2013,2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#include <getopt.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_histogram.h>

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
    float t[MAX_FLOAT_ITEMS];
    unsigned char tmp_8;
    size_t idx_p, idx_l;
    size_t n_float_items_read, n_uchar_items_read;
    int coordinates = LOCAL;
    double M[9];
    double origin[3];

    while (1) {
	int c;
	int option_index = 0;
	static struct option long_options[] = {
	    {"global", no_argument, 0, 'g'},
	    {"Version", no_argument, 0, 'V'},
	    {"help", no_argument, 0, 'h'},
	    {0, 0, 0, 0}
	};

	c = getopt_long(argc, argv, "gVh", long_options, &option_index);

	if (c == -1)
	    break;

	switch (c) {

	case 'g':
	    coordinates = GLOBAL;
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

    if (get_idx(stdin, &idx_l, &idx_p) == ERR)
	exit(EXIT_FAILURE);

    if (coordinates == LOCAL) {
	if (skip_header(stdin) == ERR)
	    return ERR;
    } else
	read_transformation(stdin, M, origin);

    /*
     * read data. (x,y,[z,]power,lambda)
     */
    n_float_items_read = fread(t, sizeof(float), idx_l + 1, stdin);

    if (!n_float_items_read) {
	fprintf(stderr, "No data found\n");
	exit(EXIT_FAILURE);
    }
    n_uchar_items_read = fread(&tmp_8, sizeof(unsigned char), 1, stdin);

    while (n_float_items_read && n_uchar_items_read) {

	if ((n_float_items_read < idx_l + 1) || !n_uchar_items_read) {	/* insufficient data read */
	    fprintf(stderr,
		    "Incomplete data set read (%d floats instead of %d and %u chars instead of 1)\n",
		    n_float_items_read, idx_l + 1, n_uchar_items_read);
	    exit(EXIT_FAILURE);
	}

	if (idx_l == MAX_FLOAT_ITEMS - 1) {
	    if (coordinates == GLOBAL) {
		double g_xyz[3];
		double l_xyz[3];

		l_xyz[0] = t[0];
		l_xyz[1] = t[1];
		l_xyz[2] = t[2];

		l2g(M, origin, l_xyz, g_xyz);
		fprintf(stdout, "%e\t%e\t%e\t%e\t%e\t%u\n", g_xyz[0],
			g_xyz[1], g_xyz[2], t[idx_p], t[idx_l], tmp_8);

	    } else
		fprintf(stdout, "%e\t%e\t%e\t%e\t%e\t%u\n", t[0], t[1],
			t[2], t[idx_p], t[idx_l], tmp_8);

	} else {
	    if (coordinates == GLOBAL) {
		double g_xyz[3];
		double l_xyz[3];

		l_xyz[0] = t[0];
		l_xyz[1] = t[1];
		l_xyz[2] = 0.0;

		l2g(M, origin, l_xyz, g_xyz);
		fprintf(stdout, "%e\t%e\t%e\t%e\t%e\t%u\n", g_xyz[0],
			g_xyz[1], g_xyz[2], t[idx_p], t[idx_l], tmp_8);

	    } else
		fprintf(stdout, "%e\t%e\t%e\t%e\t%u\n", t[0], t[1],
			t[idx_p], t[idx_l], tmp_8);
	}

	n_float_items_read = fread(t, sizeof(float), idx_l + 1, stdin);
	n_uchar_items_read =
	    fread(&tmp_8, sizeof(unsigned char), 1, stdin);

    }

    exit(EXIT_SUCCESS);
}

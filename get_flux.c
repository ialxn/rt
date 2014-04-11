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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_histogram.h>

#include "io_utils.h"
#include "version.h"



static void help(void)
{
    fprintf(stdout, "\nget_flux Version %s(%s) AI52\n\n", RELEASE,
	    RELEASE_DATE);
    fprintf(stdout, "Usage: get_flux\n");
    fprintf(stdout, "       --help, -h        Print this help message\n");
    fprintf(stdout, "       --Version, -V     Print version number\n");
    fprintf(stdout, "\n");
}


int main(int argc, char **argv)
{
    float t[MAX_ITEMS];
    size_t idx_p, idx_l;
    size_t n_items_read;

    while (1) {
	int c;
	int option_index = 0;
	static struct option long_options[] = {
	    {"Version", no_argument, 0, 'V'},
	    {"help", no_argument, 0, 'h'},
	    {0, 0, 0, 0}
	};

	c = getopt_long(argc, argv, "Vh", long_options, &option_index);

	if (c == -1)
	    break;

	switch (c) {

	case 'V':
	    fprintf(stdout, " get_flux Version %s (%s) AI52\n", RELEASE,
		    RELEASE_DATE);
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

    if (skip_header(stdin) == ERR)
	return ERR;

    /*
     * read data. (x,y,[z,]power,lambda)
     */
    n_items_read = fread(t, sizeof(float), idx_l + 1, stdin);

    if (!n_items_read) {
	fprintf(stderr, "No data found\n");
	exit(EXIT_FAILURE);
    }

    while (n_items_read) {

	if (n_items_read < idx_l + 1) {	/* insufficient data read */
	    fprintf(stderr,
		    "Incomplete data set read (%d instead of %d)\n",
		    n_items_read, idx_l + 1);
	    exit(EXIT_FAILURE);
	}

	if (idx_l == MAX_ITEMS - 1)
	    fprintf(stdout, "%f\t%f\t%f\t%f\t%f\n", t[0], t[1], t[2],
		    t[idx_p], t[idx_l]);
	else
	    fprintf(stdout, "%f\t%f\t%f\t%f\n", t[0], t[1], t[idx_p],
		    t[idx_l]);

	n_items_read = fread(t, sizeof(float), idx_l + 1, stdin);

    }

    exit(EXIT_SUCCESS);
}

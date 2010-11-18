/*	main.c
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>

#include "version.h"



static void help(void)
{
    fprintf(stdout, "\nrt Version %s (AI52)\n\n", VERSION);
    fprintf(stdout, "Usage: rt\n");
    fprintf(stdout,
	    "       --Version, -V          Print version number\n");
    fprintf(stdout,
	    "       --help, -h             Print this help message\n");
    fprintf(stdout, "\n");
}



int main(int argc, char **argv)
{

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
	    fprintf(stdout, "rt Version %s (AI52)\n", VERSION);
	    exit(EXIT_SUCCESS);
	    break;

	case 'h':
	    help();
	    exit(EXIT_SUCCESS);
	    break;

	default:
	    exit(EXIT_FAILURE);

	}
    }

    if (optind < argc) {

	fprintf(stderr, "ERROR: non-option ARGV-elements: ");

	while (optind < argc)
	    fprintf(stderr, "%s ", argv[optind++]);

	putc('\n', stderr);

	exit(EXIT_FAILURE);

    }

    return EXIT_SUCCESS;
}

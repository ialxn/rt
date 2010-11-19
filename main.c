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
#include <libconfig.h>
#include <stdlib.h>
#include <stdio.h>

#include "sources.h"
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



static int parse_input(config_t * cfg)
{
    config_init(cfg);

    /*
     * parse input file
     * on error, report it and exit.
     */
    if (!config_read(cfg, stdin)) {
	const char *fname = config_error_file(cfg);
	const char *text = config_error_text(cfg);
	const int line_nr = config_error_line(cfg);

	fprintf(stderr, "%s found during parsing of ", text);

	if (fname)
	    fprintf(stderr, "%s ", fname);
	else
	    fprintf(stderr, "stdin, ");

	fprintf(stderr, "line %d\n", line_nr);

	config_destroy(cfg);
	return (EXIT_FAILURE);
    }

    if (!check_sources(cfg))
	return EXIT_SUCCESS;
    else {
	config_destroy(cfg);
	return EXIT_FAILURE;
    }

}



int main(int argc, char **argv)
{
    config_t cfg;

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

    if (parse_input(&cfg))
	exit(EXIT_FAILURE);

    config_write(&cfg, stdout);

    config_destroy(&cfg);

    return EXIT_SUCCESS;
}

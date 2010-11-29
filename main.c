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
#include <gsl/gsl_rng.h>
#include <libconfig.h>
#include <stdlib.h>
#include <stdio.h>

#include "obj_lists.h"
#include "sources.h"
#include "targets.h"
#include "version.h"

#define CHECK_CONFIG 0
#define PRINT_GEOMETRY 1
#define RUN 2



static void output_geometry(config_t * cfg)
{
}

static void run_simulation(source_list_t * source_list,
			   target_list_t * target_list, int seed)
{
    const gsl_rng_type *T;
    gsl_rng *r;

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, (unsigned long int) abs(seed));

    gsl_rng_free(r);
}

static void help(void)
{
    fprintf(stdout, "\nrt Version %s (AI52)\n\n", VERSION);
    fprintf(stdout, "Usage: rt\n");
    fprintf(stdout, "       --mode, -m        select run mode [0].\n");
    fprintf(stdout,
	    "                         0: check and print input.\n");
    fprintf(stdout,
	    "                         1: output geometry (OFF files).\n");
    fprintf(stdout, "                         2: run simulation.\n");
    fprintf(stdout, "       --Version, -V     Print version number\n");
    fprintf(stdout, "       --help, -h        Print this help message\n");
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

    if (check_sources(cfg)) {
	config_destroy(cfg);
	return EXIT_FAILURE;
    }
    if (check_targets(cfg)) {
	config_destroy(cfg);
	return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;


}



int main(int argc, char **argv)
{
    config_t cfg;
    int mode = CHECK_CONFIG;

    while (1) {
	int c;
	int option_index = 0;
	static struct option long_options[] = {
	    {"mode", required_argument, 0, 'm'},
	    {"Version", no_argument, 0, 'V'},
	    {"help", no_argument, 0, 'h'},
	    {0, 0, 0, 0}
	};

	c = getopt_long(argc, argv, "m:Vh", long_options, &option_index);

	if (c == -1)
	    break;

	switch (c) {

	case 'm':
	    mode = atoi(optarg);
	    if ((mode > RUN) || (mode < CHECK_CONFIG)) {
		fprintf(stderr, "unknown 'mode' (%d) given.\n", mode);
		exit(EXIT_FAILURE);
	    }
	    break;
	case 'V':
	    fprintf(stdout, " rt Version %s (AI52)\n", VERSION);
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

	fprintf(stderr, "ERROR: non-option ARGV-elements:");

	while (optind < argc)
	    fprintf(stderr, "%s", argv[optind++]);

	putc('\n', stderr);

	exit(EXIT_FAILURE);

    }

    /*
     * parse and initialize input. input is checked for syntax and completness.
     * errors are reported on stderr.
     */
    if (parse_input(&cfg))
	exit(EXIT_FAILURE);

    switch (mode) {
	int seed;
	target_list_t *target_list;
	source_list_t *source_list;

    case CHECK_CONFIG:		/* print parsed input */
	config_write(&cfg, stdout);
	config_destroy(&cfg);
	exit(EXIT_SUCCESS);
	break;

    case PRINT_GEOMETRY:	/* print geometry files */
	output_geometry(&cfg);
	config_destroy(&cfg);
	exit(EXIT_SUCCESS);
	break;

    case RUN:			/* do the simulation */
	source_list = init_sources(&cfg);
	target_list = init_targets(&cfg);

	config_lookup_int(&cfg, "seed", &seed);
	config_destroy(&cfg);

	run_simulation(source_list, target_list, seed);

	source_list_free(source_list);
	target_list_free(target_list);
	free(source_list);
	free(target_list);
    default:
	break;

    }


    return EXIT_SUCCESS;
}

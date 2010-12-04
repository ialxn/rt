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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <libconfig.h>

#include <gsl/gsl_machine.h>
#include <gsl/gsl_rng.h>

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
    struct list_head *s_pos;
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    int dump_flag = 0;

    gsl_rng_set(r, (unsigned long int) abs(seed));

    list_for_each(s_pos, &(source_list->list)) {
	source_list_t *this_s = list_entry(s_pos, source_list_t, list);
	source_t *current_source = this_s->s;
	ray_t *ray = new_ray(current_source, r);	/* get first ray */

	while (ray) {		/* loop until source is exhausted */
	    while (ray) {	/* loop until ray is absorbed or leaves system */
		struct list_head *t_pos;
		target_t *nearest_target;
		double *nearest_intercept = NULL;
		double min_dist = GSL_DBL_MAX;
		int hits_target = 0;

		list_for_each(t_pos, &(target_list->list)) {	/* find closest target intercepted by ray */
		    target_list_t *this_t =
			list_entry(t_pos, target_list_t, list);
		    target_t *current_target = this_t->t;
		    double *current_intercept =
			interception(current_target, ray, &dump_flag);

		    if (current_intercept) {	/* interception found */
			int i;
			double dist = 0;

			hits_target = 1;
			for (i = 0; i < 3; i++) {
			    const double t =
				current_intercept[i] - ray->origin[i];
			    dist += t * t;
			}
			dist = sqrt(dist);	/* absolute distance origin of ray to intercept */

			if (dist < min_dist) {	/* new nearest target identified */

			    nearest_target = current_target;
			    if (nearest_intercept)
				free(nearest_intercept);

			    nearest_intercept = current_intercept;
			    min_dist = dist;

			} else	/* hit on far target. not used */
			    free(current_intercept);

		    }

		}		/* all targets tried */

		if (hits_target) {
		    ray =
			out_ray(nearest_target, ray, nearest_intercept,
				&dump_flag);
		    /* ray=NULL if it is absorbed by target */
		    free(nearest_intercept);
		    nearest_intercept = NULL;
		} else {	/* no target hit, ray is lost */
		    free(ray);
		    ray = NULL;	/* mark as absorbed */
		}

	    }			/* ray absorbed or lost */

	    ray = new_ray(current_source, r);	/* start next ray */

	}

    }

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
     * parse and initialize input. input is checked for syntax and completeness.
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

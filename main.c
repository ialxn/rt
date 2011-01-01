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

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_rng.h>

#include "obj_lists.h"
#include "off.h"
#include "sources.h"
#include "targets.h"
#include "version.h"

#define CHECK_CONFIG 0
#define PRINT_GEOMETRY 1
#define RUN 2

#define DZ 0.005
#define AXES_LENGTH 5.0

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
	    int j;
	    double P[3], N[3];
	    double X[3], Y[3];
	    double norm;
	    config_setting_t *this;

	    this = config_setting_get_member(this_t, "point");
	    for (j = 0; j < 3; j++)
		P[j] = config_setting_get_float_elem(this, j);

	    this = config_setting_get_member(this_t, "normal");
	    for (j = 0; j < 3; j++)
		N[j] = config_setting_get_float_elem(this, j);

	    /* normalize N */
	    norm = cblas_dnrm2(3, N, 1);
	    cblas_dscal(3, 1.0 / norm, N, 1);

	    this = config_setting_get_member(this_t, "x");
	    for (j = 0; j < 3; j++)
		X[j] = config_setting_get_float_elem(this, j);

	    /* normalize X */
	    norm = cblas_dnrm2(3, X, 1);
	    cblas_dscal(3, 1.0 / norm, X, 1);

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
	    int j;
	    double P[3], N[3];
	    double X[3], Y[3];
	    double norm;
	    config_setting_t *this;

	    this = config_setting_get_member(this_t, "point");
	    for (j = 0; j < 3; j++)
		P[j] = config_setting_get_float_elem(this, j);

	    this = config_setting_get_member(this_t, "normal");
	    for (j = 0; j < 3; j++)
		N[j] = config_setting_get_float_elem(this, j);

	    /* normalize N */
	    norm = cblas_dnrm2(3, N, 1);
	    cblas_dscal(3, 1.0 / norm, N, 1);

	    this = config_setting_get_member(this_t, "x");
	    for (j = 0; j < 3; j++)
		X[j] = config_setting_get_float_elem(this, j);

	    /* normalize X */
	    norm = cblas_dnrm2(3, X, 1);
	    cblas_dscal(3, 1.0 / norm, X, 1);

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

	} else if (!strcmp(type, "square")) {
	    int j;
	    double P[3], N[3];
	    double X[3], Y[3];
	    double norm;
	    double nX, nY;
	    config_setting_t *this;

	    this = config_setting_get_member(this_t, "point");
	    for (j = 0; j < 3; j++)
		P[j] = config_setting_get_float_elem(this, j);

	    this = config_setting_get_member(this_t, "x");
	    for (j = 0; j < 3; j++)
		X[j] = config_setting_get_float_elem(this, j);

	    this = config_setting_get_member(this_t, "y");
	    for (j = 0; j < 3; j++)
		Y[j] = config_setting_get_float_elem(this, j);

	    /* N = X cross Y */
	    cross_product(X, Y, N);

	    /* normalize N */
	    norm = cblas_dnrm2(3, N, 1);
	    cblas_dscal(3, 1.0 / norm, N, 1);

	    /* make 'P' point to center of square */
	    for (j = 0; j < 3; j++)
		P[j] += (X[j] + Y[j]) / 2.0;

	    nX = cblas_dnrm2(3, X, 1);
	    nY = cblas_dnrm2(3, Y, 1);
	    /*
	     * draw plane:
	     *   front (white, mirror) at 'P'+'DZ'
	     *   rear side (black, absorbs) at 'P'
	     */
	    off_plane(name, P, N, nX, nY, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
		      DZ);

	    cblas_dscal(3, 1.0 / nX, X, 1);
	    cblas_dscal(3, 1.0 / nY, Y, 1);
	    off_axes(name, P, X, Y, N);	/*local system */

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
	    int j;
	    double O[3];
	    const config_setting_t *origin =
		config_setting_get_member(this_s, "origin");

	    for (j = 0; j < 3; j++)
		O[j] = config_setting_get_float_elem(origin, j);

	    /*
	     * draw yellow (rgb=1.0,1.0,0.0) octahedron
	     * with size 1.2
	     * at origin 'O'
	     */
	    off_sphere(name, O, 1.2, 1.0, 1.0, 0.0);
	} /* end 'uniform point source' */
	else if (!strcmp(type, "spot source")) {
	    int j;
	    double O[3], dir[3];
	    config_setting_t *this_group;

	    this_group = config_setting_get_member(this_s, "origin");
	    for (j = 0; j < 3; j++)
		O[j] = config_setting_get_float_elem(this_group, j);

	    this_group = config_setting_get_member(this_s, "direction");
	    for (j = 0; j < 3; j++)
		dir[j] = config_setting_get_float_elem(this_group, j);

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

static void output_geometry(config_t * cfg)
{
    const double O[] = { 0.0, 0.0, 0.0 };
    const double X[] = { AXES_LENGTH, 0.0, 0.0 };
    const double Y[] = { 0.0, AXES_LENGTH, 0.0 };
    const double Z[] = { 0.0, 0.0, AXES_LENGTH };

    off_axes(NULL, O, X, Y, Z);
    output_sources(cfg);
    output_targets(cfg);
}

static void run_simulation(source_list_t * source_list,
			   target_list_t * target_list, const int seed,
			   const int n_targets)
{
    struct list_head *s_pos;
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    int dump_flag = 0;

    gsl_rng_set(r, (unsigned long int) abs(seed));

    fprintf(stdout,
	    "    using random number generator %s from Gnu Scientif Library\n",
	    gsl_rng_name(r));

    list_for_each(s_pos, &(source_list->list)) {
	source_list_t *this_s = list_entry(s_pos, source_list_t, list);
	source_t *current_source = this_s->s;
	const double ppr = get_ppr(current_source);	/* power per ray of 'current_source' */
	ray_t *ray = new_ray(current_source, r);	/* get first 'ray' of 'current_source' */

	fprintf(stdout, "        %s %s ... ",
		get_source_type(current_source),
		get_source_name(current_source));
	fflush(stdout);

	while (ray) {		/* loop until 'current_source' is exhausted */

	    while (ray) {	/* loop until 'ray' is absorbed or leaves system */
		struct list_head *t_pos;
		/*
		 * keep track of nearest target being hit by 'ray'.
		 * 'nearest_intercept' is used in calculation of
		 * 'ray' leaving 'nearest_target'.
		 */
		target_t *nearest_target;
		double *nearest_intercept = NULL;
		double min_dist = GSL_DBL_MAX;
		int hits_target = 0;	/* flag indicating that 'ray' hits any target */

		list_for_each(t_pos, &(target_list->list)) {	/* find closest target intercepted by 'ray' */
		    target_list_t *this_t =
			list_entry(t_pos, target_list_t, list);
		    target_t *current_target = this_t->t;
		    double *current_intercept =
			interception(current_target, ray, &dump_flag);

		    if (current_intercept) {	/* 'ray' hits 'current_target' */
			int i;
			double dist = 0;

			hits_target = 1;
			for (i = 0; i < 3; i++) {
			    const double t =
				current_intercept[i] - ray->origin[i];
			    dist += t * t;
			}
			dist = sqrt(dist);	/* absolute distance origin of ray to intercept */

			if (dist < min_dist) {	/* 'current targets' is closest target found until now */

			    nearest_target = current_target;
			    if (nearest_intercept)
				free(nearest_intercept);

			    nearest_intercept = current_intercept;
			    min_dist = dist;

			} else	/* hit on far target. not used */
			    free(current_intercept);

		    }
		    /* end 'if(current_intercept)' */
		}		/* all targets tried */

		if (hits_target) {	/* 'ray' hits 'nearest_target' */
		    ray =	/* 'out_ray' returns NULL if 'ray' is absorbed by target */
			out_ray(nearest_target, ray, ppr,
				nearest_intercept, &dump_flag, n_targets);
		    free(nearest_intercept);
		    nearest_intercept = NULL;
		} else {	/* no target hit, 'ray' is lost */
		    free(ray);
		    ray = NULL;	/* mark as absorbed */
		}
	    }			/* 'ray' absorbed or lost */
	    ray = new_ray(current_source, r);	/* start next ray */
	}			/* 'current_source' is exhausted */
	fprintf(stdout, "exhausted\n");
    }				/* all sources exhausted */
    gsl_rng_free(r);
}

static void help(void)
{
    fprintf(stdout, "\nrt Version %s (AI52)\n\n", VERSION);
    fprintf(stdout, "Usage: rt\n");
    fprintf(stdout,
	    "       --append, -a      append to output files. new seed must be given.\n");
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
    if (!config_read(cfg, stdin)) {	/* grammatically correct? */
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

    /*
     * check that sources and targets are correctly specified.
     * additional (unknown) configuration options are silently
     * ignored.
     */
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
    int seed;			/* seed for rng */
    char file_mode[2] = "w";	/* default file mode for output per target */

    while (1) {
	int c;
	int option_index = 0;
	static struct option long_options[] = {
	    {"append", required_argument, 0, 'a'},
	    {"mode", required_argument, 0, 'm'},
	    {"Version", no_argument, 0, 'V'},
	    {"help", no_argument, 0, 'h'},
	    {0, 0, 0, 0}
	};

	c = getopt_long(argc, argv, "a:m:Vh", long_options, &option_index);

	if (c == -1)
	    break;

	switch (c) {
	case 'a':
	    snprintf(file_mode, 2, "a");
	    seed = atoi(optarg);
	    break;
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
	}			/* end 'switch (c)' */
    }				/* end 'while(1) */

    if (optind < argc) {
	fprintf(stderr, "ERROR: non-option ARGV-elements:");

	while (optind < argc)
	    fprintf(stderr, "%s", argv[optind++]);
	putc('\n', stderr);

	exit(EXIT_FAILURE);
    }

    /*
     * parse and initialize configuration. input is checked for syntax
     * and completeness. errors are reported on stderr.
     */
    if (parse_input(&cfg))
	exit(EXIT_FAILURE);

    switch (mode) {
	int n_targets;		/* needed for dump cycle */
	int n_sources;
	target_list_t *target_list;	/* list of all sources */
	source_list_t *source_list;	/* list of all targets */

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
	fprintf(stdout, "rt version %s running ...\n", VERSION);

	source_list = init_sources(&cfg, &n_sources);
	fprintf(stdout, "    %d sources initialized\n", n_sources);
	target_list = init_targets(&cfg, &n_targets, file_mode);
	fprintf(stdout, "    %d targets initialized\n", n_targets);

	if (file_mode[0] == 'w') {	/* use seed from cfg, otherwise from command line */
	    config_lookup_int(&cfg, "seed", &seed);
	    fprintf(stdout,
		    "    using %d as seed for random number generator from config file\n",
		    seed);
	} else {
	    fprintf(stdout, "    appending to existing files\n");
	    fprintf(stdout,
		    "    using %d as seed for random number generator from command line\n",
		    seed);
	}

	config_destroy(&cfg);

	run_simulation(source_list, target_list, seed, n_targets);

	source_list_free(source_list);
	target_list_free(target_list);
	free(source_list);
	free(target_list);
    default:
	break;

    }
    return EXIT_SUCCESS;
}

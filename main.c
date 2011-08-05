/*	main.c
 *
 * Copyright (C) 2010, 2011 Ivo Alxneit
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

#include "io_util.h"
#include "obj_lists.h"
#include "off.h"
#include "sources.h"
#include "targets.h"
#include "vector_math.h"
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
	    double P[3], N[3];
	    double X[3], Y[3];

	    read_vector(this_t, "point", P);
	    read_vector_normalize(this_t, "normal", N);
	    read_vector_normalize(this_t, "x", X);

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
	    double P[3], N[3];
	    double X[3], Y[3];

	    read_vector(this_t, "point", P);
	    read_vector_normalize(this_t, "normal", N);
	    read_vector_normalize(this_t, "x", X);

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

	} else if (!strcmp(type, "rectangle")) {
	    int j;
	    double P[3], N[3];
	    double X[3], Y[3];
	    double nX, nY;

	    /*
	     * read three corner points:
	     * P1:
	     * P2: -> P2-P1 defines x axis
	     * P3: -> P3-P1 defines y axis
	     */
	    read_vector(this_t, "P1", P);
	    read_vector(this_t, "P2", X);
	    read_vector(this_t, "P3", Y);
	    for (j = 0; j < 3; j++) {
		X[j] -= P[j];
		Y[j] -= P[j];
	    }

	    /*
	     * draw plane:
	     *   front (white, mirror) at 'P'+'DZ'
	     *   rear side (black, absorbs) at 'P'
	     */
	    off_rectangle(name, P, X, Y, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, DZ);

	    /* make 'P1' point to center of rectangle */
	    for (j = 0; j < 3; j++)
		P[j] += (X[j] + Y[j]) / 2.0;

	    nX = normalize(X);
	    nY = normalize(Y);
	    /* surface normal N = X cross Y */
	    cross_product(X, Y, N);

	    off_axes(name, P, X, Y, N);	/*local system */

	} else if (!strcmp(type, "triangle")) {
	    int j;
	    double P1[3];
	    double P2[3];	/* redefined 'P2'='P2'-'P1' */
	    double P3[3];	/* redefined 'P3'='P3'-'P1' */
	    double N[3];	/* 'P2' cross 'P3' */
	    double Y[3];	/* local system. 'X' parallel 'P2' */
	    config_setting_t *this;

	    read_vector(this_t, "P1", P1);

	    this = config_setting_get_member(this_t, "P2");
	    for (j = 0; j < 3; j++)
		P2[j] = config_setting_get_float_elem(this, j) - P1[j];

	    this = config_setting_get_member(this_t, "P3");
	    for (j = 0; j < 3; j++)
		P3[j] = config_setting_get_float_elem(this, j) - P1[j];


	    /* N = P2 cross P3 */
	    cross_product(P2, P3, N);
	    normalize(N);

	    /*
	     * draw triangle:
	     *   front (white, mirror) at 'P1' + 'N' * 'DZ'
	     *   rear side (black, absorbs) at 'P1'
	     */
	    off_triangle(name, P1, P2, P3, N, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
			 DZ);

	    /* 'Y' = 'N' cross 'P2' */
	    normalize(P2);
	    cross_product(N, P2, Y);

	    off_axes(name, P1, P2, Y, N);	/*local system */

	} else if (!strcmp(type, "ellipsoid")) {
	    double O[3], X[3], Y[3], Z[3];
	    double axes[3];
	    double z_min, z_max;

	    read_vector(this_t, "center", O);
	    read_vector_normalize(this_t, "x", X);
	    read_vector_normalize(this_t, "z", Z);
	    read_vector(this_t, "axes", axes);

	    config_setting_lookup_float(this_t, "z_min", &z_min);
	    config_setting_lookup_float(this_t, "z_max", &z_max);

	    orthonormalize(X, Y, Z);

	    off_axes(name, O, X, Y, Z);	/*local system */

	    /*
	     * draw ellipsoid:
	     *   inside: (white, reflecting) scaled by 1-'DZ'
	     *   outside (black, absorbing)
	     */
	    off_ellipsoid(name, O, Z, axes, z_min, z_max, 1.0, 1.0, 1.0,
			  0.0, 0.0, 0.0, 1.0 - DZ);

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
	    double O[3];

	    read_vector(this_s, "origin", O);

	    /*
	     * draw yellow (rgb=1.0,1.0,0.0) octahedron
	     * with size 1.2
	     * at origin 'O'
	     */
	    off_sphere(name, O, 1.2, 1.0, 1.0, 0.0);
	} /* end 'uniform point source' */
	else if (!strcmp(type, "sphere")) {
	    double O[3];
	    double radius;

	    read_vector(this_s, "origin", O);
	    config_setting_lookup_float(this_s, "radius", &radius);

	    /*
	     * draw yellow (rgb=1.0,1.0,0.0) octahedron
	     * with size 'radius'
	     * at origin 'O'
	     */
	    off_sphere(name, O, radius, 1.0, 1.0, 0.0);
	} /* end 'sphere' */
	else if (!strcmp(type, "spot source")) {
	    double O[3], dir[3];

	    read_vector(this_s, "origin", O);
	    read_vector(this_s, "direction", dir);

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

static double distance(const double a[3], const double b[3])
{
    size_t i;
    double dist = 0.0;
    for (i = 0; i < 3; i++) {
	const double t = a[i] - b[i];
	dist += t * t;
    }
    return (sqrt(dist));

}

static void run_simulation(source_list_t * source_list,
			   target_list_t * target_list, const int seed,
			   const int n_targets)
{
    struct list_head *s_pos;
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    int dump_flag = 0;
    unsigned int n_lost = 0;

    gsl_rng_set(r, (unsigned long int) abs(seed));

    fprintf(stdout,
	    "    using random number generator %s from Gnu Scientif Library\n",
	    gsl_rng_name(r));

    list_for_each(s_pos, &(source_list->list)) {
	source_list_t *this_s = list_entry(s_pos, source_list_t, list);
	source_t *current_source = this_s->s;
	ray_t *ray;

	fprintf(stdout, "        %s %s ... ",
		get_source_type(current_source),
		get_source_name(current_source));
	fflush(stdout);

	while ((ray = new_ray(current_source, r))) {

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
			const double dist =
			    distance(current_intercept, ray->origin);

			hits_target = 1;

			if (dist < min_dist) {	/* 'current targets' is closest target found until now */

			    free(nearest_intercept);
			    nearest_target = current_target;
			    nearest_intercept = current_intercept;
			    min_dist = dist;

			} else	/* hit on far target. not used */
			    free(current_intercept);

		    }
		    /* end 'if(current_intercept)' */
		}		/* all targets tried */

		if (hits_target) {	/* 'ray' hits 'nearest_target' */
		    ray =	/* 'out_ray' returns NULL if 'ray' is absorbed by target */
			out_ray(nearest_target, ray, nearest_intercept, r,
				&dump_flag, n_targets);
		    free(nearest_intercept);
		    nearest_intercept = NULL;
		} else {	/* no target hit, 'ray' is lost */
		    n_lost++;
		    free(ray);
		    ray = NULL;	/* mark as absorbed */
		}
	    }			/* 'ray' absorbed or lost */
	}			/* 'current_source' is exhausted */
	fprintf(stdout, "exhausted\n");
    }				/* all sources exhausted */
    gsl_rng_free(r);

    fprintf(stdout, "%u rays lost\n", n_lost);
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

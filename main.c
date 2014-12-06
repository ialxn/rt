/*	main.c
 *
 * Copyright (C) 2010,2011,2012,2013,2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#define __STDC_FORMAT_MACROS

#include <getopt.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#include <libconfig.h>

#include <gsl/gsl_rng.h>

#include "obj_lists.h"
#include "off.h"

#define CHECK_CONFIG 0
#define PRINT_GEOMETRY 1
#define RUN 2

struct run_simulation_args {
    source_list_t *source_list;	/* list of all sources */
    target_list_t *target_list;	/* list of all targets */
    int seed_base;
    int seed_incr;
};

struct worker_retval {
    int64_t n_total;		/* total number of rays traced */
    int64_t n_lost;		/* number of rays lost i.e. rays that
				   are not absorbed anywhere. they may
				   have hit some targets and have been
				   reflected several times, however. */
    double p_total;		/* total power of all rays */
    double p_lost;		/* total power of all lost rays */
};

static pthread_mutex_t mutex_seed_incr = PTHREAD_MUTEX_INITIALIZER;

static pthread_mutex_t mutex_print1 = PTHREAD_MUTEX_INITIALIZER;
static int print1 = 0;

static pthread_mutex_t mutex_print2 = PTHREAD_MUTEX_INITIALIZER;
static int print2 = 0;

static int n_threads = 1;	/* default single threaded */
static int n_log = 0;		/* default do not log path of rays */
static int log_OFF;		/* TRUE: log OFF file (FALSE: raw) */


static void init_PTD(source_list_t * source_list,
		     target_list_t * target_list)
{
    struct list_head *s_pos;
    struct list_head *t_pos;

    list_for_each(t_pos, &(target_list->list)) {
	target_list_t *this_t = list_entry(t_pos, target_list_t, list);
	target_t *current_target = this_t->t;

	init_PTDT(current_target);
    }

    list_for_each(s_pos, &(source_list->list)) {
	source_list_t *this_s = list_entry(s_pos, source_list_t, list);
	source_t *current_source = this_s->s;

	init_rays_remain(current_source);
    }
}

static void flush_outbufs(target_list_t * target_list)
{
    struct list_head *t_pos;

    list_for_each(t_pos, &(target_list->list)) {
	target_list_t *this_t = list_entry(t_pos, target_list_t, list);
	target_t *current_target = this_t->t;

	flush_PTDT_outbuf(current_target);
    }
}

static void print1_once(const char *rng_name)
{
    pthread_mutex_lock(&mutex_print1);

    if (!print1)		/* only true during first call */
	fprintf(stdout,
		"    using random number generator %s from Gnu Scientific Library\n",
		rng_name);

    print1++;
    pthread_mutex_unlock(&mutex_print1);
}

static void print2_once(const char *s_type, const char *s_name, int64_t n,
			double power)
{
    pthread_mutex_lock(&mutex_print2);

    if (print2 % n_threads == 0) {	/* only true during first call */
	fprintf(stdout, "        %s (%s) started\n", s_name, s_type);
	fprintf(stdout, "            %" PRId64 " rays to trace\n", n);
	fprintf(stdout, "            with total power of %e\n", power);
	fflush(stdout);
    }

    print2++;
    pthread_mutex_unlock(&mutex_print2);
}


static void *run_simulation(void *args)
{
    struct run_simulation_args *a = (struct run_simulation_args *) args;

    int this_source = 0;
    struct list_head *s_pos;
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    struct worker_retval *retval =
	(struct worker_retval *) malloc(sizeof(struct worker_retval));

    retval->n_total = 0;
    retval->n_lost = 0;
    retval->p_total = 0.0;
    retval->p_lost = 0.0;

    pthread_mutex_lock(&mutex_seed_incr);	/* begin critical section */

    gsl_rng_set(r, (unsigned long int) abs(a->seed_base + a->seed_incr));
    a->seed_incr++;

    pthread_mutex_unlock(&mutex_seed_incr);	/* end critical section */

    init_PTD(a->source_list, a->target_list);	/* initialize per thread data */

    print1_once(gsl_rng_name(r));

    list_for_each(s_pos, &(a->source_list->list)) {
	/*
	 * iterate through all sources
	 */
	source_list_t *this_s = list_entry(s_pos, source_list_t, list);
	source_t *current_source = this_s->s;
	ray_t *ray;
	int this_source_n_log = n_log;
	FILE *log_file;

	if (n_log) {		/* log ray path enabled */
	    char fname[256];

	    if (log_OFF) {
		snprintf(fname, 250, "%03d", this_source);
		log_file = open_off(fname);
		fprintf(log_file, "LIST\n");
	    } else {
		snprintf(fname, 250, "%03d.dat", this_source);
		log_file = fopen(fname, "w");
	    }
	}

	print2_once(get_source_type(current_source),
		    get_source_name(current_source),
		    get_source_n_rays(current_source),
		    get_source_power(current_source));

	while ((ray = emit_ray(current_source, r))) {
	    /*
	     * loop until 'current_source' is exhausted indicated
	     * by 'new_ray()' returning 'NULL'.
	     */
	    if (n_log)		/* log ray path enabled */
		if (!this_source_n_log) {
		    /*
		     * sufficient ray logged, stop this source
		     */
		    free(ray);
		    break;
		}

	    retval->n_total++;
	    retval->p_total += ray->power;

	    while (ray) {
		/*
		 * loop until 'ray' is absorbed or leaves system
		 */
		struct list_head *t_pos;
		/*
		 * keep track of target closest to origin of the
		 * current ray that is intercepted by it.
		 * 'closest_icpt' is later used in the calculation
		 * of 'ray' leaving 'closest_target' in 'out_ray()' as
		 * origin of the new ray.
		 */
		target_t *closest_target = NULL;
		double *closest_icpt = NULL;
		double min_d_sqrd = GSL_DBL_MAX;

		list_for_each(t_pos, &(a->target_list->list)) {
		    /*
		     * iterate through all targets.
		     * calculate point of interception on each target
		     * to identify target closest to origin of 'ray'.
		     * this is the actual target that is hit by 'ray'.
		     */
		    target_list_t *this_t =
			list_entry(t_pos, target_list_t, list);
		    target_t *current_target = this_t->t;
		    double *current_icpt = icpt(current_target, ray);

		    if (current_icpt) {
			/*
			 * 'ray' is intercepted by 'current_target'
			 */
			const double d_sqrd =
			    d_sqr(current_icpt, ray->orig);

			if (d_sqrd < min_d_sqrd) {
			    /*
			     * 'current_target' is target closest
			     * to origin of 'ray' found so far.
			     * no need to take sqrt:
			     *   if d1 > d2 then also d1^2 > d2^2
			     */

			    free(closest_icpt);
			    closest_target = current_target;
			    closest_icpt = current_icpt;
			    min_d_sqrd = d_sqrd;

			} else
			    /*
			     * 'ray' has been intercepted by a target
			     * placed at a larger distance to origin of
			     * 'ray' that 'closest_target'.
			     */
			    free(current_icpt);

		    }		/* end 'if(current_icpt)' */
		}		/* all targets tried */

		if (closest_icpt) {
		    /*
		     * 'ray' has been intercepted by a target.
		     * update 'ray'.
		     * Note: 'out_ray' returns NULL if 'ray' is
		     *        absorbed by target. this will terminate
		     *        the "while (ray) {}" loop and a new
		     *        ray will be emitted by the current source.
		     */
#define LOG_RAY(FLAG,FILE,O,P,R,G,B) do { \
if(FLAG) \
    write_ray(FILE,O,P,R,G,B); \
else \
    fprintf(FILE,"%e\t%e\t%e\t%e\t%e\t%e\n",O[0],O[1],O[2],P[0],P[1],P[2]); \
} while(0);
		    if (n_log)	/* log ray path enabled */
			LOG_RAY(log_OFF, log_file, ray->orig, closest_icpt,
				0, 0, 0);

		    ray = out_ray(closest_target, ray, closest_icpt, r);
		    free(closest_icpt);
		    closest_icpt = NULL;

		} else {
		    /*
		     * 'ray' has not been intercepted by any target.
		     * thus, 'ray' is lost and terminates.
		     * trigger emmision of new ray from current source
		     * by assigning NULL to 'ray'
		     */
		    if (n_log) {	/* log ray path enabled */
			double end[3];

			a_plus_cb(end, ray->orig, 1.0, ray->dir);
			/*
			 * this is a ray to infinity. mark red.
			 */
			LOG_RAY(log_OFF, log_file, ray->orig, end, 1, 0,
				0);
		    }

		    retval->n_lost++;
		    retval->p_lost += ray->power;
		    free(ray);
		    ray = NULL;

		}
	    }			/* 'ray' absorbed or lost */

	    if (n_log) {	/* log ray path enabled */
		this_source_n_log--;
		if (!log_OFF)
		    fprintf(log_file, "\n");
	    }

	}			/* 'current_source' is exhausted */

	if (n_log) {		/* log ray path enabled */
	    fclose(log_file);
	    this_source++;
	}

    }				/* all sources exhausted */
    gsl_rng_free(r);
    flush_outbufs(a->target_list);

    return retval;
}

static void help(void)
{
    fprintf(stdout, "\nrt Version: %s  %s\n\n", RELEASE, RELEASE_INFO);
    fprintf(stdout, "Usage: rt\n");
    fprintf(stdout,
	    "       --append, -a      append to output files. new seed must be given.\n");
    fprintf(stdout,
	    "       --Log, -L         log path of first n rays per source [0].\n");
    fprintf(stdout, "                         Raw format used.\n");
    fprintf(stdout,
	    "       --log, -l         log path of first l rays per source [0].\n");
    fprintf(stdout, "                         OFF format used.\n");
    fprintf(stdout, "       --mode, -m        select run mode [0].\n");
    fprintf(stdout,
	    "                         0: check and print input.\n");
    fprintf(stdout,
	    "                         1: output geometry (OFF files).\n");
    fprintf(stdout, "                         2: run simulation.\n");
    fprintf(stdout,
	    "       --threads, -t     number of threads to use [1]\n");
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
    int file_mode = O_TRUNC;	/* default file mode for output per target */

    struct run_simulation_args *rs_args;	/* arguments for worker threads */
    pthread_t *tids;		/* vector with thread ids */
    pthread_attr_t attr;

    int i;

    while (1) {
	int c;
	int option_index = 0;
	static struct option long_options[] = {
	    {"append", required_argument, 0, 'a'},
	    {"Log", required_argument, 0, 'L'},
	    {"log", required_argument, 0, 'l'},
	    {"mode", required_argument, 0, 'm'},
	    {"threads", required_argument, 0, 't'},
	    {"Version", no_argument, 0, 'V'},
	    {"help", no_argument, 0, 'h'},
	    {0, 0, 0, 0}
	};

	c = getopt_long(argc, argv, "a:L:l:m:t:Vh", long_options,
			&option_index);

	if (c == -1)
	    break;

	switch (c) {

	case 'a':
	    file_mode = O_APPEND;
	    seed = atoi(optarg);
	    break;

	case 'L':
	    n_log = atoi(optarg);
	    log_OFF = 0;
	    break;

	case 'l':
	    n_log = atoi(optarg);
	    log_OFF = 1;
	    break;

	case 'm':
	    mode = atoi(optarg);
	    if ((mode > RUN) || (mode < CHECK_CONFIG)) {
		fprintf(stderr, "unknown 'mode' (%d) given.\n", mode);
		exit(EXIT_FAILURE);
	    }
	    break;

	case 't':
	    n_threads = atoi(optarg);
	    break;

	case 'V':
	    fprintf(stdout, "rt Version: %s  %s\n", RELEASE, RELEASE_INFO);
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
	struct worker_retval *retval;
	int64_t n_total;
	int64_t n_lost;
	double p_total;
	double p_lost;

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
	fprintf(stdout, "rt version %s running ...\n", RELEASE);

	source_list = init_sources(&cfg, &n_sources);
	fprintf(stdout, "    %d sources initialized\n", n_sources);
	target_list = init_targets(&cfg, &n_targets, file_mode);
	fprintf(stdout, "    %d targets initialized\n", n_targets);

	if (file_mode == O_TRUNC) {	/* use seed from cfg, otherwise from command line */
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

	if (n_threads > 1) {
	    if (!n_log)
		fprintf(stdout,
			"    starting %d threads (seed %d, %d, ...) \n",
			n_threads, seed, seed + 1);
	    else {
		n_threads = 1;
		fprintf(stdout,
			"    logging path of rays requires single threaded execution.\n");
	    }
	}

	config_destroy(&cfg);

	rs_args = (struct run_simulation_args *)
	    malloc(sizeof(struct run_simulation_args));
	rs_args->source_list = source_list;
	rs_args->target_list = target_list;
	rs_args->seed_base = seed;
	rs_args->seed_incr = 0;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	tids =
	    (pthread_t *) malloc((size_t) n_threads * sizeof(pthread_t *));

	for (i = 0; i < n_threads; i++)
	    pthread_create(&(tids[i]), &attr, run_simulation,
			   (void *) rs_args);

	/* Wait for the other threads to finish */
	for (i = 0, n_total = 0, n_lost = 0, p_total = 0.0, p_lost = 0.0;
	     i < n_threads; i++) {
	    pthread_join(tids[i], (void **) &retval);
	    n_total += retval->n_total;
	    n_lost += retval->n_lost;
	    p_total += retval->p_total;
	    p_lost += retval->p_lost;
	    free(retval);
	}

	fprintf(stdout, "  total number of rays traced: %" PRId64 "\n",
		n_total);
	fprintf(stdout, "   total power of all sources: %e\n", p_total);
	fprintf(stdout, "    total number of rays lost: %" PRId64 "\n",
		n_lost);
	fprintf(stdout, "             total power lost: %e\n", p_lost);
	fprintf(stdout, "total number of rays absorbed: %" PRId64 "\n",
		n_total - n_lost);
	fprintf(stdout, "         total power absorbed: %e\n",
		p_total - p_lost);

	free(tids);
	free(rs_args);
	source_list_free(source_list);
	target_list_free(target_list);
	free(source_list);
	free(target_list);

    default:
	break;

    }
    return EXIT_SUCCESS;
}

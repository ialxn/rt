/*	spectrum.c
 *
 * Copyright (C) 2011 Ivo Alxneit
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

#include "io_util.h"
#include "version.h"



static int get_idx(FILE * f_in, int *idx_lambda, int *idx_power)
{
    char line[LINE_LEN];
    int status = NO_ERR;

    fgets(line, LINE_LEN, f_in);

    /* check file type and set indices of columns to be read */
    if (strcmp(line, "(one-sided plane screen)")) {
	*idx_lambda = 3;
	*idx_power = 2;
    } else if (strcmp(line, "(two-sided plane screen)")) {
	*idx_lambda = 3;
	*idx_power = 2;
    } else if (strcmp(line, "(rectangle)")) {
	*idx_lambda = 3;
	*idx_power = 2;
    } else if (strcmp(line, "(triangle)")) {
	*idx_lambda = 3;
	*idx_power = 2;
    } else if (strcmp(line, "(ellipsoid)")) {
	*idx_lambda = 4;
	*idx_power = 3;
    } else {
	fprintf(stderr, "Unknown target type (%s) found\n", line);
	status = ERR;
    }
    return status;
}


static int read_hist(FILE * f_in, gsl_histogram * h, int *n_missed,
		     double *p_missed, const int idx_lambda,
		     const int idx_power)
{
    char line[LINE_LEN];

    while (fgets(line, LINE_LEN, f_in)) {	/* read all lines */
	if (strncmp(line, "#", 1) != 0) {	/* not a comment */
#define MAX_ITEMS 5		/* maximum number of items per line.
				   if changed, also change sscanf() statement below */

	    double t[MAX_ITEMS];
	    int n_items;

	    n_items =
		sscanf(line, "%lf%lf%lf%lf%lf", &t[0], &t[1], &t[2], &t[3],
		       &t[4]);

	    if (n_items < idx_lambda + 1)	/* insufficient data read */
		break;

	    if (gsl_histogram_accumulate(h, t[idx_lambda], t[idx_power])) {
		/* data lies outside of allowed range */
		(*n_missed)++;
		*p_missed += t[idx_power];
	    }
	}
    }
    return NO_ERR;
}


static void output_hist(FILE * f_out, gsl_histogram * h,
			const int n_missed, const double p_missed)
{
    size_t i;
    double t;
    double min, max;
    const size_t n = gsl_histogram_bins(h);

    /*
     * spectrum is power per wavelength interval.
     * divide by bin width. all bins are identical so use bin 0.
     */
    gsl_histogram_get_range(h, 0, &min, &max);
    t = max - min;
    gsl_histogram_scale(h, 1.0 / t);

    /* print header with statistics */
    fprintf(f_out, "#   histogram definition\n");
    t = gsl_histogram_min(h);
    fprintf(f_out, "#      minimum x-value: %e\n", t);
    t = gsl_histogram_max(h);
    fprintf(f_out, "#      maximum x-value: %e\n", t);
    i = gsl_histogram_bins(h);
    fprintf(f_out, "#       number of bins: %d\n", i);

    fprintf(stdout, "#\n#   histogram statistics\n");

    fprintf(stdout, "#       number of data point not included: %d\n",
	    n_missed);
    fprintf(stdout, "#                      total power missed: %e\n#\n",
	    p_missed);


    t = gsl_histogram_min_val(h);
    i = gsl_histogram_min_bin(h);
    gsl_histogram_get_range(h, i, &min, &max);
    fprintf(stdout, "#        minimum y-value: %e\n", t);
    fprintf(stdout, "#            at bin number (range): %d (%e -- %e)\n",
	    i, min, max);

    t = gsl_histogram_max_val(h);
    i = gsl_histogram_max_bin(h);
    gsl_histogram_get_range(h, i, &min, &max);
    fprintf(stdout, "#        maximum y-value: %e\n", t);
    fprintf(stdout, "#            at bin number (range): %d (%e -- %e)\n",
	    i, min, max);

    fprintf(stdout, "# mean y-value (+-sigma): %e (+-%e)\n",
	    gsl_histogram_mean(h), gsl_histogram_sigma(h));
    fprintf(stdout, "#\n");

    fprintf(stdout, "# x\t\ty\t\tbin_min\t\tbin_max\n");

    for (i = 0; i < n; i++) {
	double t_x;

	gsl_histogram_get_range(h, i, &min, &max);
	t_x = (min + max) / 2.0;
	t = gsl_histogram_get(h, i);

	fprintf(f_out, "%e\t%e\t%e\t%e\n", t_x, t, min, max);
    }
}


static gsl_histogram *init_hist(const double start_wl,
				const double stop_wl, const size_t n_bins)
{
    gsl_histogram *h;
    double *range;
    double width;
    size_t i;

    width = (stop_wl - start_wl) / n_bins;

    range = (double *) malloc((n_bins + 1) * sizeof(double));

    range[0] = start_wl;
    for (i = 1; i <= n_bins; i++)
	range[i] = range[i - 1] + width;

    h = gsl_histogram_alloc(n_bins);
    gsl_histogram_set_ranges(h, range, n_bins + 1);

    free(range);

    return h;
}


static void help(void)
{
    fprintf(stdout, "\nspectrum Version %s (AI52)\n\n", VERSION);
    fprintf(stdout, "Usage: spectrum\n");
    fprintf(stdout, "       --num, -n         number of bins [10]\n");
    fprintf(stdout, "       --start, -a       Start wavelength [0.0]\n");
    fprintf(stdout, "       --stop, -o        Stop wavelength [1000.0]\n");
    fprintf(stdout, "       --help, -h        Print this help message\n");
    fprintf(stdout, "       --Version, -V     Print version number\n");
    fprintf(stdout, "\n");
}

int main(int argc, char **argv)
{
    size_t n_bins = 10;
    double start_wl = 0.0;
    double stop_wl = 1000.0;
    gsl_histogram *h;
    int idx_p, idx_l;
    int n_missed = 0;
    double p_missed = 0.0;

    while (1) {
	int c;
	int option_index = 0;
	static struct option long_options[] = {
	    {"num", required_argument, 0, 'n'},
	    {"start", required_argument, 0, 'a'},
	    {"stop", required_argument, 0, 'o'},
	    {"Version", no_argument, 0, 'V'},
	    {"help", no_argument, 0, 'h'},
	    {0, 0, 0, 0}
	};

	c = getopt_long(argc, argv, "a:n:o:Vh", long_options,
			&option_index);

	if (c == -1)
	    break;

	switch (c) {

	case 'a':
	    start_wl = atof(optarg);
	    break;

	case 'n':
	    n_bins = (size_t) atoi(optarg);
	    break;

	case 'o':
	    stop_wl = atof(optarg);
	    break;

	case 'V':
	    fprintf(stdout, " spectrum version %s (AI52)\n", VERSION);
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

    if (get_idx(stdin, &idx_l, &idx_p))
	exit(EXIT_FAILURE);

    h = init_hist(start_wl, stop_wl, n_bins);

    read_hist(stdin, h, &n_missed, &p_missed, idx_l, idx_p);
    output_hist(stdout, h, n_missed, p_missed);

    gsl_histogram_free(h);

    exit(EXIT_SUCCESS);
}

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


static int read_hist(FILE * f_in, gsl_histogram * h, int *n_inc,
		     double *p_inc, int *n_missed,
		     double *p_missed, const size_t idx_lambda,
		     const size_t idx_power)
{
    float t[MAX_ITEMS];
    size_t n_items_read;

    if (skip_header(f_in) == ERR)
	return ERR;

    /*
     * read data. (x,y,[z,]power,lambda)
     */
    n_items_read = fread(t, sizeof(float), idx_lambda + 1, f_in);

    if (!n_items_read) {
	fprintf(stderr, "No data found\n");
	return (ERR);
    }

    while (n_items_read) {

	if (n_items_read < idx_lambda + 1) {	/* insufficient data read */
	    fprintf(stderr,
		    "Incomplete data set read (%d instead of %d)\n",
		    n_items_read, idx_lambda + 1);
	    return (ERR);
	}

	if (gsl_histogram_accumulate(h, t[idx_lambda], t[idx_power])) {
	    /* data lies outside of range of histogram */
	    (*n_missed)++;
	    *p_missed += t[idx_power];
	} else {
	    (*n_inc)++;
	    *p_inc += t[idx_power];
	}

	n_items_read = fread(t, sizeof(float), idx_lambda + 1, f_in);

    }

    return NO_ERR;
}


static void output_hist(FILE * f_out, gsl_histogram * h, const int n_inc,
			const double p_inc, const int n_missed,
			const double p_missed)
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
    fprintf(f_out, "#      minimum x-value: %.0f\n", t);
    t = gsl_histogram_max(h);
    fprintf(f_out, "#      maximum x-value: %.0f\n", t);
    i = gsl_histogram_bins(h);
    fprintf(f_out, "#       number of bins: %d\n", i);

    fprintf(f_out, "#\n#   histogram statistics\n");
    fprintf(f_out, "#      number of data points not included: %d\n",
	    n_missed);
    fprintf(f_out, "#                      total power missed: %e\n",
	    p_missed);
    fprintf(f_out, "#          number of data points included: %d\n",
	    n_inc);
    fprintf(f_out, "#               total power accounted for: %e\n#\n",
	    p_inc);

    t = gsl_histogram_min_val(h);
    i = gsl_histogram_min_bin(h);
    gsl_histogram_get_range(h, i, &min, &max);
    fprintf(f_out, "#        minimum y-value: %e\n", t);
    fprintf(f_out,
	    "#            at bin number (range): %d (%.0f - %.0f)\n", i,
	    min, max);

    t = gsl_histogram_max_val(h);
    i = gsl_histogram_max_bin(h);
    gsl_histogram_get_range(h, i, &min, &max);
    fprintf(f_out, "#        maximum y-value: %e\n", t);
    fprintf(f_out,
	    "#            at bin number (range): %d (%.0f - %.0f)\n", i,
	    min, max);

    fprintf(f_out, "# mean y-value (+-sigma): %e (+-%e)\n",
	    gsl_histogram_mean(h), gsl_histogram_sigma(h));
    fprintf(f_out, "#\n");

    fprintf(f_out, "# x\t\ty\tbin_min\tbin_max\n");

    for (i = 0; i < n; i++) {
	double t_x;

	gsl_histogram_get_range(h, i, &min, &max);
	t_x = (min + max) / 2.0;
	t = gsl_histogram_get(h, i);

	fprintf(f_out, "%.1f\t%e\t%.0f\t%.0f\n", t_x, t, min, max);
    }
}


static gsl_histogram *init_hist(const double start_wl,
				const double stop_wl, const size_t n_bins)
{
    gsl_histogram *h;
    double *wl_range = range(start_wl, stop_wl, n_bins);

    h = gsl_histogram_alloc(n_bins);
    gsl_histogram_set_ranges(h, wl_range, n_bins + 1);

    free(wl_range);

    return h;
}


static void help(void)
{
    fprintf(stdout, "\nrt Version: %s  %s\n\n", RELEASE, RELEASE_INFO);
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
    size_t idx_p, idx_l;
    int n_missed = 0;		/* data not included in histogram */
    double p_missed = 0.0;
    int n_inc = 0;		/* data included in histogram */
    double p_inc = 0.0;

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
	    fprintf(stdout, "spectrum Version: %s  %s\n", RELEASE, RELEASE_INFO);	    exit(EXIT_SUCCESS);
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

    h = init_hist(start_wl, stop_wl, n_bins);

    if (read_hist
	(stdin, h, &n_inc, &p_inc, &n_missed, &p_missed, idx_l,
	 idx_p) == NO_ERR)
	output_hist(stdout, h, n_inc, p_inc, n_missed, p_missed);

    gsl_histogram_free(h);

    exit(EXIT_SUCCESS);
}

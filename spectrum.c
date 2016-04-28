/*	spectrum.c
 *
 * Copyright (C) 2011,2012,2013,2014,2015,2016 Ivo Alxneit
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
		     double *p_inc, int *n_missed, double *p_missed,
		     const size_t idx_lambda)
{
    float t[MAX_FLOAT_ITEMS];
    unsigned char tmp_8;
    size_t n_float_items_read, n_uchar_items_read;
    double P_factor;

    P_factor = get_P_factor(f_in);
    if (P_factor < 0.0)
	return ERR;

    if (skip_N_comments(f_in, REST_HEADER_LINES) == ERR)
	return ERR;

    /*
     * read data. (x,y,[z,]lambda)
     */
    n_float_items_read = fread(t, sizeof(float), idx_lambda + 1, f_in);
    if (!n_float_items_read) {
	fprintf(stderr, "No data found\n");
	return (ERR);
    }
    n_uchar_items_read = fread(&tmp_8, sizeof(unsigned char), 1, f_in);

    while (n_float_items_read && n_uchar_items_read) {

	if ((n_float_items_read < idx_lambda + 1) || !n_uchar_items_read) {	/* insufficient data read */

#if __STDC_VERSION__ >= 199901L
	    fprintf(stderr,
		    "Incomplete data set read (%zu floats instead of %zu and %zu chars instead of 1)\n",
		    n_float_items_read, idx_lambda + 1,
		    n_uchar_items_read);
#else
	    fprintf(stderr,
		    "Incomplete data set read (%lu floats instead of %lu and %lu chars instead of 1)\n",
		    (unsigned long) n_float_items_read,
		    (unsigned long) idx_lambda + 1,
		    (unsigned long) n_uchar_items_read);
#endif
	    return (ERR);
	}

	if (gsl_histogram_accumulate(h, t[idx_lambda], P_factor)) {
	    /* data lies outside of range of histogram */
	    (*n_missed)++;
	    *p_missed += P_factor;
	} else {
	    (*n_inc)++;
	    *p_inc += P_factor;
	}

	n_float_items_read = fread(t, sizeof(float), idx_lambda + 1, f_in);
	n_uchar_items_read = fread(&tmp_8, sizeof(unsigned char), 1, f_in);

    }

    return NO_ERR;
}


static void output_hist(FILE * f_out, gsl_histogram * h, const int n_inc,
			const double p_inc, const int n_missed,
			const double p_missed)
{
    size_t i;
    double min, max;
    const size_t n = gsl_histogram_bins(h);

    /*
     * spectrum is power per wavelength interval.
     * divide by bin width. all bins are identical so use bin 0.
     */
    gsl_histogram_get_range(h, 0, &min, &max);
    gsl_histogram_scale(h, 1.0 / (max - min));

    /* print header with statistics */
    fprintf(f_out, "#   histogram definition\n");
    fprintf(f_out, "#      minimum x-value: %.0f\n", gsl_histogram_min(h));
    fprintf(f_out, "#      maximum x-value: %.0f\n", gsl_histogram_max(h));
#if __STDC_VERSION__ >= 199901L
    fprintf(f_out, "#       number of bins: %zu\n", gsl_histogram_bins(h));
#else
    fprintf(f_out, "#       number of bins: %lu\n",
	    (unsigned long) gsl_histogram_bins(h));
#endif

    fprintf(f_out, "#\n#   histogram statistics\n");
    fprintf(f_out, "#      number of data points not included: %d\n",
	    n_missed);
    fprintf(f_out, "#                      total power missed: %e\n",
	    p_missed);
    fprintf(f_out, "#          number of data points included: %d\n",
	    n_inc);
    fprintf(f_out, "#               total power accounted for: %e\n",
	    p_inc);

    i = gsl_histogram_min_bin(h);
    gsl_histogram_get_range(h, i, &min, &max);
    fprintf(f_out, "#        minimum y-value: %e\n",
	    gsl_histogram_min_val(h));
#if __STDC_VERSION__ >= 199901L
    fprintf(f_out,
	    "#            at bin number (range): %zu (%.0f - %.0f)\n", i,
	    min, max);
#else
    fprintf(f_out,
	    "#            at bin number (range): %lu (%.0f - %.0f)\n",
	    (unsigned long) i, min, max);
#endif

    i = gsl_histogram_max_bin(h);
    gsl_histogram_get_range(h, i, &min, &max);
    fprintf(f_out, "#        maximum y-value: %e\n",
	    gsl_histogram_max_val(h));
#if __STDC_VERSION__ >= 199901L
    fprintf(f_out,
	    "#            at bin number (range): %zu (%.0f - %.0f)\n", i,
	    min, max);
#else
    fprintf(f_out,
	    "#            at bin number (range): %lu (%.0f - %.0f)\n",
	    (unsigned long) i, min, max);
#endif

    fprintf(f_out, "# mean y-value (+-sigma): %e (+-%e)\n",
	    gsl_histogram_mean(h), gsl_histogram_sigma(h));
    fprintf(f_out, "#\n");

    fprintf(f_out, "# x\t\ty\tbin_min\tbin_max\n");

    for (i = 0; i < n; i++) {
	double t_x;

	gsl_histogram_get_range(h, i, &min, &max);
	t_x = (min + max) / 2.0;

	fprintf(f_out, "%.1f\t%e\t%.0f\t%.0f\n", t_x,
		gsl_histogram_get(h, i), min, max);
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
    size_t idx_l;
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
	    fprintf(stdout, "spectrum Version: %s  %s\n", RELEASE,
		    RELEASE_INFO);
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

    if (get_idx(stdin, &idx_l) == ERR)
	exit(EXIT_FAILURE);

    h = init_hist(start_wl, stop_wl, n_bins);

    if (read_hist(stdin, h, &n_inc, &p_inc, &n_missed, &p_missed, idx_l) ==
	NO_ERR)
	output_hist(stdout, h, n_inc, p_inc, n_missed, p_missed);

    gsl_histogram_free(h);

    exit(EXIT_SUCCESS);
}

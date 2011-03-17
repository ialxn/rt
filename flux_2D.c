/*	flux_2D.c
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

#include <gsl/gsl_histogram2d.h>

#include "io_util.h"
#include "version.h"

static int wrong_target_type(FILE * f_in)
{
    char line[LINE_LEN];
    int status = ERR;

    fgets(line, LINE_LEN, f_in);

    /*
     * check file type and set status
     * only planar targets are supported, so list them all in the
     * if-then-else ladder.
     */
    if (strstr(line, "(one-sided plane screen)"))
	status = NO_ERR;
    else if (strstr(line, "(two-sided plane screen)"))
	status = NO_ERR;
    else if (strstr(line, "(rectangle)"))
	status = NO_ERR;
    else if (strstr(line, "(triangle)"))
	status = NO_ERR;
    else
	fprintf(stderr, "Unknown target type %s\n", line);

    return status;
}


static double *range(const double min, const double max, const size_t n)
{
    size_t i;
    const double width = (max - min) / n;
    double *r = (double *) malloc((n + 1) * sizeof(double));

    r[0] = min;
    for (i = 1; i <= n; i++)
	r[i] = r[i - 1] + width;

    return r;
}

static gsl_histogram2d *init_hist(const double xmax,
				  const double xmin, const size_t x_bins,
				  const double ymax,
				  const double ymin, const size_t y_bins)
{
    gsl_histogram2d *h = gsl_histogram2d_alloc(x_bins, y_bins);
    double *xrange = range(xmin, xmax, x_bins);
    double *yrange = range(ymin, ymax, y_bins);

    gsl_histogram2d_set_ranges(h, xrange, x_bins + 1, yrange, y_bins + 1);

    free(xrange);
    free(yrange);

    return h;
}

static void read_hist(FILE * f_in, gsl_histogram2d * h, int *n_missed,
		      double *p_missed)
{
    char line[LINE_LEN];

    *n_missed = 0;
    *p_missed = 0.0;

    while (fgets(line, LINE_LEN, f_in)) {	/* read all lines */
	if (strncmp(line, "#", 1) != 0) {	/* not a comment */
#define INDX 0
#define INDY 1
#define INDP 2
#define MAX_ITEMS 4		/* maximum number of items per line.
				   if changed, also change sscanf() statement below */

	    double t[MAX_ITEMS];
	    int n_items;

	    n_items =
		sscanf(line, "%lf%lf%lf%lf", &t[0], &t[1], &t[2], &t[3]);

	    if (n_items < INDP + 1)	/* insufficient data read */
		break;

	    if (gsl_histogram2d_accumulate(h, t[INDX], t[INDY], t[INDP])) {
		/* data lies outside of allowed range */
		(*n_missed)++;
		*p_missed += t[INDP];
	    }
	}
    }
}

static void output_hist(FILE * f_out, gsl_histogram2d * h,
			const int n_missed, const double p_missed)
{
    size_t i, j;
    size_t nx, ny;
    double t;
    double A;
    double xmin, xmax;
    double ymin, ymax;

    /*
     * histogram should report flux i.e. power per unit surface area.
     * thus scale the histogram by the area A of one bin:
     *     A=a_width x b_width
     * all bins are identical, thus use (0,0)
     */
    gsl_histogram2d_get_xrange(h, 0, &xmin, &xmax);
    gsl_histogram2d_get_yrange(h, 0, &ymin, &ymax);
    A = (xmax - xmin) * (ymax - ymin);
    gsl_histogram2d_scale(h, 1.0 / A);

    fprintf(f_out, "#   histogram definition\n");
    t = gsl_histogram2d_xmin(h);
    fprintf(f_out, "#      minimum x-value: %e\n", t);
    t = gsl_histogram2d_xmax(h);
    fprintf(f_out, "#      maximum x-value: %e\n", t);
    nx = gsl_histogram2d_nx(h);
    fprintf(f_out, "#     number of x-bins: %d\n", nx);
    t = gsl_histogram2d_ymin(h);
    fprintf(f_out, "#      minimum y-value: %e\n", t);
    t = gsl_histogram2d_ymax(h);
    fprintf(f_out, "#      maximum y-value: %e\n", t);
    ny = gsl_histogram2d_ny(h);
    fprintf(f_out, "#     number of y-bins: %d\n", ny);


    fprintf(f_out, "#\n#   histogram statistics\n");

    fprintf(f_out, "#        number of data point not included: %d\n",
	    n_missed);
    fprintf(f_out, "#                       total power missed: %e\n#\n",
	    p_missed);
    fprintf(f_out, "#        area represented by one bin: %e\n#\n", A);

    t = gsl_histogram2d_min_val(h);
    gsl_histogram2d_min_bin(h, &i, &j);
    gsl_histogram2d_get_xrange(h, i, &xmin, &xmax);
    gsl_histogram2d_get_yrange(h, j, &ymin, &ymax);
    fprintf(f_out, "#        minimum value: %e\n", t);
    fprintf(f_out, "#            at bin (xrange): %d (%e -- %e)\n", i,
	    xmin, xmax);
    fprintf(f_out, "#            at bin (yrange): %d (%e -- %e)\n", j,
	    ymin, ymax);
    fprintf(f_out, "#\n");

    t = gsl_histogram2d_max_val(h);
    gsl_histogram2d_max_bin(h, &i, &j);
    gsl_histogram2d_get_xrange(h, i, &xmin, &xmax);
    gsl_histogram2d_get_yrange(h, j, &ymin, &ymax);
    fprintf(f_out, "#        maximum value: %e\n", t);
    fprintf(f_out, "#            at bin (xrange): %d (%e -- %e)\n", i,
	    xmin, xmax);
    fprintf(f_out, "#            at bin (yrange): %d (%e -- %e)\n", j,
	    ymin, ymax);
    fprintf(f_out, "#\n");

    fprintf(f_out, "# x mean value (+-sigma): %e (+-%e)\n",
	    gsl_histogram2d_xmean(h), gsl_histogram2d_xsigma(h));
    fprintf(f_out, "# y mean value (+-sigma): %e (+-%e)\n",
	    gsl_histogram2d_ymean(h), gsl_histogram2d_ysigma(h));
    fprintf(f_out, "#             covariance: %e\n",
	    gsl_histogram2d_cov(h));

    fprintf(f_out, "#\n");

    fprintf(f_out,
	    "#     x    \t      y    \t    Flux    \t  xbin_min  \t  xbin_max  \t  ybin_min  \t  ybin_max\n");

    for (i = 0; i < nx; i++) {
	double x, x_lo, x_hi;

	gsl_histogram2d_get_xrange(h, i, &x_lo, &x_hi);
	x = (x_hi + x_lo) / 2.0;

	for (j = 0; j < ny; j++) {
	    double y, y_lo, y_hi;

	    gsl_histogram2d_get_yrange(h, j, &y_lo, &y_hi);
	    y = (y_hi + y_lo) / 2.0;

	    t = gsl_histogram2d_get(h, i, j);
	    fprintf(f_out, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", x, y, t, x_lo,
		    x_hi, y_lo, y_hi);
	}
    }
}

static void help(void)
{
    fprintf(stdout, "\nflux_2D Version %s (AI52)\n\n", VERSION);
    fprintf(stdout, "Usage: flux_2D\n");
    fprintf(stdout, "       --nx, -a          number of x bins  [10]\n");
    fprintf(stdout, "       --ny, -b          number of y bins  [10]\n");
    fprintf(stdout, "       --minx, -x        Minimum x value  [-10.0]\n");
    fprintf(stdout, "       --maxx, -X        Maximum x value  [10.0]\n");
    fprintf(stdout, "       --miny, -y        Minimum y value  [-10.0]\n");
    fprintf(stdout, "       --maxy, -Y        Maximum y value  [10.0]\n");
    fprintf(stdout, "       --help, -h        Print this help message\n");
    fprintf(stdout, "       --Version, -V     Print version number\n");
    fprintf(stdout, "\n");
}

int main(int argc, char **argv)
{
    double x_max = 10.0;
    double x_min = -10.0;
    size_t x_bins = 10;
    double y_max = 10.0;
    double y_min = -10.0;
    size_t y_bins = 10;

    double p_missed;
    int n_missed;

    gsl_histogram2d *h;

    while (1) {
	int c;
	int option_index = 0;
	static struct option long_options[] = {
	    {"minx", required_argument, 0, 'x'},
	    {"maxx", required_argument, 0, 'X'},
	    {"miny", required_argument, 0, 'y'},
	    {"maxy", required_argument, 0, 'Y'},
	    {"nx", required_argument, 0, 'a'},
	    {"ny", required_argument, 0, 'b'},
	    {"Version", no_argument, 0, 'V'},
	    {"help", no_argument, 0, 'h'},
	    {0, 0, 0, 0}
	};

	c = getopt_long(argc, argv, "a:b:x:X:y:Y:Vh", long_options,
			&option_index);

	if (c == -1)
	    break;

	switch (c) {

	case 'a':
	    x_bins = (size_t) atoi(optarg);
	    break;

	case 'b':
	    y_bins = (size_t) atoi(optarg);
	    break;

	case 'x':
	    x_min = atof(optarg);
	    break;

	case 'X':
	    x_max = atof(optarg);
	    break;

	case 'y':
	    y_min = atof(optarg);
	    break;

	case 'Y':
	    y_max = atof(optarg);
	    break;

	case 'V':
	    fprintf(stdout, " flux_2D version %s (AI52)\n", VERSION);
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

    if (wrong_target_type(stdin))
	exit(EXIT_FAILURE);

    h = init_hist(x_max, x_min, x_bins, y_max, y_min, y_bins);

    read_hist(stdin, h, &n_missed, &p_missed);
    output_hist(stdout, h, n_missed, p_missed);

    gsl_histogram2d_free(h);
    exit(EXIT_SUCCESS);
}

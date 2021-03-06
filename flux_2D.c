/*	flux_2D.c
 *
 * Copyright (C) 2011 - 2018 Ivo Alxneit, Paul Scherrer Institute
 *
 * This file is part of rt
 *
 * rt is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rt is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rt. If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <getopt.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram2d.h>

#include "io_utils.h"
#include "math_utils.h"

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
    else if (strstr(line, "(annulus)"))
	status = NO_ERR;
    else if (strstr(line, "(disk)"))
	status = NO_ERR;
    else if (strstr(line, "(window)"))
	status = NO_ERR;
    else
	fprintf(stderr, "Unknown target type %s\n", line);

    return status;
}


static gsl_histogram2d *init_hist(const double xmax, const double xmin,
				  const size_t x_bins, const double ymax,
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

static int read_hist(FILE * f_in, gsl_histogram2d * h, const double l_min,
		     const double l_max, int *n_inc, double *p_inc,
		     int *n_missed, double *p_missed,
		     const int coordinates, double M[], double origin[])
{
    float t[MAX_FLOAT_ITEMS];
    unsigned char tmp_8;
    size_t n_float_items_read, n_uchar_items_read;
    double P_factor;

    P_factor = get_P_factor(f_in);
    if (P_factor < 0.0)
	return ERR;

    if (coordinates == LOCAL) {
	if (skip_N_comments(f_in, REST_HEADER_LINES) == ERR)
	    return ERR;
    } else
	read_transformation(f_in, M, origin);

    /*
     * read data. (x,y,lambda)
     */
    n_float_items_read =
	fread(t, sizeof(float), MAX_FLOAT_ITEMS - 1, f_in);

    if (!n_float_items_read) {
	fprintf(stderr, "No data found\n");
	return (ERR);
    }
    n_uchar_items_read = fread(&tmp_8, sizeof(unsigned char), 1, f_in);

    while (n_float_items_read && n_uchar_items_read) {
	if ((n_float_items_read < MAX_FLOAT_ITEMS - 1)
	    || !n_uchar_items_read) {
	    /*
	     * insufficient data read
	     */
#if __STDC_VERSION__ >= 199901L
	    fprintf(stderr,
		    "Incomplete data set read (%zu floats instead of %d and %zu chars instead of 1)\n",
		    n_float_items_read, MAX_FLOAT_ITEMS - 1,
		    n_uchar_items_read);
#else
	    fprintf(stderr,
		    "Incomplete data set read (%lu floats instead of %d and %lu chars instead of 1)\n",
		    (unsigned long) n_float_items_read,
		    MAX_FLOAT_ITEMS - 1,
		    (unsigned long) n_uchar_items_read);
#endif
	    return (ERR);
	}

	if ((t[2] >= l_min) && (t[2] <= l_max)) {
	    /*
	     * lambda is within selected intervall
	     */
	    if (gsl_histogram2d_accumulate(h, t[0], t[1], P_factor)) {
		/* data lies outside of range of histogram */
		(*n_missed)++;
		*p_missed += P_factor;
	    } else {
		(*n_inc)++;
		*p_inc += P_factor;
	    }
	} else {		/* lambda outside interval */
	    (*n_missed)++;
	    *p_missed += P_factor;
	}

	n_float_items_read =
	    fread(t, sizeof(float), MAX_FLOAT_ITEMS - 1, f_in);
	n_uchar_items_read = fread(&tmp_8, sizeof(unsigned char), 1, f_in);
    }

    return NO_ERR;

}

static void output_hist(FILE * f_out, gsl_histogram2d * h,
			const double l_min, const double l_max,
			const int n_inc, const double p_inc,
			const int n_missed, const double p_missed,
			const int coordinates, const double M[],
			const double origin[])
{
    size_t i, j;
    size_t nx, ny;
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

    fprintf(f_out, "#   histogram definition (local system)\n");
    fprintf(f_out,
	    "#      wavelength interval: %12.6g <= lambda <= %12.6g\n",
	    l_min, l_max);
    fprintf(f_out, "#          minimum x-value: %12.6g\n",
	    gsl_histogram2d_xmin(h));
    fprintf(f_out, "#          maximum x-value: %12.6g\n",
	    gsl_histogram2d_xmax(h));
    nx = gsl_histogram2d_nx(h);
#if __STDC_VERSION__ >= 199901L
    fprintf(f_out, "#         number of x-bins: %zu\n", nx);
#else
    fprintf(f_out, "#         number of x-bins: %lu\n",
	    (unsigned long) nx);
#endif
    fprintf(f_out, "#          minimum y-value: %12.6g\n",
	    gsl_histogram2d_ymin(h));
    fprintf(f_out, "#          maximum y-value: %12.6g\n",
	    gsl_histogram2d_ymax(h));
    ny = gsl_histogram2d_ny(h);
#if __STDC_VERSION__ >= 199901L
    fprintf(f_out, "#         number of y-bins: %zu\n", ny);
#else
    fprintf(f_out, "#         number of y-bins: %lu\n",
	    (unsigned long) ny);
#endif

    fprintf(f_out, "#\n#   histogram statistics\n");
    fprintf(f_out, "#        number of data point not included: %d\n",
	    n_missed);
    fprintf(f_out, "#                       total power missed: %12.6g\n",
	    p_missed);
    fprintf(f_out, "#           number of data points included: %d\n",
	    n_inc);
    fprintf(f_out,
	    "#                total power accounted for: %12.6g\n#\n",
	    p_inc);

    fprintf(f_out, "#        area represented by one bin: %12.6g\n#\n", A);

    gsl_histogram2d_min_bin(h, &i, &j);
    gsl_histogram2d_get_xrange(h, i, &xmin, &xmax);
    gsl_histogram2d_get_yrange(h, j, &ymin, &ymax);
    fprintf(f_out, "#        minimum value: %12.6g\n",
	    gsl_histogram2d_min_val(h));
#if __STDC_VERSION__ >= 199901L
    fprintf(f_out, "#            at bin (xrange): %zu (%12.6g - %12.6g)\n",
	    i, xmin, xmax);
    fprintf(f_out, "#            at bin (yrange): %zu (%12.6g - %12.6g)\n",
	    j, ymin, ymax);
#else
    fprintf(f_out, "#            at bin (xrange): %lu (%12.6g - %12.6g)\n",
	    (unsigned long) i, xmin, xmax);
    fprintf(f_out, "#            at bin (yrange): %lu (%12.6g - %12.6g)\n",
	    (unsigned long) j, ymin, ymax);
#endif
    fprintf(f_out, "#\n");

    gsl_histogram2d_max_bin(h, &i, &j);
    gsl_histogram2d_get_xrange(h, i, &xmin, &xmax);
    gsl_histogram2d_get_yrange(h, j, &ymin, &ymax);
    fprintf(f_out, "#        maximum value: %12.6g\n",
	    gsl_histogram2d_max_val(h));
#if __STDC_VERSION__ >= 199901L
    fprintf(f_out, "#            at bin (xrange): %zu (%12.6g - %12.6g)\n",
	    i, xmin, xmax);
    fprintf(f_out, "#            at bin (yrange): %zu (%12.6g - %12.6g)\n",
	    j, ymin, ymax);
#else
    fprintf(f_out, "#            at bin (xrange): %lu (%12.6g - %12.6g)\n",
	    (unsigned long) i, xmin, xmax);
    fprintf(f_out, "#            at bin (yrange): %lu (%12.6g - %12.6g)\n",
	    (unsigned long) j, ymin, ymax);
#endif
    fprintf(f_out, "#\n");

    fprintf(f_out, "# x mean value (+-sigma): %12.6g (+-%12.6g)\n",
	    gsl_histogram2d_xmean(h), gsl_histogram2d_xsigma(h));
    fprintf(f_out, "# y mean value (+-sigma): %12.6g (+-%12.6g)\n",
	    gsl_histogram2d_ymean(h), gsl_histogram2d_ysigma(h));
    fprintf(f_out, "#             covariance: %12.6g\n",
	    gsl_histogram2d_cov(h));

    fprintf(f_out, "#\n");

    if (coordinates == GLOBAL) {
	fprintf(f_out, "# coordinates are reported in global system\n");
	fprintf(f_out,
		"#     x    \t      y    \t      z    \t    Flux\n");
    } else {
	fprintf(f_out, "# coordinates are reported in local system\n");
	fprintf(f_out, "#     x    \t      y    \t    Flux    \t  ");
	fprintf(f_out,
		"xbin_min  \t  xbin_max  \t  ybin_min  \t  ybin_max\n");
    }

    for (i = 0; i < nx; i++) {
	double x, x_lo, x_hi;

	gsl_histogram2d_get_xrange(h, i, &x_lo, &x_hi);
	x = (x_hi + x_lo) / 2.0;

	for (j = 0; j < ny; j++) {
	    double y, y_lo, y_hi, val;

	    gsl_histogram2d_get_yrange(h, j, &y_lo, &y_hi);
	    y = (y_hi + y_lo) / 2.0;

	    val = gsl_histogram2d_get(h, i, j);

	    if (coordinates == GLOBAL) {
		double l_xyz[3];
		double g_xyz[3];

		l_xyz[0] = x;
		l_xyz[1] = y;
		l_xyz[2] = 0.0;

		l2g(M, origin, l_xyz, g_xyz);
		fprintf(f_out, "%12.6g\t%12.6g\t%12.6g\t%12.6g\n",
			g_xyz[0], g_xyz[1], g_xyz[2], val);

	    } else
		fprintf(f_out,
			"%12.6g\t%12.6g\t%12.6g\t%12.6g\t%12.6g\t%12.6g\t%12.6g\n",
			x, y, val, x_lo, x_hi, y_lo, y_hi);
	}
    }
}

static void help(void)
{
    fprintf(stdout, "\nflux_2D Version: %s  %s\n", RELEASE, RELEASE_INFO);
    fprintf(stdout, "Usage: flux_2D\n");
    fprintf(stdout, "       --nx, -a          number of x bins  [10]\n");
    fprintf(stdout, "       --ny, -b          number of y bins  [10]\n");
    fprintf(stdout,
	    "       --global, -g      report flux distribution in global\n");
    fprintf(stdout,
	    "                         coordinate system [local]\n");
    fprintf(stdout,
	    "       --minl, -l        Minimum wavelength  [0.0]\n");
    fprintf(stdout,
	    "       --maxl, -L        Maximum wavelength  [1E308]\n");
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
    double lambda_max = GSL_DBL_MAX;
    double lambda_min = 0.0;
    double x_max = 10.0;
    double x_min = -10.0;
    size_t x_bins = 10;
    double y_max = 10.0;
    double y_min = -10.0;
    size_t y_bins = 10;

    int n_inc = 0;		/* data included in histogram */
    double p_inc = 0.0;
    int n_missed = 0;		/* data not included in histogram */
    double p_missed = 0.0;
    int coordinates = LOCAL;
    double M[9];
    double origin[3];

    gsl_histogram2d *h;

    while (1) {
	int c;
	int option_index = 0;
	static struct option long_options[] = {
	    {"minl", required_argument, 0, 'l'},
	    {"maxl", required_argument, 0, 'L'},
	    {"minx", required_argument, 0, 'x'},
	    {"maxx", required_argument, 0, 'X'},
	    {"miny", required_argument, 0, 'y'},
	    {"maxy", required_argument, 0, 'Y'},
	    {"nx", required_argument, 0, 'a'},
	    {"ny", required_argument, 0, 'b'},
	    {"global", no_argument, 0, 'g'},
	    {"Version", no_argument, 0, 'V'},
	    {"help", no_argument, 0, 'h'},
	    {0, 0, 0, 0}
	};

	c = getopt_long(argc, argv, "a:b:gl:L:x:X:y:Y:Vh", long_options,
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

	case 'g':
	    coordinates = GLOBAL;
	    break;

	case 'l':
	    lambda_min = atof(optarg);
	    break;

	case 'L':
	    lambda_max = atof(optarg);
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
	    fprintf(stdout, "flux_2D Version: %s  %s\n", RELEASE,
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

    if (wrong_target_type(stdin))
	exit(EXIT_FAILURE);

    h = init_hist(x_max, x_min, x_bins, y_max, y_min, y_bins);

    if (!read_hist
	(stdin, h, lambda_min, lambda_max, &n_inc, &p_inc, &n_missed,
	 &p_missed, coordinates, M, origin))
	output_hist(stdout, h, lambda_min, lambda_max, n_inc, p_inc,
		    n_missed, p_missed, coordinates, M, origin);

    gsl_histogram2d_free(h);
    exit(EXIT_SUCCESS);
}

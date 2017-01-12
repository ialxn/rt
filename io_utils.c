/*	io_utils.c
 *
 * Copyright (C) 2011,2012,2013,2014,2015,2016,2017 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "io_utils.h"
#include "math_utils.h"

static int is_present(const char *section, const config_setting_t * s,
		      const char *name, const int nr)
{

    if (!config_setting_get_member(s, name)) {
	fprintf(stderr,
		"missing keyword '%s' in '%s' section %d\n", name, section,
		nr + 1);
	return ERR;
    } else
	return NO_ERR;
}

int check_array(const char *section, const config_setting_t * s,
		const char *name, const int nr)
{
    int status = NO_ERR;
    config_setting_t *m;

    if ((m = config_setting_get_member(s, name)) == NULL) {
	fprintf(stderr,
		"missing '%s' array in '%s' section %d\n", name, section,
		nr + 1);
	status = ERR;
    } else if (config_setting_is_array(m) == CONFIG_FALSE
	       || config_setting_length(m) != 3) {
	fprintf(stderr,
		"setting '%s' in '%s' section %d is not array with 3 coordinates\n",
		name, section, nr + 1);
	status = ERR;
    }

    return status;
}

int check_string(const char *section, const config_setting_t * s,
		 const char *name, const int nr)
{
    int status = NO_ERR;
    const char *S;

    if ((status = is_present(section, s, name, nr)) != ERR
	&& config_setting_lookup_string(s, name, &S) != CONFIG_TRUE) {
	fprintf(stderr,
		"'%s' keyword in '%s' section %d does not define string\n",
		name, section, nr + 1);
	status = ERR;
    }

    return status;
}

int check_return_string(const char *section, const config_setting_t * s,
			const char *name, const int nr,
			const char **string)
{
    int status = NO_ERR;

    if ((status = is_present(section, s, name, nr)) != ERR
	&& config_setting_lookup_string(s, name, string) != CONFIG_TRUE) {
	fprintf(stderr,
		"'%s' keyword in '%s' section %d does not define string\n",
		name, section, nr + 1);
	status = ERR;
    }

    return status;
}

int check_bool(const char *section, const config_setting_t * s,
	       const char *name, const int nr)
{
    int status = NO_ERR;
    int B;

    if ((status = is_present(section, s, name, nr)) != ERR
	&& config_setting_lookup_bool(s, name, &B) != CONFIG_TRUE) {
	fprintf(stderr,
		"'%s' keyword in '%s' section %d does not define boolean\n",
		name, section, nr + 1);
	status = ERR;
    }

    return status;
}

int check_float(const char *section, const config_setting_t * s,
		const char *name, const int nr)
{
    int status = NO_ERR;
    double F;

    if ((status = is_present(section, s, name, nr)) != ERR
	&& config_setting_lookup_float(s, name, &F) != CONFIG_TRUE) {
	fprintf(stderr,
		"'%s' keyword in '%s' section %d does not define float\n",
		name, section, nr + 1);
	status = ERR;
    }

    return status;
}

int check_int(const char *section, const config_setting_t * s,
	      const char *name, const int nr)
{
    int status = NO_ERR;
    int I;

    if ((status = is_present(section, s, name, nr)) != ERR
	&& config_setting_lookup_int(s, name, &I) != CONFIG_TRUE) {
	fprintf(stderr,
		"'%s' keyword in '%s' section %d does not define int\n",
		name, section, nr + 1);
	status = ERR;
    }

    return status;
}

int check_reflectivity_model(const char *section,
			     const config_setting_t * s, const char *name,
			     const int nr)
{
    int status = ERR;
    const char *S;

    if (config_setting_lookup_string(s, name, &S) == CONFIG_TRUE) {
	if (!strcmp(S, "specular"))
	    status = NO_ERR;
	else if (!strcmp(S, "lambertian"))
	    status = NO_ERR;
	else if (!strcmp(S, "microfacet_gaussian")) {
	    double sigma;

	    if (config_setting_lookup_float
		(s, "microfacet_gaussian_sigma", &sigma) == CONFIG_FALSE) {
		fprintf(stderr,
			"missing model parameter 'microfacet_gaussian_sigma' in reflectivity model '%s' found in '%s' section %d\n",
			S, section, nr + 1);

		status = ERR;
	    } else
		status = NO_ERR;

	} else {
	    fprintf(stderr,
		    "unknow reflectivity model '%s' found in '%s' section %d\n",
		    S, section, nr + 1);
	    status = ERR;

	}
    }

    return status;
}

int check_file(const char *section, const config_setting_t * s,
	       const char *name, const int nr)
{
    int status = NO_ERR;
    const char *S;
    FILE *f;

    if (config_setting_lookup_string(s, name, &S) == CONFIG_FALSE)
	return ERR;

    if ((f = fopen(S, "r")) == NULL) {
	fprintf(stderr,
		"could not open file '%s' defined by keyword '%s' in section '%s' (%d)\n",
		S, name, section, nr + 1);
	status = ERR;
    } else {			/* sufficient data ? */

#define N_DATA_MIN 3

	char line[LINE_LEN];
	int line_number = 1;
	int n_data = 0;

	fgets(line, LINE_LEN, f);
	while (!feof(f) && (n_data < N_DATA_MIN)) {

	    if (strncmp(line, "#", 1) != 0) {	/* not a comment */
		int n_read;
		double tx, ty;

		if ((n_read = sscanf(line, "%lg%lg", &tx, &ty)) > 0) {	/* not a blank line */
		    if (n_read != 2) {	/* insufficient data on line */
			fprintf(stderr,
				"insufficient data (%d items instead of 2) read on line %d of file '%s'\n",
				n_read, line_number, S);
			status = ERR;
		    } else	/* line with correct data */
			n_data++;
		}
	    }
	    fgets(line, LINE_LEN, f);
	    line_number++;
	}

	if (n_data != N_DATA_MIN) {	/* insufficient data  */
	    fprintf(stderr,
		    "insufficient data (%d items instead of %d) in file '%s'\n",
		    n_data, N_DATA_MIN, S);
	    status = ERR;
	}

	fclose(f);
    }
    return status;
}

void read_vector(const config_setting_t * s, const char *name,
		 double *const vec)
{
    int i;
    const config_setting_t *setting = config_setting_get_member(s, name);

    for (i = 0; i < 3; i++)
	vec[i] = config_setting_get_float_elem(setting, i);
}

void read_vector_normalize(const config_setting_t * s, const char *name,
			   double *const vec)
{
    read_vector(s, name, vec);
    normalize(vec);
}

void init_M_from_z(config_setting_t * this, const char *kw, double M[9])
{
    read_vector(this, kw, &M[6]);

    M[0] = 1.0;
    M[1] = 0.0;
    M[2] = 0.0;

    orthonormalize(&M[0], &M[3], &M[6]);
}

double get_A_disk(config_setting_t * this_s, const char *kw,
		  const double radius)
{
    int ans;

    config_setting_lookup_bool(this_s, kw, &ans);
    if (ans)			/* surface emits */
	return M_PI * radius * radius;
    else			/* surface does not emit */
	return 0.0;
}


static int inc_buffers(size_t * n_alloc, double **x, double **y)
{
    double *t;

    (*n_alloc) += BSIZE;

    if ((t = (double *) realloc(*x, (*n_alloc) * sizeof(double))) == NULL) {
	fprintf(stderr, "could not (re)allocate buffer\n");
	return ERR;
    }
    *x = t;

    if ((t = (double *) realloc(*y, (*n_alloc) * sizeof(double))) == NULL) {
	fprintf(stderr, "could not (re)allocate buffer\n");
	return ERR;
    }
    *y = t;

    return NO_ERR;
}


int read_data(FILE * f, double **x, double **y, size_t * n)
{
    size_t n_alloc = 0;
    char line[LINE_LEN];
    int line_number = 1;

    *n = 0;
    *x = NULL;
    *y = NULL;

    fgets(line, LINE_LEN, f);
    while (!feof(f)) {

	if (strncmp(line, "#", 1) != 0) {	/* not a comment */
	    int n_read;
	    double tx, ty;

	    if ((n_read = sscanf(line, "%lg%lg", &tx, &ty)) > 0) {	/* not a blank line */

		if (n_read != 2) {	/* insufficient data */
		    fprintf(stderr,
			    "insufficient data (%d items instead of 2) read on line %d\n",
			    n_read, line_number);
		    free(*x);
		    free(*y);
		    return ERR;
		}

		if (*n >= n_alloc)
		    if (inc_buffers(&n_alloc, x, y)) {
			free(*x);
			free(*y);
			return ERR;
		    }

		(*x)[(*n)] = tx;
		(*y)[(*n)] = ty;

		(*n)++;
	    }

	}

	fgets(line, LINE_LEN, f);
	line_number++;

    }

    return NO_ERR;

}

int skip_N_comments(FILE * f, const int N)
{
    char line[LINE_LEN];
    int n;

    for (n = 0; n < N; n++) {
	fgets(line, LINE_LEN, f);
	if (line[0] != '#') {	/* not a header line */
	    fprintf(stderr,
		    "Wrong file type (unexpected size of header)\n");
	    return (ERR);
	}
    }
    return NO_ERR;
}

double get_P_factor(FILE * f)
{
    char line[LINE_LEN];
    double P_factor;

    /* first line read in 'get_idx()' */
    if (skip_N_comments(f, 1))
	return -1.0;

    fgets(line, LINE_LEN, f);
    sscanf(&line[30], "%lf", &P_factor);

    return P_factor;
}

void read_transformation(FILE * f_in, double M[], double origin[])
{
    char line[LINE_LEN];
    size_t i;

    skip_N_comments(f_in, 7);

    for (i = 0; i < 3; i++) {
	fgets(line, LINE_LEN, f_in);
	sscanf(&line[1], "%lf%lf%lf", &M[3 * i], &M[3 * i + 1],
	       &M[3 * i + 2]);
    }

    fgets(line, LINE_LEN, f_in);	/* skip comment */

    fgets(line, LINE_LEN, f_in);
    sscanf(&line[1], "%lf%lf%lf", &origin[0], &origin[1], &origin[2]);

    skip_N_comments(f_in, 4);
}

double *range(const double min, const double max, const size_t n)
{
    size_t i;
    const double width = (max - min) / (double) n;
    double *r = (double *) malloc((n + 1) * sizeof(double));

    for (i = 1, r[0] = min; i <= n; i++)
	r[i] = r[i - 1] + width;

    return r;
}

int get_idx(FILE * f_in, size_t * idx_lambda)
{
    char line[LINE_LEN];
    int status = NO_ERR;

    fgets(line, LINE_LEN, f_in);

    /* check file type and set indices of columns to be read */
    if (strstr(line, "(annulus)")) {
	*idx_lambda = 3;
    } else if (strstr(line, "(cone)")) {
	*idx_lambda = 4;
    } else if (strstr(line, "(cpc)")) {
	*idx_lambda = 4;
    } else if (strstr(line, "(cylinder)")) {
	*idx_lambda = 4;
    } else if (strstr(line, "(disk)")) {
	*idx_lambda = 3;
    } else if (strstr(line, "(ellipsoid)")) {
	*idx_lambda = 4;
    } else if (strstr(line, "(paraboloid)")) {
	*idx_lambda = 4;
    } else if (strstr(line, "(one-sided plane screen)")) {
	*idx_lambda = 3;
    } else if (strstr(line, "(two-sided plane screen)")) {
	*idx_lambda = 3;
    } else if (strstr(line, "(rectangle)")) {
	*idx_lambda = 3;
    } else if (strstr(line, "(sphere)")) {
	*idx_lambda = 4;
    } else if (strstr(line, "(triangle)")) {
	*idx_lambda = 3;
    } else if (strstr(line, "(window)")) {
	*idx_lambda = 3;
    } else {
	fprintf(stderr, "Unknown target type (%s) found\n", line);
	status = ERR;
    }
    return status;
}

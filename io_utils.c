/*	io_utils.c
 *
 * Copyright (C) 2011,2012,2013,2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_cblas.h>

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
    int status;
    const char *S = ERR;

    if (config_setting_lookup_string(s, name, &S) == CONFIG_TRUE) {
	if (!strcmp(S, "specular"))
	    status = NO_ERR;
	else {
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

    config_setting_lookup_string(s, name, &S);

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

int skip_header(FILE * f)
{
    /* skip header (first line is has already been read */
    char line[LINE_LEN];
    int n;

    for (n = 0; n < T_HEADER_LINES - 1; n++) {
	fgets(line, LINE_LEN, f);
	if (line[0] != '#') {	/* not a header line */
	    fprintf(stderr,
		    "Wrong file type (unexpected size of header)\n");
	    return (ERR);
	}
    }

    return (NO_ERR);
}

double *range(const double min, const double max, const size_t n)
{
    size_t i;
    const double width = (max - min) / n;
    double *r = (double *) malloc((n + 1) * sizeof(double));

    for (i = 1, r[0] = min; i <= n; i++)
	r[i] = r[i - 1] + width;

    return r;
}

int get_idx(FILE * f_in, size_t * idx_lambda, size_t * idx_power)
{
    char line[LINE_LEN];
    int status = NO_ERR;

    fgets(line, LINE_LEN, f_in);

    /* check file type and set indices of columns to be read */
    if (strstr(line, "(one-sided plane screen)")) {
	*idx_lambda = 3;
	*idx_power = 2;
    } else if (strstr(line, "(two-sided plane screen)")) {
	*idx_lambda = 3;
	*idx_power = 2;
    } else if (strstr(line, "(rectangle)")) {
	*idx_lambda = 3;
	*idx_power = 2;
    } else if (strstr(line, "(triangle)")) {
	*idx_lambda = 3;
	*idx_power = 2;
    } else if (strstr(line, "(ellipsoid)")) {
	*idx_lambda = 4;
	*idx_power = 3;
    } else if (strstr(line, "(annulus)")) {
	*idx_lambda = 3;
	*idx_power = 2;
    } else if (strstr(line, "(disk)")) {
	*idx_lambda = 3;
	*idx_power = 2;
    } else {
	fprintf(stderr, "Unknown target type (%s) found\n", line);
	status = ERR;
    }
    return status;
}

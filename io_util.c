/*	io_util.c
 *
 * Copyright (C) 2011 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <gsl/gsl_cblas.h>

#include "io_util.h"

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
    config_setting_t *m;

    if ((m = config_setting_get_member(s, name)) == NULL) {
	fprintf(stderr,
		"missing keyword '%s' in '%s' section %d\n", name, section,
		nr + 1);
	status = ERR;
    } else if (config_setting_lookup_string(s, name, &S) != CONFIG_TRUE) {
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
    config_setting_t *m;

    if ((m = config_setting_get_member(s, name)) == NULL) {
	fprintf(stderr,
		"missing keyword '%s' in '%s' section %d\n", name, section,
		nr + 1);
	status = ERR;
    } else if (config_setting_lookup_string(s, name, string) !=
	       CONFIG_TRUE) {
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
    config_setting_t *m;

    if ((m = config_setting_get_member(s, name)) == NULL) {
	fprintf(stderr,
		"missing keyword '%s' in '%s' section %d\n", name, section,
		nr + 1);
	status = ERR;
    } else if (config_setting_lookup_float(s, name, &F) != CONFIG_TRUE) {
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
    config_setting_t *m;

    if ((m = config_setting_get_member(s, name)) == NULL) {
	fprintf(stderr,
		"missing keyword '%s' in '%s' section %d\n", name, section,
		nr + 1);
	status = ERR;
    } else if (config_setting_lookup_int(s, name, &I) != CONFIG_TRUE) {
	fprintf(stderr,
		"'%s' keyword in '%s' section %d does not define int\n",
		name, section, nr + 1);
	status = ERR;
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
    double norm;

    read_vector(s, name, vec);
    norm = cblas_dnrm2(3, vec, 1);
    cblas_dscal(3, 1.0 / norm, vec, 1);
}

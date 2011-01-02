/*	io_util.c
 *
 * Copyright (C) 2011 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

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

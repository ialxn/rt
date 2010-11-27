/*	sources.c
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <string.h>

#include "sources.h"

int check_sources(config_t * cfg)
{
    int status = NO_ERR;

    const config_setting_t *s = config_lookup(cfg, "sources");

    if (s != NULL) {
	const int count = config_setting_length(s);
	int i;

	if (count == 0) {
	    fprintf(stderr, "empty 'sources' section\n");
	    status = ERR;
	}

	for (i = 0; i < count; ++i) {
	    config_setting_t *this_s =
		config_setting_get_elem(s, (unsigned int) i);

	    const char *S, *type;
	    double F;
	    int I;

	    /*
	     * keywords common to all sources
	     * 'name':  identifier / string
	     * 'type':  type of source / string
	     *          - "point_source": uniform point source
	     * 'power': power [W] of source / double
	     * 'n_rays': number of rays used for this source / int
	     */
	    if (config_setting_lookup_string(this_s, "name", &S) !=
		CONFIG_TRUE) {
		fprintf(stderr,
			"missing 'name' keyword in 'sources' section %u\n",
			i + 1);
		status = ERR;
	    }
	    if (config_setting_lookup_string(this_s, "type", &type) !=
		CONFIG_TRUE) {
		fprintf(stderr,
			"missing 'type' keyword in 'sources' section %u\n",
			i + 1);
		status = ERR;
	    }
	    if (config_setting_lookup_float(this_s, "power", &F) !=
		CONFIG_TRUE) {
		fprintf(stderr,
			"missing 'power' keyword in 'sources' section %u\n",
			i + 1);
		status = ERR;
	    }
	    if (config_setting_lookup_int(this_s, "n_rays", &I) !=
		CONFIG_TRUE) {
		fprintf(stderr,
			"missing 'n_rays' keyword in 'sources' section %u\n",
			i + 1);
		status = ERR;
	    }

	    /*
	     * check source specific settings
	     */

	    /*
	     * uniform point source:
	     *  - group 'origin'
	     *          'x', 'y', 'z': coordinates / double
	     */
	    if (strstr(type, "point_source") == type) {
		config_setting_t *origin;

		if ((origin =
		     config_setting_get_member(this_s,
					       "origin")) == NULL) {
		    fprintf(stderr,
			    "missing 'origin' group in 'sources' section %u\n",
			    i + 1);
		    status = ERR;
		} else {

		    if (config_setting_lookup_float(origin, "x", &F)
			!= CONFIG_TRUE) {
			fprintf(stderr,
				"missing keyword x in group 'origin' in 'sources' section %u\n",
				i + 1);
			status = ERR;
		    }

		    if (config_setting_lookup_float(origin, "y", &F)
			!= CONFIG_TRUE) {
			fprintf(stderr,
				"missing keyword y in group 'origin' in 'sources' section %u\n",
				i + 1);
			status = ERR;
		    }

		    if (config_setting_lookup_float(origin, "z", &F)
			!= CONFIG_TRUE) {
			fprintf(stderr,
				"missing keyword y in group 'origin' in 'sources' section %u\n",
				i + 1);
			status = ERR;
		    }

		}

	    }
	    /* end 'uniform point source' */
	}
    }

    else {
	fprintf(stderr, "missing 'sources' keyword\n");
	status = ERR;
    }

    return status;
}

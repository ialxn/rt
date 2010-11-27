/*	targets.c
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <string.h>

#include "targets.h"

int check_targets(config_t * cfg)
{
    int status = NO_ERR;

    const config_setting_t *t = config_lookup(cfg, "targets");

    if (t != NULL) {
	const int count = config_setting_length(t);
	int i;

	if (count == 0) {
	    fprintf(stderr, "empty 'targets' section\n");
	    status = ERR;
	}

	for (i = 0; i < count; ++i) {
	    config_setting_t *this_t =
		config_setting_get_elem(t, (unsigned int) i);

	    const char *S, *type;
	    double F;

	    /*
	     * keywords common to all targets
	     * 'name':  identifier / string
	     * 'type':  type of source / string
	     *          - "plane_screen": non-absorbing counter plane
	     */
	    if (config_setting_lookup_string(this_t, "name", &S) !=
		CONFIG_TRUE) {
		fprintf(stderr,
			"missing 'name' keyword in 'targets' section %u\n",
			i + 1);
		status = ERR;
	    }
	    if (config_setting_lookup_string(this_t, "type", &type) !=
		CONFIG_TRUE) {
		fprintf(stderr,
			"missing 'type' keyword in 'targets' section %u\n",
			i + 1);
		status = ERR;
	    }

	    /*
	     * check target specific settings
	     */

	    /*
	     * plane_screen:
	     *  - group 'point' (point on plane)
	     *          'x', 'y', 'z': coordinates / double
	     *  - group 'normal' (normal vector of plane)
	     *          'x', 'y', 'z': coordinates / double
	     */
	    if (strstr(type, "plane_screen") == type) {
		config_setting_t *point, *normal;

		if ((point =
		     config_setting_get_member(this_t, "point")) == NULL) {
		    fprintf(stderr,
			    "missing 'point' group in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else {

		    if (config_setting_lookup_float(point, "x", &F)
			!= CONFIG_TRUE) {
			fprintf(stderr,
				"missing keyword x in group 'point' in 'targets' section %u\n",
				i + 1);
			status = ERR;
		    }

		    if (config_setting_lookup_float(point, "y", &F)
			!= CONFIG_TRUE) {
			fprintf(stderr,
				"missing keyword y in group 'point' in 'targets' section %u\n",
				i + 1);
			status = ERR;
		    }

		    if (config_setting_lookup_float(point, "z", &F)
			!= CONFIG_TRUE) {
			fprintf(stderr,
				"missing keyword y in group 'point' in 'targets' section %u\n",
				i + 1);
			status = ERR;
		    }

		}		/* end 'point' */

		if ((normal =
		     config_setting_get_member(this_t,
					       "normal")) == NULL) {
		    fprintf(stderr,
			    "missing 'normal' group in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else {

		    if (config_setting_lookup_float(normal, "x", &F)
			!= CONFIG_TRUE) {
			fprintf(stderr,
				"missing keyword x in group 'normal' in 'targets' section %u\n",
				i + 1);
			status = ERR;
		    }

		    if (config_setting_lookup_float(normal, "y", &F)
			!= CONFIG_TRUE) {
			fprintf(stderr,
				"missing keyword y in group 'normal' in 'targets' section %u\n",
				i + 1);
			status = ERR;
		    }

		    if (config_setting_lookup_float(normal, "z", &F)
			!= CONFIG_TRUE) {
			fprintf(stderr,
				"missing keyword y in group 'normal' in 'targets' section %u\n",
				i + 1);
			status = ERR;
		    }

		}		/* end 'normal' */

	    }			/* end 'plane_screen' */
	}
    }

    else {
	fprintf(stderr, "missing 'targets' keyword\n");
	status = ERR;
    }

    return status;
}

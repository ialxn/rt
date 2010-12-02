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
#include <stdlib.h>

#include "targets.h"

target_t *target_alloc(const target_type_t * type, config_t * cfg,
		       const char *name)
{
/*
 * generic part
 */

    target_t *T;

    if ((T = (target_t *) malloc(sizeof(target_t))) == NULL) {
	fprintf(stderr, "failed to allocate space for target\n");
	return NULL;
    }

    if ((T->state = malloc(type->size)) == NULL) {
	free(T);
	fprintf(stderr,
		"failed to allocate space for target internal state\n");
	return NULL;
    };

    T->type = type;

    /*
     * specific part
     */
    if ((T->type->alloc_state) (T->state) == ERR) {
	fprintf(stderr,
		"failed to allocate space for target state content\n");
	free(T->state);
	free(T);
	return NULL;
    }

    (T->type->init_state) (T->state, cfg, name);

    return T;

}

void target_free(target_t * T)
{
/*
 * free specific internal data (prevent leaks).
 */

    (T->type->free_state) (T->state);
    free(T->state);
    free(T);

}

double *interception(const target_t * T, ray_t * in_ray, int *dump_flag)
{
    return (T->type->get_intercept) (T->state, in_ray, dump_flag);
}

ray_t *out_ray(const target_t * T, ray_t * in_ray, double *hit,
	       int *dump_flag)
{
    return (T->type->get_out_ray) (T->state, in_ray, hit, dump_flag);
}

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
	     * plane screen:
	     *  - group 'point' (point on plane)
	     *          'x', 'y', 'z': coordinates / double
	     *  - group 'normal' (normal vector of plane)
	     *          'x', 'y', 'z': coordinates / double
	     */
	    if (strstr(type, "plane screen") == type) {
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

void dump_data(FILE * f, double *data, const size_t n_data,
	       const size_t n_items)
{
    size_t i, j;

    for (i = 0; i < n_data; i++) {
	const size_t N = i * n_items;

	for (j = 0; j < n_items; j++) {

	    fprintf(f, "%g", data[N + j]);

	    if (j == n_items - 1)
		fprintf(f, "\n");
	    else
		fprintf(f, "\t");

	}
    }
}

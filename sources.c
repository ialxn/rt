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

source_t *source_alloc(const source_type_t * T, config_t * cfg,
		       const char *name)
{
/*
 * generic part
 */

    source_t *S;

    if ((S = (source_t *) malloc(sizeof(source_t))) == NULL) {
	fprintf(stderr, "failed to allocate space for property struct\n");
	return NULL;
    }

    if ((S->state = malloc(T->size)) == NULL) {
	free(S);
	fprintf(stderr, "failed to allocate space for source state\n");
	return NULL;
    };

    S->type = T;

    /*
     * specific part
     */
    if ((S->type->alloc_state) (S->state) == ERR) {
	fprintf(stderr,
		"failed to allocate space for source state content\n");
	free(S->state);
	free(S);
	return NULL;
    }

    (S->type->init_state) (S->state, cfg, name);

    return S;

}

void source_free(source_t * S)
{
/*
 * free specific internal data (prevent leaks).
 */

    (S->type->free_state) (S->state);
    free(S->state);
    free(S);

}

ray_t *new_ray(const source_t * S, const gsl_rng * r)
{
    ray_t *ray;

    ray = (S->type->get_new_ray) (S->state, r);
    return ray;
}



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

/*	sources.c
 *
 * Copyright (C) 2010, 2011 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <string.h>

#include "io_util.h"
#include "sources.h"

source_t *source_alloc(const source_type_t * T, config_t * cfg,
		       const char *name)
{
    source_t *S;

    /* generic part */
    if ((S = (source_t *) malloc(sizeof(source_t))) == NULL) {
	fprintf(stderr, "failed to allocate space for 'source_t'\n");
	return NULL;
    }

    if ((S->state = malloc(T->size)) == NULL) {
	free(S);
	fprintf(stderr,
		"failed to allocate space for generic state of source\n");
	return NULL;
    };

    S->type = T;

    (S->type->init_state) (S->state, cfg, name);	/* initialize data structures */

    return S;
}

void source_free(source_t * S)
{
    (S->type->free_state) (S->state);
    free(S->state);
    free(S);
}

ray_t *new_ray(const source_t * S, const gsl_rng * r)
{
    return (S->type->get_new_ray) (S->state, r);
}

double get_ppr(const source_t * S)
{
    return (S->type->get_ppr) (S->state);
}

const char *get_source_type(const source_t * S)
{
    return S->type->type;
}

const char *get_source_name(const source_t * S)
{
    return (S->type->get_source_name) (S->state);
}



int check_sources(config_t * cfg)
{
    int status = NO_ERR;
    const config_setting_t *s = config_lookup(cfg, "sources");

    if (s == NULL) {
	fprintf(stderr, "missing 'sources' keyword\n");
	status += ERR;
    } else {			/* 'sources' section present */
	int i;
	const int n_sources = config_setting_length(s);

	if (n_sources == 0) {
	    fprintf(stderr, "empty 'sources' section\n");
	    status += ERR;
	}

	for (i = 0; i < n_sources; ++i) {
	    const char *type = NULL;
	    config_setting_t *this_s =
		config_setting_get_elem(s, (unsigned int) i);

	    /*
	     * keywords common to all sources
	     * 'name':  identifier / string
	     * 'type':  type of source / string
	     *          - "uniform point_source": uniform point source
	     * 'power': power [W] of source / double
	     * 'n_rays': number of rays used for this source / int
	     */
	    status += check_string("sources", this_s, "name", i);
	    status += check_float("sources", this_s, "power", i);
	    status += check_int("sources", this_s, "n_rays", i);

	    status +=
		check_return_string("sources", this_s, "type", i, &type);
	    if (!type)
		continue;

	    /* check source specific settings */
	    if (!strcmp(type, "uniform point source")) {
		/*
		 * uniform point source:
		 *  - array 'origin' [x,y,z] / double            */
		status += check_array("sources", this_s, "origin", i);

	    } /* end 'uniform point source' */
	    else if (!strcmp(type, "spot source")) {
		/*
		 * spot source:
		 *  - array 'origin' [x,y,z] / double
		 *  - array 'direction' [x,y,z] / double
		 *  - 'theta' / double
		 */
		status += check_array("sources", this_s, "origin", i);
		status += check_array("sources", this_s, "direction", i);
		status += check_float("sources", this_s, "theta", i);

	    }			/* end 'uniform point source' */
	}			/* end 'this_s', check next source */
    }

    return status;
}

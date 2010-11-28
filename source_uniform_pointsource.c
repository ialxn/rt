/*	source_uniform_pointsource.c
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

typedef struct ups_state_t {
    const char *name;
    const char *type;
    int n_rays;
    double power;
    double ppr;
} ups_state_t;


static int ups_alloc_state(void *vstate)
{
    return NO_ERR;
}


static int ups_init_state(void *vstate, config_t * cfg, const char *name)
{
    ups_state_t *state = (ups_state_t *) vstate;
    state->name = strdup(name);

    return NO_ERR;
}

static void ups_free_state(void *vstate)
{
    ups_state_t *state = (ups_state_t *) vstate;

    free(state->name);
    free(state->type);
}


static int ups_get_new_ray(void *vstate, const gsl_rng * r)
{
    ups_state_t *state = (ups_state_t *) vstate;

    if (state->n_rays) {
	state->n_rays--;
	return 1;
    }
    return 0;
}



static const source_type_t ups_type_t = {
    sizeof(struct ups_state_t),
    &ups_alloc_state,
    &ups_init_state,
    &ups_free_state, &ups_get_new_ray
};

const source_type_t *source_uniform_pointsource = &ups_type_t;

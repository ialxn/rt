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

static void ups_init_state(void *vstate, config_t * cfg, const char *name)
{
    ups_state_t *state = (ups_state_t *) vstate;

    unsigned int i = 0;
    const char *S;
    int I;
    double F;

    config_setting_t *this_s;
    const config_setting_t *s = config_lookup(cfg, "sources");

    state->name = strdup(name);

    while (1) {			/* find setting for source 'name' */
	this_s = config_setting_get_elem(s, i);

	config_setting_lookup_string(this_s, "name", &S);
	if (strstr(S, name))
	    break;

	i++;
    }

    config_setting_lookup_string(this_s, "type", &S);
    state->type = strdup(S);

    config_setting_lookup_int(this_s, "n_rays", &I);
    state->n_rays = I;

    config_setting_lookup_float(this_s, "power", &F);
    state->power = F;

    state->ppr = F / I;

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



static const source_type_t ups_t = {
    sizeof(struct ups_state_t),
    &ups_alloc_state,
    &ups_init_state,
    &ups_free_state,
    &ups_get_new_ray
};

const source_type_t *source_uniform_point_source = &ups_t;

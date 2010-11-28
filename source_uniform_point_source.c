/*	source_uniform_pointsource.c
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

#include "ray.h"
#include "sources.h"

typedef struct ups_state_t {
    const char *name;
    const char *type;
    double origin[3];
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

static ray_t *ups_get_new_ray(void *vstate, const gsl_rng * r)
{
    ups_state_t *state = (ups_state_t *) vstate;
    ray_t *ray = NULL;

    if (state->n_rays) {
	double t;
	double theta, sin_theta, phi;

	ray = (ray_t *) malloc(sizeof(ray_t));

	t = gsl_rng_uniform(r);
	theta = acos(1.0 - 2.0 * t);
	sin_theta = sin(theta);

	t = gsl_rng_uniform(r);
	phi = 2.0 * M_PI * t;

	ray->direction.x = sin_theta * cos(phi);
	ray->direction.y = sin_theta * sin(phi);
	ray->direction.z = cos(theta);

	ray->origin.x = state->origin[0];
	ray->origin.y = state->origin[1];
	ray->origin.z = state->origin[2];

	ray->power = state->ppr;

	state->n_rays--;
    }

    return ray;
}



static const source_type_t ups_t = {
    sizeof(struct ups_state_t),
    &ups_alloc_state,
    &ups_init_state,
    &ups_free_state,
    &ups_get_new_ray
};

const source_type_t *source_uniform_point_source = &ups_t;

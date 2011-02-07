/*	source_uniform_point_source.c
 *
 * Copyright (C) 2010, 2011 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <math.h>
#include <string.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>

#include "io_util.h"
#include "ray.h"
#include "sources.h"

typedef struct ups_state_t {
    char *name;			/* name (identifier) of uniform point source */
    double origin[3];
    int n_rays;			/* number of rays remaining until source is exhausted */
    double power;		/* power of source */
    double ppr;			/* power allotted to one ray */
    gsl_spline *spline;		/* spline holding cdf of source spectrum */
    gsl_interp_accel *acc;	/* cached data for spline interpolation */
} ups_state_t;


static void ups_init_state(void *vstate, config_t * cfg, const char *name)
{
    ups_state_t *state = (ups_state_t *) vstate;

    unsigned int i = 0;
    const char *S;
    config_setting_t *this_s;
    const config_setting_t *s = config_lookup(cfg, "sources");

    state->name = strdup(name);

    while (1) {			/* find setting for source 'name' */
	this_s = config_setting_get_elem(s, i);

	config_setting_lookup_string(this_s, "name", &S);
	if (!strcmp(S, name))
	    break;

	i++;
    }

    config_setting_lookup_int(this_s, "n_rays", &state->n_rays);
    config_setting_lookup_float(this_s, "power", &state->power);
    state->ppr = state->power / state->n_rays;

    read_vector(this_s, "origin", state->origin);

    /* initialize source spectrum */
    config_setting_lookup_string(this_s, "spectrum", &S);
    init_spectrum(S, &state->spline, &state->acc);
}

static void ups_free_state(void *vstate)
{
    ups_state_t *state = (ups_state_t *) vstate;

    free(state->name);
    gsl_spline_free(state->spline);
    gsl_interp_accel_free(state->acc);
}

static ray_t *ups_get_new_ray(void *vstate, const gsl_rng * r)
{
    ups_state_t *state = (ups_state_t *) vstate;
    ray_t *ray = NULL;

    if (state->n_rays) {
	double t;
	double cos_theta, sin_theta, phi;

	ray = (ray_t *) malloc(sizeof(ray_t));

	t = gsl_rng_uniform(r);
	cos_theta = 1.0 - 2.0 * t;
	sin_theta = sin(acos(cos_theta));

	t = gsl_rng_uniform(r);
	phi = 2.0 * M_PI * t;

	ray->direction[0] = sin_theta * cos(phi);
	ray->direction[1] = sin_theta * sin(phi);
	ray->direction[2] = cos_theta;

	memcpy(ray->origin, state->origin, 3 * sizeof(double));

	ray->power = state->ppr;
	t = gsl_rng_uniform(r);
	ray->lambda = gsl_spline_eval(state->spline, t, state->acc);


	state->n_rays--;
    }

    return ray;
}

static const char *ups_get_source_name(void *vstate)
{
    ups_state_t *state = (ups_state_t *) vstate;

    return state->name;
}


static const source_type_t ups_t = {
    "uniform point source",
    sizeof(struct ups_state_t),
    &ups_init_state,
    &ups_free_state,
    &ups_get_new_ray,
    &ups_get_source_name
};

const source_type_t *source_uniform_point_source = &ups_t;

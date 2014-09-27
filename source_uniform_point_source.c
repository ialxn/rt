/*	source_uniform_point_source.c
 *
 * Copyright (C) 2010,2011,2012,2013,2014 Ivo Alxneit
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

#include "io_utils.h"
#include "math_utils.h"
#include "ray.h"
#include "sources.h"



typedef struct ups_state_t {
    char *name;			/* name (identifier) of uniform point source */
    double orig[3];
    int64_t n_rays;		/* number of rays remaining until source is exhausted */
    pthread_mutex_t mutex_n_rays;	/* protect n_rays */
    pthread_key_t rays_remain_key;	/* no of ray remain in group (PTD) */
    double power;		/* power of source */
    double ppr;			/* power allotted to one ray */
    gsl_spline *spectrum;	/* spline holding cdf of source spectrum */
} ups_state_t;


static void ups_init_state(void *vstate, config_setting_t * this_s)
{
    ups_state_t *state = (ups_state_t *) vstate;

    const char *S;

    config_setting_lookup_string(this_s, "name", &S);
    state->name = strdup(S);
    config_setting_lookup_int64(this_s, "n_rays", &state->n_rays);
    pthread_mutex_init(&state->mutex_n_rays, NULL);
    config_setting_lookup_float(this_s, "power", &state->power);
    state->ppr = state->power / (double) state->n_rays;

    read_vector(this_s, "origin", state->orig);

    /* initialize source spectrum */
    config_setting_lookup_string(this_s, "spectrum", &S);
    init_spectrum(S, &state->spectrum);

    pthread_key_create(&state->rays_remain_key, free);
}

static void ups_free_state(void *vstate)
{
    ups_state_t *state = (ups_state_t *) vstate;

    free(state->name);
    gsl_spline_free(state->spectrum);
}

static ray_t *ups_emit_ray(void *vstate, const gsl_rng * r)
{
    ups_state_t *state = (ups_state_t *) vstate;
    ray_t *ray = NULL;
    int64_t *rays_remain = pthread_getspecific(state->rays_remain_key);

    if (!*rays_remain)
	*rays_remain =
	    per_thread_get_new_raygroup(&state->mutex_n_rays,
					&state->n_rays);

    if (*rays_remain > 0) {	/* rays still available in group */

	(*rays_remain)--;
	pthread_setspecific(state->rays_remain_key, rays_remain);

	ray = (ray_t *) malloc(sizeof(ray_t));

	get_uniform_random_vector(ray->dir, 1.0, r);
	memcpy(ray->orig, state->orig, 3 * sizeof(double));
	ray->power = state->ppr;
	ray->lambda =
	    gsl_spline_eval(state->spectrum, gsl_rng_uniform(r), NULL);
    }

    return ray;
}

static const char *ups_get_source_name(void *vstate)
{
    return ((ups_state_t *) vstate)->name;
}

static int64_t ups_get_source_n_rays(void *vstate)
{
    ups_state_t *state = (ups_state_t *) vstate;

    return per_thread_get_source_n_rays(&state->mutex_n_rays,
					&state->n_rays);
}

static double ups_get_source_power(void *vstate)
{
    return ((ups_state_t *) vstate)->power;
}

static void ups_init_rays_remain(void *vstate)
{
    per_thread_init_rays_remain(((ups_state_t *) vstate)->rays_remain_key);
}


static const source_type_t ups_t = {
    "uniform point source",
    sizeof(struct ups_state_t),
    &ups_init_state,
    &ups_free_state,
    &ups_emit_ray,
    &ups_get_source_name,
    &ups_get_source_n_rays,
    &ups_get_source_power,
    &ups_init_rays_remain
};

const source_type_t *source_uniform_point_source = &ups_t;

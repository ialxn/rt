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
    pthread_mutex_t mutex_n_rays;	/* protect n_rays */
    pthread_key_t rays_remain_key;	/* no of ray remain in group (PTD) */
    double power;		/* power of source */
    double ppr;			/* power allotted to one ray */
    gsl_spline *spline;		/* spline holding cdf of source spectrum */
    double lambda_min;		/* minimum wavelength in spectrum */
} ups_state_t;


static void ups_init_state(void *vstate, config_setting_t * this_s,
			   config_t * cfg)
{
    ups_state_t *state = (ups_state_t *) vstate;

    const char *S;

    (void) cfg;			/* avoid warning: unused parameter 'cfg' */

    config_setting_lookup_string(this_s, "name", &S);
    state->name = strdup(S);
    config_setting_lookup_int(this_s, "n_rays", &state->n_rays);
    pthread_mutex_init(&state->mutex_n_rays, NULL);
    config_setting_lookup_float(this_s, "power", &state->power);
    state->ppr = state->power / state->n_rays;

    read_vector(this_s, "origin", state->origin);

    /* initialize source spectrum */
    config_setting_lookup_string(this_s, "spectrum", &S);
    init_spectrum(S, &state->spline, &state->lambda_min);

    pthread_key_create(&state->rays_remain_key, free);
}

static void ups_free_state(void *vstate)
{
    ups_state_t *state = (ups_state_t *) vstate;

    free(state->name);
    gsl_spline_free(state->spline);
}

static ray_t *ups_get_new_ray(void *vstate, const gsl_rng * r)
{
    ups_state_t *state = (ups_state_t *) vstate;
    ray_t *ray = NULL;
    int *rays_remain = pthread_getspecific(state->rays_remain_key);

    if (!*rays_remain) {
	/*
	 * group of rays has been consumed. check if source is not
	 * yet exhausted
	 */
	int work_needed;

	pthread_mutex_lock(&state->mutex_n_rays);
	work_needed = state->n_rays;

	if (work_needed >= RAYS_PER_GROUP) {	/* get new group */
	    state->n_rays -= RAYS_PER_GROUP;
	    *rays_remain = RAYS_PER_GROUP;
	} else {		/* make source empty */
	    state->n_rays = 0;
	    *rays_remain = work_needed;
	    /*
	     * if source was already exhausted, work_needed is zero
	     * and no ray will be emitted
	     */
	}
	pthread_mutex_unlock(&state->mutex_n_rays);
    }

    if (*rays_remain > 0) {	/* rays still available in group */
	double t;
	double cos_theta, sin_theta, phi;

	(*rays_remain)--;
	pthread_setspecific(state->rays_remain_key, rays_remain);

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

	/* choose random wavelength */
	ray->lambda =
	    state->lambda_min + gsl_spline_eval(state->spline,
						gsl_rng_uniform(r), NULL);
    }

    return ray;
}

static const char *ups_get_source_name(void *vstate)
{
    ups_state_t *state = (ups_state_t *) vstate;

    return state->name;
}

static double ups_get_source_ppr(void *vstate)
{
    ups_state_t *state = (ups_state_t *) vstate;

    return state->ppr;
}

static void ups_init_rays_remain(void *vstate)
{
    ups_state_t *state = (ups_state_t *) vstate;

    int *rays_remain = (int *) malloc(sizeof(int));

    *rays_remain = 0;
    pthread_setspecific(state->rays_remain_key, rays_remain);
}


static const source_type_t ups_t = {
    "uniform point source",
    sizeof(struct ups_state_t),
    &ups_init_state,
    &ups_free_state,
    &ups_get_new_ray,
    &ups_get_source_name,
    &ups_get_source_ppr,
    &ups_init_rays_remain
};

const source_type_t *source_uniform_point_source = &ups_t;

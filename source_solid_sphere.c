/*	source_sphere.c
 *
 * Copyright (C) 2011,2012,2013,2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#include <string.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include "io_utils.h"
#include "math_utils.h"
#include "ray.h"
#include "sources.h"

typedef struct ssp_state_t {
    char *name;			/* name (identifier) of uniform point source */
    double orig[3];
    double radius;
    int n_rays;			/* number of rays remaining until source is exhausted */
    pthread_mutex_t mutex_n_rays;	/* protect n_rays */
    pthread_key_t rays_remain_key;	/* no of ray remain in group (PTD) */
    double power;		/* power of source */
    double ppr;			/* power allotted to one ray */
    gsl_spline *spectrum;	/* spline holding cdf of source spectrum */
} ssp_state_t;


static void ssp_init_state(void *vstate, config_setting_t * this_s)
{
    ssp_state_t *state = (ssp_state_t *) vstate;

    const char *S;

    config_setting_lookup_string(this_s, "name", &S);
    state->name = strdup(S);
    config_setting_lookup_int(this_s, "n_rays", &state->n_rays);
    pthread_mutex_init(&state->mutex_n_rays, NULL);
    config_setting_lookup_float(this_s, "power", &state->power);
    state->ppr = state->power / state->n_rays;

    config_setting_lookup_float(this_s, "radius", &state->radius);
    read_vector(this_s, "origin", state->orig);

    /* initialize source spectrum */
    config_setting_lookup_string(this_s, "spectrum", &S);
    init_spectrum(S, &state->spectrum);

    pthread_key_create(&state->rays_remain_key, free);
}

static void ssp_free_state(void *vstate)
{
    ssp_state_t *state = (ssp_state_t *) vstate;

    free(state->name);
    gsl_spline_free(state->spectrum);
}

static ray_t *ssp_emit_ray(void *vstate, const gsl_rng * r)
{
    ssp_state_t *state = (ssp_state_t *) vstate;
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
	double point[3];

	(*rays_remain)--;
	pthread_setspecific(state->rays_remain_key, rays_remain);

	ray = (ray_t *) malloc(sizeof(ray_t));

	/*
	 * select random point on sphere with radius 'state->radius'
	 * and initialize 'ray->origin' with it.
	 */
	get_uniform_random_vector(point, state->radius, r);

	ray->orig[0] = state->orig[0] + point[0];
	ray->orig[1] = state->orig[1] + point[1];
	ray->orig[2] = state->orig[2] + point[2];

	/*
	 * choose random direction
	 * make shure ray does not enter sphere i.e.
	 * 'point' dot 'ray->dir' is positive. 'point', the
	 * vector from the origin of the sphere to the point
	 * on the sphere's surface is parallel to the normal
	 * vector (but not normalized, which is no problem here).
	 */
	do {
	    get_uniform_random_vector(ray->dir, 1.0, r);
	} while (cblas_ddot(3, point, 1, ray->dir, 1) < 0.0);

	ray->power = state->ppr;

	/* choose random wavelength */
	ray->lambda =
	    gsl_spline_eval(state->spectrum, gsl_rng_uniform(r), NULL);
    }
    return ray;
}

static const char *ssp_get_source_name(void *vstate)
{
    ssp_state_t *state = (ssp_state_t *) vstate;

    return state->name;
}

static void ssp_init_rays_remain(void *vstate)
{
    ssp_state_t *state = (ssp_state_t *) vstate;

    int *rays_remain = (int *) malloc(sizeof(int));

    *rays_remain = 0;
    pthread_setspecific(state->rays_remain_key, rays_remain);
}


static const source_type_t ssp_t = {
    "solid sphere",
    sizeof(struct ssp_state_t),
    &ssp_init_state,
    &ssp_free_state,
    &ssp_emit_ray,
    &ssp_get_source_name,
    &ssp_init_rays_remain
};

const source_type_t *source_solid_sphere = &ssp_t;

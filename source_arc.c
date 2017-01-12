/*	source_arc.c
 *
 * Copyright (C) 2015,2016,2017 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#define _ISOC99_SOURCE		/* because of llrint() */
#define _GNU_SOURCE		/* for sincos() */

#include <math.h>
#include <string.h>

#include <gsl/gsl_machine.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include "io_utils.h"
#include "math_utils.h"
#include "off.h"
#include "ray.h"
#include "sources.h"

typedef struct sarc_state_t {
    char *name;			/* name (identifier) of uniform point source */
    double orig[3];		/* origin (center of base face) */
    double dir[3];		/* direction of cylinder axis */
    double radius;		/* radius of cylinder */
    double length;		/* length of cylinder */
    double alpha, beta;		/* convert local to global coordinates */
    int64_t n_rays;		/* number of rays remaining until source is exhausted */
    double power;		/* power of source */
    gsl_spline *spectrum;	/* spline holding cdf of source spectrum */
    pthread_mutex_t mutex_n_rays;	/* protect n_rays */
    pthread_key_t rays_remain_key;	/* no of ray remain in group (PTD) */
} sarc_state_t;


static void sarc_init_state(void *vstate, config_setting_t * this_s,
			    const double P_factor)
{
    sarc_state_t *state = (sarc_state_t *) vstate;
    const char *S;
    double t[3];

    config_setting_lookup_string(this_s, "name", &S);
    state->name = strdup(S);

    config_setting_lookup_float(this_s, "power", &state->power);
    state->n_rays = llrint(state->power / P_factor);
    pthread_mutex_init(&state->mutex_n_rays, NULL);

    read_vector(this_s, "origin", state->orig);
    read_vector(this_s, "direction", state->dir);
    config_setting_lookup_float(this_s, "radius", &state->radius);
    config_setting_lookup_float(this_s, "length", &state->length);

    /* determine alpha, beta. discard t */
    g2l_off_rot(state->dir, t, &state->alpha, &state->beta);

    /* initialize source spectrum */
    init_source_spectrum(this_s, "spectrum", &state->spectrum);

    pthread_key_create(&state->rays_remain_key, free);
}

static void sarc_free_state(void *vstate)
{
    sarc_state_t *state = (sarc_state_t *) vstate;

    free(state->name);
    gsl_spline_free(state->spectrum);
}

static ray_t *sarc_emit_ray(void *vstate, const gsl_rng * r)
{
    sarc_state_t *state = (sarc_state_t *) vstate;
    ray_t *ray = NULL;
    int64_t *rays_remain = pthread_getspecific(state->rays_remain_key);

    if (!*rays_remain)
	*rays_remain =
	    per_thread_get_new_raygroup(&state->mutex_n_rays,
					&state->n_rays);

    if (*rays_remain > 0) {	/* rays still available in group */
	double point[3];
	double t;
	double sin_phi, cos_phi;

	(*rays_remain)--;
	pthread_setspecific(state->rays_remain_key, rays_remain);

	ray = (ray_t *) malloc(sizeof(ray_t));

	/*
	 * select random point in cylinder.
	 */
	t = 2.0 * M_PI * gsl_rng_uniform(r);
	sincos(t, &sin_phi, &cos_phi);

	t = state->radius * sqrt(gsl_rng_uniform(r));

	point[0] = t * cos_phi;
	point[1] = t * sin_phi;
	point[2] = state->length * gsl_rng_uniform(r);

	l2g_off(state->orig, point, ray->orig, state->alpha, state->beta);

	/*
	 * choose random direction
	 */
	get_uniform_random_vector(point, 1.0, r);
	l2g_off_rot(point, ray->dir, state->alpha, state->beta);

	ray->lambda =
	    gsl_spline_eval(state->spectrum, gsl_rng_uniform(r), NULL);
	ray->n_refl = 0;
    }

    return ray;
}

static const char *sarc_get_source_name(void *vstate)
{
    return ((sarc_state_t *) vstate)->name;
}

static int64_t sarc_get_source_n_rays(void *vstate)
{
    sarc_state_t *state = (sarc_state_t *) vstate;

    return per_thread_get_source_n_rays(&state->mutex_n_rays,
					&state->n_rays);
}

static double sarc_get_source_power(void *vstate)
{
    return ((sarc_state_t *) vstate)->power;
}

static void sarc_init_rays_remain(void *vstate)
{
    per_thread_init_rays_remain(((sarc_state_t *)
				 vstate)->rays_remain_key);
}


static const source_type_t sarc_t = {
    "arc",
    sizeof(struct sarc_state_t),
    &sarc_init_state,
    &sarc_free_state,
    &sarc_emit_ray,
    &sarc_get_source_name,
    &sarc_get_source_n_rays,
    &sarc_get_source_power,
    &sarc_init_rays_remain
};

const source_type_t *source_arc = &sarc_t;

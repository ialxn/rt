/*	source_sphere.c
 *
 * Copyright (C) 2011 - 2018 Ivo Alxneit, Paul Scherrer Institute
 *
 * This file is part of rt
 *
 * rt is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rt is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rt. If not, see <http://www.gnu.org/licenses/>.
 *
 */
#define _ISOC99_SOURCE		/* because of llrint() */

#include <math.h>
#include <string.h>

#include <gsl/gsl_machine.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include "io_utils.h"
#include "math_utils.h"
#include "ray.h"
#include "sources.h"

typedef struct sp_state_t {
    char *name;			/* name (identifier) of uniform point source */
    double orig[3];
    double radius;
    int64_t n_rays;		/* number of rays remaining until source is exhausted */
    double power;		/* power of source */
    gsl_spline *spectrum;	/* spline holding cdf of source spectrum */
    pthread_mutex_t mutex_n_rays;	/* protect n_rays */
    pthread_key_t rays_remain_key;	/* no of ray remain in group (PTD) */
} sp_state_t;


static void sp_init_state(void *vstate, config_setting_t * this_s,
			  const double P_factor)
{
    sp_state_t *state = (sp_state_t *) vstate;
    const char *S;

    config_setting_lookup_string(this_s, "name", &S);
    state->name = strdup(S);

    config_setting_lookup_float(this_s, "power", &state->power);
    state->n_rays = llrint(state->power / P_factor);
    pthread_mutex_init(&state->mutex_n_rays, NULL);

    config_setting_lookup_float(this_s, "radius", &state->radius);
    read_vector(this_s, "origin", state->orig);

    /* initialize source spectrum */
    init_source_spectrum(this_s, "spectrum", &state->spectrum);

    pthread_key_create(&state->rays_remain_key, free);
}

static void sp_free_state(void *vstate)
{
    sp_state_t *state = (sp_state_t *) vstate;

    free(state->name);
    gsl_spline_free(state->spectrum);
}

static ray_t *sp_emit_ray(void *vstate, const gsl_rng * r)
{
    sp_state_t *state = (sp_state_t *) vstate;
    ray_t *ray = NULL;
    int64_t *rays_remain = pthread_getspecific(state->rays_remain_key);

    if (!*rays_remain)
	*rays_remain =
	    per_thread_get_new_raygroup(&state->mutex_n_rays,
					&state->n_rays);

    if (*rays_remain > 0) {	/* rays still available in group */
	double R, point[3];

	(*rays_remain)--;
	pthread_setspecific(state->rays_remain_key, rays_remain);

	ray = (ray_t *) malloc(sizeof(ray_t));

	/*
	 * select random point inside unit sphere.
	 * there are two possible algorithms:
	 * 1) choose random 'r', 'theta', 'phi' and then
	 *    convert these to cartesian coordinates.
	 *    this requires the calculation of four
	 *    trigonometric functions (as below for
	 *    the direction vector) and one cube root
	 *    in addition to the three random numbers.
	 * 2) select random cartesian coordinates in
	 *    the range +-1 and discard numbers that
	 *    lie outside of the sphere. taking the
	 *    square root of the result is not needed
	 *    for the comparison with radius 1. to
	 *    further speed up things, restart if
	 *    x and y coordinates lie outside of
	 *    unit circle.
	 * 1) is generally faster (for optimization
	 *    better than O1 on an atom 450 processor)
	 */

	R = state->radius * sqrt(gsl_rng_uniform(r));
	get_uniform_random_vector(point, R, r);

	ray->orig[0] = state->orig[0] + point[0];
	ray->orig[1] = state->orig[1] + point[1];
	ray->orig[2] = state->orig[2] + point[2];

	/* choose random direction */
	get_uniform_random_vector(ray->dir, 1.0, r);

	ray->lambda =
	    gsl_spline_eval(state->spectrum, gsl_rng_uniform(r), NULL);
	ray->n_refl = 0;
    }

    return ray;
}

static const char *sp_get_source_name(void *vstate)
{
    return ((sp_state_t *) vstate)->name;
}

static int64_t sp_get_source_n_rays(void *vstate)
{
    sp_state_t *state = (sp_state_t *) vstate;

    return per_thread_get_source_n_rays(&state->mutex_n_rays,
					&state->n_rays);
}

static double sp_get_source_power(void *vstate)
{
    return ((sp_state_t *) vstate)->power;
}

static void sp_init_rays_remain(void *vstate)
{
    per_thread_init_rays_remain(((sp_state_t *) vstate)->rays_remain_key);
}


static const source_type_t sp_t = {
    "sphere",
    sizeof(struct sp_state_t),
    &sp_init_state,
    &sp_free_state,
    &sp_emit_ray,
    &sp_get_source_name,
    &sp_get_source_n_rays,
    &sp_get_source_power,
    &sp_init_rays_remain
};

const source_type_t *source_sphere = &sp_t;

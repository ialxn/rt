/*	source_solid_rod.c
 *
 * Copyright (C) 2015 Ivo Alxneit
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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include "io_utils.h"
#include "math_utils.h"
#include "off.h"
#include "ray.h"
#include "sources.h"

typedef struct srod_state_t {
    char *name;			/* name (identifier) of uniform point source */
    double orig[3];		/* origin (center of base face) */
    double radius;		/* radius of cylinder */
    double length;		/* length of cylinder */
    double barrier1;		/* barriers to select emiting surface */
    double barrier2;
    double M[9];		/* convert local <-> global coordinates */
    int64_t n_rays;		/* number of rays remaining until source is exhausted */
    double power;		/* power of source */
    gsl_spline *spectrum;	/* spline holding cdf of source spectrum */
    pthread_mutex_t mutex_n_rays;	/* protect n_rays */
    pthread_key_t rays_remain_key;	/* no of ray remain in group (PTD) */
} srod_state_t;

static double get_A_disk(srod_state_t * state, config_setting_t * this_s,
			 const char *kw)
{
    int ans;

    config_setting_lookup_bool(this_s, kw, &ans);
    if (ans)			/* surface emits */
	return M_PI * state->radius * state->radius;
    else			/* surface does not emit */
	return 0.0;
}

static void init_barriers(config_setting_t * this_s, srod_state_t * state)
{
    double sum = 0.0;
    double A_wall, A_base, A_top;

    A_wall = 2.0 * M_PI * state->radius * state->length;
    sum += A_wall;
    A_base = get_A_disk(state, this_s, "base_face_emits");
    sum += A_base;
    A_top = get_A_disk(state, this_s, "base_face_emits");
    sum += A_top;

    /*
     * a random number will be used to select face that emits ray:
     * 
     *  0 <= r < A_wall -> from wall
     *  A_wall <= r < (A_wall + A_base) -> from base disk
     *  r >= (A_wall + A_base)  -> from top disk
     */
    state->barrier1 = A_wall / sum;
    state->barrier2 = (A_wall + A_base) / sum;
}

static void srod_init_state(void *vstate, config_setting_t * this_s,
			    const double P_factor)
{
    srod_state_t *state = (srod_state_t *) vstate;
    const char *S;

    config_setting_lookup_string(this_s, "name", &S);
    state->name = strdup(S);

    config_setting_lookup_float(this_s, "power", &state->power);
    state->n_rays = llrint(state->power / P_factor);
    pthread_mutex_init(&state->mutex_n_rays, NULL);

    read_vector(this_s, "origin", state->orig);
    config_setting_lookup_float(this_s, "radius", &state->radius);
    config_setting_lookup_float(this_s, "length", &state->length);

    init_barriers(this_s, state);

    read_vector(this_s, "direction", &state->M[6]);
    state->M[0] = 1.0;
    state->M[1] = 0.0;
    state->M[2] = 0.0;
    orthonormalize(&state->M[0], &state->M[3], &state->M[6]);

    /* initialize source spectrum */
    init_source_spectrum(this_s, "spectrum", &state->spectrum);

    pthread_key_create(&state->rays_remain_key, free);
}

static void srod_free_state(void *vstate)
{
    srod_state_t *state = (srod_state_t *) vstate;

    free(state->name);
    gsl_spline_free(state->spectrum);
}

static ray_t *srod_emit_ray(void *vstate, const gsl_rng * r)
{
    srod_state_t *state = (srod_state_t *) vstate;
    ray_t *ray = NULL;
    int64_t *rays_remain = pthread_getspecific(state->rays_remain_key);

    if (!*rays_remain)
	*rays_remain =
	    per_thread_get_new_raygroup(&state->mutex_n_rays,
					&state->n_rays);

    if (*rays_remain > 0) {	/* rays still available in group */
	double point[3];
	double normal[3];
	double phi, sin_phi, cos_phi;
	double selector = gsl_rng_uniform(r);

	(*rays_remain)--;
	pthread_setspecific(state->rays_remain_key, rays_remain);

	ray = (ray_t *) malloc(sizeof(ray_t));

	phi = 2.0 * M_PI * gsl_rng_uniform(r);
	sincos(phi, &sin_phi, &cos_phi);

	if (selector < state->barrier1) {
	    /* ray originates from wall */
	    point[0] = state->radius * cos_phi;
	    point[1] = state->radius * sin_phi;
	    point[2] = state->length * gsl_rng_uniform(r);

	    normal[0] = point[0];
	    normal[1] = point[1];
	    normal[2] = 0.0;
	    normalize(normal);
	} else if (state->barrier2 < selector) {
	    /* ray originates from top disk */
	    double t = sqrt(gsl_rng_uniform(r)) * state->radius;

	    point[0] = t * sin_phi;
	    point[1] = t * cos_phi;
	    point[2] = state->length;

	    normal[0] = 0.0;
	    normal[1] = 0.0;
	    normal[2] = 1.0;
	} else {
	    /* ray originates from base disk */
	    double t = sqrt(gsl_rng_uniform(r)) * state->radius;

	    point[0] = t * sin_phi;
	    point[1] = t * cos_phi;
	    point[2] = 0.0;

	    normal[0] = 0.0;
	    normal[1] = 0.0;
	    normal[2] = -1.0;
	}

	l2g(state->M, state->orig, point, ray->orig);

	/*
	 * choose random direction
	 */
	get_uniform_random_vector_hemisphere(point, 1.0, normal, r);
	l2g_rot(state->M, point, ray->dir);

	ray->lambda =
	    gsl_spline_eval(state->spectrum, gsl_rng_uniform(r), NULL);
	ray->n_refl = 0;
    }

    return ray;
}

static const char *srod_get_source_name(void *vstate)
{
    return ((srod_state_t *) vstate)->name;
}

static int64_t srod_get_source_n_rays(void *vstate)
{
    srod_state_t *state = (srod_state_t *) vstate;

    return per_thread_get_source_n_rays(&state->mutex_n_rays,
					&state->n_rays);
}

static double srod_get_source_power(void *vstate)
{
    return ((srod_state_t *) vstate)->power;
}

static void srod_init_rays_remain(void *vstate)
{
    per_thread_init_rays_remain(((srod_state_t *)
				 vstate)->rays_remain_key);
}


static const source_type_t srod_t = {
    "solid rod",
    sizeof(struct srod_state_t),
    &srod_init_state,
    &srod_free_state,
    &srod_emit_ray,
    &srod_get_source_name,
    &srod_get_source_n_rays,
    &srod_get_source_power,
    &srod_init_rays_remain
};

const source_type_t *source_solid_rod = &srod_t;

/*	source_sphere.c
 *
 * Copyright (C) 2011 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <math.h>
#include <string.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include "io_util.h"
#include "ray.h"
#include "sources.h"

typedef struct sp_state_t {
    char *name;			/* name (identifier) of uniform point source */
    double origin[3];
    double radius;
    int n_rays;			/* number of rays remaining until source is exhausted */
    pthread_mutex_t mutex_n_rays;	/* protect n_rays */
    pthread_key_t rays_remain_key;	/* no of ray remain in group (PTD) */
    double power;		/* power of source */
    double ppr;			/* power allotted to one ray */
    gsl_spline *spline;		/* spline holding cdf of source spectrum */
    double lambda_min;		/* minimum wavelength in spectrum */
} sp_state_t;


static void sp_init_state(void *vstate, config_setting_t * this_s)
{
    sp_state_t *state = (sp_state_t *) vstate;

    const char *S;

    config_setting_lookup_string(this_s, "name", &S);
    state->name = strdup(S);
    config_setting_lookup_int(this_s, "n_rays", &state->n_rays);
    pthread_mutex_init(&state->mutex_n_rays, NULL);
    config_setting_lookup_float(this_s, "power", &state->power);
    state->ppr = state->power / state->n_rays;

    config_setting_lookup_float(this_s, "radius", &state->radius);
    read_vector(this_s, "origin", state->origin);

    /* initialize source spectrum */
    config_setting_lookup_string(this_s, "spectrum", &S);
    init_spectrum(S, &state->spline, &state->lambda_min);

    pthread_key_create(&state->rays_remain_key, free);
}

static void sp_free_state(void *vstate)
{
    sp_state_t *state = (sp_state_t *) vstate;

    free(state->name);
    gsl_spline_free(state->spline);
}

static ray_t *sp_get_new_ray(void *vstate, const gsl_rng * r)
{
    sp_state_t *state = (sp_state_t *) vstate;
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
	double sin_theta, cos_theta;
	double phi;
	double R, R_sin_theta;

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

	cos_theta = 1.0 - 2.0 * gsl_rng_uniform(r);
	sin_theta = sin(acos(cos_theta));
	phi = 2.0 * M_PI * gsl_rng_uniform(r);

	R = state->radius * sqrt(gsl_rng_uniform(r));
	R_sin_theta = R * sin_theta;

	ray->origin[0] = state->origin[0] + R_sin_theta * cos(phi);
	ray->origin[1] = state->origin[1] + R_sin_theta * sin(phi);
	ray->origin[2] = state->origin[2] + state->radius * cos_theta;

	/* choose random direction */
	cos_theta = 1.0 - 2.0 * gsl_rng_uniform(r);
	sin_theta = sin(acos(cos_theta));
	phi = 2.0 * M_PI * gsl_rng_uniform(r);

	ray->direction[0] = sin_theta * cos(phi);
	ray->direction[1] = sin_theta * sin(phi);
	ray->direction[2] = cos_theta;

	ray->power = state->ppr;

	/* choose random wavelength */
	ray->lambda =
	    state->lambda_min + gsl_spline_eval(state->spline,
						gsl_rng_uniform(r), NULL);
    }

    return ray;
}

static const char *sp_get_source_name(void *vstate)
{
    sp_state_t *state = (sp_state_t *) vstate;

    return state->name;
}

static double sp_get_source_ppr(void *vstate)
{
    sp_state_t *state = (sp_state_t *) vstate;

    return state->ppr;
}

static void sp_init_rays_remain(void *vstate)
{
    sp_state_t *state = (sp_state_t *) vstate;

    int *rays_remain = (int *) malloc(sizeof(int));

    *rays_remain = 0;
    pthread_setspecific(state->rays_remain_key, rays_remain);
}


static const source_type_t sp_t = {
    "uniform point source",
    sizeof(struct sp_state_t),
    &sp_init_state,
    &sp_free_state,
    &sp_get_new_ray,
    &sp_get_source_name,
    &sp_get_source_ppr,
    &sp_init_rays_remain
};

const source_type_t *source_sphere = &sp_t;

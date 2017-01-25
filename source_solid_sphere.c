/*	source_solid_sphere.c
 *
 * Copyright (C) 2011,2012,2013,2014,2015,2016,2017 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#define _ISOC99_SOURCE		/* because of llrint() */

#include <math.h>
#include <string.h>

#include <gsl/gsl_machine.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include "intercept.h"
#include "io_utils.h"
#include "likely.h"
#include "math_utils.h"
#include "ray.h"
#include "reflect.h"
#include "sources.h"
#include "targets.h"


typedef struct ssp_state_t {
    char *name;			/* name (identifier) of uniform point source */
    double orig[3];
    double radius;
    int64_t n_rays;		/* number of rays remaining until source is exhausted */
    double power;		/* power of source */
    gsl_spline *spectrum;	/* spline holding cdf of source spectrum */
    pthread_mutex_t mutex_n_rays;	/* protect n_rays */
    pthread_key_t rays_remain_key;	/* no of ray remain in group (PTD) */
} ssp_state_t;

typedef struct vtssp_state_t {
    double center[3];		/* center coordinate of sphere */
    double radius;		/* radius of sphere */
    double radius2;		/* radius^2 of sphere */
    gsl_spline *spectrum;	/* interpolated emission spectrum */
    gsl_spline *reflectivity;	/* interpolated reflectivity spectrum */
    refl_model_t *refl_model;	/* reflection model */
} vtssp_state_t;


static void random_ray_on_sphere(ray_t * ray, const double *origin,
				 const double radius,
				 const gsl_spline * spectrum,
				 const gsl_rng * r)
{
    /*
     * select random point on sphere with radius 'radius' and origin 'origin'
     * and initialize 'ray->origin' with it.
     */
    double point[3];

    get_uniform_random_vector(point, radius, r);

    ray->orig[0] = origin[0] + point[0];
    ray->orig[1] = origin[1] + point[1];
    ray->orig[2] = origin[2] + point[2];

    get_uniform_random_vector_hemisphere(ray->dir, 1.0, point, r);

    ray->lambda = gsl_spline_eval(spectrum, gsl_rng_uniform(r), NULL);
}

static void ssp_init_state(void *vstate, config_setting_t * this_s,
			   const double P_factor)
{
    ssp_state_t *state = (ssp_state_t *) vstate;
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
    int64_t *rays_remain = pthread_getspecific(state->rays_remain_key);

    if (!*rays_remain)
	*rays_remain =
	    per_thread_get_new_raygroup(&state->mutex_n_rays,
					&state->n_rays);

    if (*rays_remain > 0) {	/* rays still available in group */

	(*rays_remain)--;
	pthread_setspecific(state->rays_remain_key, rays_remain);

	ray = (ray_t *) malloc(sizeof(ray_t));
	random_ray_on_sphere(ray, state->orig, state->radius,
			     state->spectrum, r);
	ray->n_refl = 0;
    }

    return ray;
}

static const char *ssp_get_source_name(void *vstate)
{
    return ((ssp_state_t *) vstate)->name;
}

static int64_t ssp_get_source_n_rays(void *vstate)
{
    ssp_state_t *state = (ssp_state_t *) vstate;

    return per_thread_get_source_n_rays(&state->mutex_n_rays,
					&state->n_rays);
}

static double ssp_get_source_power(void *vstate)
{
    return ((ssp_state_t *) vstate)->power;
}

static void ssp_init_rays_remain(void *vstate)
{
    per_thread_init_rays_remain(((ssp_state_t *) vstate)->rays_remain_key);
}


static int vtssp_init_state(void *vstate, config_setting_t * this_target, const int
			    __attribute__ ((__unused__)) file_mode, const int
			    __attribute__ ((__unused__)) keep_closed, const double
			    __attribute__ ((__unused__)) P_factor)
{
    vtssp_state_t *state = (vtssp_state_t *) vstate;

    read_vector(this_target, "origin", state->center);
    config_setting_lookup_float(this_target, "radius", &state->radius);
    state->radius2 = state->radius * state->radius;

    init_source_spectrum(this_target, "spectrum", &state->spectrum);

    init_spectrum(this_target, "reflectivity", &state->reflectivity);
    state->refl_model = init_refl_model(this_target);

    return NO_ERR;
}

static void vtssp_free_state(void *vstate)
{
    vtssp_state_t *state = (vtssp_state_t *) vstate;

    gsl_spline_free(state->spectrum);
    gsl_spline_free(state->reflectivity);

    free(state->refl_model);
}

static double *vtssp_get_intercept(void *vstate, ray_t * ray)
{
    vtssp_state_t *state = (vtssp_state_t *) vstate;

    return intercept_sphere(ray, NULL, state->center,
			    state->radius2, 0.0, 0.0, NULL);
}

static ray_t *vtssp_get_out_ray(void *vstate, ray_t * ray, double *hit,
				const gsl_rng * r)
{
    vtssp_state_t *state = (vtssp_state_t *) vstate;

    if (gsl_rng_uniform(r) >
	gsl_spline_eval(state->reflectivity, ray->lambda, NULL)) {
	/*
	 * ray is absorbed. we emit the same power from a random point on the
	 * source but with the emissionspectrum of this source.
	 */
	random_ray_on_sphere(ray, state->center, state->radius,
			     state->spectrum, r);

	if (unlikely(ray->n_refl == UCHAR_MAX))	/* too many reflections is unlikely */
	    fprintf(stderr,
		    "                INFO: maximum path length exceeded (wrap around of counter occurs)\n");
	++ray->n_refl;

    } else {			/* reflect 'ray' */
	double normal[3];

	diff(normal, hit, state->center);
	normalize(normal);
	reflect_ray(ray, normal, hit, r, state->refl_model);
    }

    return ray;
}

static void vtssp_init_PTDT(void __attribute__ ((__unused__)) * vstate)
{
    return;
}


static const source_type_t ssp_t = {
    "solid sphere",
    sizeof(struct ssp_state_t),
    &ssp_init_state,
    &ssp_free_state,
    &ssp_emit_ray,
    &ssp_get_source_name,
    &ssp_get_source_n_rays,
    &ssp_get_source_power,
    &ssp_init_rays_remain
};

static const target_type_t vt_ssp_t = {
    NULL,
    sizeof(struct vtssp_state_t),
    &vtssp_init_state,
    &vtssp_free_state,
    &vtssp_get_intercept,
    &vtssp_get_out_ray,
    &vtssp_init_PTDT,
    NULL
};


const source_type_t *source_solid_sphere = &ssp_t;
const target_type_t *virtual_target_solid_sphere = &vt_ssp_t;

/*	virtual_target_solid_sphere.c
 *
 * Copyright (C) 2010,2011,2012,2013,2014,2015 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#include "io_utils.h"
#include "intercept.h"
#include "reflect.h"
#include "targets.h"
#include "virtual_targets.h"


typedef struct vtssp_state_t {
    double center[3];		/* center coordinate of sphere */
    double radius;		/* radius^2 of disk */
    gsl_spline *spectrum;	/* interpolated emission spectrum */
    gsl_spline *reflectivity;	/* interpolated reflectivity spectrum */
    refl_func_pointer_t refl_func;	/* reflection model */
    void *refl_func_pars;	/* model specific parameters */
} vtssp_state_t;


static int vtssp_init_state(void *vstate, config_setting_t * this_target,
			    const int file_mode, const int keep_closed,
			    const double P_factor)
{
    vtssp_state_t *state = (vtssp_state_t *) vstate;

    read_vector(this_target, "origin", state->center);
    config_setting_lookup_float(this_target, "radius", &state->radius);

    if (init_spectrum(this_target, "spectrum", &state->spectrum))
	return ERR;

    if (init_spectrum(this_target, "reflectivity", &state->reflectivity))
	return ERR;

    init_refl_model(this_target, &state->refl_func,
		    &state->refl_func_pars);

    return NO_ERR;
}

static void vtssp_free_state(void *vstate)
{
    vtssp_state_t *state = (vtssp_state_t *) vstate;

    gsl_spline_free(state->spectrum);
    gsl_spline_free(state->reflectivity);

    if (state->refl_func == reflect_microfacet_gaussian)
	free((double *) state->refl_func_pars);
}

static double *vtssp_get_intercept(void *vstate, ray_t * ray)
{
    vtssp_state_t *state = (vtssp_state_t *) vstate;

    return intercept_sphere(ray, NULL, state->center, state->radius, 0.0,
			    0.0, NULL);
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
	get_uniform_random_vector(ray->orig, state->radius, r);
	ray->orig[0] += state->center[0];
	ray->orig[1] += state->center[1];
	ray->orig[2] += state->center[2];

	get_uniform_random_vector_hemisphere(ray->dir, 1.0, ray->orig, r);

	ray->lambda =
	    gsl_spline_eval(state->spectrum, gsl_rng_uniform(r), NULL);

    } else {			/* reflect 'ray' */
	double normal[3];

	diff(normal, hit, state->center);
	normalize(normal);
	state->refl_func(ray, normal, hit, r, state->refl_func_pars);
    }

    return ray;
}

static void vtssp_init_PTDT(void *vstate)
{
    return;
}

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

const target_type_t *virtual_target_solid_sphere = &vt_ssp_t;

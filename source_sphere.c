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
    double power;		/* power of source */
    double ppr;			/* power allotted to one ray */
    gsl_spline *spline;		/* spline holding cdf of source spectrum */
    gsl_interp_accel *acc;	/* cached data for spline interpolation */
} sp_state_t;


static void sp_init_state(void *vstate, config_t * cfg, const char *name)
{
    sp_state_t *state = (sp_state_t *) vstate;

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

    config_setting_lookup_float(this_s, "radius", &state->radius);
    read_vector(this_s, "origin", state->origin);

    /* initialize source spectrum */
    config_setting_lookup_string(this_s, "spectrum", &S);
    init_spectrum(S, &state->spline, &state->acc);
}

static void sp_free_state(void *vstate)
{
    sp_state_t *state = (sp_state_t *) vstate;

    free(state->name);
    gsl_spline_free(state->spline);
    gsl_interp_accel_free(state->acc);
}

static ray_t *sp_get_new_ray(void *vstate, const gsl_rng * r)
{
    sp_state_t *state = (sp_state_t *) vstate;
    ray_t *ray = NULL;

    if (state->n_rays) {
	double sin_theta, cos_theta;
	double phi;
	double R, R_sin_theta;

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

	/* choose random wavelength */
	ray->power = state->ppr;
	ray->lambda =
	    gsl_spline_eval(state->spline, gsl_rng_uniform(r), state->acc);

	state->n_rays--;
    }

    return ray;
}

static const char *sp_get_source_name(void *vstate)
{
    sp_state_t *state = (sp_state_t *) vstate;

    return state->name;
}


static const source_type_t sp_t = {
    "uniform point source",
    sizeof(struct sp_state_t),
    &sp_init_state,
    &sp_free_state,
    &sp_get_new_ray,
    &sp_get_source_name
};

const source_type_t *source_sphere = &sp_t;

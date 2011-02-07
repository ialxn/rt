/*	source_spot.c
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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_cblas.h>

#include "io_util.h"
#include "off.h"
#include "ray.h"
#include "sources.h"

typedef struct sp_state_t {
    char *name;			/* name (identifier) of uniform point source */
    double origin[3];
    double dir[3];		/* direction of cone */
    double cos_theta;		/* cosine of opening angle of light cone */
    double alpha;		/* used to convert from local to global system */
    double beta;		/* used to convert from local to global system */
    int n_rays;			/* number of rays remaining until source is exhausted */
    double power;		/* power of source */
    double ppr;			/* power allotted to one ray */
    gsl_spline *spline;		/* spline holding cdf of source spectrum */
    gsl_interp_accel *acc;	/* cached data for spline interpolation */
} sp_state_t;


static void sp_init_state(void *vstate, config_t * cfg, const char *name)
{
    sp_state_t *state = (sp_state_t *) vstate;

    const double O[] = { 0.0, 0.0, 0.0 };
    double t[3];
    unsigned int i = 0;
    const char *S;
    FILE *spectrum;
    double *lambda;
    double *cumul_I;
    size_t n_lambda;
    config_setting_t *this_s;
    const config_setting_t *s = config_lookup(cfg, "sources");

    state->name = strdup(name);

    while (1) {			/* find setting for source 'name' */
	this_s = config_setting_get_elem(s, i);

	config_setting_lookup_string(this_s, "name", &S);
	if (strstr(S, name))
	    break;

	i++;
    }

    config_setting_lookup_int(this_s, "n_rays", &state->n_rays);
    config_setting_lookup_float(this_s, "power", &state->power);
    state->ppr = state->power / state->n_rays;
    config_setting_lookup_float(this_s, "theta", &state->cos_theta);
    state->cos_theta = cos(state->cos_theta / 180.0 * M_PI);

    read_vector(this_s, "origin", state->origin);
    read_vector_normalize(this_s, "direction", state->dir);

    /* determine alpha, beta. discard t */
    g2l_off(O, state->dir, t, &state->alpha, &state->beta);

    /* read source spectrum */
    config_setting_lookup_string(this_s, "spectrum", &S);
    spectrum = fopen(S, "r");
    read_data(spectrum, &lambda, &cumul_I, &n_lambda);
    fclose(spectrum);

    /*
     * calculate normalized cumulative spectrum.
     * this will be the cumulative distribution function
     * needed to obtain a random wavelength.
     */
    for (i = 1; i < n_lambda; i++)
	cumul_I[i] += cumul_I[i - 1];
    for (i = 0; i < n_lambda; i++)
	cumul_I[i] /= cumul_I[n_lambda - 1];

    /* cspline will be used to interpolate */
    state->acc = gsl_interp_accel_alloc();
    state->spline = gsl_spline_alloc(gsl_interp_cspline, n_lambda);

    /* 'cumul_I' -> x and 'lambda' -> y */
    gsl_spline_init(state->spline, cumul_I, lambda, n_lambda);

    free(lambda);
    free(cumul_I);

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
    double t;

    if (state->n_rays) {	/* source not exhausted */

	ray = (ray_t *) malloc(sizeof(ray_t));

	if (state->cos_theta < 1.0) {	/* theta != 0.0 */
	    const double O[] = { 0.0, 0.0, 0.0 };
	    double l_ray[3];	/* direction of ray in local system */
	    double sin_theta, cos_theta, phi;

	    t = gsl_rng_uniform(r);
	    cos_theta = t * (state->cos_theta - 1.0) + 1.0;	/* theta: random 0 .. theta */
	    sin_theta = sin(acos(cos_theta));

	    t = gsl_rng_uniform(r);
	    phi = 2.0 * M_PI * t;	/* phi: random 0 .. 2*pi */

	    l_ray[0] = sin_theta * cos(phi);
	    l_ray[1] = sin_theta * sin(phi);
	    l_ray[2] = cos_theta;

	    l2g_off(O, l_ray, ray->direction, state->alpha, state->beta);

	} else			/* all rays in direction 'dir' */
	    memcpy(ray->direction, state->dir, 3 * sizeof(double));

	/* copy / initialize rest of structure */
	memcpy(ray->origin, state->origin, 3 * sizeof(double));
	ray->power = state->ppr;
	t = gsl_rng_uniform(r);
	ray->lambda = gsl_spline_eval(state->spline, t, state->acc);

	state->n_rays--;
    }

    return ray;
}

static const char *sp_get_source_name(void
				      *vstate)
{
    sp_state_t *state = (sp_state_t *) vstate;
    return state->name;
}


static const source_type_t sp_t = {
    "spot source",
    sizeof(struct sp_state_t),
    &sp_init_state,
    &sp_free_state,
    &sp_get_new_ray,
    &sp_get_source_name
};

const source_type_t *source_spot = &sp_t;

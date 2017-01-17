/*	target_ellipsoid.c
 *
 * Copyright (C) 2011,2012,2013,2014,2015,2016,2017 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#include <cblas.h>
#include <math.h>
#include <string.h>

#include "math_utils.h"
#include "io_utils.h"
#include "intercept.h"
#include "targets.h"

#define TARGET_TYPE "ellipsoid"
#define NO_ITEMS 4


typedef struct ell_state_t {
    double center[3];		/* center coordinate, origin of local system */
    double axes[3];		/* a^2, b^2, c^2 parameters (semi axes) */
    double z_min, z_max;	/* range of valid values of 'z' in local system */
    double *M;			/* transform matrix local -> global coordinates */
    gsl_spline *refl_spectrum;	/* for interpolated reflectivity spectrum */
    refl_model_t *refl_model;	/* reflection models */
    union fh_t output;		/* output file handle or name */
    int flags;
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
} ell_state_t;


static int ell_init_state(void *vstate, config_setting_t * this_target,
			  const int file_mode, const int keep_closed,
			  const double P_factor)
{
    ell_state_t *state = (ell_state_t *) vstate;
    int i;

    read_vector(this_target, "center", state->center);

    state->M = init_M(this_target, "x", "z");

    read_vector(this_target, "axes", state->axes);
    for (i = 0; i < 3; i++)	/* we will use only 1/a^2, 1/b^2, 1/c^2 */
	state->axes[i] = 1.0 / (state->axes[i] * state->axes[i]);

    config_setting_lookup_float(this_target, "z_min", &state->z_min);
    config_setting_lookup_float(this_target, "z_max", &state->z_max);
    if (state->z_max < state->z_min)	/* safety */
	SWAP(state->z_max, state->z_min);

    state->flags = 0;
    if (keep_closed)
	state->flags |= KEEP_CLOSED;

    if (init_output
	(TARGET_TYPE, this_target, file_mode, P_factor, &state->output,
	 &state->flags, state->center, state->M) == ERR) {
	state->refl_spectrum = NULL;
	return ERR;
    }

    /* initialize reflectivity spectrum */
    init_spectrum(this_target, "reflectivity", &state->refl_spectrum);
    state->refl_model = init_refl_model(this_target);

    state->flags |= init_reflecting_surface(this_target);

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);

    return NO_ERR;
}

static void ell_free_state(void *vstate)
{
    ell_state_t *state = (ell_state_t *) vstate;

    state_free(state->output, state->flags, state->M,
	       state->refl_spectrum, state->refl_model);
}

static double *ell_get_intercept(void *vstate, ray_t * ray)
{
    ell_state_t *state = (ell_state_t *) vstate;

    double *intercept;
    int hits_outside;

    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & LAST_WAS_HIT) {
	/*
	 * ray starts on this target's convex side, no hit posible
	 */
	data->flag &= ~LAST_WAS_HIT;
	return NULL;
    }

    intercept =
	intercept_ellipsoid(ray, state->M, state->center, state->axes,
			    state->z_min, state->z_max, &hits_outside);

    if (!intercept)		/* ray does not hit target */
	return NULL;

    /*
     * mark as absorbed if non-reflecting surface is hit.
     * mark if convex (outside) surface is hit.
     */
    if ((!(state->flags & OUTSIDE) && hits_outside)
	|| (state->flags & OUTSIDE && !hits_outside))
	data->flag |= ABSORBED;

    if (hits_outside)
	data->flag |= ICPT_ON_CONVEX_SIDE;

    return intercept;
}

static ray_t *ell_get_out_ray(void *vstate, ray_t * ray, double *hit,
			      const gsl_rng * r)
{
    ell_state_t *state = (ell_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & ABSORBED
	|| (gsl_rng_uniform(r) >
	    gsl_spline_eval(state->refl_spectrum, ray->lambda, NULL))) {
	/*
	 * if ABSORBED is set we know ray has been absorbed
	 * because it was intercepted by a surface with absorptivity=1
	 * (reflectivity=0) e.g. the backside of the target. this was
	 * checked (and the flag was set) in 'xxx_get_intercept()'
	 * above.
	 * then we check if ray is absorbed because the reflectivity of
	 * the mirror surface is less than 1.0 (absorptivity > 0.0).
	 */

	if (state->flags & OUTPUT_REQUIRED)
	    store_xyz(state->output, state->flags, ray, hit, state->M,
		      state->center, data, &state->mutex_writefd);

	data->flag &= ~(LAST_WAS_HIT | ABSORBED | ICPT_ON_CONVEX_SIDE);	/* clear flags */

	free(ray);
	return NULL;

    } else {			/* reflect 'ray' */
	double l_N[3], N[3];
	double hit_local[3];

	g2l(state->M, state->center, hit, hit_local);	/* transform to local coordinates */
	ell_surf_normal(hit_local, state->axes, l_N);	/* normal vector local system */
	l2g_rot(state->M, l_N, N);	/* normal vector global system */

	if (!(state->flags & OUTSIDE))
	    cblas_dscal(3, -1.0, N, 1);	/* make normal point inwards */

	reflect_ray(ray, N, hit, r, state->refl_model);

	if (data->flag & ICPT_ON_CONVEX_SIDE) {
	    data->flag |= LAST_WAS_HIT;	/* mark as hit */
	    data->flag &= ~ICPT_ON_CONVEX_SIDE;
	}

	return ray;
    }
}

static void ell_init_PTDT(void *vstate)
{
    per_thread_init(((ell_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void ell_flush_PTDT_outbuf(void *vstate)
{
    ell_state_t *state = (ell_state_t *) vstate;

    per_thread_flush(state->output, state->flags, state->PTDT_key,
		     &state->mutex_writefd);
}


static const target_type_t ell_t = {
    TARGET_TYPE,
    sizeof(struct ell_state_t),
    &ell_init_state,
    &ell_free_state,
    &ell_get_intercept,
    &ell_get_out_ray,
    &ell_init_PTDT,
    &ell_flush_PTDT_outbuf
};

const target_type_t *target_ellipsoid = &ell_t;

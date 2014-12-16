/*	target_cone.c
 *
 * Copyright (C) 2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <math.h>
#include <string.h>

#include "io_utils.h"
#include "intercept.h"
#include "targets.h"

#define TARGET_TYPE "cone"
#define NO_ITEMS 5


typedef struct cone_state_t {
    double origin[3];		/* center of base disk, origin of local system */
    double H;			/* (full) height of cone */
    double tan2_a;		/* a=opening angle of cone */
    double z_max;		/* height at radius 'r' */
    double *M;			/* transform matrix local -> global coordinates */
    gsl_spline *refl_spectrum;	/* for interpolated reflectivity spectrum */
    void *refl_model_params;
    union fh_t output;		/* output file handle or name */
    int flags;
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
} cone_state_t;


static int cone_init_state(void *vstate, config_setting_t * this_target,
			   const int file_mode, const int keep_closed)
{
    cone_state_t *state = (cone_state_t *) vstate;

    const char *S;
    double R, r;

    read_vector(this_target, "origin", state->origin);

    state->M = init_M(this_target, "x", "axis");

    config_setting_lookup_float(this_target, "R", &R);
    config_setting_lookup_float(this_target, "r", &r);
    config_setting_lookup_float(this_target, "h", &state->z_max);

    state->H = R * state->z_max / (R - r);
    state->tan2_a = (R * R) / (state->H * state->H);

    state->flags = 0;
    if (keep_closed)
	state->flags |= KEEP_CLOSED;

    if (init_output
	(TARGET_TYPE, this_target, file_mode, &state->output,
	 &state->flags, state->origin, state->M) == ERR) {
	state->refl_spectrum = NULL;
	state->flags |= MODEL_NONE;
	return ERR;
    }

    /* initialize reflectivity spectrum */
    config_setting_lookup_string(this_target, "reflectivity", &S);
    if (init_spectrum(S, &state->refl_spectrum)) {
	state->flags |= MODEL_NONE;
	return ERR;
    }
    init_refl_model(this_target, &state->flags, &state->refl_model_params);

    state->flags |= init_reflecting_surface(this_target);

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);

    return NO_ERR;
}

static void cone_free_state(void *vstate)
{
    cone_state_t *state = (cone_state_t *) vstate;

    state_free(state->output, state->flags, state->M,
	       state->refl_spectrum, state->refl_model_params);
}

static double *cone_get_intercept(void *vstate, ray_t * ray)
{
    cone_state_t *state = (cone_state_t *) vstate;

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
	intercept_cone(ray, state->M, state->origin, state->tan2_a,
		       state->H, state->z_max, &hits_outside);

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

static ray_t *cone_get_out_ray(void *vstate, ray_t * ray, double *hit,
			       const gsl_rng * r)
{
    cone_state_t *state = (cone_state_t *) vstate;
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
		      state->origin, data, &state->mutex_writefd);

	data->flag &= ~(LAST_WAS_HIT | ABSORBED | ICPT_ON_CONVEX_SIDE);	/* clear flags */

	free(ray);
	return NULL;

    } else {			/* reflect 'in_ray' */
	double l_N[3], N[3];
	double hit_local[3];

	/*
	 * calculate normal vector 'N' @ 'hit':
	 * - convert 'hit' to local coordinates
	 * - calculate normal vector in local coordinates
	 * - transform normal vector to global coordinates
	 */
	g2l(state->M, state->origin, hit, hit_local);
	cone_surf_normal(hit_local, state->tan2_a, state->H, l_N);
	l2g_rot(state->M, l_N, N);

	if (!(state->flags & OUTSIDE))
	    cblas_dscal(3, -1.0, N, 1);	/* make normal point inwards */

	reflect(ray, N, hit, state->flags, r, state->refl_model_params);

	if (data->flag & ICPT_ON_CONVEX_SIDE) {
	    data->flag |= LAST_WAS_HIT;	/* mark as hit */
	    data->flag &= ~ICPT_ON_CONVEX_SIDE;
	}

	return ray;
    }
}

static void cone_init_PTDT(void *vstate)
{
    per_thread_init(((cone_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void cone_flush_PTDT_outbuf(void *vstate)
{
    cone_state_t *state = (cone_state_t *) vstate;

    per_thread_flush(state->output, state->flags, state->PTDT_key,
		     &state->mutex_writefd);
}

static const target_type_t cone_t = {
    TARGET_TYPE,
    sizeof(struct cone_state_t),
    &cone_init_state,
    &cone_free_state,
    &cone_get_intercept,
    &cone_get_out_ray,
    &cone_init_PTDT,
    &cone_flush_PTDT_outbuf
};

const target_type_t *target_cone = &cone_t;

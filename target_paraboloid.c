/*	target_paraboloid.c
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

#include "math_utils.h"
#include "intercept.h"
#include "io_utils.h"
#include "targets.h"

#define TARGET_TYPE "paraboloid"
#define NO_ITEMS 5


typedef struct par_state_t {
    double vertex[3];		/* vertex coordinate, origin of local system */
    double foc2;		/* 1/(2*focal length of paraboloid) */
    double foc4;		/* 1/(4*focal length of paraboloid) */
    double z_min, z_max;	/* range of valid values of 'z' in local system */
    double *M;			/* transform matrix local -> global coordinates */
    gsl_spline *refl_spectrum;	/* for interpolated reflectivity spectrum */
    char reflectivity_model;	/* reflectivity model used for this target */
    char reflecting_surface;
    void *refl_model_params;
    union fh_t output;		/* output file handle or name */
    int flags;
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
} par_state_t;


static int par_init_state(void *vstate, config_setting_t * this_target,
			  const int file_mode, const int keep_closed)
{
    par_state_t *state = (par_state_t *) vstate;

    const char *S;

    read_vector(this_target, "vertex", state->vertex);

    state->M = init_M(this_target, "x", "z");

    config_setting_lookup_float(this_target, "focal_length", &state->foc4);
    state->foc4 = 1.0 / (4.0 * state->foc4);
    state->foc2 = 2.0 * state->foc4;

    config_setting_lookup_float(this_target, "z_min", &state->z_min);
    config_setting_lookup_float(this_target, "z_max", &state->z_max);
    if (state->z_max < state->z_min)	/* safety */
	SWAP(state->z_max, state->z_min);

    state->flags = keep_closed;
    if (init_output
	(TARGET_TYPE, this_target, file_mode, &state->output,
	 &state->flags, state->vertex, state->M) == ERR) {
	state->refl_spectrum = NULL;
	state->reflectivity_model = MODEL_NONE;
	return ERR;
    }

    /* initialize reflectivity spectrum */
    config_setting_lookup_string(this_target, "reflectivity", &S);
    if (init_spectrum(S, &state->refl_spectrum)) {
	state->reflectivity_model = MODEL_NONE;
	return ERR;
    }
    init_refl_model(this_target, &state->reflectivity_model,
		    &state->refl_model_params);

    state->reflecting_surface = init_reflecting_surface(this_target);

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);

    return NO_ERR;
}

static void par_free_state(void *vstate)
{
    par_state_t *state = (par_state_t *) vstate;

    state_free(state->output, state->flags, state->M,
	       state->refl_spectrum, state->reflectivity_model,
	       state->refl_model_params);
}

static double *par_get_intercept(void *vstate, ray_t * ray)
{
    par_state_t *state = (par_state_t *) vstate;

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
	intercept_paraboloid(ray, state->M, state->vertex, state->foc2,
			     state->foc4, state->z_min, state->z_max,
			     &hits_outside);

    if (!intercept)		/* ray does not hit target */
	return NULL;

    /*
     * mark as absorbed if non-reflecting surface is hit.
     * mark if convex (outside) surface is hit.
     */
    if ((state->reflecting_surface == INSIDE && hits_outside)
	|| (state->reflecting_surface == OUTSIDE && !hits_outside))
	data->flag |= ABSORBED;

    if (hits_outside)
	data->flag |= ICPT_ON_CONVEX_SIDE;

    return intercept;
}

static ray_t *par_get_out_ray(void *vstate, ray_t * ray, double *hit,
			      const gsl_rng * r)
{
    par_state_t *state = (par_state_t *) vstate;
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
		      state->vertex, data, &state->mutex_writefd);

	data->flag &= ~(LAST_WAS_HIT | ABSORBED | ICPT_ON_CONVEX_SIDE);	/* clear flags */

	free(ray);
	return NULL;

    } else {			/* reflect 'ray' */
	double l_N[3], N[3];
	double hit_local[3];

	g2l(state->M, state->vertex, hit, hit_local);	/* transform to local coordinates */
	par_surf_normal(hit_local, state->foc2, l_N);	/* normal vector local system */
	l2g_rot(state->M, l_N, N);	/* normal vector global system */

	if (state->reflecting_surface == INSIDE)
	    cblas_dscal(3, -1.0, N, 1);	/* make normal point inwards */

	reflect(ray, N, hit, state->reflectivity_model, r,
		state->refl_model_params);

	if (data->flag & ICPT_ON_CONVEX_SIDE) {
	    data->flag |= LAST_WAS_HIT;	/* mark as hit */
	    data->flag &= ~ICPT_ON_CONVEX_SIDE;
	}

	return ray;
    }
}

static void par_init_PTDT(void *vstate)
{
    per_thread_init(((par_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void par_flush_PTDT_outbuf(void *vstate)
{
    par_state_t *state = (par_state_t *) vstate;

    per_thread_flush(state->output, state->flags, state->PTDT_key,
		     &state->mutex_writefd);
}


static const target_type_t par_t = {
    TARGET_TYPE,
    sizeof(struct par_state_t),
    &par_init_state,
    &par_free_state,
    &par_get_intercept,
    &par_get_out_ray,
    &par_init_PTDT,
    &par_flush_PTDT_outbuf
};

const target_type_t *target_paraboloid = &par_t;

/*	target_annulus.c
 *
 * Copyright (C) 2010,2011,2012,2013,2014 Ivo Alxneit
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

#define TARGET_TYPE "annulus"
#define NO_ITEMS 4


typedef struct ann_state_t {
    double point[3];		/* center coordinate of annulus */
    double R2;			/* inner radius^2 of annulus */
    double r2;			/* inner radius^2 of annulus */
    double *M;			/* transform matrix local -> global coordinates */
    gsl_spline *refl_spectrum;	/* for interpolated reflectivity spectrum */
    char reflectivity_model;	/* reflectivity model used for this target */
    void *refl_model_params;
    union fh_t output;		/* output file handle or name */
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
} ann_state_t;


static int ann_init_state(void *vstate, config_setting_t * this_target,
			  const int file_mode, const int keep_closed)
{
    ann_state_t *state = (ann_state_t *) vstate;

    const char *S;
    double t;

    read_vector(this_target, "P", state->point);

    state->M = init_M(this_target, "x", "N");

    if (init_output
	(TARGET_TYPE, this_target, file_mode, &state->output, state->point,
	 state->M) == ERR) {
	state->refl_spectrum = NULL;
	state->reflectivity_model = MODEL_NONE;
	return ERR;
    }

    /* initialize reflectivity spectrum */
    config_setting_lookup_string(this_target, "reflectivity", &S);
    if (init_refl_spectrum(S, &state->refl_spectrum)) {
	state->reflectivity_model = MODEL_NONE;
	return ERR;
    }
    init_refl_model(this_target, &state->reflectivity_model,
		    &state->refl_model_params);

    config_setting_lookup_float(this_target, "R", &t);
    state->R2 = t * t;
    config_setting_lookup_float(this_target, "r", &t);
    state->r2 = t * t;

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);

    return NO_ERR;
}

static void ann_free_state(void *vstate)
{
    ann_state_t *state = (ann_state_t *) vstate;

    state_free(state->output.fh, state->M, state->refl_spectrum,
	       state->reflectivity_model, state->refl_model_params);
}

static double *ann_get_intercept(void *vstate, ray_t * ray)
{
    ann_state_t *state = (ann_state_t *) vstate;

    double *intercept;
    double l_intercept[3];
    double r2_intercept;

    int hits_front;

    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & LAST_WAS_HIT) {	/* ray starts on this target, no hit posible */
	data->flag &= ~LAST_WAS_HIT;
	return NULL;
    }

    intercept =
	intercept_plane(ray, &state->M[6], state->point, &hits_front);

    if (!intercept)		/* ray does not hit target */
	return NULL;

    /* convert to local coordinates, origin is 'state->point' */
    g2l(state->M, state->point, intercept, l_intercept);

    /*
     * r2_intercep is squared distance from center of annulus to intercept
     * in the plane of the annulus. we are in local system.
     * compare r^2 to avoid sqrt()
     */
    r2_intercept =
	l_intercept[0] * l_intercept[0] + l_intercept[1] * l_intercept[1];
    if ((r2_intercept > state->R2)
	|| (r2_intercept < state->r2)) {

	/* hit not within boundaries */
	free(intercept);
	return NULL;

    } else {			/* hits disk between radius 'r' and 'R' */

	if (!hits_front)	/* hits rear side, absorbed */
	    data->flag |= ABSORBED;

	return intercept;

    }
}

static ray_t *ann_get_out_ray(void *vstate, ray_t * ray, double *hit,
			      const gsl_rng * r)
{
    ann_state_t *state = (ann_state_t *) vstate;
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

	if (state->output.fh != -1)
	    store_xy(state->output.fh, ray, hit, state->M, state->point,
		     data, &state->mutex_writefd);

	data->flag &= ~(LAST_WAS_HIT | ABSORBED);	/* clear flags */

	free(ray);
	return NULL;

    } else {			/* reflect 'ray' */
	reflect(ray, &state->M[6], hit, state->reflectivity_model, r,
		state->refl_model_params);

	data->flag |= LAST_WAS_HIT;	/* mark as hit */

	return ray;
    }
}

static void ann_init_PTDT(void *vstate)
{
    per_thread_init(((ann_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void ann_flush_PTDT_outbuf(void *vstate)
{
    ann_state_t *state = (ann_state_t *) vstate;

    per_thread_flush(state->output.fh, state->PTDT_key,
		     &state->mutex_writefd);
}

static const target_type_t ann_t = {
    TARGET_TYPE,
    sizeof(struct ann_state_t),
    &ann_init_state,
    &ann_free_state,
    &ann_get_intercept,
    &ann_get_out_ray,
    &ann_init_PTDT,
    &ann_flush_PTDT_outbuf
};

const target_type_t *target_annulus = &ann_t;

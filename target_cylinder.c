/*	target_cylinder.c
 *
 * Copyright (C) 2014 - 2018 Ivo Alxneit, Paul Scherrer Institute
 *
 * This file is part of rt
 *
 * rt is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rt is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rt. If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <cblas.h>
#include <math.h>
#include <string.h>

#include "io_utils.h"
#include "intercept.h"
#include "targets.h"

#define TARGET_TYPE "cylinder"
#define NO_ITEMS 4


typedef struct cyl_state_t {
    double C[3];		/* center of face 1, origin of local system */
    double r;			/* radius of cylinder */
    double l;			/* length of cylinder */
    double *M;			/* transform matrix local -> global coordinates */
    gsl_spline *refl_spectrum;	/* for interpolated reflectivity spectrum */
    refl_model_t *refl_model;	/* reflection models */
    union fh_t output;		/* output file handle or name */
    int flags;
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
} cyl_state_t;


static int cyl_init_state(void *vstate, config_setting_t * this_target,
			  const int file_mode, const int keep_closed,
			  const double P_factor)
{
    cyl_state_t *state = (cyl_state_t *) vstate;

    read_vector(this_target, "C", state->C);

    state->M = init_M(this_target, "x", "a");

    config_setting_lookup_float(this_target, "r", &state->r);
    config_setting_lookup_float(this_target, "l", &state->l);

    state->flags = 0;
    if (keep_closed)
	state->flags |= KEEP_CLOSED;

    if (init_output
	(TARGET_TYPE, this_target, file_mode, P_factor, &state->output,
	 &state->flags, state->C, state->M) == ERR) {
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

static void cyl_free_state(void *vstate)
{
    cyl_state_t *state = (cyl_state_t *) vstate;

    state_free(state->output, state->flags, state->M,
	       state->refl_spectrum, state->refl_model);
}

static double *cyl_get_intercept(void *vstate, ray_t * ray)
{
    cyl_state_t *state = (cyl_state_t *) vstate;

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
	intercept_cylinder(ray, state->C, &state->M[6], state->r, state->l,
			   &hits_outside);

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

static ray_t *cyl_get_out_ray(void *vstate, ray_t * ray, double *hit,
			      const gsl_rng * r)
{
    cyl_state_t *state = (cyl_state_t *) vstate;
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
		      state->C, data, &state->mutex_writefd);

	data->flag &= ~(LAST_WAS_HIT | ABSORBED | ICPT_ON_CONVEX_SIDE);	/* clear flags */

	free(ray);
	return NULL;

    } else {			/* reflect 'in_ray' */
	double N[3];		/* normal vector at 'hit' */

	cyl_surf_normal(hit, state->C, &state->M[6], state->r, N);

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

static void cyl_init_PTDT(void *vstate)
{
    per_thread_init(((cyl_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void cyl_flush_PTDT_outbuf(void *vstate)
{
    cyl_state_t *state = (cyl_state_t *) vstate;

    per_thread_flush(state->output, state->flags, state->PTDT_key,
		     &state->mutex_writefd);
}

static const target_type_t cyl_t = {
    TARGET_TYPE,
    sizeof(struct cyl_state_t),
    &cyl_init_state,
    &cyl_free_state,
    &cyl_get_intercept,
    &cyl_get_out_ray,
    &cyl_init_PTDT,
    &cyl_flush_PTDT_outbuf
};

const target_type_t *target_cylinder = &cyl_t;

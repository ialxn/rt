/*	target_disk.c
 *
 * Copyright (C) 2010 - 2018 Ivo Alxneit, Paul Scherrer Institute
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

#include <math.h>
#include <string.h>

#include "io_utils.h"
#include "intercept.h"
#include "targets.h"

#define TARGET_TYPE "disk"
#define NO_ITEMS 3


typedef struct disk_state_t {
    double point[3];		/* center coordinate of disk */
    double r2;			/* radius^2 of disk */
    double *M;			/* transform matrix local -> global coordinates */
    gsl_spline *refl_spectrum;	/* for interpolated reflectivity spectrum */
    refl_model_t *refl_model;	/* reflection models */
    union fh_t output;		/* output file handle or name */
    int flags;
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
} disk_state_t;

static int disk_init_state(void *vstate, config_setting_t * this_target,
			   const int file_mode, const int keep_closed,
			   const double P_factor)
{
    disk_state_t *state = (disk_state_t *) vstate;
    double t;

    read_vector(this_target, "P", state->point);

    state->M = init_M(this_target, "x", "N");

    state->flags = 0;
    if (keep_closed)
	state->flags |= KEEP_CLOSED;

    if (init_output
	(TARGET_TYPE, this_target, file_mode, P_factor, &state->output,
	 &state->flags, state->point, state->M) == ERR) {
	state->refl_spectrum = NULL;
	return ERR;
    }

    /* initialize reflectivity spectrum */
    init_spectrum(this_target, "reflectivity", &state->refl_spectrum);
    state->refl_model = init_refl_model(this_target);

    config_setting_lookup_float(this_target, "r", &t);
    state->r2 = t * t;

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);

    return NO_ERR;
}

static void disk_free_state(void *vstate)
{
    disk_state_t *state = (disk_state_t *) vstate;

    state_free(state->output, state->flags, state->M,
	       state->refl_spectrum, state->refl_model);
}

static double *disk_get_intercept(void *vstate, ray_t * ray)
{
    disk_state_t *state = (disk_state_t *) vstate;

    double *intercept;

    int hits_front;

    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & LAST_WAS_HIT) {	/* ray starts on this target, no hit posible */
	data->flag &= ~LAST_WAS_HIT;
	return NULL;
    }

    intercept =
	intercept_disk(ray, state->point, state->M, state->r2,
		       &hits_front);

    if (intercept && !hits_front)	/* hits rear side, absorbed */
	data->flag |= ABSORBED;

    return intercept;

}

static ray_t *disk_get_out_ray(void *vstate, ray_t * ray, double *hit,
			       const gsl_rng * r)
{
    disk_state_t *state = (disk_state_t *) vstate;
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
	    store_xy(state->output, state->flags, ray, hit, state->M,
		     state->point, data, &state->mutex_writefd);

	data->flag &= ~(LAST_WAS_HIT | ABSORBED);	/* clear flags */

	free(ray);
	return NULL;

    } else {			/* reflect 'ray' */
	reflect_ray(ray, &state->M[6], hit, r, state->refl_model);

	data->flag |= LAST_WAS_HIT;	/* mark as hit */

	return ray;
    }
}

static void disk_init_PTDT(void *vstate)
{
    per_thread_init(((disk_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void disk_flush_PTDT_outbuf(void *vstate)
{
    disk_state_t *state = (disk_state_t *) vstate;

    per_thread_flush(state->output, state->flags, state->PTDT_key,
		     &state->mutex_writefd);
}


static const target_type_t disk_t = {
    TARGET_TYPE,
    sizeof(struct disk_state_t),
    &disk_init_state,
    &disk_free_state,
    &disk_get_intercept,
    &disk_get_out_ray,
    &disk_init_PTDT,
    &disk_flush_PTDT_outbuf
};

const target_type_t *target_disk = &disk_t;

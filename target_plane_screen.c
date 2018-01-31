/*	target_plane_screen.c
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

#define TARGET_TYPE_A "one-sided plane screen"
#define TARGET_TYPE_B "two-sided plane screen"
#define NO_ITEMS 3


typedef struct ps_state_t {
    double point[3];		/* point on plane */
    char one_sided;		/* flag [one-sided|two-sided] */
    double *M;			/* transform matrix local -> global coordinates */
    union fh_t output;		/* output file handle or name */
    int flags;
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
} ps_state_t;


static int ps_init_state(void *vstate, config_setting_t * this_target,
			 const int file_mode, const int keep_closed,
			 const double P_factor)
{
    ps_state_t *state = (ps_state_t *) vstate;

    read_vector(this_target, "point", state->point);

    state->M = init_M(this_target, "x", "normal");

    state->flags = 0;
    if (keep_closed)
	state->flags |= KEEP_CLOSED;

    if (state->one_sided) {
	if (init_output
	    (TARGET_TYPE_A, this_target, file_mode, P_factor,
	     &state->output, &state->flags, state->point, state->M) == ERR)
	    return ERR;
    } else {
	if (init_output
	    (TARGET_TYPE_B, this_target, file_mode, P_factor,
	     &state->output, &state->flags, state->point, state->M) == ERR)
	    return ERR;
    }

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);

    return NO_ERR;
}

static int ps1_init_state(void *vstate, config_setting_t * this_target,
			  const int file_name, const int keep_closed,
			  const double P_factor)
{
    ps_state_t *state = (ps_state_t *) vstate;

    state->one_sided = 1;
    return ps_init_state(vstate, this_target, file_name, keep_closed,
			 P_factor);
}

static int ps2_init_state(void *vstate, config_setting_t * this_target,
			  const int file_name, const int keep_closed,
			  const double P_factor)
{
    ps_state_t *state = (ps_state_t *) vstate;

    state->one_sided = 0;
    return ps_init_state(vstate, this_target, file_name, keep_closed,
			 P_factor);
}

static void ps_free_state(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    state_free(state->output, state->flags, state->M, NULL, NULL);
}

static double *ps_get_intercept(void *vstate, ray_t * ray)
{
    ps_state_t *state = (ps_state_t *) vstate;

    double *intercept;

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

    /*
     * only accept intercepts if front side is hit
     */
    if (state->one_sided && !hits_front) {
	free(intercept);
	return NULL;
    }

    return intercept;
}

static ray_t *ps_get_out_ray(void *vstate, ray_t * ray, double *hit,
			     const gsl_rng
			     __attribute__ ((__unused__)) * r)
{
    ps_state_t *state = (ps_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (state->flags & OUTPUT_REQUIRED)
	store_xy(state->output, state->flags, ray, hit, state->M,
		 state->point, data, &state->mutex_writefd);

    data->flag &= ~LAST_WAS_HIT;	/* clear flag */

    memcpy(ray->orig, hit, 3 * sizeof(double));	/* update origin */
    return ray;
}

static void ps_init_PTDT(void *vstate)
{
    per_thread_init(((ps_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void ps_flush_PTDT_outbuf(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    per_thread_flush(state->output, state->flags, state->PTDT_key,
		     &state->mutex_writefd);
}


static const target_type_t ps1_t = {
    TARGET_TYPE_A,
    sizeof(struct ps_state_t),
    &ps1_init_state,
    &ps_free_state,
    &ps_get_intercept,
    &ps_get_out_ray,
    &ps_init_PTDT,
    &ps_flush_PTDT_outbuf
};

static const target_type_t ps2_t = {
    TARGET_TYPE_B,
    sizeof(struct ps_state_t),
    &ps2_init_state,
    &ps_free_state,
    &ps_get_intercept,
    &ps_get_out_ray,
    &ps_init_PTDT,
    &ps_flush_PTDT_outbuf
};

const target_type_t *target_plane_screen_one_sided = &ps1_t;
const target_type_t *target_plane_screen_two_sided = &ps2_t;

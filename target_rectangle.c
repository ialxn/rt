/*	target_rectangle.c
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

#define TARGET_TYPE "rectangle"
#define NO_ITEMS 4


typedef struct sq_state_t {
    double point[3];		/* center coordinate */
    double dx;			/* rectangle is '2*dx' times '2*dy' local coordinates */
    double dy;
    double *M;			/* transform matrix local -> global coordinates */
    gsl_spline *refl_spectrum;	/* for interpolated reflectivity spectrum */
    int reflectivity_model;	/* reflectivity model used for this target */
    void *refl_model_params;
    int dump_file;
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
} sq_state_t;


static void sq_init_state(void *vstate, config_setting_t * this_target,
			  const int file_mode)
{
    sq_state_t *state = (sq_state_t *) vstate;

    int i;
    const char *S;

    read_vector(this_target, "P1", state->point);
    /*
     * generate transform matrix M to convert
     * between local and global coordinates
     * l2g:   g(x, y, z) = MT l(x, y, z) + o(x, y, z)
     * g2l:   l(x, y, z) = M (g(x, y, z) - o(x, y, z))
     */
    /*
     * get the other two corner points that define the 'x' and 'y'
     * axis of the plane. note for the axes 'x'='P2'-'P1',
     * 'y'='P3'-'P1'.
     */
    state->M = (double *) malloc(9 * sizeof(double));
    read_vector(this_target, "P2", state->M);
    read_vector(this_target, "P3", &state->M[3]);

    for (i = 0; i < 3; i++) {
	state->M[i] -= state->point[i];
	state->M[3 + i] -= state->point[i];
    }

    /* make 'point' point to center of rectangle */
    for (i = 0; i < 3; i++)
	state->point[i] += (state->M[i] + state->M[3 + i]) / 2.0;

    state->dx = normalize(state->M) / 2.0;
    state->dy = normalize(&state->M[3]) / 2.0;

    cross_product(state->M, &state->M[3], &state->M[6]);

    /* initialize reflectivity spectrum */
    config_setting_lookup_string(this_target, "reflectivity", &S);
    init_spectrum(S, &state->refl_spectrum);
    init_refl_model(this_target, &state->reflectivity_model,
		    &state->refl_model_params);

    state->dump_file =
	init_output(file_mode, TARGET_TYPE, this_target, state->point,
		    state->M);

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);
}

static void sq_free_state(void *vstate)
{
    sq_state_t *state = (sq_state_t *) vstate;

    state_free(state->dump_file, state->M, state->refl_spectrum,
	       state->reflectivity_model, state->refl_model_params);
}

static double *sq_get_intercept(void *vstate, ray_t * ray)
{
    sq_state_t *state = (sq_state_t *) vstate;

    double *intercept;
    double l_intercept[3];
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

    if ((l_intercept[0] <= -state->dx) || (l_intercept[0] >= state->dx)
	|| (l_intercept[1] <= -state->dy)
	|| (l_intercept[1] >= state->dy)) {

	/* hit not within boundaries */
	free(intercept);
	return NULL;

    } else {			/* hits within target dimensions 'dx' times 'dy' */

	if (!hits_front)	/* hits rear side, absorbed */
	    data->flag |= ABSORBED;

	return intercept;

    }
}

static ray_t *sq_get_out_ray(void *vstate, ray_t * ray, double *hit,
			     const gsl_rng * r)
{
    sq_state_t *state = (sq_state_t *) vstate;
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

	if (state->dump_file != -1)
	    store_xy(state->dump_file, ray, hit, state->M, state->point,
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

static void sq_init_PTDT(void *vstate)
{
    per_thread_init(((sq_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void sq_flush_PTDT_outbuf(void *vstate)
{
    sq_state_t *state = (sq_state_t *) vstate;

    per_thread_flush(state->dump_file, state->PTDT_key,
		     &state->mutex_writefd);
}


static const target_type_t sq_t = {
    TARGET_TYPE,
    sizeof(struct sq_state_t),
    &sq_init_state,
    &sq_free_state,
    &sq_get_intercept,
    &sq_get_out_ray,
    &sq_init_PTDT,
    &sq_flush_PTDT_outbuf
};

const target_type_t *target_rectangle = &sq_t;

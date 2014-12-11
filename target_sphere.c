/*	target_sphere.c
 *
 * Copyright (C) 2010,2011,2012,2013,2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#include <string.h>

#include "io_utils.h"
#include "intercept.h"
#include "targets.h"

#define TARGET_TYPE "sphere"
#define NO_ITEMS 5


typedef struct sph_state_t {
    double origin[3];		/* center coordinate of sphere */
    double radius;		/* radius^2 of disk */
    double z_min, z_max;	/* range of valid values of 'z' in local system */
    double *M;			/* transform matrix local -> global coordinates */
    gsl_spline *refl_spectrum;	/* for interpolated reflectivity spectrum */
    void *refl_model_params;
    union fh_t output;		/* output file handle or name */
    int flags;
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
} sph_state_t;


static int sph_init_state(void *vstate, config_setting_t * this_target,
			  const int file_mode, const int keep_closed)
{
    sph_state_t *state = (sph_state_t *) vstate;

    const char *S;

    read_vector(this_target, "origin", state->origin);

    state->M = init_M(this_target, "x", "z");

    config_setting_lookup_float(this_target, "radius", &state->radius);

    config_setting_lookup_float(this_target, "z_min", &state->z_min);
    config_setting_lookup_float(this_target, "z_max", &state->z_max);
    if (state->z_max < state->z_min)	/* safety */
	SWAP(state->z_max, state->z_min);

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

static void sph_free_state(void *vstate)
{
    sph_state_t *state = (sph_state_t *) vstate;

    state_free(state->output, state->flags, state->M,
	       state->refl_spectrum, state->refl_model_params);
}


static double *sph_get_intercept(void *vstate, ray_t * ray)
{
    sph_state_t *state = (sph_state_t *) vstate;

    double *intercept;
    int hits_outside;

    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & LAST_WAS_HIT) {
	/*
	 * ray starts on this target, no hit possible.
	 * this test fails for ray initially emitted by
	 * this solid source (flag cannot be set by source)
	 */
	data->flag &= ~LAST_WAS_HIT;
	return NULL;
    }

    intercept =
	intercept_sphere(ray, state->M, state->origin, state->radius,
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

static ray_t *sph_get_out_ray(void *vstate, ray_t * ray, double *hit,
			      const gsl_rng * r)
{
    sph_state_t *state = (sph_state_t *) vstate;
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
	double O[] = { 0.0, 0.0, 0.0 };
	double hit_local[3];

	g2l(state->M, state->origin, hit, hit_local);	/* transform to local coordinates */
	sph_surf_normal(hit_local, l_N);	/* normal vector local system */
	l2g(state->M, O, l_N, N);	/* normal vector global system */

	if (!(state->flags & OUTSIDE))
	    a_times_const(N, N, -1.0);

	reflect(ray, N, hit, state->flags, r, state->refl_model_params);

	if (data->flag & ICPT_ON_CONVEX_SIDE) {
	    data->flag |= LAST_WAS_HIT;	/* mark as hit */
	    data->flag &= ~ICPT_ON_CONVEX_SIDE;
	}

	return ray;
    }

}

static void sph_init_PTDT(void *vstate)
{
    per_thread_init(((sph_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void sph_flush_PTDT_outbuf(void *vstate)
{
    sph_state_t *state = (sph_state_t *) vstate;

    per_thread_flush(state->output, state->flags, state->PTDT_key,
		     &state->mutex_writefd);
}

static const target_type_t sph_t = {
    TARGET_TYPE,
    sizeof(struct sph_state_t),
    &sph_init_state,
    &sph_free_state,
    &sph_get_intercept,
    &sph_get_out_ray,
    &sph_init_PTDT,
    &sph_flush_PTDT_outbuf
};

const target_type_t *target_sphere = &sph_t;

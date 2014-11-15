/*	target_disk.c
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
#include "targets.h"

#define TARGET_TYPE "disk"
#define NO_ITEMS 4


typedef struct disk_state_t {
    char reflectivity_model;	/* reflectivity model used for this target */
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
    int dump_file;
    double point[3];		/* center coordinate of disk */
    double r2;			/* radius^2 of disk */
    gsl_spline *refl_spectrum;	/* for interpolated reflectivity spectrum */
    double M[9];		/* transform matrix local -> global coordinates */
    void *refl_model_params;
} disk_state_t;

static void disk_init_state(void *vstate, config_setting_t * this_target,
			    const int file_mode)
{
    disk_state_t *state = (disk_state_t *) vstate;

    const char *S;
    double t;

    read_vector(this_target, "P", state->point);
    /*
     * generate transform matrix M to convert
     * between local and global coordinates
     * l2g:   g(x, y, z) = MT l(x, y, z) + o(x, y, z)
     * g2l:   l(x, y, z) = M (g(x, y, z) - o(x, y, z))
     */
    /* get normal vector of plane (serving also as basis vector z) */
    read_vector_normalize(this_target, "N", &state->M[6]);

    /* get basis vector x */
    read_vector_normalize(this_target, "x", state->M);

    /* state->M[3-5] = y = z cross x */
    cross_product(&state->M[6], state->M, &state->M[3]);

    /* initialize reflectivity spectrum */
    config_setting_lookup_string(this_target, "reflectivity", &S);
    init_refl_spectrum(S, &state->refl_spectrum);
    init_refl_model(this_target, &state->reflectivity_model,
		    &state->refl_model_params);

    config_setting_lookup_float(this_target, "r", &t);
    state->r2 = t * t;

    state->dump_file =
	init_output(file_mode, TARGET_TYPE, this_target, state->point,
		    state->M);

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);
}

static void disk_free_state(void *vstate)
{
    disk_state_t *state = (disk_state_t *) vstate;

    state_free(state->dump_file, state->refl_spectrum,
	       state->reflectivity_model, state->refl_model_params);
}

static double *disk_get_intercept(void *vstate, ray_t * ray)
{
    disk_state_t *state = (disk_state_t *) vstate;

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
     * r2_intercep is squared distance from center of disk to intercept
     * in the plane of the disk. we are in local system.
     * compare r^2 to avoid sqrt()
     */
    r2_intercept =
	l_intercept[0] * l_intercept[0] + l_intercept[1] * l_intercept[1];
    if (r2_intercept > state->r2) {

	/* hit not within boundaries */
	free(intercept);
	return NULL;

    } else {			/* hits within disk radius 'r' */

	if (!hits_front)	/* hits rear side, absorbed */
	    data->flag |= ABSORBED;

	return intercept;

    }
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

static void disk_init_PTDT(void *vstate)
{
    per_thread_init(((disk_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void disk_flush_PTDT_outbuf(void *vstate)
{
    disk_state_t *state = (disk_state_t *) vstate;

    per_thread_flush(state->dump_file, state->PTDT_key,
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

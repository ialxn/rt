/*	target_cylinder.c
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
#include "math_utils.h"
#include "targets.h"

#define TARGET_TYPE "cylinder"
#define NO_ITEMS 5

#define INSIDE 0
#define OUTSIDE 1


typedef struct cyl_state_t {
    char reflectivity_model;	/* reflectivity model used for this target */
    char reflecting_surface;
    void *refl_model_params;
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
    int dump_file;
    double C[3];		/* center of face 1, origin of local system */
    double a[3];		/* cylinder axis (local z) pointing to face 2 */
    double r;			/* radius of cylinder */
    double l;			/* length of cylinder */
    gsl_spline *refl_spectrum;	/* for interpolated reflectivity spectrum */
    double M[9];		/* transform matrix local -> global coordinates */
} cyl_state_t;


static void cyl_init_state(void *vstate, config_setting_t * this_target,
			   const int file_mode)
{
    cyl_state_t *state = (cyl_state_t *) vstate;

    const char *S;

    read_vector(this_target, "C", state->C);
    read_vector_normalize(this_target, "a", state->a);
    /*
     * generate transform matrix M to convert
     * between local and global coordinates
     * l2g:   g(x, y, z) = M l(x, y, z) + o(x, y, z))
     * g2l:   l(x, y, z) = MT (g(x, y, z) - o(x, y, z))
     */
    /*
     * get 'x' vector,
     * save normalized vector in the transformation matrix 'M'
     */
    read_vector_normalize(this_target, "x", state->M);
    /*
     * get 'z' vector,
     * save normalized vector in the transformation matrix 'M'
     */
    read_vector_normalize(this_target, "a", &state->M[6]);

    orthonormalize(state->M, &state->M[3], &state->M[6]);

    config_setting_lookup_float(this_target, "r", &state->r);
    config_setting_lookup_float(this_target, "l", &state->l);

    /* initialize reflectivity spectrum */
    config_setting_lookup_string(this_target, "reflectivity", &S);
    init_refl_spectrum(S, &state->refl_spectrum);
    init_refl_model(this_target, &state->reflectivity_model,
		    &state->refl_model_params);

    config_setting_lookup_string(this_target, "reflecting_surface", &S);
    if (!strcmp(S, "inside"))
	state->reflecting_surface = INSIDE;
    else
	state->reflecting_surface = OUTSIDE;

    state->dump_file =
	init_output(file_mode, TARGET_TYPE, this_target, state->C,
		    state->M);

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);
}

static void cyl_free_state(void *vstate)
{
    cyl_state_t *state = (cyl_state_t *) vstate;

    state_free(state->dump_file, state->refl_spectrum,
	       state->reflectivity_model, state->refl_model_params);
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
	intercept_cylinder(ray, state->C, state->a, state->r, state->l,
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
	store_xyz(state->dump_file, ray, hit, state->M, state->C,
		  data, &state->mutex_writefd);

	data->flag &= ~(LAST_WAS_HIT | ABSORBED | ICPT_ON_CONVEX_SIDE);	/* clear flags */

	free(ray);
	return NULL;

    } else {			/* reflect 'in_ray' */
	double l_N[3], N[3];
	double O[] = { 0.0, 0.0, 0.0 };
	double hit_local[3];

	g2l(state->M, state->C, hit, hit_local);	/* transform to local coordinates */
	cyl_surf_normal(hit_local, state->C, state->a, state->r, l_N);	/* normal vector local system */
	l2g(state->M, O, l_N, N);	/* normal vector global system */

	if (state->reflecting_surface == INSIDE) {
	    int i;

	    for (i = 0; i < 3; i++)
		N[i] = -N[i];
	}
	reflect(ray, N, hit, state->reflectivity_model, r,
		state->refl_model_params);

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

    per_thread_flush(state->dump_file, state->PTDT_key,
		     state->mutex_writefd);
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

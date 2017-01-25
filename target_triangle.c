/*	target_triangle.c
 *
 * Copyright (C) 2010,2011,2012,2013,2014,2015,2016,2017 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <cblas.h>
#include <math.h>
#include <string.h>

#include "io_utils.h"
#include "intercept.h"
#include "targets.h"

#define TARGET_TYPE "triangle"
#define NO_ITEMS 4


typedef struct tr_state_t {
    double P1[3];		/* corner point of triangle */
    double E2[3];		/* edge 'P2' - 'P1' */
    double E3[3];		/* edge 'P3' - 'P1' */
    double *M;			/* transform matrix local -> global coordinates */
    gsl_spline *refl_spectrum;	/* for interpolated reflectivity spectrum */
    refl_model_t *refl_model;	/* reflection models */
    union fh_t output;		/* output file handle or name */
    int flags;
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
} tr_state_t;


static int tr_init_state(void *vstate, config_setting_t * this_target,
			 const int file_mode, const int keep_closed,
			 const double P_factor)
{
    tr_state_t *state = (tr_state_t *) vstate;

    int i;
    config_setting_t *point;

    read_vector(this_target, "P1", state->P1);

    point = config_setting_get_member(this_target, "P2");
    for (i = 0; i < 3; i++)
	state->E2[i] =
	    config_setting_get_float_elem(point, i) - state->P1[i];

    point = config_setting_get_member(this_target, "P3");
    for (i = 0; i < 3; i++)
	state->E3[i] =
	    config_setting_get_float_elem(point, i) - state->P1[i];

    /*
     * generate transform matrix M to convert
     * between local and global coordinates
     * l2g:   g(x, y, z) = MT l(x, y, z) + o(x, y, z)
     * g2l:   l(x, y, z) = M (g(x, y, z) - o(x, y, z))
     */
    state->M = (double *) malloc(9 * sizeof(double));
    /* x = 'E2' */
    memcpy(state->M, state->E2, 3 * sizeof(double));
    normalize(state->M);

    /* z = 'E2' cross 'E3' */
    cross_product(state->E2, state->E3, &state->M[6]);
    normalize(&state->M[6]);

    /* y = state->M[3-5] = z cross x */
    cross_product(&state->M[6], state->M, &state->M[3]);

    state->flags = 0;
    if (keep_closed)
	state->flags |= KEEP_CLOSED;

    if (init_output
	(TARGET_TYPE, this_target, file_mode, P_factor, &state->output,
	 &state->flags, state->P1, state->M) == ERR) {
	state->refl_spectrum = NULL;
	return ERR;
    }

    /* initialize reflectivity spectrum */
    init_spectrum(this_target, "reflectivity", &state->refl_spectrum);
    state->refl_model = init_refl_model(this_target);

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);

    return NO_ERR;
}

static void tr_free_state(void *vstate)
{
    tr_state_t *state = (tr_state_t *) vstate;

    state_free(state->output, state->flags, state->M,
	       state->refl_spectrum, state->refl_model);
}

static double *tr_get_intercept(void *vstate, ray_t * ray)
{
    tr_state_t *state = (tr_state_t *) vstate;

    double *intercept;
    double P[3], Q[3], T[3];
    double u, v;
    double t;
    double det;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & LAST_WAS_HIT) {	/* ray starts on this target, no hit posible */
	data->flag &= ~LAST_WAS_HIT;
	return NULL;
    }

    /*
     * use barycentric coordinates (u,v) to determine whether ray is
     * intercepted by triangle.
     * implementation according to Möller and Trumbore
     */
    cross_product(ray->dir, state->E3, P);

    det = cblas_ddot(3, state->E2, 1, P, 1);

    if (fabs(det) < GSL_SQRT_DBL_EPSILON)	/* parallel to triangle */
	return NULL;

    diff(T, ray->orig, state->P1);
    u = cblas_ddot(3, T, 1, P, 1);
    if (u < 0.0 || u > det)	/* outside */
	return NULL;

    cross_product(T, state->E2, Q);
    v = cblas_ddot(3, ray->dir, 1, Q, 1);

    if (v < 0.0 || u + v > det)	/* outside */
	return NULL;

    t = cblas_ddot(3, state->E3, 1, Q, 1) / det;

/*
 * if (t<0)	ray points away from target
 *   return NULL;
 *
 * is this test needed?
 */
    intercept = (double *) malloc(3 * sizeof(double));

    a_plus_cb(intercept, ray->orig, t, ray->dir);

    if (det < 0.0)		/* hits rear side (parallel to surface normal) */
	data->flag |= ABSORBED;

    return intercept;

}

static ray_t *tr_get_out_ray(void *vstate, ray_t * ray, double *hit,
			     const gsl_rng * r)
{
    tr_state_t *state = (tr_state_t *) vstate;
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
		     state->P1, data, &state->mutex_writefd);

	data->flag &= ~(LAST_WAS_HIT | ABSORBED);	/* clear flags */

	free(ray);
	return NULL;

    } else {			/* reflect 'ray' */
	reflect_ray(ray, &state->M[6], hit, r, state->refl_model);

	data->flag |= LAST_WAS_HIT;	/* mark as hit */

	return ray;
    }
}

static void tr_init_PTDT(void *vstate)
{
    per_thread_init(((tr_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void tr_flush_PTDT_outbuf(void *vstate)
{
    tr_state_t *state = (tr_state_t *) vstate;

    per_thread_flush(state->output, state->flags, state->PTDT_key,
		     &state->mutex_writefd);
}


static const target_type_t tr_t = {
    TARGET_TYPE,
    sizeof(struct tr_state_t),
    &tr_init_state,
    &tr_free_state,
    &tr_get_intercept,
    &tr_get_out_ray,
    &tr_init_PTDT,
    &tr_flush_PTDT_outbuf
};

const target_type_t *target_triangle = &tr_t;

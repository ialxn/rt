/*	target_triangle.c
 *
 * Copyright (C) 2010,2011,2012,2013 Ivo Alxneit
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

#define NO_ITEMS 4


typedef struct tr_state_t {
    char *name;			/* name (identifier) of target */
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    int dump_file;
    double P1[3];		/* corner point of triangle */
    double E2[3];		/* edge 'P2' - 'P1' */
    double E3[3];		/* edge 'P3' - 'P1' */
    double normal[3];		/* normal vector of plane */
    gsl_spline *spline;		/* for interpolated reflectivity spectrum */
    double M[9];		/* transform matrix local -> global coordinates */
} tr_state_t;


static void tr_init_state(void *vstate, config_setting_t * this_target,
			  const int file_mode)
{
    tr_state_t *state = (tr_state_t *) vstate;

    int i;
    const char *S;
    char f_name[256];
    config_setting_t *point;

    config_setting_lookup_string(this_target, "name", &S);
    state->name = strdup(S);

    snprintf(f_name, 256, "%s.dat", state->name);
    state->dump_file =
	open(f_name, O_CREAT | O_WRONLY | file_mode, S_IRUSR | S_IWUSR);

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
     * l2g:   g(x, y, z) = M l(x, y, z) + o(x, y, z))
     * g2l:   l(x, y, z) = MT (g(x, y, z) - o(x, y, z))
     */
    /* x = 'E2' */
    for (i = 0; i < 3; i++)
	state->M[i] = state->E2[i];
    normalize(state->M);

    /* z = 'E2' cross 'E3' */
    cross_product(state->E2, state->E3, &state->M[6]);
    normalize(&state->M[6]);

    /* y = state->M[3-5] = z cross x */
    cross_product(&state->M[6], state->M, &state->M[3]);

    /* copy normal vector of plane */
    memcpy(state->normal, &state->M[6], 3 * sizeof(double));

    /* initialize reflectivity spectrum */
    config_setting_lookup_string(this_target, "reflectivity", &S);
    init_refl_spectrum(S, &state->spline);

    pthread_key_create(&state->PTDT_key, free_PTDT);
}

static void tr_free_state(void *vstate)
{
    tr_state_t *state = (tr_state_t *) vstate;

    close(state->dump_file);

    free(state->name);
    gsl_spline_free(state->spline);
}

static double *tr_get_intercept(void *vstate, ray_t * in_ray)
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
    cross_product(in_ray->direction, state->E3, P);

    det = cblas_ddot(3, state->E2, 1, P, 1);

    if (fabs(det) < GSL_SQRT_DBL_EPSILON)	/* parallel to triangle */
	return NULL;

    diff(T, in_ray->origin, state->P1);
    u = cblas_ddot(3, T, 1, P, 1);
    if (u < 0.0 || u > det)	/* outside */
	return NULL;

    cross_product(T, state->E2, Q);
    v = cblas_ddot(3, in_ray->direction, 1, Q, 1);

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

    a_plus_cb(intercept, in_ray->origin, t, in_ray->direction);

    if (det < 0.0)		/* hits rear side (parallel to surface normal) */
	data->flag |= ABSORBED;

    return intercept;

}

static ray_t *tr_get_out_ray(void *vstate, ray_t * in_ray, double *hit,
			     const gsl_rng * r)
{
    tr_state_t *state = (tr_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & ABSORBED
	|| (gsl_rng_uniform(r) >
	    gsl_spline_eval(state->spline, in_ray->lambda, NULL))) {
	/*
	 * if ABSORBED is set we know ray has been absorbed
	 * because it was intercepted by a surface with absorptivity=1
	 * (reflectivity=0) e.g. the backside of the target. this was
	 * checked (and the flag was set) in 'xxx_get_intercept()'
	 * above.
	 * then we check if ray is absorbed because the reflectivity of
	 * the mirror surface is less than 1.0 (absorptivity > 0.0).
	 */

	store_xy(state->dump_file, in_ray, hit, state->M, state->P1, data);

	data->flag &= ~(LAST_WAS_HIT | ABSORBED);	/* clear flags */

	free(in_ray);
	return NULL;

    } else {			/* reflect 'in_ray' */
	reflect(in_ray, state->normal, hit);

	data->flag |= LAST_WAS_HIT;	/* mark as hit */

	return in_ray;
    }
}

static const char *tr_get_target_name(void *vstate)
{
    tr_state_t *state = (tr_state_t *) vstate;

    return state->name;
}

static void tr_dump_string(void *vstate, const char *str)
{
    tr_state_t *state = (tr_state_t *) vstate;

    write(state->dump_file, str, strlen(str));
}

static double *tr_M(void *vstate)
{
    tr_state_t *state = (tr_state_t *) vstate;

    return state->M;
}

static void tr_init_PTDT(void *vstate)
{
    tr_state_t *state = (tr_state_t *) vstate;
    PTDT_t *data = (PTDT_t *) malloc(sizeof(PTDT_t));

    data->buf = (float *) malloc(BUF_SIZE * NO_ITEMS * sizeof(float));
    data->i = 0;
    data->flag = 0;

    pthread_setspecific(state->PTDT_key, data);
}

static void tr_flush_PTDT_outbuf(void *vstate)
{
    tr_state_t *state = (tr_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->i != 0)		/* write rest of buffer to file. */
	write(state->dump_file, data->buf, sizeof(float) * data->i);
}


static const target_type_t tr_t = {
    "triangle",
    sizeof(struct tr_state_t),
    &tr_init_state,
    &tr_free_state,
    &tr_get_intercept,
    &tr_get_out_ray,
    &tr_get_target_name,
    &tr_dump_string,
    &tr_M,
    &tr_init_PTDT,
    &tr_flush_PTDT_outbuf
};

const target_type_t *target_triangle = &tr_t;

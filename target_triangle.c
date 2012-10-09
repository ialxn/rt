/*	target_triangle.c
 *
 * Copyright (C) 2010, 2011 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <math.h>
#include <string.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include "io_util.h"
#include "ray.h"
#include "targets.h"
#include "vector_math.h"


#define N_COORDINATES 2		/* store only x,y */

typedef struct tr_state_t {
    char *name;			/* name (identifier) of target */
    char last_was_hit;		/* flag */
    char absorbed;
    FILE *dump_file;
    double P1[3];		/* corner point of triangle */
    double E2[3];		/* edge 'P2' - 'P1' */
    double E3[3];		/* edge 'P3' - 'P1' */
    double normal[3];		/* normal vector of plane */
    gsl_spline *spline;		/* for interpolated reflectivity spectrum */
    gsl_interp_accel *acc;	/* cache for spline */
    double M[9];		/* transform matrix local -> global coordinates */
    size_t n_alloc;		/* buffer 'data' can hold 'n_alloc' data sets */
    size_t n_data;		/* buffer 'data' currently holds 'n_data' sets */
    double *data;		/* buffer to store hits */
} tr_state_t;


static int tr_alloc_state(void *vstate)
{
    tr_state_t *state = (tr_state_t *) vstate;

    /* 4 items per data set (x,y,ppr,lambda) */
    if (!
	(state->data =
	 (double *) malloc((N_COORDINATES + 2) * BLOCK_SIZE *
			   sizeof(double))))
	return ERR;

    state->n_alloc = BLOCK_SIZE;

    return NO_ERR;
}

static void tr_init_state(void *vstate, config_setting_t * this_target,
			  config_t * cfg, const char *file_mode)
{
    tr_state_t *state = (tr_state_t *) vstate;

    int i;
    const char *S;
    char f_name[256];
    config_setting_t *point;

    (void) cfg;			/* avoid warning: unused parameter 'cfg' */

    config_setting_lookup_string(this_target, "name", &S);
    state->name = strdup(S);

    snprintf(f_name, 256, "%s.dat", state->name);
    state->dump_file = fopen(f_name, file_mode);

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
    init_refl_spectrum(S, &state->spline, &state->acc);

    state->last_was_hit = 0;
    state->absorbed = 0;
    state->n_data = 0;
}

static void tr_free_state(void *vstate)
{
    tr_state_t *state = (tr_state_t *) vstate;

    fclose(state->dump_file);

    free(state->name);
    free(state->data);
    gsl_spline_free(state->spline);
    gsl_interp_accel_free(state->acc);
}

static double *tr_get_intercept(void *vstate, ray_t * in_ray,
				int *dump_flag)
{
    tr_state_t *state = (tr_state_t *) vstate;

    double *intercept;
    double P[3], Q[3], T[3];
    double u, v;
    double t;
    double det;

    if (state->last_was_hit) {	/* ray starts on this target, definitely no hit */
	state->last_was_hit = 0;
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

    v_diff(T, in_ray->origin, state->P1);
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

    v_a_plus_cb(intercept, in_ray->origin, t, in_ray->direction);

    if (det < 0.0)		/* hits rear side (parallel to surface normal) */
	state->absorbed = 1;

    return intercept;

}

static ray_t *tr_get_out_ray(void *vstate, ray_t * in_ray, double *hit,
			     const gsl_rng * r)
{
    tr_state_t *state = (tr_state_t *) vstate;

    if (state->absorbed
	|| (gsl_rng_uniform(r) >
	    gsl_spline_eval(state->spline, in_ray->lambda, state->acc))) {
	/*
	 * if 'state->absorbed'is true we know ray has been absorbed
	 * because it was intercepted by a surface with absorptivity=1
	 * (reflectivity=0) e.g. the backside of the target. this was
	 * checked (and 'state->absorbed' was set) in 'xxx_get_intercept()'
	 * above.
	 * then we check if ray is absorbed because the reflectivity of
	 * the mirror surface is less than 1.0 (absorptivity > 0.0).
	 */
	double hit_local[3];

	/* transform to local coordinates */
	memcpy(hit_local, hit, 3 * sizeof(double));
	g2l(state->M, state->P1, hit, hit_local);

	/*
	 * store 4 items per data set (x,y,ppr,lambda)
	 * first x,y then ppr,lambda
	 */
	fprintf(state->dump_file, "%g\t%g\t%g\t%g\n", hit_local[0],
		hit_local[1], in_ray->power, in_ray->lambda);

	state->absorbed = 0;	/* reset flags */
	state->last_was_hit = 0;

	free(in_ray);
	return NULL;

    } else {			/* reflect 'in_ray' */
	reflect(in_ray, state->normal, hit);

	state->last_was_hit = 1;	/* mark as hit */

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

    fprintf(state->dump_file, "%s", str);
}

static double *tr_M(void *vstate)
{
    tr_state_t *state = (tr_state_t *) vstate;

    return state->M;
}


static const target_type_t tr_t = {
    "triangle",
    sizeof(struct tr_state_t),
    &tr_alloc_state,
    &tr_init_state,
    &tr_free_state,
    &tr_get_intercept,
    &tr_get_out_ray,
    &tr_get_target_name,
    &tr_dump_string,
    &tr_M
};

const target_type_t *target_triangle = &tr_t;

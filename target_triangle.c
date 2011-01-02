/*	target_triangle.c
 *
 * Copyright (C) 2010 Ivo Alxneit
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

#include "io_util.h"
#include "ray.h"
#include "targets.h"

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
    double M[9];		/* transform matrix local -> global coordinates */
    size_t n_alloc;		/* buffer 'data' can hold 'n_alloc' data sets */
    size_t n_data;		/* buffer 'data' currently holds 'n_data' sets */
    double *data;		/* buffer to store hits */
} tr_state_t;


static int tr_alloc_state(void *vstate)
{
    tr_state_t *state = (tr_state_t *) vstate;

    /* 3 items per data set (x,y,ppr) */
    if (!
	(state->data =
	 (double *) malloc((N_COORDINATES + 1) * BLOCK_SIZE *
			   sizeof(double))))
	return ERR;

    state->n_alloc = BLOCK_SIZE;

    return NO_ERR;
}

static void tr_init_state(void *vstate, config_t * cfg, const char *name,
			  const char *file_mode)
{
    tr_state_t *state = (tr_state_t *) vstate;

    unsigned int i = 0;
    const char *S;
    double norm;
    char f_name[256];
    int j;

    config_setting_t *this_target, *point;
    const config_setting_t *targets = config_lookup(cfg, "targets");

    state->name = strdup(name);

    state->last_was_hit = 0;

    snprintf(f_name, 256, "%s.dat", name);
    state->dump_file = fopen(f_name, file_mode);

    while (1) {			/* find setting for target 'name' */
	this_target = config_setting_get_elem(targets, i);

	config_setting_lookup_string(this_target, "name", &S);
	if (!strcmp(S, name))
	    break;

	i++;
    }

    point = config_setting_get_member(this_target, "P1");
    for (j = 0; j < 3; j++)
	state->P1[j] = config_setting_get_float_elem(point, j);

    point = config_setting_get_member(this_target, "P2");
    for (j = 0; j < 3; j++)
	state->E2[j] =
	    config_setting_get_float_elem(point, j) - state->P1[j];

    point = config_setting_get_member(this_target, "P3");
    for (j = 0; j < 3; j++)
	state->E3[j] =
	    config_setting_get_float_elem(point, j) - state->P1[j];

    /*
     * generate transform matrix M to convert
     * between local and global coordinates
     * l2g:   g(x, y, z) = M l(x, y, z) + o(x, y, z))
     * g2l:   l(x, y, z) = MT (g(x, y, z) - o(x, y, z))
     */
    /* x = 'E2' */
    for (j = 0; j < 3; j++)
	state->M[j] = state->E2[j];
    norm = cblas_dnrm2(3, state->M, 1);
    cblas_dscal(3, 1.0 / norm, state->M, 1);


    /* z = 'E2' cross 'E3' */
    cross_product(state->E2, state->E3, &state->M[6]);

    /* normalize z */
    norm = cblas_dnrm2(3, &state->M[6], 1);
    cblas_dscal(3, 1.0 / norm, &state->M[6], 1);

    /* y = state->M[3-5] = z cross x */
    cross_product(&state->M[6], state->M, &state->M[3]);

    /* copy normal vector of plane */
    memcpy(state->normal, &state->M[6], 3 * sizeof(double));

    state->last_was_hit = 0;
    state->absorbed = 0;
    state->n_data = 0;
}

static void tr_free_state(void *vstate)
{
    tr_state_t *state = (tr_state_t *) vstate;

    /* first write remaining data to file. 3 items per data (x,y,ppr) */
    dump_data(state->dump_file, state->data, state->n_data,
	      N_COORDINATES + 1);
    fclose(state->dump_file);

    free(state->name);
    free(state->data);
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
    int i;

    if (*dump_flag) {		/* we are in a dump cycle and have not yet written data */
	dump_data(state->dump_file, state->data, state->n_data,
		  N_COORDINATES + 1);
	shrink_memory(&(state->data), &(state->n_data), &(state->n_alloc),
		      N_COORDINATES + 1);
	(*dump_flag)--;
    }

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

    for (i = 0; i < 3; i++)
	T[i] = in_ray->origin[i] - state->P1[i];

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

    intercept[0] = in_ray->origin[0] + t * in_ray->direction[0];
    intercept[1] = in_ray->origin[1] + t * in_ray->direction[1];
    intercept[2] = in_ray->origin[2] + t * in_ray->direction[2];

    if (det < 0.0)		/* hits rear side (parallel to surface normal) */
	state->absorbed = 1;

    return intercept;

}

static ray_t *tr_get_out_ray(void *vstate, ray_t * in_ray,
			     const double ppr, double *hit,
			     int *dump_flag, const int n_targets)
{
    tr_state_t *state = (tr_state_t *) vstate;

    if (state->absorbed) {	/* rear surface was hit */
	double hit_copy[3];

	/* transform to local coordinates */
	memcpy(hit_copy, hit, 3 * sizeof(double));
	g2l(state->M, state->P1, hit, hit_copy);

	/*
	 * store 3 items per data set (x,y,ppr)
	 * first x,y then ppr
	 */
	memcpy(&(state->data[(N_COORDINATES + 1) * state->n_data]),
	       hit_copy, N_COORDINATES * sizeof(double));
	state->data[(N_COORDINATES + 1) * state->n_data + N_COORDINATES] =
	    ppr;
	state->n_data++;
	state->last_was_hit = 1;	/* mark as hit */

	/*
	 * increase data size for next interception
	 * or
	 * initiate dump cycle
	 */
	if (state->n_data == state->n_alloc)	/* buffer full */
	    try_increase_memory(&(state->data), &(state->n_data),
				&(state->n_alloc), N_COORDINATES + 1,
				state->dump_file, dump_flag, n_targets);

	state->absorbed = 0;	/* reset flag */
	free(in_ray);
	return NULL;

    } else {			/* reflect 'in_ray' */
	const double t = cblas_ddot(3, state->normal, 1, in_ray->direction, 1);	/* 'N' dot 'in_ray' */

	/* 'in_ray' - 2 * 'N' dot 'in_ray' * 'N' */
	cblas_daxpy(3, -2.0 * t, state->normal, 1, in_ray->direction, 1);
	memcpy(in_ray->origin, hit, 3 * sizeof(double));	/* update origin */

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

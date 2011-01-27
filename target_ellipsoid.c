/*	target_ellipsoid.c
 *
 * Copyright (C) 2011 Ivo Alxneit
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
#include "vector_math.h"

#define N_COORDINATES 3		/* store x,y,z */

typedef struct ell_state_t {
    char *name;			/* name (identifier) of target */
    FILE *dump_file;
    double center[3];		/* center coordinate, origin of local system */
    double axes[3];		/* a^2, b^2, c^2 parameters (semi axes) */
    double z_min, z_max;	/* range of valid values of 'z' in local system */
    double reflectivity;	/* reflectivity of target */
    int absorbed;		/* flag to indicated that ray was absorbed */
    double M[9];		/* transform matrix local -> global coordinates */
    size_t n_alloc;		/* buffer 'data' can hold 'n_alloc' data sets */
    size_t n_data;		/* buffer 'data' currently holds 'n_data' sets */
    double *data;		/* buffer to store hits */
} ell_state_t;


static void ell_surf_normal(const double *point, const double *axes,
			    double *const normal)
{
    int i;
    double norm = 0.0;

    for (i = 0; i < 3; i++) {
	normal[i] = 2.0 * point[i] / axes[i];
	norm += normal[i] * normal[i];
    }
    norm = sqrt(norm);
    cblas_dscal(3, 1.0 / norm, normal, 1);

}

static int ell_alloc_state(void *vstate)
{
    ell_state_t *state = (ell_state_t *) vstate;

    /* 4 items per data set (x,y,z,ppr) */
    if (!
	(state->data =
	 (double *) malloc((N_COORDINATES + 1) * BLOCK_SIZE *
			   sizeof(double))))
	return ERR;

    state->n_alloc = BLOCK_SIZE;

    return NO_ERR;
}

static void ell_init_state(void *vstate, config_t * cfg, const char *name,
			   const char *file_mode)
{
    ell_state_t *state = (ell_state_t *) vstate;

    unsigned int i = 0;
    const char *S;
    char f_name[256];
    config_setting_t *this_target;
    const config_setting_t *targets = config_lookup(cfg, "targets");

    state->name = strdup(name);

    snprintf(f_name, 256, "%s.dat", name);
    state->dump_file = fopen(f_name, file_mode);

    while (1) {			/* find setting for target 'name' */
	this_target = config_setting_get_elem(targets, i);

	config_setting_lookup_string(this_target, "name", &S);
	if (!strcmp(S, name))
	    break;

	i++;
    }

    read_vector(this_target, "center", state->center);
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
    read_vector_normalize(this_target, "z", &state->M[6]);

    orthonormalize(state->M, &state->M[3], &state->M[6]);

    read_vector(this_target, "axes", state->axes);
    for (i = 0; i < 3; i++)	/* we will use only a^2, b^2, c^2 */
	state->axes[i] *= state->axes[i];

    config_setting_lookup_float(this_target, "z_min", &state->z_min);
    config_setting_lookup_float(this_target, "z_max", &state->z_max);
    config_setting_lookup_float(this_target, "reflectivity",
				&state->reflectivity);
    state->absorbed = 0;
    state->n_data = 0;
}

static void ell_free_state(void *vstate)
{
    ell_state_t *state = (ell_state_t *) vstate;

    /* first write remaining data to file. 3 items per data (x,y,ppr) */
    dump_data(state->dump_file, state->data, state->n_data,
	      N_COORDINATES + 1);
    fclose(state->dump_file);

    free(state->name);
    free(state->data);
}

static double *ell_get_intercept(void *vstate, ray_t * in_ray,
				 int *dump_flag, const gsl_rng * r)
{
    ell_state_t *state = (ell_state_t *) vstate;

    int i;
    double r_O[3], r_N[3];	/* origin, direction of ray in local system */
    double O[] = { 0.0, 0.0, 0.0 };
    double A = 0.0, B = 0.0, C = -1.0;
    double D;

    if (*dump_flag) {		/* we are in a dump cycle and have not yet written data */
	dump_data(state->dump_file, state->data, state->n_data,
		  N_COORDINATES + 1);
	shrink_memory(&(state->data), &(state->n_data), &(state->n_alloc),
		      N_COORDINATES + 1);
	(*dump_flag)--;
    }

    /*
     * calculate point of interception D
     */
    /*
     * transform 'in_ray' from global to local system
     * origin 'in_ray': rotate / translate by origin of local system
     * dir 'in_ray': rotate only
     */
    g2l(state->M, state->center, in_ray->origin, r_O);
    g2l(state->M, O, in_ray->direction, r_N);

    /*
     * solve quadratic equation
     */
    for (i = 0; i < 3; i++) {
	A += r_N[i] * r_N[i] / state->axes[i];
	B += 2.0 * r_O[i] * r_N[i] / state->axes[i];
	C += r_O[i] * r_O[i] / state->axes[i];
    }

    D = B * B - 4.0 * A * C;

    if (D < GSL_SQRT_DBL_EPSILON)	/* no or one (tangent ray) interception */
	return NULL;
    else {			/* two interceptions, h- and h+ */
	const double t = sqrt(D);
	const double A2 = 2.0 * A;
	double h, z;
	/*
	 * - 'A' must be positive as it is the sum of squares.
	 * - the solution that contains the term '-t' must be smaller
	 *   i.e. further towards -inf.
	 * - the ray travels in forward direction thus only positive
	 *   solutions are of interest.
	 * thus:
	 * 1) check if the solution (h-) with '-t' is positive and
	 *    its z component 'z' is inside the allowed range, i.e.
	 *    'z_min' <= 'z' <= 'z_max'.
	 * 2) only if 1) does not produce a valid solution repeat with
	 *    '+t' (h+)
	 * 3) if 2) produces no valid solution return 'NULL'
	 */

	h = (-B - t) / A2;
	if (h <= GSL_SQRT_DBL_EPSILON)	/* h- not positive / valid */
	    h = (-B + t) / A2;
	if (h <= GSL_SQRT_DBL_EPSILON)	/* h+ not positive / valid */
	    return NULL;

	z = r_O[2] + h * r_N[2];	/* z component of h */
	if (z < state->z_min || z > state->z_max)
	    return NULL;	/* z component outside range */
	else {
	    double dot;
	    double l_intercept[3];
	    double normal[3];
	    double *intercept = (double *) malloc(3 * sizeof(double));

	    l_intercept[0] = r_O[0] + h * r_N[0];
	    l_intercept[1] = r_O[1] + h * r_N[1];
	    l_intercept[2] = z;	/* use precomputed value */

	    ell_surf_normal(l_intercept, state->axes, normal);
	    dot = cblas_ddot(3, normal, 1, r_N, 1);
	    if (dot < 0.0)	/* anti-parallel. hits outside, absorbed */
		state->absorbed = 1;
	    else if (gsl_rng_uniform(r) > state->reflectivity)	/* hits inside */
		state->absorbed = 1;

	    /* convert to global coordinates, origin is 'state->center' */
	    l2g(state->M, state->center, l_intercept, intercept);

	    return intercept;
	}
    }
}

static ray_t *ell_get_out_ray(void *vstate, ray_t * in_ray,
			      const double ppr, double *hit,
			      int *dump_flag, const int n_targets)
{
    ell_state_t *state = (ell_state_t *) vstate;
    double l_hit[3];

    /* transform to local coordinates */
    g2l(state->M, state->center, hit, l_hit);

    if (state->absorbed) {
	/*
	 * store 4 items per data set (x,y,ppr)
	 * first x,y,z then ppr
	 */
	memcpy(&(state->data[(N_COORDINATES + 1) * state->n_data]),
	       l_hit, N_COORDINATES * sizeof(double));
	state->data[(N_COORDINATES + 1) * state->n_data + N_COORDINATES] =
	    ppr;
	state->n_data++;

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
	double l_N[3], N[3];
	double O[] = { 0.0, 0.0, 0.0 };

	ell_surf_normal(l_hit, state->axes, l_N);	/* normal vector local system */
	l2g(state->M, O, l_N, N);	/* normal vector global system */

	reflect(in_ray, N, hit);

	return in_ray;

    }
}

static const char *ell_get_target_name(void *vstate)
{
    ell_state_t *state = (ell_state_t *) vstate;

    return state->name;
}

static void ell_dump_string(void *vstate, const char *str)
{
    ell_state_t *state = (ell_state_t *) vstate;

    fprintf(state->dump_file, "%s", str);
}

static double *ell_M(void *vstate)
{
    ell_state_t *state = (ell_state_t *) vstate;

    return state->M;
}


static const target_type_t ell_t = {
    "ellipsoid",
    sizeof(struct ell_state_t),
    &ell_alloc_state,
    &ell_init_state,
    &ell_free_state,
    &ell_get_intercept,
    &ell_get_out_ray,
    &ell_get_target_name,
    &ell_dump_string,
    &ell_M
};

const target_type_t *target_ellipsoid = &ell_t;

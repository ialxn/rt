/*	target_rectangle.c
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

typedef struct sq_state_t {
    char *name;			/* name (identifier) of target */
    char last_was_hit;		/* flag */
    FILE *dump_file;
    double point[3];		/* center coordinate */
    double dx;			/* rectangle is '2*dx' times '2*dy' local coordinates */
    double dy;
    gsl_spline *spline;		/* for interpolated reflectivity spectrum */
    gsl_interp_accel *acc;	/* cache for spline */
    int absorbed;		/* flag to indicated hit on rear surface == absorbed */
    double normal[3];		/* normal vector of plane */
    double M[9];		/* transform matrix local -> global coordinates */
    size_t n_alloc;		/* buffer 'data' can hold 'n_alloc' data sets */
    size_t n_data;		/* buffer 'data' currently holds 'n_data' sets */
    double *data;		/* buffer to store hits */
} sq_state_t;


static int sq_alloc_state(void *vstate)
{
    sq_state_t *state = (sq_state_t *) vstate;

    /* 4 items per data set (x,y,ppr,lambda) */
    if (!
	(state->data =
	 (double *) malloc((N_COORDINATES + 2) * BLOCK_SIZE *
			   sizeof(double))))
	return ERR;

    state->n_alloc = BLOCK_SIZE;

    return NO_ERR;
}

static void sq_init_state(void *vstate, config_setting_t * this_target,
			  config_t * cfg, const char *file_mode)
{
    sq_state_t *state = (sq_state_t *) vstate;

    int i;
    const char *S;
    char f_name[256];

    config_setting_lookup_string(this_target, "name", &S);
    state->name = strdup(S);

    snprintf(f_name, 256, "%s.dat", state->name);
    state->dump_file = fopen(f_name, file_mode);

    read_vector(this_target, "P1", state->point);
    /*
     * generate transform matrix M to convert
     * between local and global coordinates
     * l2g:   g(x, y, z) = M l(x, y, z) + o(x, y, z))
     * g2l:   l(x, y, z) = MT (g(x, y, z) - o(x, y, z))
     */
    /*
     * get the other two corner points that define the 'x' and 'y'
     * axis of the plane. note for the axes 'x'='P2'-'P1',
     * 'y'='P3'-'P1'.
     */
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

    /* state->normal = state->M[6,7,8] = z = x cross y */
    cross_product(state->M, &state->M[3], state->normal);
    memcpy(&state->M[6], state->normal, 3 * sizeof(double));

    /* initialize reflectivity spectrum */
    config_setting_lookup_string(this_target, "reflectivity", &S);
    init_refl_spectrum(S, &state->spline, &state->acc);

    state->last_was_hit = 0;
    state->absorbed = 0;
    state->n_data = 0;
}

static void sq_free_state(void *vstate)
{
    sq_state_t *state = (sq_state_t *) vstate;

    /* first write remaining data to file. 4 items per data (x,y,ppr,lambda) */
    dump_data(state->dump_file, state->data, state->n_data,
	      N_COORDINATES + 2);
    fclose(state->dump_file);

    free(state->name);
    free(state->data);
    gsl_spline_free(state->spline);
    gsl_interp_accel_free(state->acc);
}

static double *sq_get_intercept(void *vstate, ray_t * in_ray,
				int *dump_flag, const gsl_rng * r)
{
    sq_state_t *state = (sq_state_t *) vstate;

    double t1, t2[3], t3;
    double d;

    if (*dump_flag) {		/* we are in a dump cycle and have not yet written data */
	dump_data(state->dump_file, state->data, state->n_data,
		  N_COORDINATES + 2);
	shrink_memory(&(state->data), &(state->n_data), &(state->n_alloc),
		      N_COORDINATES + 2);
	(*dump_flag)--;
    }

    if (state->last_was_hit) {	/* ray starts on this target, definitely no hit */
	state->last_was_hit = 0;
	return NULL;
    }

    /*
     * calculate point of interception d
     *
     * d = {(\mathbf{p_0}-\mathbf{l_0})\cdot\mathbf{n} \over \mathbf{l}\cdot\mathbf{n}}
     *
     * with
     *       p_0: point on the plane
     *         n: normal vector of the plane (|n|=1)
     *       l_0: origin of the line
     *         l: unit vector in direction of the line
     *
     * If the line starts outside the plane and is parallel to the plane, there is no intersection.
     * In this case, the above denominator will be zero and the numerator will be non-zero. If the
     * line starts inside the plane and is parallel to the plane, the line intersects the plane
     * everywhere. In this case, both the numerator and denominator above will be zero. In all other
     * cases, the line intersects the plane once and d represents the intersection as the distance
     * along the line from \mathbf{l_0}.
     */

    t1 = cblas_ddot(3, in_ray->direction, 1, state->normal, 1);	/* l dot n */
    if (fabs(t1) < GSL_SQRT_DBL_EPSILON)	/* line is parallel to target, no hit possible */
	return NULL;

    t2[0] = state->point[0] - in_ray->origin[0];	/* p_0 - l_0 */
    t2[1] = state->point[1] - in_ray->origin[1];
    t2[2] = state->point[2] - in_ray->origin[2];

    t3 = cblas_ddot(3, t2, 1, state->normal, 1);	/* (p_0 - l_0) dot N */
    if (fabs(t3) < GSL_SQRT_DBL_EPSILON)	/* line does start in target, conservative */
	return NULL;

    d = t3 / t1;
    if (d < 0.0)		/* intercepted target is not in front */
	return NULL;
    else {			/* 'in_ray' intercepts target plane */
	double l_intercept[3];
	double *intercept = (double *) malloc(3 * sizeof(double));

	intercept[0] = in_ray->origin[0] + d * in_ray->direction[0];
	intercept[1] = in_ray->origin[1] + d * in_ray->direction[1];
	intercept[2] = in_ray->origin[2] + d * in_ray->direction[2];

	/* convert to local coordinates, origin is 'state->point' */
	g2l(state->M, state->point, intercept, l_intercept);

	if ((l_intercept[0] <= -state->dx) || (l_intercept[0] >= state->dx)
	    || (l_intercept[1] <= -state->dy)
	    || (l_intercept[1] >= state->dy)) {

	    /* hit not within boundaries */
	    free(intercept);
	    return NULL;

	} else {		/* hits within target dimensions 'dx' times 'dy' */

	    if (t1 > 0.0)	/* hits rear side, absorbed */
		state->absorbed = 1;
	    else if (gsl_rng_uniform(r) >
		     gsl_spline_eval(state->spline, in_ray->lambda,
				     state->acc))
		state->absorbed = 1;

	    return intercept;

	}
    }

}

static ray_t *sq_get_out_ray(void *vstate, ray_t * in_ray, double *hit,
			     int *dump_flag, const int n_targets)
{
    sq_state_t *state = (sq_state_t *) vstate;

    if (state->absorbed) {	/* rear surface was hit */
	double hit_copy[3];

	/* transform to local coordinates */
	memcpy(hit_copy, hit, 3 * sizeof(double));
	g2l(state->M, state->point, hit, hit_copy);

	/*
	 * store 4 items per data set (x,y,ppr,lambda)
	 * first x,y then ppr,lambda
	 */
	memcpy(&(state->data[(N_COORDINATES + 2) * state->n_data]),
	       hit_copy, N_COORDINATES * sizeof(double));
	state->data[(N_COORDINATES + 2) * state->n_data + N_COORDINATES] =
	    in_ray->power;
	state->data[(N_COORDINATES + 2) * state->n_data + N_COORDINATES +
		    1] = in_ray->lambda;
	state->n_data++;

	/*
	 * increase data size for next interception
	 * or
	 * initiate dump cycle
	 */
	if (state->n_data == state->n_alloc)	/* buffer full */
	    try_increase_memory(&(state->data), &(state->n_data),
				&(state->n_alloc), N_COORDINATES + 2,
				state->dump_file, dump_flag, n_targets);

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

static const char *sq_get_target_name(void *vstate)
{
    sq_state_t *state = (sq_state_t *) vstate;

    return state->name;
}

static void sq_dump_string(void *vstate, const char *str)
{
    sq_state_t *state = (sq_state_t *) vstate;

    fprintf(state->dump_file, "%s", str);
}

static double *sq_M(void *vstate)
{
    sq_state_t *state = (sq_state_t *) vstate;

    return state->M;
}


static const target_type_t sq_t = {
    "rectangle",
    sizeof(struct sq_state_t),
    &sq_alloc_state,
    &sq_init_state,
    &sq_free_state,
    &sq_get_intercept,
    &sq_get_out_ray,
    &sq_get_target_name,
    &sq_dump_string,
    &sq_M
};

const target_type_t *target_rectangle = &sq_t;

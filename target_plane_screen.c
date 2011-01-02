/*	target_plane_screen.c
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

#include "io_util.h"
#include "ray.h"
#include "targets.h"

#define N_COORDINATES 2		/* store only x,y */

typedef struct ps_state_t {
    char *name;			/* name (identifier) of target */
    char last_was_hit;		/* flag */
    char one_sided;		/* flag [one-sided|two-sided] */
    FILE *dump_file;
    double point[3];		/* point on plane */
    double normal[3];		/* normal vector of plane */
    double M[9];		/* transform matrix local -> global coordinates */
    size_t n_alloc;		/* buffer 'data' can hold 'n_alloc' data sets */
    size_t n_data;		/* buffer 'data' currently holds 'n_data' sets */
    double *data;		/* buffer to store hits */
} ps_state_t;


static int ps_alloc_state(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    /* 3 items per data set (x,y,ppr) */
    if (!
	(state->data =
	 (double *) malloc((N_COORDINATES + 1) * BLOCK_SIZE *
			   sizeof(double))))
	return ERR;

    state->n_alloc = BLOCK_SIZE;

    return NO_ERR;
}

static void ps_init_state(void *vstate, config_t * cfg, const char *name,
			  const char *file_mode)
{
    ps_state_t *state = (ps_state_t *) vstate;

    unsigned int i = 0;
    const char *S;
    double norm;
    char f_name[256];
    int j;

    config_setting_t *this_target, *point, *normal, *x;
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

    point = config_setting_get_member(this_target, "point");
    for (j = 0; j < 3; j++)
	state->point[j] = config_setting_get_float_elem(point, j);

    /*
     * generate transform matrix M to convert
     * between local and global coordinates
     * l2g:   g(x, y, z) = M l(x, y, z) + o(x, y, z))
     * g2l:   l(x, y, z) = MT (g(x, y, z) - o(x, y, z))
     */
    /* get normal vector of plane (serving also as basis vector z) */
    normal = config_setting_get_member(this_target, "normal");
    for (j = 0; j < 3; j++)
	state->normal[j] = config_setting_get_float_elem(normal, j);

    /* normalize normal vector */
    norm = cblas_dnrm2(3, state->normal, 1);
    cblas_dscal(3, 1.0 / norm, state->normal, 1);

    memcpy(&state->M[6], state->normal, 3 * sizeof(double));

    /* get basis vector x */
    x = config_setting_get_member(this_target, "x");
    for (j = 0; j < 3; j++)
	state->M[j] = config_setting_get_float_elem(x, j);
    /* normalize basis vector x */
    norm = cblas_dnrm2(3, state->M, 1);
    cblas_dscal(3, 1.0 / norm, state->M, 1);

    cross_product(&state->M[6], state->M, &state->M[3]);
    /* state->M[3-5] = z cross x */
    state->n_data = 0;
}

static void ps1_init_state(void *vstate, config_t * cfg, const char *name,
			   const char *file_name)
{
    ps_state_t *state = (ps_state_t *) vstate;

    ps_init_state(vstate, cfg, name, file_name);
    state->one_sided = 1;
}

static void ps2_init_state(void *vstate, config_t * cfg, const char *name,
			   const char *file_name)
{
    ps_state_t *state = (ps_state_t *) vstate;

    ps_init_state(vstate, cfg, name, file_name);
    state->one_sided = 0;
}

static void ps_free_state(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    /* first write remaining data to file. 3 items per data (x,y,ppr) */
    dump_data(state->dump_file, state->data, state->n_data,
	      N_COORDINATES + 1);
    fclose(state->dump_file);

    free(state->name);
    free(state->data);
}

static double *ps_get_intercept(void *vstate, ray_t * in_ray,
				int *dump_flag, const gsl_rng * r)
{
    ps_state_t *state = (ps_state_t *) vstate;

    double t1, t2[3], t3;
    double d;
    double *intercept;

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
    /*
     * in case of a one-sided target, check that t1 is negative.
     * if l and n are parallel (dot product is positive),
     * the intersection does not count
     */
    if (state->one_sided && (t1 > 0.0))
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

    intercept = (double *) malloc(3 * sizeof(double));

    intercept[0] = in_ray->origin[0] + d * in_ray->direction[0];
    intercept[1] = in_ray->origin[1] + d * in_ray->direction[1];
    intercept[2] = in_ray->origin[2] + d * in_ray->direction[2];

    return intercept;

}

static ray_t *ps_get_out_ray(void *vstate, ray_t * in_ray,
			     const double ppr, double *hit,
			     int *dump_flag, const int n_targets)
{
    ps_state_t *state = (ps_state_t *) vstate;

    ray_t *out;
    double hit_copy[3];

    /* transform to local coordinates */
    memcpy(hit_copy, hit, 3 * sizeof(double));
    g2l(state->M, state->point, hit, hit_copy);

    /*
     * store 3 items per data set (x,y,ppr)
     * first x,y then ppr
     */
    memcpy(&(state->data[(N_COORDINATES + 1) * state->n_data]), hit_copy,
	   N_COORDINATES * sizeof(double));
    state->data[(N_COORDINATES + 1) * state->n_data + N_COORDINATES] = ppr;
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

    out = (ray_t *) malloc(sizeof(ray_t));

    memcpy(out->origin, hit, 3 * sizeof(double));	/* update origin */
    memcpy(out->direction, in_ray->direction, 3 * sizeof(double));	/* copy direction */
    out->power = in_ray->power;

    free(in_ray);

    return out;
}

static const char *ps_get_target_name(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    return state->name;
}

static void ps_dump_string(void *vstate, const char *str)
{
    ps_state_t *state = (ps_state_t *) vstate;

    fprintf(state->dump_file, "%s", str);
}

static double *ps_M(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    return state->M;
}


static const target_type_t ps1_t = {
    "one-sided plane screen",
    sizeof(struct ps_state_t),
    &ps_alloc_state,
    &ps1_init_state,
    &ps_free_state,
    &ps_get_intercept,
    &ps_get_out_ray,
    &ps_get_target_name,
    &ps_dump_string,
    &ps_M
};

static const target_type_t ps2_t = {
    "two-sided plane screen",
    sizeof(struct ps_state_t),
    &ps_alloc_state,
    &ps2_init_state,
    &ps_free_state,
    &ps_get_intercept,
    &ps_get_out_ray,
    &ps_get_target_name,
    &ps_dump_string,
    &ps_M
};

const target_type_t *target_plane_screen_one_sided = &ps1_t;
const target_type_t *target_plane_screen_two_sided = &ps2_t;

/*	target_plane_screen.c
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <string.h>

#include "ray.h"
#include "targets.h"

typedef struct ps_state_t {
    char *name;
    FILE *dump_file;
    double point[3];
    double normal[3];
    size_t n_alloc;
    size_t n_data;
    double *data;
} ps_state_t;


static int ps_alloc_state(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    if (!
	(state->data = (double *) malloc(3 * BLOCK_SIZE * sizeof(double))))
	return ERR;

    state->n_alloc = BLOCK_SIZE;

    return NO_ERR;
}

static void ps_init_state(void *vstate, config_t * cfg, const char *name)
{
    ps_state_t *state = (ps_state_t *) vstate;

    unsigned int i = 0;
    const char *S;
    double F, norm;
    char f_name[256];
    gsl_vector_view N;

    config_setting_t *this_target, *point, *normal;
    const config_setting_t *targets = config_lookup(cfg, "targets");

    state->name = strdup(name);

    snprintf(f_name, 256, "%s.dat", name);
    state->dump_file = fopen(f_name, "w");

    while (1) {			/* find setting for target 'name' */
	this_target = config_setting_get_elem(targets, i);

	config_setting_lookup_string(this_target, "name", &S);
	if (strstr(S, name))
	    break;

	i++;
    }

    point = config_setting_get_member(this_target, "point");
    config_setting_lookup_float(point, "x", &F);
    state->point[0] = F;
    config_setting_lookup_float(point, "y", &F);
    state->point[1] = F;
    config_setting_lookup_float(point, "z", &F);
    state->point[2] = F;

    normal = config_setting_get_member(this_target, "normal");
    config_setting_lookup_float(normal, "x", &F);
    state->normal[0] = F;
    config_setting_lookup_float(normal, "y", &F);
    state->normal[1] = F;
    config_setting_lookup_float(normal, "z", &F);
    state->normal[2] = F;

    /* normalize normal vector */
    N = gsl_vector_view_array(state->normal, 3);
    norm = gsl_blas_dnrm2(&N.vector);
    gsl_vector_scale(&N.vector, 1.0 / norm);

    state->n_data = 0;
}

static void ps_free_state(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    size_t i;

    const size_t n = state->n_data;
    double *data = state->data;

    /* first write data to file */
    for (i = 0; i < n; i++) {
	const size_t idx = 3 * i;

	fprintf(state->dump_file, "%g\t%g\t%g\n", data[idx],
		data[idx + 1], data[idx + 2]);
    }

    fclose(state->dump_file);
    free(state->name);
    free(state->data);
}

static double *ps_get_intercept(void *vstate, ray_t * in_ray,
				int *dump_flag)
{
    ps_state_t *state = (ps_state_t *) vstate;

    double t1[3], t2;
    gsl_vector_view N, T1;

    if (*dump_flag) {
	if (state->n_data) {	/* we have not yet dumped our data */
	    double *t;
	    size_t i;
	    const size_t n = state->n_data;
	    double *data = state->data;

	    for (i = 0; i < n; i++) {
		const size_t idx = 3 * i;

		fprintf(state->dump_file, "%g\t%g\t%g\n", data[idx],
			data[idx + 1], data[idx + 2]);
	    }

	    /* shrink memory to minimum (BLOCK) */
	    t = (double *) realloc(state->data,
				   3 * BLOCK_SIZE * sizeof(double));
	    state->data = t;
	    state->n_data = 0;
	    state->n_alloc = BLOCK_SIZE;

	} else			/* we have allready dumped our data. mark cycle complete */
	    *dump_flag = 0;

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

    t1[0] = state->point[0] - in_ray->origin[0];	/* p_0 - l_0 */
    t1[1] = state->point[1] - in_ray->origin[1];
    t1[2] = state->point[2] - in_ray->origin[2];

    T1 = gsl_vector_view_array(t1, 3);
    N = gsl_vector_view_array(state->normal, 3);

    gsl_blas_ddot(&T1.vector, &N.vector, &t2);	/* (p_0 - l_0) dot N */

    if (t2) {			/* line does not start in plane */
	gsl_vector_view L = gsl_vector_view_array(in_ray->direction, 3);
	double t3;

	gsl_blas_ddot(&L.vector, &N.vector, &t3);	/* l dot n */

	if (t3) {		/* line not parallel to plane */
	    double d;
	    double *intercept = (double *) malloc(3 * sizeof(double));

	    d = t2 / t3;

	    intercept[0] = in_ray->origin[0] + d * in_ray->direction[0];
	    intercept[1] = in_ray->origin[1] + d * in_ray->direction[1];
	    intercept[2] = in_ray->origin[2] + d * in_ray->direction[2];

	    return intercept;
	}

    }

    return NULL;		/* no interception found */
}

static ray_t *ps_get_out_ray(void *vstate, ray_t * in_ray,
			     double *hit, int *dump_flag)
{
    ps_state_t *state = (ps_state_t *) vstate;

    ray_t *out = (ray_t *) malloc(sizeof(ray_t));

    memcpy(&(state->data[3 * state->n_data]), hit, 3 * sizeof(double));	/* store hit */

    if (state->n_data == state->n_alloc) {	/* inc data size for next hit */
	const size_t n = state->n_data + BLOCK_SIZE;
	double *t;

	t = (double *) realloc(state->data, 3 * n * sizeof(double));
	if (t) {		/* success, update state */
	    state->data = t;
	    state->n_data = n;
	} else {		/* memory exhausted, dump data to file and shrink memory to default size */
	    const size_t m = state->n_data;
	    double *data = state->data;

	    size_t i;

	    for (i = 0; i < m; i++) {
		const size_t idx = 3 * i;

		fprintf(state->dump_file, "%g\t%g\t%g\n", data[idx],
			data[idx + 1], data[idx + 2]);
	    }

	    /* shrink memory to minimum (BLOCK_SIZE) */
	    t = (double *) realloc(state->data,
				   3 * BLOCK_SIZE * sizeof(double));
	    state->data = t;
	    state->n_data = 0;
	    state->n_alloc = BLOCK_SIZE;

	    *dump_flag = 1;
	}
    } else
	state->n_data++;

    memcpy(out->origin, hit, 3 * sizeof(double));	/* update origin */
    memcpy(out->direction, in_ray->direction, 3 * sizeof(double));	/* copy direction */
    out->power = in_ray->power;

    return out;

}


static const target_type_t ps_t = {
    "plane screen",
    sizeof(struct ps_state_t),
    &ps_alloc_state,
    &ps_init_state,
    &ps_free_state,
    &ps_get_intercept,
    &ps_get_out_ray
};

const target_type_t *target_plane_screen = &ps_t;

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


typedef struct sq_state_t {
    char *name;			/* name (identifier) of target */
    pthread_key_t flags_key;	/* flags see target.h */
    int dump_file;
    double point[3];		/* center coordinate */
    double dx;			/* rectangle is '2*dx' times '2*dy' local coordinates */
    double dy;
    gsl_spline *spline;		/* for interpolated reflectivity spectrum */
    double normal[3];		/* normal vector of plane */
    double M[9];		/* transform matrix local -> global coordinates */
} sq_state_t;


static void sq_init_state(void *vstate, config_setting_t * this_target,
			  config_t * cfg, const int file_mode)
{
    sq_state_t *state = (sq_state_t *) vstate;

    int i;
    const char *S;
    char f_name[256];

    (void) cfg;			/* avoid warning: unused parameter 'cfg' */

    config_setting_lookup_string(this_target, "name", &S);
    state->name = strdup(S);

    snprintf(f_name, 256, "%s.dat", state->name);
    state->dump_file =
	open(f_name, O_CREAT | O_WRONLY | file_mode, S_IRUSR | S_IWUSR);

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
    init_refl_spectrum(S, &state->spline);

    pthread_key_create(&state->flags_key, free);
}

static void sq_free_state(void *vstate)
{
    sq_state_t *state = (sq_state_t *) vstate;

    close(state->dump_file);

    free(state->name);
    gsl_spline_free(state->spline);
}

static double *sq_get_intercept(void *vstate, ray_t * in_ray)
{
    sq_state_t *state = (sq_state_t *) vstate;

    double t1, t2[3], t3;
    double d;
    int *flag = pthread_getspecific(state->flags_key);

    if (*flag & LAST_WAS_HIT) {	/* ray starts on this target, no hit posible */
	*flag &= ~LAST_WAS_HIT;
	pthread_setspecific(state->flags_key, flag);
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

    v_diff(t2, state->point, in_ray->origin);	/* p_0 - l_0 */
    t3 = cblas_ddot(3, t2, 1, state->normal, 1);	/* (p_0 - l_0) dot N */
    if (fabs(t3) < GSL_SQRT_DBL_EPSILON)	/* line does start in target, conservative */
	return NULL;

    d = t3 / t1;
    if (d < 0.0)		/* intercepted target is not in front */
	return NULL;
    else {			/* 'in_ray' intercepts target plane */
	double l_intercept[3];
	double *intercept = (double *) malloc(3 * sizeof(double));

	v_a_plus_cb(intercept, in_ray->origin, d, in_ray->direction);
	/* convert to local coordinates, origin is 'state->point' */
	g2l(state->M, state->point, intercept, l_intercept);

	if ((l_intercept[0] <= -state->dx) || (l_intercept[0] >= state->dx)
	    || (l_intercept[1] <= -state->dy)
	    || (l_intercept[1] >= state->dy)) {

	    /* hit not within boundaries */
	    free(intercept);
	    return NULL;

	} else {		/* hits within target dimensions 'dx' times 'dy' */

	    if (t1 > 0.0) {	/* hits rear side, absorbed */
		*flag |= ABSORBED;
		pthread_setspecific(state->flags_key, flag);
	    }

	    return intercept;

	}
    }

}

static ray_t *sq_get_out_ray(void *vstate, ray_t * in_ray, double *hit,
			     const gsl_rng * r)
{
    sq_state_t *state = (sq_state_t *) vstate;
    int *flag = pthread_getspecific(state->flags_key);

    if (*flag & ABSORBED
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
	double hit_local[3];

	/* transform to local coordinates */
	memcpy(hit_local, hit, 3 * sizeof(double));
	g2l(state->M, state->point, hit, hit_local);

	/*
	 * store 4 items per data set (x,y,ppr,lambda)
	 * first x,y then ppr,lambda
	 */
	write(state->dump_file, hit_local, sizeof(double) * 2);
	write(state->dump_file, &in_ray->power, sizeof(double));
	write(state->dump_file, &in_ray->lambda, sizeof(double));

	*flag &= ~(LAST_WAS_HIT | ABSORBED);	/* clear flags */
	pthread_setspecific(state->flags_key, flag);

	free(in_ray);
	return NULL;

    } else {			/* reflect 'in_ray' */
	reflect(in_ray, state->normal, hit);

	*flag |= LAST_WAS_HIT;	/* mark as hit */
	pthread_setspecific(state->flags_key, flag);

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

    write(state->dump_file, str, strlen(str));
}

static double *sq_M(void *vstate)
{
    sq_state_t *state = (sq_state_t *) vstate;

    return state->M;
}

static void sq_init_flags(void *vstate)
{
    sq_state_t *state = (sq_state_t *) vstate;

    int *flag = (int *) malloc(sizeof(int));

    *flag = 0;
    pthread_setspecific(state->flags_key, flag);
}


static const target_type_t sq_t = {
    "rectangle",
    sizeof(struct sq_state_t),
    &sq_init_state,
    &sq_free_state,
    &sq_get_intercept,
    &sq_get_out_ray,
    &sq_get_target_name,
    &sq_dump_string,
    &sq_M,
    &sq_init_flags
};

const target_type_t *target_rectangle = &sq_t;

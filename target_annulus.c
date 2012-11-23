/*	target_annulus.c
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

#define NO_ITEMS 4


typedef struct ann_state_t {
    char *name;			/* name (identifier) of target */
    pthread_key_t flags_key;	/* flags see target.h */
    pthread_key_t outbuf_key;	/* access to output buffer for each target */
    int dump_file;
    double point[3];		/* center coordinate of annulus */
    double normal[3];		/* normal vector of annulus */
    double R;			/* inner radius of annulus */
    double r;			/* inner radius of annulus */
    gsl_spline *spline;		/* for interpolated reflectivity spectrum */
    int absorbed;		/* flag to indicated hit on rear surface == absorbed */
    double M[9];		/* transform matrix local -> global coordinates */
} ann_state_t;


static void ann_init_state(void *vstate, config_setting_t * this_target,
			   config_t * cfg, const int file_mode)
{
    ann_state_t *state = (ann_state_t *) vstate;

    const char *S;
    char f_name[256];

    (void) cfg;			/* avoid warning: unused parameter 'cfg' */

    config_setting_lookup_string(this_target, "name", &S);
    state->name = strdup(S);

    snprintf(f_name, 256, "%s.dat", state->name);
    state->dump_file =
	open(f_name, O_CREAT | O_WRONLY | file_mode, S_IRUSR | S_IWUSR);


    read_vector(this_target, "P", state->point);
    /*
     * generate transform matrix M to convert
     * between local and global coordinates
     * l2g:   g(x, y, z) = M l(x, y, z) + o(x, y, z))
     * g2l:   l(x, y, z) = MT (g(x, y, z) - o(x, y, z))
     */
    /* get normal vector of plane (serving also as basis vector z) */
    read_vector_normalize(this_target, "N", state->normal);
    memcpy(&state->M[6], state->normal, 3 * sizeof(double));

    /* get basis vector x */
    read_vector_normalize(this_target, "x", state->M);

    /* state->M[3-5] = y = z cross x */
    cross_product(&state->M[6], state->M, &state->M[3]);

    /* initialize reflectivity spectrum */
    config_setting_lookup_string(this_target, "reflectivity", &S);
    init_refl_spectrum(S, &state->spline);

    config_setting_lookup_float(this_target, "R", &state->R);
    config_setting_lookup_float(this_target, "r", &state->r);

    pthread_key_create(&state->flags_key, free);
    pthread_key_create(&state->outbuf_key, free_outbuf);
}

static void ann_free_state(void *vstate)
{
    ann_state_t *state = (ann_state_t *) vstate;

    close(state->dump_file);

    free(state->name);
    gsl_spline_free(state->spline);
}

static double *ann_get_intercept(void *vstate, ray_t * in_ray)
{
    ann_state_t *state = (ann_state_t *) vstate;

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
	double r_intercept;

	v_a_plus_cb(intercept, in_ray->origin, d, in_ray->direction);
	/* convert to local coordinates, origin is 'state->point' */
	g2l(state->M, state->point, intercept, l_intercept);

	/*
	 * r_intercep is distance from center of annulus to intercept
	 * in the plane of the annulus. we are in local system that
	 * is offset by intercept[2], so leave the latter out
	 */
	r_intercept =
	    sqrt(intercept[0] * intercept[0] +
		 intercept[1] * intercept[1]);
	if ((r_intercept > state->R)
	    || (r_intercept < state->r)) {

	    /* hit not within boundaries */
	    free(intercept);
	    return NULL;

	} else {		/* hits disk between radius 'r' and 'R' */

	    if (t1 > 0.0) {	/* hits rear side, absorbed */
		*flag |= ABSORBED;
		pthread_setspecific(state->flags_key, flag);
	    }

	    return intercept;

	}
    }

}

static ray_t *ann_get_out_ray(void *vstate, ray_t * in_ray, double *hit,
			      const gsl_rng * r)
{
    ann_state_t *state = (ann_state_t *) vstate;
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
	outbuf_t *out = pthread_getspecific(state->outbuf_key);

	/* transform to local coordinates */
	memcpy(hit_local, hit, 3 * sizeof(double));
	g2l(state->M, state->point, hit, hit_local);

	/*
	 * store 4 items per data set (x,y,ppr,lambda)
	 * first x,y then ppr,lambda
	 */
	if (out->i == BUF_SIZE * NO_ITEMS) {
	    write(state->dump_file, out->buf, sizeof(float) * out->i);
	    out->i = 0;
	}

	out->buf[out->i++] = (float) hit_local[0];
	out->buf[out->i++] = (float) hit_local[1];
	out->buf[out->i++] = (float) in_ray->power;
	out->buf[out->i++] = (float) in_ray->lambda;

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

static const char *ann_get_target_name(void *vstate)
{
    ann_state_t *state = (ann_state_t *) vstate;

    return state->name;
}

static void ann_dump_string(void *vstate, const char *str)
{
    ann_state_t *state = (ann_state_t *) vstate;

    write(state->dump_file, str, strlen(str));
}

static double *ann_M(void *vstate)
{
    ann_state_t *state = (ann_state_t *) vstate;

    return state->M;
}

static void ann_init_flags(void *vstate)
{
    ann_state_t *state = (ann_state_t *) vstate;

    int *flag = (int *) malloc(sizeof(int));

    *flag = 0;
    pthread_setspecific(state->flags_key, flag);
}

static void ann_init_outbuf(void *vstate)
{
    ann_state_t *state = (ann_state_t *) vstate;
    outbuf_t *out = (outbuf_t *) malloc(sizeof(outbuf_t));

    out->buf = (float *) malloc(BUF_SIZE * NO_ITEMS * sizeof(float));
    out->i = 0;

    pthread_setspecific(state->outbuf_key, out);
}

static void ann_flush_outbuf(void *vstate)
{
    ann_state_t *state = (ann_state_t *) vstate;
    outbuf_t *out = pthread_getspecific(state->outbuf_key);

    if (out->i != 0)		/* write rest of buffer to file. */
	write(state->dump_file, out->buf, sizeof(float) * out->i);
}

static const target_type_t ann_t = {
    "annulus",
    sizeof(struct ann_state_t),
    &ann_init_state,
    &ann_free_state,
    &ann_get_intercept,
    &ann_get_out_ray,
    &ann_get_target_name,
    &ann_dump_string,
    &ann_M,
    &ann_init_flags,
    &ann_init_outbuf,
    &ann_flush_outbuf
};

const target_type_t *target_annulus = &ann_t;

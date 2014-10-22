/*	target_ellipsoid.c
 *
 * Copyright (C) 2011,2012,2013,2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <math.h>
#include <string.h>

#include "math_utils.h"
#include "io_utils.h"
#include "targets.h"

#define TARGET_TYPE "ellipsoid"
#define NO_ITEMS 5


typedef struct ell_state_t {
    char reflectivity_model;	/* reflectivity model used for this target */
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
    int dump_file;
    double center[3];		/* center coordinate, origin of local system */
    double axes[3];		/* a^2, b^2, c^2 parameters (semi axes) */
    double z_min, z_max;	/* range of valid values of 'z' in local system */
    gsl_spline *refl_spectrum;	/* for interpolated reflectivity spectrum */
    double M[9];		/* transform matrix local -> global coordinates */
    void *refl_model_params;
} ell_state_t;


static void ell_surf_normal(const double *point, const double *axes,
			    double *const normal)
{
/*
 * from
 * http://www.freemathhelp.com/forum/threads/80712-Question-about-Normal-to-an-Ellipsoid
 *
 * this normal is for ellisoid cetered at (0,0,0) i.e.
 *
 *      x^2/a^2 + y^2/b^2 + z^2/c^2 - 1 = 0
 *
 * and points inwards! (note -2.0 * .....)
 */
    int i;
    double norm;

    for (i = 0, norm = 0.0; i < 3; i++) {
	normal[i] = -2.0 * point[i] / axes[i];
	norm += normal[i] * normal[i];
    }
    norm = sqrt(norm);
    cblas_dscal(3, 1.0 / norm, normal, 1);

}

static void ell_init_state(void *vstate, config_setting_t * this_target,
			   const int file_mode)
{
    ell_state_t *state = (ell_state_t *) vstate;

    const char *S;
    int i;

    read_vector(this_target, "center", state->center);
    /*
     * generate transform matrix M to convert
     * between local and global coordinates
     * l2g:   g(x, y, z) = MT l(x, y, z) + o(x, y, z)
     * g2l:   l(x, y, z) = M (g(x, y, z) - o(x, y, z))
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
    if (state->z_max < state->z_min)	/* safety */
	SWAP(state->z_max, state->z_min);

    /* initialize reflectivity spectrum */
    config_setting_lookup_string(this_target, "reflectivity", &S);
    init_refl_spectrum(S, &state->refl_spectrum);
    init_refl_model(this_target, &state->reflectivity_model,
		    &state->refl_model_params);

    state->dump_file =
	init_output(file_mode, TARGET_TYPE, this_target, state->center,
		    state->M);

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);
}

static void ell_free_state(void *vstate)
{
    ell_state_t *state = (ell_state_t *) vstate;

    state_free(state->dump_file, state->refl_spectrum,
	       state->reflectivity_model, state->refl_model_params);
}

static double *ell_get_intercept(void *vstate, ray_t * ray)
{
    ell_state_t *state = (ell_state_t *) vstate;

    int i;
    double r_O[3], r_N[3];	/* origin, direction of ray in local system */
    double O[] = { 0.0, 0.0, 0.0 };
    double A = 0.0, B = 0.0, C = -1.0;
    double D;

    /*
     * calculate point of interception D
     */
    /*
     * transform 'ray' from global to local system
     * origin 'ray': rotate / translate by origin of local system
     * dir 'ray': rotate only
     */
    g2l(state->M, state->center, ray->orig, r_O);
    g2l(state->M, O, ray->dir, r_N);

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
	    PTDT_t *data = pthread_getspecific(state->PTDT_key);
	    double dot;
	    double l_intercept[3];
	    double normal[3];
	    double *intercept = (double *) malloc(3 * sizeof(double));

	    l_intercept[0] = r_O[0] + h * r_N[0];
	    l_intercept[1] = r_O[1] + h * r_N[1];
	    l_intercept[2] = z;	/* use precomputed value */

	    ell_surf_normal(l_intercept, state->axes, normal);
	    dot = cblas_ddot(3, normal, 1, r_N, 1);

	    if (dot > 0.0)	/* parallel. hits outside, absorbed */
		data->flag |= ABSORBED;

	    /* convert to global coordinates, origin is 'state->center' */
	    l2g(state->M, state->center, l_intercept, intercept);

	    return intercept;
	}
    }
}

static ray_t *ell_get_out_ray(void *vstate, ray_t * ray, double *hit,
			      const gsl_rng * r)
{
    ell_state_t *state = (ell_state_t *) vstate;
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

	if (state->dump_file != -1)
	    store_xyz(state->dump_file, ray, hit, state->M, state->center,
		      data, &state->mutex_writefd);

	data->flag &= ~ABSORBED;	/* clear flag */

	free(ray);
	return NULL;

    } else {			/* reflect 'ray' */
	double l_N[3], N[3];
	double O[] = { 0.0, 0.0, 0.0 };
	double hit_local[3];

	g2l(state->M, state->center, hit, hit_local);	/* transform to local coordinates */
	ell_surf_normal(hit_local, state->axes, l_N);	/* normal vector local system */
	l2g(state->M, O, l_N, N);	/* normal vector global system */

	reflect(ray, N, hit, state->reflectivity_model, r,
		state->refl_model_params);

	return ray;
    }
}

static void ell_init_PTDT(void *vstate)
{
    per_thread_init(((ell_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void ell_flush_PTDT_outbuf(void *vstate)
{
    ell_state_t *state = (ell_state_t *) vstate;

    per_thread_flush(state->dump_file, state->PTDT_key,
		     &state->mutex_writefd);
}


static const target_type_t ell_t = {
    TARGET_TYPE,
    sizeof(struct ell_state_t),
    &ell_init_state,
    &ell_free_state,
    &ell_get_intercept,
    &ell_get_out_ray,
    &ell_init_PTDT,
    &ell_flush_PTDT_outbuf
};

const target_type_t *target_ellipsoid = &ell_t;

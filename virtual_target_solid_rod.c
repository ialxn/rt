/*	virtual_target_solid_rod.c
 *
 * Copyright (C) 2015 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#define _GNU_SOURCE		/* for sincos() */

#include<string.h>

#include "io_utils.h"
#include "intercept.h"
#include "off.h"
#include "likely.h"
#include "reflect.h"
#include "targets.h"
#include "virtual_targets.h"

#define NO_ITEMS 0


typedef struct vtsrod_state_t {
    double origin[3];		/* origin of cylinder */
    double radius;		/* radius of cylinder */
    double length;		/* length of cylinder */
    double barrier1;		/* barriers to select emitting surface */
    double barrier2;
    double M[9];		/* convert local <-> global coordinates */
    gsl_spline *spectrum;	/* interpolated emission spectrum */
    gsl_spline *reflectivity;	/* interpolated reflectivity spectrum */
    refl_func_pointer_t refl_func;	/* reflection model */
    void *refl_func_pars;	/* model specific parameters */
    pthread_key_t PTDT_key;	/* access to per thread flags */
} vtsrod_state_t;


static double get_A_disk(vtsrod_state_t * state, config_setting_t * this_s,
			 const char *kw)
{
    int ans;

    config_setting_lookup_bool(this_s, kw, &ans);
    if (ans)			/* surface emits */
	return M_PI * state->radius * state->radius;
    else			/* surface does not emit */
	return 0.0;
}

static void init_barriers(config_setting_t * this_s,
			  vtsrod_state_t * state)
{
    double sum = 0.0;
    double A_wall, A_base, A_top;

    A_wall = 2.0 * M_PI * state->radius * state->length;
    sum += A_wall;
    A_base = get_A_disk(state, this_s, "base_face_emits");
    sum += A_base;
    A_top = get_A_disk(state, this_s, "base_face_emits");
    sum += A_top;

    /*
     * a random number will be used to select face that emits ray:
     * 
     *  0 <= r < A_wall -> from wall
     *  A_wall <= r < (A_wall + A_base) -> from base disk
     *  r >= (A_wall + A_base)  -> from top disk
     */
    state->barrier1 = A_wall / sum;
    state->barrier2 = (A_wall + A_base) / sum;
}

static double *intercept_face(const ray_t * in_ray,
			      const vtsrod_state_t * state,
			      const double *center)
{
    double *intercept;
    int dummy_i;

    if ((intercept =
	 intercept_plane(in_ray, &state->M[6], center,
			 &dummy_i)) != NULL) {
	double l_intercept[3];

	g2l(state->M, center, intercept, l_intercept);

	if ((l_intercept[0] * l_intercept[0] +
	     l_intercept[1] * l_intercept[1]) >
	    (state->radius * state->radius)) {
	    /* hit not within boundaries */
	    free(intercept);
	    intercept = NULL;	/* mark as not valid */
	}
    }

    return intercept;
}


static int vtsrod_init_state(void *vstate, config_setting_t * this_target,
			     const int file_mode, const int keep_closed,
			     const double P_factor)
{
    vtsrod_state_t *state = (vtsrod_state_t *) vstate;

    read_vector(this_target, "origin", state->origin);
    config_setting_lookup_float(this_target, "radius", &state->radius);
    config_setting_lookup_float(this_target, "length", &state->length);

    init_barriers(this_target, state);

    read_vector(this_target, "direction", &state->M[6]);
    state->M[0] = 1.0;
    state->M[1] = 0.0;
    state->M[2] = 0.0;
    orthonormalize(&state->M[0], &state->M[3], &state->M[6]);

    if (init_spectrum(this_target, "spectrum", &state->spectrum))
	return ERR;

    if (init_spectrum(this_target, "reflectivity", &state->reflectivity))
	return ERR;

    init_refl_model(this_target, &state->refl_func,
		    &state->refl_func_pars);

    pthread_key_create(&state->PTDT_key, free_PTDT);

    return NO_ERR;
}

static void vtsrod_free_state(void *vstate)
{
    vtsrod_state_t *state = (vtsrod_state_t *) vstate;

    gsl_spline_free(state->spectrum);
    gsl_spline_free(state->reflectivity);

    if (state->refl_func == reflect_microfacet_gaussian)
	free((double *) state->refl_func_pars);
}

static double *vtsrod_get_intercept(void *vstate, ray_t * ray)
{
    vtsrod_state_t *state = (vtsrod_state_t *) vstate;
    double *intercept;
    int side;

    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & LAST_WAS_HIT) {
	data->flag &= ~LAST_WAS_HIT;
	return NULL;
    }

    intercept = intercept_cylinder(ray, state->origin, &state->M[6],
				   state->radius, state->length, &side);

    if (intercept && (side == OUTSIDE))
	return intercept;

    if (intercept)
	free(intercept);

    if (my_ddot(&state->M[6], ray->dir) > 0)
	intercept = intercept_face(ray, state, state->origin);
    else {
	double center_face2[3];

	a_plus_cb(center_face2, state->origin, state->length,
		  &state->M[6]);
	intercept = intercept_face(ray, state, center_face2);
    }

    return intercept;		/* might be NULL */
}

static ray_t *vtsrod_get_out_ray(void *vstate, ray_t * ray, double *hit,
				 const gsl_rng * r)
{
    vtsrod_state_t *state = (vtsrod_state_t *) vstate;

    PTDT_t *data = pthread_getspecific(state->PTDT_key);


    if (gsl_rng_uniform(r) >
	gsl_spline_eval(state->reflectivity, ray->lambda, NULL)) {
	/*
	 * ray is absorbed. we emit the same power from a random point on the
	 * source but with the emissionspectrum of this source.
	 */
	double point[3];
	double normal[3];
	double phi, sin_phi, cos_phi;
	double selector = gsl_rng_uniform(r);

	phi = 2.0 * M_PI * gsl_rng_uniform(r);
	sincos(phi, &sin_phi, &cos_phi);

	if (selector < state->barrier1) {
	    /* ray originates from wall */
	    point[0] = state->radius * cos_phi;
	    point[1] = state->radius * sin_phi;
	    point[2] = state->length * gsl_rng_uniform(r);

	    normal[0] = point[0];
	    normal[1] = point[1];
	    normal[2] = 0.0;
	    normalize(normal);

	} else if (state->barrier2 < selector) {
	    /* ray originates from top disk */
	    double t = sqrt(gsl_rng_uniform(r)) * state->radius;

	    point[0] = t * sin_phi;
	    point[1] = t * cos_phi;
	    point[2] = state->length;

	    normal[0] = 0.0;
	    normal[1] = 0.0;
	    normal[2] = 1.0;

	} else {
	    /* ray originates from base disk */
	    double t = sqrt(gsl_rng_uniform(r)) * state->radius;

	    point[0] = t * sin_phi;
	    point[1] = t * cos_phi;
	    point[2] = 0.0;

	    normal[0] = 0.0;
	    normal[1] = 0.0;
	    normal[2] = -1.0;

	}

	l2g(state->M, state->origin, point, ray->orig);

	/*
	 * choose random direction
	 */
	get_uniform_random_vector_hemisphere(point, 1.0, normal, r);
	l2g_rot(state->M, point, ray->dir);

	ray->lambda =
	    gsl_spline_eval(state->spectrum, gsl_rng_uniform(r), NULL);

    } else {			/* reflect 'ray' */
	double normal[3];
	double tmp[3];

	/*
	 * to calculate the normal vector at 'hit' we need to know which part
	 * of the cylinder was hit.
	 *
	 * base face:
	 *      + hit - state->origin is perpendicular to state->direction
	 *      + normal is -state->direction
	 * top face:
	 *      + hit - (state->origin + state->length*state->direction) is
	 *        perpendicular to state->direction
	 *      + normal is state->direction
	 * cylinder wall:
	 *      + neither face hit
	 *      + normal is
	 */

	diff(tmp, hit, state->origin);
	/*
	 * normalize tmp to avoid problems with intercepts close to origin
	 */
	normalize(tmp);

	if (fabs(my_ddot(tmp, &state->M[6])) < GSL_SQRT_DBL_EPSILON) {	/* base */

	    normal[0] = -state->M[6];
	    normal[1] = -state->M[7];
	    normal[2] = -state->M[8];

	} else {
	    diff(tmp, hit, state->origin);
	    my_daxpy(-state->length, &state->M[6], tmp);
	    normalize(tmp);

	    if (fabs(my_ddot(tmp, &state->M[6])) < GSL_SQRT_DBL_EPSILON)	/* top */
		memcpy(normal, &state->M[6], 3 * sizeof(double));
	    else		/* wall */
		cyl_surf_normal(hit, state->origin, &state->M[6],
				state->radius, normal);
	}
	state->refl_func(ray, normal, hit, r, state->refl_func_pars);
    }

    if (unlikely(ray->n_refl == UCHAR_MAX))	/* too many reflections is unlikely */
	fprintf(stderr,
		"                INFO: maximum path length exceeded (wrap around of counter occurs)\n");
    ++ray->n_refl;

    data->flag |= LAST_WAS_HIT;	/* mark as hit */
    return ray;
}

static void vtsrod_init_PTDT(void *vstate)
{
    per_thread_init(((vtsrod_state_t *) vstate)->PTDT_key, 0);
}

static const target_type_t vt_srod_t = {
    NULL,
    sizeof(struct vtsrod_state_t),
    &vtsrod_init_state,
    &vtsrod_free_state,
    &vtsrod_get_intercept,
    &vtsrod_get_out_ray,
    &vtsrod_init_PTDT,
    NULL
};

const target_type_t *virtual_target_solid_rod = &vt_srod_t;

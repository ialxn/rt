/*	source_solid_cylinder.c
 *
 * Copyright (C) 2015,2016 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#define _ISOC99_SOURCE		/* because of llrint() */
#define _GNU_SOURCE		/* for sincos() */

#include <cblas.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include "io_utils.h"
#include "intercept.h"
#include "likely.h"
#include "math_utils.h"
#include "off.h"
#include "ray.h"
#include "reflect.h"
#include "targets.h"
#include "sources.h"


typedef struct scylinder_state_t {
    char *name;			/* name (identifier) of uniform point source */
    double orig[3];		/* origin (center of base face) */
    double radius;		/* radius of cylinder */
    double length;		/* length of cylinder */
    double barrier1;		/* barriers to select emiting surface */
    double barrier2;
    double M[9];		/* convert local <-> global coordinates */
    int64_t n_rays;		/* number of rays remaining until source is exhausted */
    double power;		/* power of source */
    gsl_spline *spectrum;	/* spline holding cdf of source spectrum */
    pthread_mutex_t mutex_n_rays;	/* protect n_rays */
    pthread_key_t rays_remain_key;	/* no of ray remain in group (PTD) */
} scylinder_state_t;

typedef struct vtscylinder_state_t {
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
} vtscylinder_state_t;


static void init_barriers(config_setting_t * this_s, double *barrier1,
			  double *barrier2, const double radius,
			  const double length)
{
    double sum = 0.0;
    double A_wall, A_base, A_top;

    A_wall = 2.0 * M_PI * radius * length;
    sum += A_wall;
    A_base = get_A_disk(this_s, "base_face_emits", radius);
    sum += A_base;
    A_top = get_A_disk(this_s, "top_face_emits", radius);
    sum += A_top;

    /*
     * a random number will be used to select face that emits ray:
     * 
     *  0 <= r < A_wall -> from wall
     *  A_wall <= r < (A_wall + A_base) -> from base disk
     *  r >= (A_wall + A_base)  -> from top disk
     */
    *barrier1 = A_wall / sum;
    *barrier2 = (A_wall + A_base) / sum;
}

static void init_ray(ray_t * ray, const double radius,
		     const double length, const double barrier1,
		     const double barrier2, const double origin[3],
		     const double *M, gsl_spline * spectrum,
		     const gsl_rng * r)
{
    double point[3];
    double normal[3];
    double phi, sin_phi, cos_phi;
    double selector = gsl_rng_uniform(r);

    phi = 2.0 * M_PI * gsl_rng_uniform(r);
    sincos(phi, &sin_phi, &cos_phi);

    if (selector < barrier1) {
	/* ray originates from wall */
	point[0] = radius * cos_phi;
	point[1] = radius * sin_phi;
	point[2] = length * gsl_rng_uniform(r);

	normal[0] = point[0];
	normal[1] = point[1];
	normal[2] = 0.0;
	normalize(normal);
    } else if (barrier2 < selector) {
	/* ray originates from top disk */
	double t = sqrt(gsl_rng_uniform(r)) * radius;

	point[0] = t * sin_phi;
	point[1] = t * cos_phi;
	point[2] = length;

	normal[0] = 0.0;
	normal[1] = 0.0;
	normal[2] = 1.0;
    } else {
	/* ray originates from base disk */
	double t = sqrt(gsl_rng_uniform(r)) * radius;

	point[0] = t * sin_phi;
	point[1] = t * cos_phi;
	point[2] = 0.0;

	normal[0] = 0.0;
	normal[1] = 0.0;
	normal[2] = -1.0;
    }

    l2g(M, origin, point, ray->orig);

    /*
     * choose random direction
     */
    get_uniform_random_vector_hemisphere(point, 1.0, normal, r);
    l2g_rot(M, point, ray->dir);

    ray->lambda = gsl_spline_eval(spectrum, gsl_rng_uniform(r), NULL);
}

static ray_t *new_ray(const double radius, const double length,
		      const double barrier1, const double barrier2,
		      const double origin[3], const double *M,
		      gsl_spline * spectrum, const gsl_rng * r)
{
    ray_t *ray = (ray_t *) malloc(sizeof(ray_t));

    init_ray(ray, radius, length, barrier1, barrier2, origin, M,
	     spectrum, r);
    ray->n_refl = 0;
    return ray;
}


static void scylinder_init_state(void *vstate, config_setting_t * this_s,
				 const double P_factor)
{
    scylinder_state_t *state = (scylinder_state_t *) vstate;
    const char *S;

    config_setting_lookup_string(this_s, "name", &S);
    state->name = strdup(S);

    config_setting_lookup_float(this_s, "power", &state->power);
    state->n_rays = llrint(state->power / P_factor);
    pthread_mutex_init(&state->mutex_n_rays, NULL);

    read_vector(this_s, "origin", state->orig);
    config_setting_lookup_float(this_s, "radius", &state->radius);
    config_setting_lookup_float(this_s, "length", &state->length);

    init_barriers(this_s, &state->barrier1, &state->barrier2,
		  state->radius, state->length);

    init_M_from_z(this_s, "direction", state->M);

    /* initialize source spectrum */
    init_source_spectrum(this_s, "spectrum", &state->spectrum);

    pthread_key_create(&state->rays_remain_key, free);
}

static void scylinder_free_state(void *vstate)
{
    scylinder_state_t *state = (scylinder_state_t *) vstate;

    free(state->name);
    gsl_spline_free(state->spectrum);
}

static ray_t *scylinder_emit_ray(void *vstate, const gsl_rng * r)
{
    scylinder_state_t *state = (scylinder_state_t *) vstate;
    ray_t *ray = NULL;
    int64_t *rays_remain = pthread_getspecific(state->rays_remain_key);

    if (!*rays_remain)
	*rays_remain =
	    per_thread_get_new_raygroup(&state->mutex_n_rays,
					&state->n_rays);

    if (*rays_remain > 0) {	/* rays still available in group */

	(*rays_remain)--;
	pthread_setspecific(state->rays_remain_key, rays_remain);

	ray = new_ray(state->radius, state->length, state->barrier1,
		      state->barrier2, state->orig, state->M,
		      state->spectrum, r);
    }

    return ray;
}

static const char *scylinder_get_source_name(void *vstate)
{
    return ((scylinder_state_t *) vstate)->name;
}

static int64_t scylinder_get_source_n_rays(void *vstate)
{
    scylinder_state_t *state = (scylinder_state_t *) vstate;

    return per_thread_get_source_n_rays(&state->mutex_n_rays,
					&state->n_rays);
}

static double scylinder_get_source_power(void *vstate)
{
    return ((scylinder_state_t *) vstate)->power;
}

static void scylinder_init_rays_remain(void *vstate)
{
    per_thread_init_rays_remain(((scylinder_state_t *)
				 vstate)->rays_remain_key);
}


static int vtscylinder_init_state(void *vstate,
				  config_setting_t * this_target, const int
				  __attribute__ ((__unused__)) file_mode, const int
				  __attribute__ ((__unused__)) keep_closed, const double
				  __attribute__ ((__unused__)) P_factor)
{
    vtscylinder_state_t *state = (vtscylinder_state_t *) vstate;

    read_vector(this_target, "origin", state->origin);
    config_setting_lookup_float(this_target, "radius", &state->radius);
    config_setting_lookup_float(this_target, "length", &state->length);

    init_barriers(this_target, &state->barrier1, &state->barrier2,
		  state->radius, state->length);

    init_M_from_z(this_target, "direction", state->M);

    init_source_spectrum(this_target, "spectrum", &state->spectrum);

    init_spectrum(this_target, "reflectivity", &state->reflectivity);
    init_refl_model(this_target, &state->refl_func,
		    &state->refl_func_pars);

    pthread_key_create(&state->PTDT_key, free_PTDT);

    return NO_ERR;
}

static void vtscylinder_free_state(void *vstate)
{
    vtscylinder_state_t *state = (vtscylinder_state_t *) vstate;

    gsl_spline_free(state->spectrum);
    gsl_spline_free(state->reflectivity);

    if (state->refl_func == reflect_microfacet_gaussian)
	free((double *) state->refl_func_pars);
}

static double *vtscylinder_get_intercept(void *vstate, ray_t * ray)
{
    vtscylinder_state_t *state = (vtscylinder_state_t *) vstate;
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

    if (cblas_ddot(3, &state->M[6], 1, ray->dir, 1) > 0)
	intercept =
	    intercept_disk(ray, state->origin, state->M,
			   state->radius * state->radius, &side);
    else {
	double center_face2[3];

	a_plus_cb(center_face2, state->origin, state->length,
		  &state->M[6]);
	intercept =
	    intercept_disk(ray, center_face2, state->M,
			   state->radius * state->radius, &side);
    }

    return intercept;		/* might be NULL */
}

static ray_t *vtscylinder_get_out_ray(void *vstate, ray_t * ray,
				      double *hit, const gsl_rng * r)
{
    vtscylinder_state_t *state = (vtscylinder_state_t *) vstate;

    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (gsl_rng_uniform(r) >
	gsl_spline_eval(state->reflectivity, ray->lambda, NULL)) {
	/*
	 * ray is absorbed. we emit the same power from a random point on the
	 * source but with the emissionspectrum of this source.
	 */
	init_ray(ray, state->radius, state->length, state->barrier1,
		 state->barrier2, state->origin, state->M,
		 state->spectrum, r);
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

	if (fabs(cblas_ddot(3, tmp, 1, &state->M[6], 1)) < GSL_SQRT_DBL_EPSILON) {	/* base */

	    normal[0] = -state->M[6];
	    normal[1] = -state->M[7];
	    normal[2] = -state->M[8];

	} else {
	    diff(tmp, hit, state->origin);
	    my_daxpy(-state->length, &state->M[6], tmp);
	    normalize(tmp);

	    if (fabs(cblas_ddot(3, tmp, 1, &state->M[6], 1)) < GSL_SQRT_DBL_EPSILON)	/* top */
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

static void vtscylinder_init_PTDT(void *vstate)
{
    per_thread_init(((vtscylinder_state_t *) vstate)->PTDT_key, 0);
}


static const source_type_t scylinder_t = {
    "solid cylinder",
    sizeof(struct scylinder_state_t),
    &scylinder_init_state,
    &scylinder_free_state,
    &scylinder_emit_ray,
    &scylinder_get_source_name,
    &scylinder_get_source_n_rays,
    &scylinder_get_source_power,
    &scylinder_init_rays_remain
};

static const target_type_t vt_scylinder_t = {
    NULL,
    sizeof(struct vtscylinder_state_t),
    &vtscylinder_init_state,
    &vtscylinder_free_state,
    &vtscylinder_get_intercept,
    &vtscylinder_get_out_ray,
    &vtscylinder_init_PTDT,
    NULL
};


const source_type_t *source_solid_cylinder = &scylinder_t;
const target_type_t *virtual_target_solid_cylinder = &vt_scylinder_t;

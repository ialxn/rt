/*	target_cone.c
 *
 * Copyright (C) 2014,2015,2016 Ivo Alxneit
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

#include "io_utils.h"
#include "intercept.h"
#include "likely.h"
#include "sources.h"
#include "targets.h"

#define TARGET_TYPE "cone"
#define NO_ITEMS 4


typedef struct scone_state_t {
    char *name;			/* name (identifier) of uniform point source */
    double orig[3];		/* center of base disk, origin of local system */
    double H;			/* (full) height of cone */
    double tan2_a;		/* a=opening angle of cone */
    double h;			/* height at radius 'r' */
    double R;			/* radius at base */
    double r;			/* radius at top */
    double barrier1;		/* barriers to select emiting surface */
    double barrier2;
    double M[9];		/* transform matrix local -> global coordinates */
    int64_t n_rays;		/* number of rays remaining until source is exhausted */
    double power;		/* power of source */
    gsl_spline *spectrum;	/* spline holding cdf of source spectrum */
    pthread_mutex_t mutex_n_rays;	/* protect n_rays */
    pthread_key_t rays_remain_key;	/* no of ray remain in group (PTD) */
} scone_state_t;

typedef struct vtscone_state_t {
    double origin[3];		/* center of base disk, origin of local system */
    double H;			/* (full) height of cone */
    double tan2_a;		/* a=opening angle of cone */
    double h;			/* height at radius 'r' */
    double R;			/* radius at base */
    double r;			/* radius at top */
    double barrier1;		/* barriers to select emiting surface */
    double barrier2;
    double M[9];		/* transform matrix local -> global coordinates */
    gsl_spline *spectrum;	/* interpolated emission spectrum */
    gsl_spline *reflectivity;	/* interpolated reflectivity spectrum */
    refl_func_pointer_t refl_func;	/* reflection model */
    void *refl_func_pars;	/* model specific parameters */
    pthread_key_t PTDT_key;	/* access to per thread flags */
} vtscone_state_t;


static void init_barriers(config_setting_t * this_s, double *barrier1,
			  double *barrier2, const double R, const double r,
			  const double H, const double h)
{
    double sum = 0.0;
    double S, s;
    double A_wall, A_base, A_top;

    S = sqrt(R * R + H * H);
    s = sqrt(r * r + (H - h) * (H - h));

    A_wall = M_PI * (R * S - r * s);
    sum += A_wall;
    A_base = get_A_disk(this_s, "base_face_emits", R);
    sum += A_base;
    A_top = get_A_disk(this_s, "top_face_emits", r);
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

static void init_ray(ray_t * ray, const double Radius,
		     const double radius, const double H, const double h,
		     const double barrier1, const double barrier2,
		     const double origin[3], const double tan2_a,
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
	double z, h_z;
	const double cutoff = radius / Radius;

	do {
	    z = sqrt(gsl_rng_uniform(r));
	} while (z < cutoff);

	h_z = z * H;

	point[0] = z * cos_phi;
	point[1] = z * sin_phi;
	point[2] = h_z;

	cone_surf_normal(point, tan2_a, H, normal);
    } else if (barrier2 < selector) {
	/* ray originates from top disk */
	double t = sqrt(gsl_rng_uniform(r)) * radius;

	point[0] = t * sin_phi;
	point[1] = t * cos_phi;
	point[2] = h;

	normal[0] = 0.0;
	normal[1] = 0.0;
	normal[2] = 1.0;
    } else {
	/* ray originates from base disk */
	double t = sqrt(gsl_rng_uniform(r)) * Radius;

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

static ray_t *new_ray(const double Radius, const double radius,
		      const double H, const double h,
		      const double barrier1, const double barrier2,
		      const double origin[3], const double tan2_a,
		      const double *M, gsl_spline * spectrum,
		      const gsl_rng * r)
{
    ray_t *ray = (ray_t *) malloc(sizeof(ray_t));

    init_ray(ray, Radius, radius, H, h, barrier1, barrier2, origin,
	     tan2_a, M, spectrum, r);
    ray->n_refl = 0;
    return ray;
}



static void scone_init_state(void *vstate, config_setting_t * this_s,
			     const double P_factor)
{
    scone_state_t *state = (scone_state_t *) vstate;
    const char *S;

    config_setting_lookup_string(this_s, "name", &S);
    state->name = strdup(S);

    read_vector(this_s, "origin", state->orig);

    init_M_from_z(this_s, "z", state->M);

    config_setting_lookup_float(this_s, "R", &state->R);
    config_setting_lookup_float(this_s, "r", &state->r);
    config_setting_lookup_float(this_s, "h", &state->h);

    state->H = state->R * state->h / (state->R - state->r);
    state->tan2_a = (state->R * state->R) / (state->H * state->H);

    init_barriers(this_s, &state->barrier1, &state->barrier2,
		  state->R, state->r, state->H, state->h);

    config_setting_lookup_float(this_s, "power", &state->power);
    state->n_rays = llrint(state->power / P_factor);

    /* initialize source spectrum */
    init_source_spectrum(this_s, "spectrum", &state->spectrum);

    pthread_mutex_init(&state->mutex_n_rays, NULL);
    pthread_key_create(&state->rays_remain_key, free);
}

static void scone_free_state(void *vstate)
{
    scone_state_t *state = (scone_state_t *) vstate;

    free(state->name);
    gsl_spline_free(state->spectrum);
}

static ray_t *scone_emit_ray(void *vstate, const gsl_rng * r)
{
    scone_state_t *state = (scone_state_t *) vstate;
    ray_t *ray = NULL;
    int64_t *rays_remain = pthread_getspecific(state->rays_remain_key);

    if (!*rays_remain)
	*rays_remain =
	    per_thread_get_new_raygroup(&state->mutex_n_rays,
					&state->n_rays);
    if (*rays_remain > 0) {	/* rays still available in group */

	(*rays_remain)--;
	pthread_setspecific(state->rays_remain_key, rays_remain);

	ray = new_ray(state->R, state->r, state->H, state->h,
		      state->barrier1, state->barrier2, state->orig,
		      state->tan2_a, state->M, state->spectrum, r);
    }

    return ray;
}

static const char *scone_get_source_name(void *vstate)
{
    return ((scone_state_t *) vstate)->name;
}

static int64_t scone_get_source_n_rays(void *vstate)
{
    scone_state_t *state = (scone_state_t *) vstate;

    return per_thread_get_source_n_rays(&state->mutex_n_rays,
					&state->n_rays);
}

static double scone_get_source_power(void *vstate)
{
    return ((scone_state_t *) vstate)->power;
}

static void scone_init_rays_remain(void *vstate)
{
    per_thread_init_rays_remain(((scone_state_t *)
				 vstate)->rays_remain_key);
}


static int vtscone_init_state(void *vstate, config_setting_t * this_target, const int
			      __attribute__ ((__unused__)) file_mode, const int
			      __attribute__ ((__unused__)) keep_closed, const double
			      __attribute__ ((__unused__)) P_factor)
{
    vtscone_state_t *state = (vtscone_state_t *) vstate;

    read_vector(this_target, "origin", state->origin);
    init_M_from_z(this_target, "z", state->M);

    config_setting_lookup_float(this_target, "R", &state->R);
    config_setting_lookup_float(this_target, "r", &state->r);
    config_setting_lookup_float(this_target, "h", &state->h);

    state->H = state->R * state->h / (state->R - state->r);
    state->tan2_a = (state->R * state->R) / (state->H * state->H);

    init_barriers(this_target, &state->barrier1, &state->barrier2,
		  state->R, state->r, state->H, state->h);

    init_source_spectrum(this_target, "spectrum", &state->spectrum);

    init_spectrum(this_target, "reflectivity", &state->reflectivity);
    init_refl_model(this_target, &state->refl_func,
		    &state->refl_func_pars);

    pthread_key_create(&state->PTDT_key, free_PTDT);

    return NO_ERR;
}

static void vtscone_free_state(void *vstate)
{
    vtscone_state_t *state = (vtscone_state_t *) vstate;

    gsl_spline_free(state->spectrum);
    gsl_spline_free(state->reflectivity);

    if (state->refl_func == reflect_microfacet_gaussian)
	free((double *) state->refl_func_pars);
}

static double *vtscone_get_intercept(void *vstate, ray_t * ray)
{
    vtscone_state_t *state = (vtscone_state_t *) vstate;
    double *intercept;
    int side;

    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & LAST_WAS_HIT) {
	data->flag &= ~LAST_WAS_HIT;
	return NULL;
    }

    intercept = intercept_cone(ray, state->M, state->origin,
			       state->tan2_a, state->H, state->h, &side);

    if (intercept && (side == OUTSIDE))
	return intercept;

    if (intercept)		/* hit inside of cone, discard */
	free(intercept);

    if (cblas_ddot(3, &state->M[6], 1, ray->dir, 1) > 0)	/* hit base face possible */
	intercept =
	    intercept_disk(ray, state->origin, state->M,
			   state->R * state->R, &side);
    else {
	double center_face2[3];

	a_plus_cb(center_face2, state->origin, state->h, &state->M[6]);
	intercept =
	    intercept_disk(ray, center_face2, state->M,
			   state->r * state->r, &side);
    }

    return intercept;		/* might be NULL */
}

static ray_t *vtscone_get_out_ray(void *vstate, ray_t * ray, double *hit,
				  const gsl_rng * r)
{
    vtscone_state_t *state = (vtscone_state_t *) vstate;

    PTDT_t *data = pthread_getspecific(state->PTDT_key);


    if (gsl_rng_uniform(r) >
	gsl_spline_eval(state->reflectivity, ray->lambda, NULL)) {
	/*
	 * ray is absorbed. we emit the same power from a random point on the
	 * source but with the emissionspectrum of this source.
	 */
	init_ray(ray, state->R, state->r, state->H, state->h,
		 state->barrier1, state->barrier2, state->origin,
		 state->tan2_a, state->M, state->spectrum, r);
    } else {			/* reflect 'ray' */
	double normal[3];
	double tmp[3];

	/*
	 * to calculate the normal vector at 'hit' we need to know which part
	 * of the cylinder was hit.
	 *
	 * base face:
	 *      + hit - state->origin is perpendicular to axis of cone
	 *      + normal is -axis_of_cone
	 * top face:
	 *      + hit - (state->origin + state->h*axis_of_cone) is
	 *        perpendicular to axis_of_cone
	 *      + normal is axis_of_cone
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
	    my_daxpy(-state->h, &state->M[6], tmp);
	    normalize(tmp);

	    if (fabs(cblas_ddot(3, tmp, 1, &state->M[6], 1)) < GSL_SQRT_DBL_EPSILON)	/* top */
		memcpy(normal, &state->M[6], 3 * sizeof(double));
	    else		/* wall */
		cone_surf_normal(hit, state->tan2_a, state->H, normal);
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

static void vtscone_init_PTDT(void *vstate)
{
    per_thread_init(((vtscone_state_t *) vstate)->PTDT_key, 0);
}


static const source_type_t scone_t = {
    "solid cone",
    sizeof(struct scone_state_t),
    &scone_init_state,
    &scone_free_state,
    &scone_emit_ray,
    &scone_get_source_name,
    &scone_get_source_n_rays,
    &scone_get_source_power,
    &scone_init_rays_remain
};

static const target_type_t vt_scone_t = {
    NULL,
    sizeof(struct vtscone_state_t),
    &vtscone_init_state,
    &vtscone_free_state,
    &vtscone_get_intercept,
    &vtscone_get_out_ray,
    &vtscone_init_PTDT,
    NULL
};


const source_type_t *source_solid_cone = &scone_t;
const target_type_t *virtual_target_solid_cone = &vt_scone_t;

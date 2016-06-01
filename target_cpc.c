/*	target_cpc.c
 *
 * Copyright (C) 2016 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#define _GNU_SOURCE		/* for sincos() */

#include <math.h>
#include <string.h>
#include <gsl/gsl_roots.h>

#include "io_utils.h"
#include "intercept.h"
#include "reflect.h"
#include "targets.h"

#define TARGET_TYPE "cpc"
#define NO_ITEMS 4

#define ENTRANCE_APERTURE 1
#define EXIT_APERTURE 2
#define MAX_ITER 100
#define ABS_EPS_ROOT 1.0e-6
#define REL_EPS_ROOT 1.0e-2

typedef struct cpc_state_t {
    double origin[3];		/* center of exit aperture */
    double phi_a;		/* acceptance angle */
    double theta_t;		/* truncation angle (>= phi_a) */
    double Rt;			/* radius of entrance aperture (at theta_t) */
    double r;			/* radius of exit aperture */
    double origin2[3];		/* center of entrance aperture */
    double f2;			/* 2* focal length of CPC */
    double *M;
    gsl_spline *refl_spectrum;	/* for interpolated reflectivity spectrum */
    refl_func_pointer_t refl_func;	/* reflection model */
    void *refl_func_pars;	/* model specific parameters */
    union fh_t output;		/* output file handle or name */
    int flags;
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
} cpc_state_t;



static inline double r_cpc(const double phi, cpc_state_t * state)
/*
 * Winston, R.; Miñano, J. C.; Benítez, P.; Shatz, N. and Bortz, J. C.;
 * Nonimaging Optics Elsevier Inc., 2005
 * Equation 4.7
 * Note, we use Patrik Coray's definition of phi with regard to the cpc axis
 * and not with regard to the parabola's axis
 */
{
    return state->f2 * sin(phi) / (1 - cos(phi + state->phi_a)) - state->r;
}

static inline double z_cpc(const double phi, cpc_state_t * state)
/*
 * Winston, R.; Miñano, J. C.; Benítez, P.; Shatz, N. and Bortz, J. C.;
 * Nonimaging Optics Elsevier Inc., 2005
 * Equation 4.7
 * Note, we use Patrik Coray's definition of phi with regard to the cpc axis
 * and not with regard to the parabola's axis
 */
{
    return state->f2 * cos(phi) / (1 - cos(phi + state->phi_a));
}

static inline double dr_dz(const double phi, cpc_state_t * state)
/*
 * we know only r(phi) and z(phi). thus
 *
 *     dr/dz = dr/phi * dphi/dz
 *
 *   dr/dphi = f2 * cos(phi - phi_a) / (-cos(phi) + 1)
 *             - f2 * sin(phi) * sin(phi - phi_a) / (-cos(phi) + 1)**2
 *
 *    dz/phi = -f2 * sin(phi - phi_a) / (-cos(phi) + 1)
 *             - f2 * sin(phi) * cos(phi - phi_a) / (-cos(phi) + 1)**2
 *
 *     dr/dz = dr/dphi / dz/dphi
 *           = (cos(phi_a) - cos(phi - phi_a)) / (sin(phi_a) + sin(phi - phi_a))
 */
{
/*    double sin_phi, cos_phi;
    const double term = phi - state->phi_a;
    double sin_term, cos_term;

    sincos(state->phi_a, &sin_phi, &cos_phi);
    sincos(term, &sin_term, &cos_term);

    return (cos_phi - cos_term) / (sin_phi + sin_term);
*/
    const double tan_phi = tan(phi);
    const double term = phi + state->phi_a;
    double sin_term, cos_term;

    sincos(term, &sin_term, &cos_term);

    return ((tan_phi * sin_term - (1 - cos_term)) /
	    (sin_term + tan_phi * (1 - cos_term)));
}

static void cpc_surf_normal(double const *p, cpc_state_t * state,
			    double *Nl)
/*
 * normal vector N and radius vector R (both at point p) define the
 * the tangent (dr/dz) at point p. in the following N points inwards
 * while R points outwards.
 *
 *     Pl = p - origin
 *     pz = Pl dot axis
 *
 * pz is the length of the projection of Pl onto axis.
 *
 *      R = Pl - pz axis
 *   lenR = length(R)
 *     Rn = lenR R
 *
 * with Rn the normalized radius vector R. 
 *
 * and finally
 *
 *  phi = atan2(lenR + state->r, pz)
 */
{
    double drdz;
    double phi;
    double R[3], lenR;
    double Pl[3], pz;

    diff(Pl, p, state->origin);
    pz = my_ddot(Pl, &state->M[6]);

    a_plus_cb(R, Pl, -pz, &state->M[6]);
    lenR = normalize(R);

    phi = atan2(lenR + state->r, pz);
    drdz = dr_dz(phi, state);

    a_plus_cb(Nl, R, -drdz, &state->M[6]);
    normalize(Nl);
    my_dscal(-1.0, Nl);		/* make point inwards */
}

typedef struct f_parameters_t {
    double *p1;			/* coordinates of starting point p1 */
    double r1;			/* radial coordinate at p1 */
    double *dir;		/* direction vector from p1 to p2 */
    cpc_state_t *state;
} f_parameters_t;

static double residual(double l, void *params)
/*
 * residual at p given by p1 + l * dir
 *
 * residual = r_CPC - R
 *
 * p_l is vector from origin of cpc to p
 * p_l = p - origin
 * p = p1 + l*dir
 *
 * R_l is radius vector of p_l
 * R_l = p_l - z_l*CPC_axis
 * z_l = p_l dot CPC_axis
 * R = length(R_l)
 *
 * phi_l = atan2(R_l + r, z_l)
 * r_CPC = r_cpc(phi)
 *
 * residual is positive: p is inside cpc
 * residual is negative: p is outside cpc
 */
{
    f_parameters_t *par = (f_parameters_t *) params;

    double p[3], p_z, p_l[3];
    double R_l[3], R;
    double phi_l, r_CPC;
    double res;

    a_plus_cb(p, par->p1, l, par->dir);
    diff(p_l, p, par->state->origin);
    p_z = my_ddot(p_l, &par->state->M[6]);
    a_plus_cb(R_l, p_l, -p_z, &par->state->M[6]);
    R = normalize(R_l);

    phi_l = atan2(R + par->state->r, p_z);
    r_CPC = r_cpc(phi_l, par->state);

    res = (r_CPC - R);

    return res;
}

static double *find_intercept(double *p1, double *p2, cpc_state_t * state,
			      double *dir)
{
/*
 * we are searching the distance l from point p1 in direction
 * of the ray's direction vector where the residual
 * 	residual = R - r
 * is zero.
 *	  R: radial component of point p1 + l*dir
 *        r: radial component of P
 * with
 *        P: atan2(r - r_exit, Rz)
 *
 * we return the coordinates of the intercept but guarantee that the point
 * lies on the inside of the cpc (numerical error!) i.e there is no sign
 * change of the residula function.
 */
    int status;
    int iter = 0, max_iter = MAX_ITER;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

    gsl_function F;
    f_parameters_t par;

    double l_min, l1, l2;
    double r1, r2, r_min;
    double P_ref[3], R_ref[3];
    double t[3];
    double *intercept;

    diff(t, p2, p1);
    l2 = normalize(t);
    l1 = 0;

    diff(P_ref, p1, state->origin);
    a_plus_cb(R_ref, P_ref, -my_ddot(P_ref, &state->M[6]), &state->M[6]);

    par.p1 = p1;
    par.r1 = normalize(R_ref);
    par.dir = dir;
    par.state = state;

    r1 = residual(l1, (void *) &par);
    r2 = residual(l2, (void *) &par);

#ifdef DEBUG
    fprintf(stderr, "residual at start point 1: %e\n", r1);
    fprintf(stderr, "residual at start point 2: %e\n", r2);
#endif

    if (r1 < 0) {

#ifdef DEBUG
	fprintf(stderr, "* * * (a) this cannot happen, (p1 outside cpc) returning NULL\n");
#endif

	gsl_root_fsolver_free(s);
	return NULL;
    }
    if (r2 > 0) {

#ifdef DEBUG
	fprintf(stderr, "* * * (b) this cannot happen, (p2 inside cpc) returning NULL\n");
#endif

	gsl_root_fsolver_free(s);
	return NULL;
    }


    F.function = &residual;
    F.params = &par;

    gsl_root_fsolver_set(s, &F, l1, l2);

    do {
	gsl_root_fsolver_iterate(s);

	l_min = gsl_root_fsolver_root(s);
	l1 = gsl_root_fsolver_x_lower(s);
	l2 = gsl_root_fsolver_x_upper(s);

#ifdef DEBUG
	fprintf(stderr, "l1, l2, l_min: %e\t%e\t%e\n", l1, l2, l_min);
#endif

	status = gsl_root_test_interval(l1, l2, ABS_EPS_ROOT, REL_EPS_ROOT);
	iter++;

    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);
    r_min = residual(l_min, (void *) &par);

#ifdef DEBUG
    fprintf(stderr, "number of iterations: %d, abs_error, residual: %e\t%e\n", iter,
	   l2 - l1, r_min);
#endif

    /*
     * calculate intercept as p1 + l_min*dir
     * check if intercept lies on inside of cpc i.e.
     * r_min > 0.
     * otherwise use l1 (x_lower)
     */
#ifdef DEBUG
    fprintf(stderr, "r_min (%e) at l_min (%e) -> ", r_min, l_min);
#endif

    intercept = (double *) malloc(3 * sizeof(double));
    if (r_min > 0) {
	/* l_min is on inside */

#ifdef DEBUG
	fprintf(stderr, "l_min is on inside\n(ok)\n");
#endif

	a_plus_cb(intercept, p1, l_min, dir);
    } else {
	/* too far, use l1 */

#ifdef DEBUG
	fprintf(stderr, "l_min is on outside\n");
	fprintf(stderr, "use l1 (%e): now r_min (%e)\n(ok)\n", l1,
	       residual(l1, (void *) &par));
#endif

	a_plus_cb(intercept, p1, l1, dir);
    }
    return intercept;
}

static void determine_far_aperture(double **origin, double *radius,
				   ray_t * ray, cpc_state_t * state)
{
    if (my_ddot(&state->M[6], ray->dir) > 0) {
	/* exit aperture hit. far aperture is entrance aperture */
	*radius = state->Rt;
	*origin = state->origin2;
    } else {
	/* entrance aperture hit. far aperture is exit aperture */
	*radius = state->r;
	*origin = state->origin;
    }

}


static int cpc_init_state(void *vstate, config_setting_t * this_target,
			  const int file_mode, const int keep_closed,
			  const double P_factor)
{
    cpc_state_t *state = (cpc_state_t *) vstate;
    double h;

    read_vector(this_target, "origin", state->origin);
    config_setting_lookup_float(this_target, "acceptance_angle",
				&state->phi_a);
    config_setting_lookup_float(this_target, "truncation_angle",
				&state->theta_t);
    config_setting_lookup_float(this_target, "exit_radius", &state->r);

    state->phi_a *= M_PI / 180.0;
    state->theta_t *= M_PI / 180.0;
    state->f2 = 2.0 * state->r * (1.0 + sin(state->phi_a));
    state->Rt = r_cpc(state->theta_t, state);
    h = z_cpc(state->theta_t, state);
    state->reflen =
	sqrt(h * h + (state->Rt + state->r) * (state->Rt + state->r));

    state->M = init_M(this_target, "x", "axis");
    a_plus_cb(state->origin2, state->origin, h, &state->M[6]);

    state->flags = 0;
    if (keep_closed)
	state->flags |= KEEP_CLOSED;

    if (init_output
	(TARGET_TYPE, this_target, file_mode, P_factor, &state->output,
	 &state->flags, state->origin, state->M) == ERR) {
	state->refl_spectrum = NULL;
	return ERR;
    }

    init_spectrum(this_target, "reflectivity", &state->refl_spectrum);
    init_refl_model(this_target, &state->refl_func,
		    &state->refl_func_pars);

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);

    return NO_ERR;
}

static void cpc_free_state(void *vstate)
{
    cpc_state_t *state = (cpc_state_t *) vstate;

    state_free(state->output, state->flags, state->M, NULL,
	       state->refl_func, state->refl_func_pars);
    gsl_spline_free(state->refl_spectrum);
}

static double *cpc_get_intercept(void *vstate, ray_t * ray)
{
/*
 * returns closest intercept of 'ray' with cpc.
 * Note: Only rays that pass through the closer aperture but
 *	 not through the far aperture are intercepted.
 */
    cpc_state_t *state = (cpc_state_t *) vstate;

    double *near_intercept;
    double *far_intercept;
    int aperture_hit;
    int dummy;

    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & LAST_WAS_HIT) {	/* ray starts on this target, no hit posible */
	data->flag &= ~LAST_WAS_HIT;	/* clear flag */
	return NULL;
    }

    /*
     * get intercept with closer aperture:
     * exit aperture: (center point is state->origin), if state->M[6] is
     * parallel to ray->dir
     * entrance aperture: (center point is state->origin2), if state->M[6]
     * is anti-parallel to ray->dir
     */
    if (my_ddot(&state->M[6], ray->dir) > 0) {
	aperture_hit = EXIT_APERTURE;

	near_intercept =
	    intercept_disk(ray, state->origin, state->M,
			   state->r * state->r, &dummy);
    } else {
	aperture_hit = ENTRANCE_APERTURE;

	near_intercept =
	    intercept_disk(ray, state->origin2, state->M,
			   state->Rt * state->Rt, &dummy);
    }
    if (!near_intercept)
	/*
	 * no intercept possible as (per definition) rays cannot
	 * intercept with outside of cpc wall (user has to ensure this)
	 */
	return NULL;

    /*
     * intercept found with first aperture. also check
     * if ray hits second aperture and passes through cpc
     */
    if (aperture_hit == ENTRANCE_APERTURE) {
	/* do we intercept exit aperture too? */

#ifdef DEBUG
	fprintf(stderr, "entering CPC by entrance aperture");
#endif

	far_intercept =
	    intercept_disk(ray, state->origin, state->M,
			   state->r * state->r, &dummy);
    } else {
	/* do we intercept entrance aperture too? */

#ifdef DEBUG
	fprintf(stderr, "entering CPC by exit aperture");
#endif

	far_intercept =
	    intercept_disk(ray, state->origin2, state->M,
			   state->Rt * state->Rt, &dummy);
    }
    if (far_intercept) {	/* both apertures hit. ray passes through cpc */
	free(near_intercept);
	free(far_intercept);

#ifdef DEBUG
	fprintf(stderr, " ... and exiting again\n");
#endif

	return NULL;
    }

#ifdef DEBUG
    fprintf(stderr, "\n");
#endif

    /*
     * this is a shortcut as we assume that no other targets extend
     * into the cpc. thus no need to find the real intercept here.
     * the real intercept will be found in cpc_get_out_ray() that
     * will also take care of multiple reflections inside the cpc.
     */
    return near_intercept;
}

static ray_t *cpc_get_out_ray(void *vstate, ray_t * ray, double *hit,
			      const gsl_rng * r)
{
    cpc_state_t *state = (cpc_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    double *center_far_aperture;
    double far_R;
    double *far_hit;
    double *point;
    double normal[3];
    double R_to_center[3], Rn;
    int dummy;
    /*
     * determine interception with plane of other_aperture. then we
     * know that real intercept lies between hit and far_hit and
     * we can start with the root finder.
     */
    determine_far_aperture(&center_far_aperture, &far_R, ray, state);
    far_hit =
	intercept_plane(ray, &state->M[6], center_far_aperture, &dummy);

    do {
	/* find next intercept with cpc wall (point).
	 * we know that intercept lies between the last intercept (hit)
	 * and far_hit i.e. the intercept with plane in which far aperture
	 * lies.
	 * when we enter the function hit represents the intercept of
	 * ray with the closest aperture of the cpc. while we are in the
	 * loop hit represents the previous intercept with the cpc's wall
	 * (plus a small offset)
	 */
	point = find_intercept(hit, far_hit, state, ray->dir);
	free(far_hit);

	if (point == NULL) {
	    /*
	     * this should never happen as it means that either
	     * hit is outside of cpc
	     * or far_hit is inside cpc
	     * and the root solver could not start as there is no
	     * sign change in the residual. these rays are counted
	     * as lost.
	     */
	    free(ray);
	    return NULL;
	}

	if (gsl_rng_uniform(r) >
	    gsl_spline_eval(state->refl_spectrum, ray->lambda, NULL)) {
	    /*
	     * if ray is absorbed at the cpc wall store intercept (point)
	     * in local coordinates
	     */
	    if (state->flags & OUTPUT_REQUIRED)
		store_xy(state->output, state->flags, ray, point,
			 state->M, state->origin, data,
			 &state->mutex_writefd);

	    data->flag &= ~LAST_WAS_HIT;	/* clear flag */

	    free(point);
	    free(ray);

#ifdef DEBUG
	    fprintf(stderr, "absorbed in cpc\n");
#endif

	    return NULL;
	}

	/*
	 * reflect ray at point.
	 */
	cpc_surf_normal(point, state, normal);
	state->refl_func(ray, normal, point, r, state->refl_func_pars);
	free(point);

#ifdef DEBUG
	fprintf(stderr, "reflected in cpc\n");
#endif
	/*
	 * reflected ray originates at point but might now point in direction
	 * of a different aperture than before reflection.
	 * thus determine parameters of far aperture again.
	 */
	determine_far_aperture(&center_far_aperture, &far_R, ray, state);
	/*
	 * determine intercept with plane of far aperture
	 */
	far_hit =
	    intercept_plane(ray, &state->M[6], center_far_aperture,
			    &dummy);

	diff(R_to_center, far_hit, state->origin);
	Rn = normalize(R_to_center);
	if (Rn < far_R) {
	    /*
	     * ray hits aperture and exits cpc and we are done
	     */
	    memcpy(ray->orig, far_hit, 3 * sizeof(double));
	    free(far_hit);
	    data->flag |= LAST_WAS_HIT;	/* mark as hit */

#ifdef DEBUG
	    fprintf(stderr, "exits cpc\n");
#endif

	    return ray;
	}
	/*
	 * reinitialize hit.
	 */
	memcpy(hit, ray->orig, 3 * sizeof(double));

    } while (1);
}

static void cpc_init_PTDT(void *vstate)
{
    per_thread_init(((cpc_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void cpc_flush_PTDT_outbuf(void *vstate)
{
    cpc_state_t *state = (cpc_state_t *) vstate;

    per_thread_flush(state->output, state->flags, state->PTDT_key,
		     &state->mutex_writefd);
}


static const target_type_t cpc_t = {
    TARGET_TYPE,
    sizeof(struct cpc_state_t),
    &cpc_init_state,
    &cpc_free_state,
    &cpc_get_intercept,
    &cpc_get_out_ray,
    &cpc_init_PTDT,
    &cpc_flush_PTDT_outbuf
};

const target_type_t *target_cpc = &cpc_t;

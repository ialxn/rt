/*	target_window.c
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

#include "io_utils.h"
#include "intercept.h"
#include "reflect.h"
#include "targets.h"

#define TARGET_TYPE "window"
#define NO_ITEMS 4

typedef struct window_state_t {
    double C[3];		/* center coordinate of first face */
    double r;			/* radius of window */
    double d;			/* tickness of window */
    double *M;			/* transform matrix local -> global coordinates */
    gsl_spline *abs_spectrum;	/* for interpolated absorptivity spectrum */
    gsl_spline *dispersion;	/* for interpolated dispersion curve */
    refl_func_pointer_t refl_func;	/* reflection model */
    void *refl_func_pars;	/* model specific parameters */
    union fh_t output;		/* output file handle or name */
    int flags;
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
} window_state_t;



static double R_fresnell(ray_t * ray, const double *normal,
			 const double n1, const double n2)
{
/*
 * return reflectivity (fresnell) at planar interface defined by
 * its normal. ray passes from medium with index of refraction
 * n1 to medium with index of refraction n2.
 * total internal reflection is no special case and R=1.0 is returned
 */
    double rs, rp;
    double R;
    double t1, t2;

    const double cos_phi = my_ddot(ray->dir, normal);
    const double sin_sqr_phi = 1.0 - cos_phi * cos_phi;
    double a;

    a = 1.0 - n1 * n1 / (n2 * n2) * sin_sqr_phi;
    if (a < 0.0)		/* total internal reflection */
	return 1.0;

    a = sqrt(a);

    t1 = n1 * cos_phi;
    t2 = n2 * a;
    rs = (t1 - t2) / (t1 + t2);
    rs *= rs;

    t1 = n1 * a;
    t2 = n2 * cos_phi;
    rp = (t1 - t2) / (t1 + t2);
    rp *= rp;

    R = (rs + rp) * 0.5;

    return R;
}

static int snell(ray_t * ray, const double *normal, const double n1,
		 const double n2)
/*
 * sets new 'ray->dir' of ray which passes from medium with index of
 * refraction 'n1' to medium with index of refraction 'n2' according to
 * Snell's law.
 * Heckbert's algorithm (Heckbert2008) is used.
 * Appeared in Introduction to Ray Tracing, (Andrew Glassner, ed.),
 * Academic Press, London, 1989, pp. 263-293.
 *
 * return 0 (and original 'ray->dir') if total internal reflection occurs
 * otherwise return 1.
 *
 * Note: - Heckbert uses eta=1/n
 *       - uses normalized vectors
 *       - 'normal' points towards propagation vector of ray
 *         we check for both cases
 */
{
    const double eta = n1 / n2;
    double c1, cs2;

    c1 = my_ddot(ray->dir, normal);
    cs2 = 1.0 - eta * eta * (1.0 - c1 * c1);

    if (cs2 < 0.0)
	return 1;
    cs2 = sqrt(cs2);

    my_dscal(eta, ray->dir);
    my_daxpy(-eta * c1 + cs2, normal, ray->dir);

    return 0;
}

static int window_init_state(void *vstate, config_setting_t * this_target,
			     const int file_mode, const int keep_closed,
			     const double P_factor)
{
    window_state_t *state = (window_state_t *) vstate;

    read_vector(this_target, "C", state->C);

    state->M = init_M(this_target, "x", "a");

    state->flags = 0;
    if (keep_closed)
	state->flags |= KEEP_CLOSED;

    if (init_output
	(TARGET_TYPE, this_target, file_mode, P_factor, &state->output,
	 &state->flags, state->C, state->M) == ERR) {
	state->abs_spectrum = NULL;
	state->dispersion = NULL;
	return ERR;
    }

    /* initialize absorptivity spectrum */
    if (init_spectrum(this_target, "absorptivity", &state->abs_spectrum)) {
	state->dispersion = NULL;
	return ERR;
    }

    /* initialize dispersion curve */
    if (init_spectrum(this_target, "idx_refraction", &state->dispersion)) {
	gsl_spline_free(state->abs_spectrum);
	state->abs_spectrum = NULL;
	return ERR;
    }

    init_refl_model(this_target, &state->refl_func,
		    &state->refl_func_pars);

    config_setting_lookup_float(this_target, "r", &state->r);
    config_setting_lookup_float(this_target, "d", &state->d);

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);

    return NO_ERR;
}

static void window_free_state(void *vstate)
{
    window_state_t *state = (window_state_t *) vstate;

    state_free(state->output, state->flags, state->M, NULL,
	       state->refl_func, state->refl_func_pars);
    gsl_spline_free(state->abs_spectrum);
    gsl_spline_free(state->dispersion);
}

static double *intercept_face(const ray_t * in_ray,
			      const window_state_t * state,
			      const double *center)
{
    double *intercept;
    double I[3];
    int dummy_i;

    if ((intercept =
	 intercept_plane(in_ray, &state->M[6], center,
			 &dummy_i)) != NULL) {

	g2l(state->M, center, intercept, I);

	if ((I[0] * I[0] + I[1] * I[1]) > (state->r * state->r)) {
	    /* hit not within boundaries */
	    free(intercept);
	    intercept = NULL;	/* mark as not valid */
	}
    }

    return intercept;
}

static double *window_get_intercept(void *vstate, ray_t * ray)
{
/*
 * returns closest intercept of 'ray' with window.
 * Note: we have to test for intercept with both faces
 *       AND with the (absorbing) cylinder wall.
 */
    window_state_t *state = (window_state_t *) vstate;

    double *intercept;
    double center_face2[3];
    int surf_hit;

    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & LAST_WAS_HIT) {	/* ray starts on this target, no hit posible */
	data->flag &= ~LAST_WAS_HIT;
	return NULL;
    }

    /*
     * get intercepts with closer face:
     * face1 (center point is state->C), if state->a is parallel to
     * ray->dir
     * face2 (center point is state->C + state->d*state->a), if
     * if state->a is anti-parallel to ray->dir
     */
    if (my_ddot(&state->M[6], ray->dir) > 0)
	intercept = intercept_face(ray, state, state->C);
    else {
	a_plus_cb(center_face2, state->C, state->d, &state->M[6]);
	intercept = intercept_face(ray, state, center_face2);
    }

    if (intercept)		/* intecept found with face, no need to check wall */
	return intercept;

    intercept =
	intercept_cylinder(ray, state->C, &state->M[6], state->r, state->d,
			   &surf_hit);

    if (intercept)		/* intercept with outside wall */
	data->flag |= ABSORBED;

    return intercept;
}

static ray_t *window_get_out_ray(void *vstate, ray_t * ray, double *hit,
				 const gsl_rng * r)
{
    window_state_t *state = (window_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);
    int origin_is_face1;
    double center[3];
    double normal[3];
    int inside = 1;
    double R_tot;
    double n_in;

    if (data->flag & ABSORBED) {
	/*
	 * if ABSORBED is set we know ray has been absorbed
	 * because it was intercepted by cylinder wall with
	 * absorptivity=1. this was checked (and the flag was set)
	 * in 'window_get_intercept()' above.
	 */

	if (state->flags & OUTPUT_REQUIRED)
	    store_xy(state->output, state->flags, ray, hit, state->M,
		     state->C, data, &state->mutex_writefd);

	data->flag &= ~(LAST_WAS_HIT | ABSORBED);	/* clear flags */

	free(ray);
	return NULL;

    }
    /*
     * test which face of window was hit.
     * face1 (center point is state->C), if state->a is parallel to
     * ray->dir
     * face2 (center point is state->C + state->d*state->a), if
     * if state->a is anti-parallel to ray->dir
     */
    if (my_ddot(ray->dir, &state->M[6]) > 0) {
	origin_is_face1 = 1;
	memcpy(center, state->C, 3 * sizeof(double));
	memcpy(normal, &state->M[6], 3 * sizeof(double));
    } else {
	origin_is_face1 = 0;
	a_plus_cb(center, state->C, state->d, &state->M[6]);
	a_times_const(normal, &state->M[6], -1.0);
    }

    n_in = gsl_spline_eval(state->dispersion, ray->lambda, NULL);
    R_tot = R_fresnell(ray, normal, 1.0, n_in);
    if (gsl_rng_uniform(r) <= R_tot) {	/* external reflection occurs */
	/*
	 * rays can reflect off both faces of the window. if face1 is hit
	 * the normal vector 'state->a' is parallel to 'ray->dir' and
	 * anti.parallel if face2 is hit ('state->dir points from face1 to
	 * face2). does not seem to matter for 'reflect()'.
	 */
	state->refl_func(ray, normal, hit, r, state->refl_func_pars);

	data->flag |= LAST_WAS_HIT;
	return ray;
    }

    /*
     * ray is not externally reflected and enters window.
     * update its origin (originates at 'hit')
     */
    memcpy(ray->orig, hit, 3 * sizeof(double));

    /*
     * calculate new direction after ray has entered window. we do not
     * need to check for total internal reflection as we pass from less
     * dense medium into the denser one.
     * Note: 'ray' cannot originate from inside the window.
     */
    snell(ray, normal, 1.0, n_in);

    while (inside) {
	/*
	 * we checked before whether we entered via face 1 or 2
	 * by definition:
	 *     - 'a' dot 'ray->dir' > 0 (parallel) at face 1
	 *     - 'a' dot 'ray->dir' < 0 (anti-parallel) at face 2
	 * because 'a' points from face 1 to face 2.
	 *
	 * inside the while loop we will keep track of the origin
	 * of the ray as follows:
	 *     - 'origin_is_face1' % 2 == 1 YES
	 *     - 'origin_is_face2' % 2 == 0 NO
	 * after every reflection at any of the two faces. 'origin_is_face1++'
	 * switches origin from one face to the other. here we just initialize
	 * the flag properly.
	 *
	 *
	 * outline of algorithm:
	 * hits other face?
	 *  YES:
	 *          ABSORBED inside window?
	 *                  YES:
	 *                          return NULL
	 *                  NO:
	 *                          exits window?
	 *                                  YES:
	 *                                          return ray
	 *                                  NO:
	 *                                          reflect ray
	 *  NO:
	 *          find intercept with wall
	 *          return NULL
	 */
	double *intercept;
	double center_other_face[3];
	double normal_other_face[3];
	/*
	 * calculate center of other face and calculate intercept
	 * with ray.
	 * Note: face2 is at a distance of 'd' in direction 'a'
	 *       from face1:
	 *       'center_other_face' = 'state->C' + 'state->d' * 'state->a'
	 */
	if (origin_is_face1 % 2) {	/* other face is face2 */
	    a_plus_cb(center_other_face, state->C, state->d, &state->M[6]);
	    memcpy(normal_other_face, &state->M[6], 3 * sizeof(double));
	} else {		/* other face is face1 */
	    memcpy(center_other_face, state->C, 3 * sizeof(double));
	    a_times_const(normal_other_face, &state->M[6], -1.0);
	}

	if ((intercept =
	     intercept_face(ray, state, center_other_face)) != NULL) {
	    /*
	     * ray hits other face.
	     * no need to check wall but check if ray is absorbed between the two faces
	     */
	    const double d_travelled = sqrt(d_sqr(ray->orig, intercept));
	    const double a_coeff =
		gsl_spline_eval(state->abs_spectrum, ray->lambda,
				NULL);
	    const double A = exp(-d_travelled / a_coeff);

	    if (gsl_rng_uniform(r) <= A) {	/* ray is absorbed */
		/*
		 * FIXME: ray is NOT absorbed at 'intercept' i.e. at face1 or face2
		 *        but earlier between 'ray->orig' and 'intercept'.
		 *        error introduced is small as, generally:
		 *            - windows are (relatively) thin and most rays
		 *              travel at a low angle relative to 'state->a'.
		 *              both result in a small lateral offset between
		 *              'intercept' and the position where ray is absorbed.
		 *            - only a few ray will be absorbed inside window
		 *              as windows absorptivity is small (definition).
		 */
		if (state->flags & OUTPUT_REQUIRED)
		    store_xy(state->output, state->flags, ray,
			     intercept, state->M, state->C, data,
			     &state->mutex_writefd);

		data->flag &= ~(LAST_WAS_HIT | ABSORBED);

		free(ray);
		free(intercept);

		return NULL;
	    }

	    /*
	     * ray reflected at other window?
	     * 'R_tot' = 1.0 in case of total internal reflection at face.
	     */
	    R_tot = R_fresnell(ray, normal_other_face, n_in, 1.0);
	    if (gsl_rng_uniform(r) <= R_tot) {
		/*
		 * fresnel reflection or total internal reflection occurs.
		 * reflect ray and switch face.
		 *
		 * rays can reflect off both faces of the window. if face1 is hit
		 * the normal vector 'state->a' is parallel to 'ray->dir' and
		 * anti.parallel if face2 is hit ('state->dir points from face1 to
		 * face2). does not seem to matter for 'reflect()'.
		 */
		state->refl_func(ray, normal_other_face, intercept, r,
				 state->refl_func_pars);

		origin_is_face1++;
		free(intercept);

	    } else {
		/*
		 * ray leaves window.
		 * calculated new direction and update origin of ray.
		 */
		snell(ray, normal_other_face, n_in, 1.0);
		memcpy(ray->orig, intercept, 3 * sizeof(double));

		free(intercept);
		data->flag |= LAST_WAS_HIT;	/* mark as hit */
		inside = 0;	/* terminate while loop */
	    }

	} else {
	    /*
	     * ray hits wall and is absorbed
	     */
	    int surf_hit;

	    intercept =
		intercept_cylinder(ray, state->C, &state->M[6], state->r,
				   state->d, &surf_hit);
	    if (state->flags & OUTPUT_REQUIRED)
		store_xy(state->output, state->flags, ray, hit,
			 state->M, state->C, data, &state->mutex_writefd);

	    data->flag &= ~(LAST_WAS_HIT | ABSORBED);

	    free(ray);
	    free(intercept);

	    return NULL;
	}
    }				/* end while(inside) */
    return ray;
}

static void window_init_PTDT(void *vstate)
{
    per_thread_init(((window_state_t *) vstate)->PTDT_key,
		    NO_ITEMS * sizeof(float) + sizeof(unsigned char));
}

static void window_flush_PTDT_outbuf(void *vstate)
{
    window_state_t *state = (window_state_t *) vstate;

    per_thread_flush(state->output, state->flags, state->PTDT_key,
		     &state->mutex_writefd);
}


static const target_type_t window_t = {
    TARGET_TYPE,
    sizeof(struct window_state_t),
    &window_init_state,
    &window_free_state,
    &window_get_intercept,
    &window_get_out_ray,
    &window_init_PTDT,
    &window_flush_PTDT_outbuf
};

const target_type_t *target_window = &window_t;

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

#include "io_util.h"
#include "targets.h"

#define NO_ITEMS 4


typedef struct sq_state_t {
    char *name;			/* name (identifier) of target */
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
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

    pthread_key_create(&state->PTDT_key, free_PTDT);
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


    double *intercept;
    double l_intercept[3];

    int hits_front;

    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & LAST_WAS_HIT) {	/* ray starts on this target, no hit posible */
	data->flag &= ~LAST_WAS_HIT;
	return NULL;
    }

    intercept =
	intercept_plane(in_ray, state->normal, state->point, &hits_front);

    if (!intercept)		/* ray does not hit target */
	return NULL;

    /* convert to local coordinates, origin is 'state->point' */
    g2l(state->M, state->point, intercept, l_intercept);

    if ((l_intercept[0] <= -state->dx) || (l_intercept[0] >= state->dx)
	|| (l_intercept[1] <= -state->dy)
	|| (l_intercept[1] >= state->dy)) {

	/* hit not within boundaries */
	free(intercept);
	return NULL;

    } else {			/* hits within target dimensions 'dx' times 'dy' */

	if (!hits_front)	/* hits rear side, absorbed */
	    data->flag |= ABSORBED;

	return intercept;

    }
}

static ray_t *sq_get_out_ray(void *vstate, ray_t * in_ray, double *hit,
			     const gsl_rng * r)
{
    sq_state_t *state = (sq_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & ABSORBED
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

	store_xy(state->dump_file, in_ray, hit, state->M, state->point,
		 data);

	data->flag &= ~(LAST_WAS_HIT | ABSORBED);	/* clear flags */

	free(in_ray);
	return NULL;

    } else {			/* reflect 'in_ray' */
	reflect(in_ray, state->normal, hit);

	data->flag |= LAST_WAS_HIT;	/* mark as hit */

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

static void sq_init_PTDT(void *vstate)
{
    sq_state_t *state = (sq_state_t *) vstate;
    PTDT_t *data = (PTDT_t *) malloc(sizeof(PTDT_t));

    data->buf = (float *) malloc(BUF_SIZE * NO_ITEMS * sizeof(float));
    data->i = 0;
    data->flag = 0;

    pthread_setspecific(state->PTDT_key, data);
}

static void sq_flush_PTDT_outbuf(void *vstate)
{
    sq_state_t *state = (sq_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->i != 0)		/* write rest of buffer to file. */
	write(state->dump_file, data->buf, sizeof(float) * data->i);
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
    &sq_init_PTDT,
    &sq_flush_PTDT_outbuf
};

const target_type_t *target_rectangle = &sq_t;

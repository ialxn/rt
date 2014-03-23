/*	target_disk.c
 *
 * Copyright (C) 2010,2011,2012,2013,2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <math.h>
#include <string.h>

#include "io_utils.h"
#include "reflect.h"
#include "targets.h"

#define NO_ITEMS 4


typedef struct disk_state_t {
    char *name;			/* name (identifier) of target */
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    int dump_file;
    double point[3];		/* center coordinate of disk */
    double normal[3];		/* normal vector of disk */
    double r2;			/* radius^2 of disk */
    gsl_spline *spline;		/* for interpolated reflectivity spectrum */
    double M[9];		/* transform matrix local -> global coordinates */
} disk_state_t;


static void disk_init_state(void *vstate, config_setting_t * this_target,
			    const int file_mode)
{
    disk_state_t *state = (disk_state_t *) vstate;

    const char *S;
    char f_name[256];
    double t;

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

    config_setting_lookup_float(this_target, "r", &t);
    state->r2 = t * t;

    pthread_key_create(&state->PTDT_key, free_PTDT);
}

static void disk_free_state(void *vstate)
{
    disk_state_t *state = (disk_state_t *) vstate;

    close(state->dump_file);

    free(state->name);
    gsl_spline_free(state->spline);
}

static double *disk_get_intercept(void *vstate, ray_t * ray)
{
    disk_state_t *state = (disk_state_t *) vstate;

    double *intercept;
    double l_intercept[3];
    double r2_intercept;

    int hits_front;

    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & LAST_WAS_HIT) {	/* ray starts on this target, no hit posible */
	data->flag &= ~LAST_WAS_HIT;
	return NULL;
    }

    intercept =
	intercept_plane(ray, state->normal, state->point, &hits_front);

    if (!intercept)		/* ray does not hit target */
	return NULL;

    /* convert to local coordinates, origin is 'state->point' */
    g2l(state->M, state->point, intercept, l_intercept);

    /*
     * r2_intercep is sqared distance from center of disk to intercept
     * in the plane of the disk. we are in local system that
     * is offset by intercept[2], so leave the latter out.
     * compare r^2 to avoid sqrt()
     */
    r2_intercept =
	intercept[0] * intercept[0] + intercept[1] * intercept[1];
    if (r2_intercept > state->r2) {

	/* hit not within boundaries */
	free(intercept);
	return NULL;

    } else {			/* hits within disk radius 'r' */

	if (!hits_front)	/* hits rear side, absorbed */
	    data->flag |= ABSORBED;

	return intercept;

    }
}

static ray_t *disk_get_out_ray(void *vstate, ray_t * ray, double *hit,
			       const gsl_rng * r)
{
    disk_state_t *state = (disk_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & ABSORBED
	|| (gsl_rng_uniform(r) >
	    gsl_spline_eval(state->spline, ray->lambda, NULL))) {
	/*
	 * if ABSORBED is set we know ray has been absorbed
	 * because it was intercepted by a surface with absorptivity=1
	 * (reflectivity=0) e.g. the backside of the target. this was
	 * checked (and the flag was set) in 'xxx_get_intercept()'
	 * above.
	 * then we check if ray is absorbed because the reflectivity of
	 * the mirror surface is less than 1.0 (absorptivity > 0.0).
	 */

	store_xy(state->dump_file, ray, hit, state->M, state->point, data);

	data->flag &= ~(LAST_WAS_HIT | ABSORBED);	/* clear flags */

	free(ray);
	return NULL;

    } else {			/* reflect 'ray' */
	reflect(ray, state->normal, hit);

	data->flag |= LAST_WAS_HIT;	/* mark as hit */

	return ray;
    }
}

static const char *disk_get_target_name(void *vstate)
{
    disk_state_t *state = (disk_state_t *) vstate;

    return state->name;
}

static void disk_dump_string(void *vstate, const char *str)
{
    disk_state_t *state = (disk_state_t *) vstate;

    write(state->dump_file, str, strlen(str));
}

static double *disk_M(void *vstate)
{
    disk_state_t *state = (disk_state_t *) vstate;

    return state->M;
}

static void disk_init_PTDT(void *vstate)
{
    disk_state_t *state = (disk_state_t *) vstate;
    PTDT_t *data = (PTDT_t *) malloc(sizeof(PTDT_t));

    data->buf = (float *) malloc(BUF_SIZE * NO_ITEMS * sizeof(float));
    data->i = 0;
    data->flag = 0;

    pthread_setspecific(state->PTDT_key, data);
}

static void disk_flush_PTDT_outbuf(void *vstate)
{
    disk_state_t *state = (disk_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->i != 0)		/* write rest of buffer to file. */
	write(state->dump_file, data->buf, sizeof(float) * data->i);
}


static const target_type_t disk_t = {
    "disk",
    sizeof(struct disk_state_t),
    &disk_init_state,
    &disk_free_state,
    &disk_get_intercept,
    &disk_get_out_ray,
    &disk_get_target_name,
    &disk_dump_string,
    &disk_M,
    &disk_init_PTDT,
    &disk_flush_PTDT_outbuf
};

const target_type_t *target_disk = &disk_t;

/*	target_plane_screen.c
 *
 * Copyright (C) 2010,2011,2012,2013 Ivo Alxneit
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


typedef struct ps_state_t {
    char *name;			/* name (identifier) of target */
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    char one_sided;		/* flag [one-sided|two-sided] */
    int dump_file;
    double point[3];		/* point on plane */
    double normal[3];		/* normal vector of plane */
    double M[9];		/* transform matrix local -> global coordinates */
} ps_state_t;


static void ps_init_state(void *vstate, config_setting_t * this_target,
			  const int file_mode)
{
    ps_state_t *state = (ps_state_t *) vstate;

    const char *S;
    char f_name[256];

    config_setting_lookup_string(this_target, "name", &S);
    state->name = strdup(S);

    snprintf(f_name, 256, "%s.dat", state->name);
    state->dump_file =
	open(f_name, O_CREAT | O_WRONLY | file_mode, S_IRUSR | S_IWUSR);

    read_vector(this_target, "point", state->point);
    /*
     * generate transform matrix M to convert
     * between local and global coordinates
     * l2g:   g(x, y, z) = M l(x, y, z) + o(x, y, z))
     * g2l:   l(x, y, z) = MT (g(x, y, z) - o(x, y, z))
     */
    /* get normal vector of plane (serving also as basis vector z) */
    read_vector_normalize(this_target, "normal", state->normal);
    memcpy(&state->M[6], state->normal, 3 * sizeof(double));

    /* get basis vector x */
    read_vector_normalize(this_target, "x", state->M);

    /* state->M[3-5] = y = z cross x */
    cross_product(&state->M[6], state->M, &state->M[3]);

    pthread_key_create(&state->PTDT_key, free_PTDT);
}

static void ps1_init_state(void *vstate, config_setting_t * this_target,
			   const int file_name)
{
    ps_state_t *state = (ps_state_t *) vstate;

    ps_init_state(vstate, this_target, file_name);
    state->one_sided = 1;
}

static void ps2_init_state(void *vstate, config_setting_t * this_target,
			   const int file_name)
{
    ps_state_t *state = (ps_state_t *) vstate;

    ps_init_state(vstate, this_target, file_name);
    state->one_sided = 0;
}

static void ps_free_state(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    close(state->dump_file);

    free(state->name);
}

static double *ps_get_intercept(void *vstate, ray_t * in_ray)
{
    ps_state_t *state = (ps_state_t *) vstate;

    double *intercept;

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

    /*
     * only accept intercepts if rear side is hit
     */
    if (state->one_sided && hits_front) {
	free(intercept);
	return NULL;
    }

    return intercept;
}

static ray_t *ps_get_out_ray(void *vstate, ray_t * in_ray, double *hit,
			     const gsl_rng * r)
{
    ps_state_t *state = (ps_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    (void) r;			/* avoid warning : unused parameter 'r' */

    store_xy(state->dump_file, in_ray, hit, state->M, state->point, data);

    data->flag &= ~LAST_WAS_HIT;	/* clear flag */

    memcpy(in_ray->origin, hit, 3 * sizeof(double));	/* update origin */
    return in_ray;
}

static const char *ps_get_target_name(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    return state->name;
}

static void ps_dump_string(void *vstate, const char *str)
{
    ps_state_t *state = (ps_state_t *) vstate;

    write(state->dump_file, str, strlen(str));
}

static double *ps_M(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    return state->M;
}

static void ps_init_PTDT(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;
    PTDT_t *data = (PTDT_t *) malloc(sizeof(PTDT_t));

    data->buf = (float *) malloc(BUF_SIZE * NO_ITEMS * sizeof(float));
    data->i = 0;
    data->flag = 0;

    pthread_setspecific(state->PTDT_key, data);
}

static void ps_flush_PTDT_outbuf(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->i != 0)		/* write rest of buffer to file. */
	write(state->dump_file, data->buf, sizeof(float) * data->i);
}


static const target_type_t ps1_t = {
    "one-sided plane screen",
    sizeof(struct ps_state_t),
    &ps1_init_state,
    &ps_free_state,
    &ps_get_intercept,
    &ps_get_out_ray,
    &ps_get_target_name,
    &ps_dump_string,
    &ps_M,
    &ps_init_PTDT,
    &ps_flush_PTDT_outbuf
};

static const target_type_t ps2_t = {
    "two-sided plane screen",
    sizeof(struct ps_state_t),
    &ps2_init_state,
    &ps_free_state,
    &ps_get_intercept,
    &ps_get_out_ray,
    &ps_get_target_name,
    &ps_dump_string,
    &ps_M,
    &ps_init_PTDT,
    &ps_flush_PTDT_outbuf
};

const target_type_t *target_plane_screen_one_sided = &ps1_t;
const target_type_t *target_plane_screen_two_sided = &ps2_t;

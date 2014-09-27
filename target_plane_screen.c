/*	target_plane_screen.c
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
#include "targets.h"

#define TARGET_TYPE_A "one-sided plane screen"
#define TARGET_TYPE_B "two-sided plane screen"
#define NO_ITEMS 4


typedef struct ps_state_t {
    char *name;			/* name (identifier) of target */
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
    pthread_mutex_t mutex_writefd;	/* protect write(2) */
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
    int i;

    config_setting_lookup_string(this_target, "name", &S);
    state->name = strdup(S);

    if (config_setting_lookup_bool(this_target, "no_output", &i) ==
	CONFIG_FALSE || i == 0) {
	snprintf(f_name, 256, "%s.dat", state->name);
	state->dump_file =
	    open(f_name, O_CREAT | O_WRONLY | file_mode,
		 S_IRUSR | S_IWUSR);
    } else
	state->dump_file = -1;

    read_vector(this_target, "point", state->point);
    /*
     * generate transform matrix M to convert
     * between local and global coordinates
     * l2g:   g(x, y, z) = MT l(x, y, z) + o(x, y, z)
     * g2l:   l(x, y, z) = M (g(x, y, z) - o(x, y, z))
     */
    /* get normal vector of plane (serving also as basis vector z) */
    read_vector_normalize(this_target, "normal", state->normal);
    memcpy(&state->M[6], state->normal, 3 * sizeof(double));

    /* get basis vector x */
    read_vector_normalize(this_target, "x", state->M);

    /* state->M[3-5] = y = z cross x */
    cross_product(&state->M[6], state->M, &state->M[3]);

    /* write header to dump file */
    if (state->dump_file != -1 && file_mode == O_TRUNC) {
	if (state->one_sided)
	    write_target_header(state->dump_file, state->name,
				TARGET_TYPE_A, state->point, state->M);
	else
	    write_target_header(state->dump_file, state->name,
				TARGET_TYPE_A, state->point, state->M);
    }

    pthread_key_create(&state->PTDT_key, free_PTDT);
    pthread_mutex_init(&state->mutex_writefd, NULL);
}

static void ps1_init_state(void *vstate, config_setting_t * this_target,
			   const int file_name)
{
    ps_state_t *state = (ps_state_t *) vstate;

    state->one_sided = 1;
    ps_init_state(vstate, this_target, file_name);
}

static void ps2_init_state(void *vstate, config_setting_t * this_target,
			   const int file_name)
{
    ps_state_t *state = (ps_state_t *) vstate;

    state->one_sided = 0;
    ps_init_state(vstate, this_target, file_name);
}

static void ps_free_state(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    if (state->dump_file != -1)
	close(state->dump_file);

    free(state->name);
}

static double *ps_get_intercept(void *vstate, ray_t * ray)
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
	intercept_plane(ray, state->normal, state->point, &hits_front);

    if (!intercept)		/* ray does not hit target */
	return NULL;

    /*
     * only accept intercepts if front side is hit
     */
    if (state->one_sided && !hits_front) {
	free(intercept);
	return NULL;
    }

    return intercept;
}

static ray_t *ps_get_out_ray(void *vstate, ray_t * ray, double *hit,
			     const gsl_rng * r)
{
    ps_state_t *state = (ps_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    (void) r;			/* avoid warning : unused parameter 'r' */

    if (state->dump_file != -1)
	store_xy(state->dump_file, ray, hit, state->M, state->point,
		 data, &state->mutex_writefd);

    data->flag &= ~LAST_WAS_HIT;	/* clear flag */

    memcpy(ray->orig, hit, 3 * sizeof(double));	/* update origin */
    return ray;
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
	if (state->dump_file != -1) {

	    pthread_mutex_lock(&state->mutex_writefd);
	    write(state->dump_file, data->buf, sizeof(float) * data->i);
	    pthread_mutex_unlock(&state->mutex_writefd);

	}
}


static const target_type_t ps1_t = {
    TARGET_TYPE_A,
    sizeof(struct ps_state_t),
    &ps1_init_state,
    &ps_free_state,
    &ps_get_intercept,
    &ps_get_out_ray,
    &ps_init_PTDT,
    &ps_flush_PTDT_outbuf
};

static const target_type_t ps2_t = {
    TARGET_TYPE_B,
    sizeof(struct ps_state_t),
    &ps2_init_state,
    &ps_free_state,
    &ps_get_intercept,
    &ps_get_out_ray,
    &ps_init_PTDT,
    &ps_flush_PTDT_outbuf
};

const target_type_t *target_plane_screen_one_sided = &ps1_t;
const target_type_t *target_plane_screen_two_sided = &ps2_t;

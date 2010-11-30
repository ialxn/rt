/*	target_plane_screen.c
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>

#include "ray.h"
#include "targets.h"

typedef struct ps_state_t {
    char *name;
    FILE *dump_file;
    double point[3];
    double normal[3];
    size_t n_alloc;
    size_t n_data;
    double *data;
} ps_state_t;


static int ps_alloc_state(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    if (!
	(state->data = (double *) malloc(3 * BLOCK_SIZE * sizeof(double))))
	return ERR;

    state->n_alloc = BLOCK_SIZE;

    return NO_ERR;
}

static void ps_init_state(void *vstate, config_t * cfg, const char *name)
{
    ps_state_t *state = (ps_state_t *) vstate;

    unsigned int i = 0;
    const char *S;
    double F;
    char f_name[256];

    config_setting_t *this_target, *point, *normal;
    const config_setting_t *targets = config_lookup(cfg, "targets");

    state->name = strdup(name);

    snprintf(f_name, 256, "%s.dat", name);
    state->dump_file = fopen(f_name, "w");

    while (1) {			/* find setting for target 'name' */
	this_target = config_setting_get_elem(targets, i);

	config_setting_lookup_string(this_target, "name", &S);
	if (strstr(S, name))
	    break;

	i++;
    }

    point = config_setting_get_member(this_target, "point");
    config_setting_lookup_float(point, "x", &F);
    state->point[0] = F;
    config_setting_lookup_float(point, "y", &F);
    state->point[1] = F;
    config_setting_lookup_float(point, "z", &F);
    state->point[2] = F;

    normal = config_setting_get_member(this_target, "normal");
    config_setting_lookup_float(normal, "x", &F);
    state->normal[0] = F;
    config_setting_lookup_float(normal, "y", &F);
    state->normal[1] = F;
    config_setting_lookup_float(normal, "z", &F);
    state->normal[2] = F;

    state->n_data = 0;
}

static void ps_free_state(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    size_t i;

    const size_t n = state->n_data;
    double *data = state->data;

    /* first write data to file */
    for (i = 0; i < n; i++) {
	const size_t idx = 3 * i;

	fprintf(state->dump_file, "%g\t%g\t%g\n", data[idx],
		data[idx + 1], data[idx + 2]);
    }

    fclose(state->dump_file);
    free(state->name);
    free(state->data);
}

static double *ps_get_intercept(void *vstate, ray_t * in_ray,
				int *dump_flag)
{
    ps_state_t *state = (ps_state_t *) vstate;

    if (*dump_flag) {
	if (state->n_data) {	/* we have not yet dumped our data */
	    double *t;
	    size_t i;
	    const size_t n = state->n_data;
	    double *data = state->data;

	    for (i = 0; i < n; i++) {
		const size_t idx = 3 * i;

		fprintf(state->dump_file, "%g\t%g\t%g\n", data[idx],
			data[idx + 1], data[idx + 2]);
	    }

	    /* shrink memory to minimum (BLOCK) */
	    t = (double *) realloc(&(state->data),
				   3 * BLOCK_SIZE * sizeof(double));
	    state->data = t;
	    state->n_data = 0;
	    state->n_alloc = BLOCK_SIZE;

	} else			/* we have allready dumped our data. mark cycle complete */
	    *dump_flag = 0;

    }

    /* calculate point of interception */

    return NULL;		/* no interception found */
}

static ray_t *ps_get_out_ray(void *vstate, ray_t * in_ray,
			     double *hit, int *dump_flag)
{
    ps_state_t *state = (ps_state_t *) vstate;

    ray_t *out = (ray_t *) malloc(sizeof(ray_t));

    memcpy(state->data, hit, 3);	/* store hit */

    if (state->n_data == state->n_alloc) {	/* inc data size for next hit */
	const size_t n = state->n_data + BLOCK_SIZE;
	double *t;

	t = (double *) realloc(&(state->data), 3 * n * sizeof(double));
	if (t) {		/* success, update state */
	    state->data = t;
	    state->n_data = n;
	} else {		/* memory exhausted, dump data to file and shrink memory to default size */
	    const size_t m = state->n_data;
	    double *data = state->data;

	    size_t i;

	    for (i = 0; i < m; i++) {
		const size_t idx = 3 * i;

		fprintf(state->dump_file, "%g\t%g\t%g\n", data[idx],
			data[idx + 1], data[idx + 2]);
	    }

	    /* shrink memory to minimum (BLOCK) */
	    t = (double *) realloc(&(state->data),
				   3 * BLOCK_SIZE * sizeof(double));
	    state->data = t;
	    state->n_data = 0;
	    state->n_alloc = BLOCK_SIZE;

	    *dump_flag = 1;
	}
    } else
	state->n_data++;

    memcpy(out->origin, hit, 3);	/* update origin */
    memcpy(out->direction, in_ray->direction, 3);	/* copy direction */
    out->power = in_ray->power;

    return out;

}


static const target_type_t ps_t = {
    "plane screen",
    sizeof(struct ps_state_t),
    &ps_alloc_state,
    &ps_init_state,
    &ps_free_state,
    &ps_get_intercept,
    &ps_get_out_ray
};

const target_type_t *target_plane_screen = &ps_t;

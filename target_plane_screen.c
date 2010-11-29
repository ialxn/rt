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
    const char *name;
    double point[3];
    double normal[3];
    double intercept[3];
    int n_alloc;
    int n_data;
    double *data;
} ps_state_t;


static int ps_alloc_state(void *vstate)
{
    ps_state_t *state = (ps_state_t *) vstate;

    if (!(state->data = (double *) malloc(3 * BLOCK * sizeof(double))))
	return ERR;

    state->n_alloc = BLOCK;

    return NO_ERR;
}

static void ps_init_state(void *vstate, config_t * cfg, const char *name)
{
    ps_state_t *state = (ps_state_t *) vstate;

    unsigned int i = 0;
    const char *S;
    double F;

    config_setting_t *this_target, *point, *normal;
    const config_setting_t *targets = config_lookup(cfg, "targets");

    state->name = strdup(name);

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

    free(state->name);
    free(state->data);
}

static vec_t *ps_get_intercept(void *vstate, ray_t * in_ray,
			       int *dump_flag)
{
    return NULL;
}

static ray_t *ps_get_out_ray(void *vstate, ray_t * in_ray,
			     vec_t * intercpt, int *dump_flag)
{
    return NULL;
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

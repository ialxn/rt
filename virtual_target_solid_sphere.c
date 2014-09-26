/*	virtual_target_solid_sphere.c
 *
 * Copyright (C) 2010,2011,2012,2013,2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#include <string.h>

#include "io_utils.h"
#include "reflect.h"
#include "virtual_targets.h"

#define NO_ITEMS 4


typedef struct vtssp_state_t {
    double center[3];		/* center coordinate of sphere */
    double radius;		/* radius^2 of disk */
    pthread_key_t PTDT_key;	/* access to output buffer and flags for each target */
} vtssp_state_t;


static void vtssp_init_state(void *vstate, config_setting_t * this_target,
			     const int file_mode)
{
    vtssp_state_t *state = (vtssp_state_t *) vstate;
    (void) file_mode;

    read_vector(this_target, "origin", state->center);
    config_setting_lookup_float(this_target, "radius", &state->radius);

    pthread_key_create(&state->PTDT_key, free_PTDT);
}


static double *vtssp_get_intercept(void *vstate, ray_t * ray)
{
    vtssp_state_t *state = (vtssp_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);

    if (data->flag & LAST_WAS_HIT) {
	/*
	 * ray starts on this target, no hit possible.
	 * this test fails for ray initially emitted by
	 * this solid source (flag cannot be set by source)
	 */
	data->flag &= ~LAST_WAS_HIT;
	return NULL;
    }

    return intercept_sphere(ray, state->center, state->radius);
}

static ray_t *vtssp_get_out_ray(void *vstate, ray_t * ray, double *hit,
				const gsl_rng * r)
{
    vtssp_state_t *state = (vtssp_state_t *) vstate;
    PTDT_t *data = pthread_getspecific(state->PTDT_key);
    double radius_vector[3];

    /*
     * rays that hit solid sphere are scattered uniformly
     * choose random direction but make sure ray does not enter sphere
     * i.e. 'radius_vector' dot 'ray->dir' is positive. Note that the
     * radius vector is parallel to the normal vector (but not normalized,
     * which is no problem here).
     * 
     * More correct:
     * - define reflectivity of source
     * - non-absorbed rays are scattered as required by reflectivity model
     * - ABSORBED rays transform a power of 'ray->power' to the source.
     *   a number of rays with the same power are re-emitted from random points
     *   on the source into random directions. CAREFUL 'ray->power' of absorbed
     *   is not necessary equal to the power_per_ray of the source that was
     *   hit (several sources involved). bookkeeping required.
     */
    memcpy(ray->orig, hit, 3 * sizeof(double));
    diff(radius_vector, hit, state->center);
    do {
	get_uniform_random_vector(ray->dir, 1.0, r);
    } while (cblas_ddot(3, radius_vector, 1, ray->dir, 1) < 0.0);

    data->flag |= LAST_WAS_HIT;	/* mark as hit */

    return ray;
}

static void vtssp_init_PTDT(void *vstate)
{
    per_thread_init(((vtssp_state_t *) vstate)->PTDT_key, NO_ITEMS);
}

static const target_type_t vt_ssp_t = {
    "",
    sizeof(struct vtssp_state_t),
    &vtssp_init_state,
    NULL,
    &vtssp_get_intercept,
    &vtssp_get_out_ray,
    &vtssp_init_PTDT,
    NULL
};

const target_type_t *virtual_target_solid_sphere = &vt_ssp_t;

/*	targets.h
 *
 * Copyright (C) 2010, 2011 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __TARGETS_H__
#define __TARGETS_H__

#include <pthread.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include <libconfig.h>

#include "ray.h"

#define LAST_WAS_HIT (1<<0)	/* target was hit by last ray */
#define ABSORBED     (1<<1)	/* ray was absorbed on target */


typedef struct target_type_t {
    const char *type;		/* type of target */
    size_t size;		/* internally used to allocate the state (individual,
				   type specific data) of the target. */
    void (*init_state) (void *state, config_setting_t * this_t, config_t * cfg, const char *file_mode);	/* initialize
													   internal data
													   from
													   configuration */
    void (*free_state) (void *state);	/* free */
    double *(*get_intercept) (void *state, ray_t * in_ray);	/* point of intersection */
    ray_t *(*get_out_ray) (void *state, ray_t * in_ray, double *hit,
			   const gsl_rng * r);
    const char *(*get_target_name) (void *state);
    void (*dump_string) (void *state, const char *str);	/* write 'str' to dump file */
    double *(*M) (void *state);	/* returns pointer to M matrix */
    void (*init_flags) (void *state);	/* allocate and zero per thread flags */
} target_type_t;

typedef struct target_t {
    const target_type_t *type;
    void *state;
} target_t;

/*
 * list of all defined targets found in individual files (target_*.c)
 */
const target_type_t *target_plane_screen_one_sided;
const target_type_t *target_plane_screen_two_sided;
const target_type_t *target_rectangle;
const target_type_t *target_triangle;
const target_type_t *target_ellipsoid;
const target_type_t *target_annulus;
const target_type_t *target_disk;
/*
 *  public functions to access/manipulate the targets (found in targets.c)
 */
extern target_t *target_alloc(const target_type_t * type,
			      config_setting_t * this_t, config_t * cfg,
			      const char *file_mode);
extern void target_free(target_t * T);
extern double *interception(const target_t * T, ray_t * in_ray);
extern ray_t *out_ray(const target_t * T, ray_t * in_ray, double *hit,
		      const gsl_rng * r);
extern const char *get_target_type(const target_t * T);
extern const char *get_target_name(const target_t * T);
extern void dump_string(const target_t * T, const char *str);
extern double *M(const target_t * T);
extern void init_flags(const target_t * T);


/*
 * utility functions
 */
extern int check_targets(config_t * cfg);
extern void init_refl_spectrum(const char *f_name, gsl_spline ** spline);

#endif				/* __TARGETS_H__ */

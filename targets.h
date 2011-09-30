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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include <libconfig.h>

#include "ray.h"

#define BLOCK_SIZE 32768
#define MAX_BLOCK_SIZE 262144

typedef struct target_type_t {
    const char *type;		/* type of target */
    size_t size;		/* internally used to allocate the state (individual,
				   type specific data) of the target. */
    int (*alloc_state) (void *state);	/* allocate */
    void (*init_state) (void *state, config_setting_t * this_t, config_t * cfg, const char *file_mode);	/* initialize
													   internal data
													   from
													   configuration */
    void (*free_state) (void *state);	/* free */
    double *(*get_intercept) (void *state, ray_t * in_ray, int *dump_flag);	/* point of intersection */
    ray_t *(*get_out_ray) (void *state, ray_t * in_ray, double *hit,
			   const gsl_rng * r, int *dump_flag,
			   const int n_targets);
    const char *(*get_target_name) (void *state);
    void (*dump_string) (void *state, const char *str);	/* write 'str' to dump file */
    double *(*M) (void *state);	/* returns pointer to M matrix */
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
/*
 *  public functions to access/manipulate the targets (found in targets.c)
 */
extern target_t *target_alloc(const target_type_t * type,
			      config_setting_t * this_t, config_t * cfg,
			      const char *file_mode);
extern void target_free(target_t * T);
extern double *interception(const target_t * T, ray_t * in_ray,
			    int *dump_flag);
extern ray_t *out_ray(const target_t * T, ray_t * in_ray, double *hit,
		      const gsl_rng * r, int *dump_flag,
		      const int n_targets);
extern const char *get_target_type(const target_t * T);
extern const char *get_target_name(const target_t * T);
extern void dump_string(const target_t * T, const char *str);
extern double *M(const target_t * T);


/*
 * utility functions
 */
extern int check_targets(config_t * cfg);
extern void dump_data(FILE * f, double *data, const size_t n_data,
		      const size_t n_items);
extern void shrink_memory(double **data, size_t * n_data,
			  size_t * n_alloc, const size_t N);
extern void try_increase_memory(double **data, size_t * n_data,
				size_t * n_alloc, const size_t N,
				FILE * dump_file, int *dump_flag,
				const int n_targets);
extern void init_refl_spectrum(const char *f_name, gsl_spline ** spline,
			       gsl_interp_accel ** acc);

#endif				/* __TARGETS_H__ */

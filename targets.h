/*	targets.h
 *
 * Copyright (C) 2010,2011,2012,2013,2014,2015,2016 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __TARGETS_H__
#define __TARGETS_H__

#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include <libconfig.h>

#include "ray.h"
#include "math_utils.h"
#include "reflect.h"

/*
 * per thread flags
 */
#define LAST_WAS_HIT (1<<0)	/* target was hit by last ray */
#define ABSORBED     (1<<1)	/* ray was absorbed on target */
#define ICPT_ON_CONVEX_SIDE	(1<<2)	/* ray intercepted by convex side */

/*
 * per target flags
 */
#define KEEP_CLOSED	(1<<0)	/* open/close output when flushing buffer */
#define OUTPUT_REQUIRED	(1<<1)	/* 'no_output' was set */
#define OUTSIDE		(1<<2)	/* mark reflecting surface of non-planar targets */

#define BUF_SIZE 4096

union fh_t {
    int fh;			/* file handle (output kept open) */
    char *fname;		/* file name (output kept close) */
} fh_t;

typedef struct PTDT_t {		/* per thread data of every target */
    char *buf;			/* output buffer */
    size_t i;			/* current position within output buffer */
    int flag;			/* various flags. */
} PTDT_t;

typedef struct target_type_t {
    const char *type;		/* type of target */
    size_t size;		/* internally used to allocate the state
				   (individual, type specific data) of
				   the target. */
    int (*init_state) (void *state, config_setting_t * this_t,
		       const int file_mode, const int keep_closed,
		       const double P_factor);
    /* initialize internal data from configuration */
    void (*free_state) (void *state);	/* free */
    double *(*get_intercept) (void *state, ray_t * ray);
    /* point of intersection */
    ray_t *(*get_out_ray) (void *state, ray_t * ray, double *hit,
			   const gsl_rng * r);
    void (*init_PTDT) (void *state);	/* allocate per thread buffer */
    void (*flush_PTDT_outbuf) (void *state);	/* flush per thread buffer */
} target_type_t;

typedef struct target_t {
    const target_type_t *type;
    void *state;
} target_t;

/*
 * list of all defined targets found in individual files (target_*.c)
 */
const target_type_t *target_annulus;
const target_type_t *target_cone;
const target_type_t *target_cylinder;
const target_type_t *target_disk;
const target_type_t *target_ellipsoid;
const target_type_t *target_paraboloid;
const target_type_t *target_plane_screen_one_sided;
const target_type_t *target_plane_screen_two_sided;
const target_type_t *target_rectangle;
const target_type_t *target_sphere;
const target_type_t *target_triangle;
const target_type_t *target_window;

/*
 * list of all virtual targets.
 * each solid source need to define a corresponding virtual target. 
 */
const target_type_t *virtual_target_solid_cone;
const target_type_t *virtual_target_solid_cylinder;
const target_type_t *virtual_target_solid_sphere;


/*
 *  public functions to access/manipulate the targets (found in targets.c)
 */
extern target_t *target_alloc(const target_type_t * type,
			      config_setting_t * this_t,
			      const int file_mode, const int keep_closed,
			      const double P_factor);
extern void target_free(target_t * T);
extern double *icpt(const target_t * T, ray_t * ray);
extern ray_t *out_ray(const target_t * T, ray_t * ray, double *hit,
		      const gsl_rng * r);
extern void init_PTDT(const target_t * T);
extern void flush_PTDT_outbuf(const target_t * T);
extern void free_PTDT(void *p);


/*
 * utility functions
 */
typedef void (*refl_func_pointer_t) (ray_t * r, const double N[3],
				     const double
				     P[3],
				     const gsl_rng *
				     rng, void *model_params);

extern int check_targets(config_t * cfg);
extern int init_output(const char *target_type,
		       config_setting_t * this_target, const int file_mode,
		       const double P_factor, union fh_t *output,
		       int *out_flag, double point[], double M[]);
extern void init_spectrum(config_setting_t * this_target, const char *kw,
			  gsl_spline ** spectrum);
extern void init_refl_model(const struct config_setting_t *s,
			    refl_func_pointer_t * refl_func,
			    void **refl_model_params);
extern int init_reflecting_surface(config_setting_t * this_target);
extern double *init_M(config_setting_t * this_target, const char *x,
		      const char *z);
extern void per_thread_init(pthread_key_t key, size_t n);
extern void per_thread_flush(union fh_t output, const int out_flag,
			     pthread_key_t key, pthread_mutex_t * mutex);
extern void state_free(union fh_t output, const int out_flag, double *M,
		       gsl_spline * s, refl_func_pointer_t refl_func,
		       void *p);
extern void store_xy(union fh_t output, const int out_flag, ray_t * ray,
		     const double *hit, const double *m,
		     const double *point, PTDT_t * data,
		     pthread_mutex_t * mutex_writefd);
extern void store_xyz(union fh_t output, const int out_flag, ray_t * ray,
		      const double *hit, const double *m,
		      const double *point, PTDT_t * data,
		      pthread_mutex_t * mutex_writefd);

#endif				/* __TARGETS_H__ */

/*	targets.h
 *
 * Copyright (C) 2010,2011,2012,2013,2014 Ivo Alxneit
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

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>

#include <libconfig.h>

#include "ray.h"
#include "math_utils.h"
#include "reflect.h"


#define LAST_WAS_HIT (1<<0)	/* target was hit by last ray */
#define ABSORBED     (1<<1)	/* ray was absorbed on target */
#define ICPT_ON_CONVEX_SIDE	(1<<2)	/* ray intercepted by convex side */

#define OUTSIDE 1

#define BUF_SIZE 4096

typedef struct PTDT_t {		/* per thread data of every target */
    char *buf;			/* output buffer */
    size_t i;			/* current position within output buffer */
    int flag;			/* various flags. */
} PTDT_t;

typedef struct target_type_t {
    const char *type;		/* type of target */
    size_t size;		/* internally used to allocate the state (individual,
				   type specific data) of the target. */
    void (*init_state) (void *state, config_setting_t * this_t, const int file_mode);	/* initialize
											   internal data
											   from
											   configuration */
    void (*free_state) (void *state);	/* free */
    double *(*get_intercept) (void *state, ray_t * ray);	/* point of intersection */
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
const target_type_t *target_plane_screen_one_sided;
const target_type_t *target_plane_screen_two_sided;
const target_type_t *target_rectangle;
const target_type_t *target_triangle;
const target_type_t *target_ellipsoid;
const target_type_t *target_annulus;
const target_type_t *target_disk;
const target_type_t *target_cylinder;
const target_type_t *target_window;
const target_type_t *target_paraboloid;
const target_type_t *target_sphere;
/*
 *  public functions to access/manipulate the targets (found in targets.c)
 */
extern target_t *target_alloc(const target_type_t * type,
			      config_setting_t * this_t,
			      const int file_mode);
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
extern int check_targets(config_t * cfg);
extern int init_output(const int file_mode, const char *target_type,
		       config_setting_t * this_target, double point[],
		       double M[]);
extern void init_spectrum(const char *f_name, gsl_spline ** refl_spectrum);
extern void init_refl_model(const struct config_setting_t *s, char *model,
			    void **refl_model_params);
extern char init_reflecting_surface(config_setting_t * this_target);
extern double *init_M(config_setting_t * this_target, const char *x,
		      const char *z);
extern void per_thread_init(pthread_key_t key, size_t n);
extern void per_thread_flush(int fh, pthread_key_t key,
			     pthread_mutex_t * mutex);
extern void state_free(int fh, double *M, gsl_spline * s, char model,
		       void *p);
extern void store_xy(const int fd, ray_t * ray, const double *hit,
		     const double *m, const double *point, PTDT_t * data,
		     pthread_mutex_t * mutex_writefd);
extern void store_xyz(const int fd, ray_t * ray, const double *hit,
		      const double *m, const double *point, PTDT_t * data,
		      pthread_mutex_t * mutex_writefd);

#endif				/* __TARGETS_H__ */

/*	sources.h
 *
 * Copyright (C) 2010,2011,2012,2013,2014,2015 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __SOURCES_H__
#define __SOURCES_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_cblas.h>

#include <libconfig.h>
#include <pthread.h>

#include "ray.h"

#define RAYS_PER_GROUP 1000


typedef struct source_type_t {
    const char *type;		/* type of source */
    size_t size;		/* internally used to allocate the state
				   (individual, type specific data) of
				   the source. */
    void (*init_state) (void *state, config_setting_t * this_s,
			const double P_factor);
    /* initialize internal data from configuration */
    void (*free_state) (void *state);
    /* free */
    ray_t *(*emit_ray) (void *state, const gsl_rng * r);
    /* returns a new ray, or NULL if exhausted */
    const char *(*get_source_name) (void *state);
    /* get name of source */
     int64_t(*get_source_n_rays) (void *state);
    /* get number of rays of source */
    double (*get_source_power) (void *state);
    /* get power of source */
    void (*init_rays_remain) (void *state);
    /* init PTD variable */
} source_type_t;

typedef struct source_t {
    const source_type_t *type;
    void *state;
} source_t;

/*
 * list of all defined sources found in individual files (source_*.c)
 */
const source_type_t *source_solid_rod;
const source_type_t *source_solid_sphere;
const source_type_t *source_sphere;
const source_type_t *source_spot;
const source_type_t *source_uniform_point_source;

/*
 *  public functions to access / manipulate the sources (found in sources.c)
 */
extern source_t *source_alloc(const source_type_t * T,
			      config_setting_t * this_s,
			      const double P_factor);
extern void source_free(source_t * S);
extern ray_t *emit_ray(const source_t * S, const gsl_rng * r);
extern const char *get_source_type(const source_t * S);
extern const char *get_source_name(const source_t * S);
extern int64_t get_source_n_rays(const source_t * S);
extern double get_source_power(const source_t * S);
extern void init_rays_remain(const source_t * S);

/*
 * utility functions
 */
extern int check_sources(config_t * cfg);
extern void init_source_spectrum(config_setting_t * this_s, const char *kw,
				 gsl_spline ** spectrum);
extern void per_thread_init_rays_remain(pthread_key_t key);
extern int64_t per_thread_get_source_n_rays(pthread_mutex_t * mutex,
					    int64_t * n_rays);
extern int64_t per_thread_get_new_raygroup(pthread_mutex_t * mutex,
					   int64_t * n_rays);

#endif				/* __SOURCES_H__ */

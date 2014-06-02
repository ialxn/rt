/*	sources.h
 *
 * Copyright (C) 2010,2011,2012,2013,2014 Ivo Alxneit
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
    size_t size;		/* internally used to allocate the state (individual,
				   type specific data) of the source. */
    void (*init_state) (void *state, config_setting_t * this_s);	/* initialize internal data
									   from configuration */
    void (*free_state) (void *state);	/* free */
    ray_t *(*emit_ray) (void *state, const gsl_rng * r);	/* returns a new ray, or NULL if exhausted */
    const char *(*get_source_name) (void *state);	/* get name of source */
    void (*init_rays_remain) (void *state);	/* init PTD variable */
} source_type_t;

typedef struct source_t {
    const source_type_t *type;
    void *state;
} source_t;

/*
 * list of all defined sources found in individual files (source_*.c)
 */
const source_type_t *source_uniform_point_source;
const source_type_t *source_spot;
const source_type_t *source_sphere;
const source_type_t *source_solid_sphere;

/*
 *  public functions to access / manipulate the sources (found in sources.c)
 */
extern source_t *source_alloc(const source_type_t * T,
			      config_setting_t * this_s);
extern void source_free(source_t * S);
extern ray_t *emit_ray(const source_t * S, const gsl_rng * r);
extern const char *get_source_type(const source_t * S);
extern const char *get_source_name(const source_t * S);
extern void init_rays_remain(const source_t * S);

/*
 * utility functions
 */
extern int check_sources(config_t * cfg);
extern void init_spectrum(const char *f_name, gsl_spline ** spectrum);

#endif				/* __SOURCES_H__ */

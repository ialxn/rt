/*	sources.h
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __SOURCES_H__
#define __SOURCES_H__

#include <gsl/gsl_rng.h>
#include <libconfig.h>

#define NO_ERR 0
#define ERR 1

typedef struct source_type_t {
    const char *name;		/* name of the source */
    const char *type;		/* type of the source */
    size_t size;		/* internally used to allocate the state (individual,
				   type specific data) of the source. */
    int (*alloc_state) (void *state);	/* allocate */
    int (*init_state) (void *state, config_t * cfg, const char *name);	/* initialize internal data from                                                                                   configuration */
    void (*free_state) (void *state);	/* free */
    int (*get_new_ray) (void *state, const gsl_rng * r);	/* returns a new ray, or NULL if exhausted */
} source_type_t;




typedef struct source_t {
    const source_type_t *type;
    void *state;
} source_t;




/*
 * list of all define sources found in individual files (source_*.c)
 */
const struct source_type_t *source_uniform_point_source;

/*
 * functions to access/manipulate the sources (found in sources.c)
 */
extern source_t *source_alloc(const source_type_t * T, config_t * cfg,
			      const char *name);
extern void source_free(source_t * S);
extern int new_ray(const source_t * S, const gsl_rng * r);

extern int check_sources(config_t * cfg);

#endif				/* __SOURCES_H__ */

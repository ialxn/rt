/*	targets.h
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __TARGETS_H__
#define __TARGETS_H__

#include <libconfig.h>

#include "ray.h"

#define NO_ERR 0
#define ERR 1
//#define BLOCK_SIZE 32768
//#define MAX_BLOCK_SIZE 262144
#define BLOCK_SIZE 2
#define MAX_BLOCK_SIZE 4

typedef struct target_type_t {
    const char *type;		/* type of target */
    size_t size;		/* internally used to allocate the state (individual,
				   type specific data) of the target. */
    int (*alloc_state) (void *state);	/* allocate */
    void (*init_state) (void *state, config_t * cfg, const char *name);	/* initialize internal data
									   from configuration */
    void (*free_state) (void *state);	/* free */
    double *(*get_intercept) (void *state, ray_t * in_ray, int *dump_flag);
    ray_t *(*get_out_ray) (void *state, ray_t * in_ray, double *hit,
			   int *dump_flag);
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
/*
 *  publicfunctions to access/manipulate the sources (found in sources.c)
 */
extern target_t *target_alloc(const target_type_t * type, config_t * cfg,
			      const char *name);
extern void target_free(target_t * T);
extern double *interception(const target_t * T, ray_t * in_ray,
			    int *dump_flag);
extern ray_t *out_ray(const target_t * T, ray_t * in_ray,
		      double *hit, int *dump_flag);

/*
 * utility functions
 */
extern int check_targets(config_t * cfg);
extern void dump_data(FILE * f, double *data, const size_t n_data,
		      const size_t n_items);

#endif				/* __TARGETS_H__ */

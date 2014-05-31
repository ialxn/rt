/*	obj_lists.h
 *
 * Copyright (C) 2010,2011,2012,2013,2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __OBJ_LISTS_H__
#define __OBJ_LISTS_H__

#include <libconfig.h>

#include "list.h"
#include "sources.h"
#include "targets.h"
#include "virtual_targets.h"

typedef struct source_list_t {
    struct list_head list;
    source_t *s;
} source_list_t;

typedef struct target_list_t {
    struct list_head list;
    target_t *t;
} target_list_t;

extern source_list_t *init_sources(config_t * cfg, int *n_sources);
extern target_list_t *init_targets(config_t * cfg, int *n_targets,
				   const int file_mode);

extern void source_list_free(source_list_t * s);
extern void target_list_free(target_list_t * t);

#endif				/* __OBJ_LISTS_H__ */

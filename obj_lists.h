/*	obj_lists.h
 *
 * Copyright (C) 2010 - 2018 Ivo Alxneit, Paul Scherrer Institute
 *
 * This file is part of rt
 *
 * rt is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rt is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rt. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef __OBJ_LISTS_H__
#define __OBJ_LISTS_H__

#include <libconfig.h>

#include "list.h"
#include "sources.h"
#include "targets.h"

typedef struct source_list_t {
    struct list_head list;
    source_t *s;
} source_list_t;

typedef struct target_list_t {
    struct list_head list;
    target_t *t;
} target_list_t;

extern source_list_t *init_sources(config_t * cfg, int *n_sources,
				   const double P_factor);
extern target_list_t *init_targets(config_t * cfg, int *n_targets,
				   const int file_mode,
				   const int keep_closed,
				   const double P_factor);

extern void source_list_free(source_list_t * s);
extern void target_list_free(target_list_t * t);

#endif				/* __OBJ_LISTS_H__ */

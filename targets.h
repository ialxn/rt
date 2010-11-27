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

#define NO_ERR 0
#define ERR 1

extern int check_targets(config_t * cfg);

#endif				/* __TARGETS_H__ */

/*	virtual_targets.h
 *
 * Copyright (C) 2010,2011,2012,2013,2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __VIRTUAL_TARGETS_H__
#define __VIRTUAL_TARGETS_H__

#include "targets.h"

/*
 * list of all virtual targets (see individual files 'virtual_target_*.c')
 * each solid source need to define a corresponding virtual target. 
 */
const target_type_t *virtual_target_solid_sphere;


#endif				/* __VIRTUAL_TARGETS_H__ */

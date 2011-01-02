/*	io_lists.h
 *
 * Copyright (C) 2011 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __IO_LISTS_H__
#define __IO_LISTS_H__

#include <libconfig.h>

#define NO_ERR 0
#define ERR 1


extern int check_array(const char *section, const config_setting_t * s,
		       const char *name, const int nr);
extern int check_string(const char *section, const config_setting_t * s,
			const char *name, const int nr);
extern int check_float(const char *section, const config_setting_t * s,
		       const char *name, const int nr);


#endif				/* __IO_LISTS_H__ */

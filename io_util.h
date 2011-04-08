/*	io_util.h
 *
 * Copyright (C) 2011 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __IO_UTIL_H__
#define __IO_UTIL_H__

#include <libconfig.h>

#define NO_ERR 0
#define ERR 1

#define LINE_LEN 256
#define BSIZE 256


extern int check_array(const char *section, const config_setting_t * s,
		       const char *name, const int nr);
extern int check_string(const char *section, const config_setting_t * s,
			const char *name, const int nr);
extern int check_return_string(const char *section,
			       const config_setting_t * s,
			       const char *name, const int nr,
			       const char **string);
extern int check_float(const char *section, const config_setting_t * s,
		       const char *name, const int nr);
extern int check_int(const char *section, const config_setting_t * s,
		     const char *name, const int nr);
extern int check_file(const char *section, const config_setting_t * s,
		      const char *name, const int nr);
extern void read_vector(const config_setting_t * s, const char *name,
			double *const vec);
extern void read_vector_normalize(const config_setting_t * s,
				  const char *name, double *const vec);
extern int read_data(FILE * f, double **x, double **y, size_t * n);


#endif				/* __IO_UTIL_H__ */

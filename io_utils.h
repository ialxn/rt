/*	io_utils.h
 *
 * Copyright (C) 2011,2012,2013,2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __IO_UTILS_H__
#define __IO_UTILS_H__

#include <libconfig.h>

#define NO_ERR 0
#define ERR 1

#define LINE_LEN 256
#define BSIZE 256
#define T_HEADER_LINES 10	/* number of header lines in target output */
#define MAX_ITEMS 5		/* maximum number of items per data set
				   (x,y,z,power,lambda) */


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
extern int skip_header(FILE * f);
extern double *range(const double min, const double max, const size_t n);
extern int get_idx(FILE * f_in, size_t * idx_lambda, size_t * idx_power);


#endif				/* __IO_UTIL_H__ */

/*	io_utils.h
 *
 * Copyright (C) 2011,2012,2013,2014,2015,2016,2017 Ivo Alxneit
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

#define LOCAL 0
#define GLOBAL 1

#define LINE_LEN 256
#define BSIZE 256
#define HEADER_LINES 16		/* number of header lines after 'get_id()' */
#define REST_HEADER_LINES 14	/* number of header lines after 'P_factor' */
#define MAX_FLOAT_ITEMS 4	/* maximum number of items per data set
				   (x,y,z,lambda) */


extern int check_array(const char *section, const config_setting_t * s,
		       const char *name, const int nr);
extern int check_string(const char *section, const config_setting_t * s,
			const char *name, const int nr);
extern int check_return_string(const char *section,
			       const config_setting_t * s,
			       const char *name, const int nr,
			       const char **string);
extern int check_bool(const char *section, const config_setting_t * s,
		      const char *name, const int nr);
extern int check_float(const char *section, const config_setting_t * s,
		       const char *name, const int nr);
extern int check_int(const char *section, const config_setting_t * s,
		     const char *name, const int nr);
extern int check_reflectivity_model(const char *section,
				    const config_setting_t * s,
				    const char *name, const int nr);
extern int check_file(const char *section, const config_setting_t * s,
		      const char *name, const int nr);
extern void read_vector(const config_setting_t * s, const char *name,
			double *const vec);
extern void read_vector_normalize(const config_setting_t * s,
				  const char *name, double *const vec);
extern void init_M_from_z(config_setting_t * this, const char *kw,
			  double M[9]);
extern double get_A_disk(config_setting_t * this_s, const char *kw,
			 const double radius);
extern int read_data(FILE * f, double **x, double **y, size_t * n);
extern double get_P_factor(FILE * f);
extern int skip_N_comments(FILE * f, const int N);
extern void read_transformation(FILE * f_in, double M[], double origin[]);
extern double *range(const double min, const double max, const size_t n);
extern int get_idx(FILE * f_in, size_t * idx_lambda);


#endif				/* __IO_UTIL_H__ */

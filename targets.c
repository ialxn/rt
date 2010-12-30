/*	targets.c
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_cblas.h>

#include "targets.h"

target_t *target_alloc(const target_type_t * type, config_t * cfg,
		       const char *name, const char *file_mode)
{
    target_t *T;

    /* generic part */
    if ((T = (target_t *) malloc(sizeof(target_t))) == NULL) {
	fprintf(stderr, "failed to allocate space for 'target_t'\n");
	return NULL;
    }

    if ((T->state = malloc(type->size)) == NULL) {
	free(T);
	fprintf(stderr,
		"failed to allocate space for generic state of target\n");
	return NULL;
    };

    T->type = type;

    /* specific part for target 'type' */
    if ((T->type->alloc_state) (T->state) == ERR) {
	fprintf(stderr,
		"failed to allocate space for specific internal state of target\n");
	free(T->state);
	free(T);
	return NULL;
    }

    (T->type->init_state) (T->state, cfg, name, file_mode);	/* initialize data structures */

    return T;
}

void target_free(target_t * T)
{
    (T->type->free_state) (T->state);
    free(T->state);
    free(T);
}

double *interception(const target_t * T, ray_t * in_ray, int *dump_flag)
{
    return (T->type->get_intercept) (T->state, in_ray, dump_flag);
}

ray_t *out_ray(const target_t * T, ray_t * in_ray, const double ppr,
	       double *hit, int *dump_flag, const int n_targets)
{
    return (T->type->get_out_ray) (T->state, in_ray, ppr, hit,
				   dump_flag, n_targets);
}

const char *get_target_type(const target_t * T)
{
    return T->type->type;
}

const char *get_target_name(const target_t * T)
{
    return (T->type->get_target_name) (T->state);
}

void dump_string(const target_t * T, const char *str)
{
    (T->type->dump_string) (T->state, str);
}

double *M(const target_t * T)
{
    return (T->type->M) (T->state);
}


int check_targets(config_t * cfg)
{
    int status = NO_ERR;
    const config_setting_t *t = config_lookup(cfg, "targets");

    if (t == NULL) {
	fprintf(stderr, "missing 'targets' keyword\n");
	status = ERR;
    } else {			/* 'targets' section present */
	const int n_targets = config_setting_length(t);
	int i;

	if (n_targets == 0) {
	    fprintf(stderr, "empty 'targets' section\n");
	    status = ERR;
	}

	for (i = 0; i < n_targets; ++i) {
	    config_setting_t *this_t =
		config_setting_get_elem(t, (unsigned int) i);

	    const char *S, *type = NULL;

	    /*
	     * keywords common to all targets
	     * 'name':  identifier / string
	     * 'type':  type of source / string
	     *          - "one-sided plane_screen": non-absorbing counter plane
	     *                                      only rays intersecting anti-parallel
	     *                                      to plane's normal vector are counted
	     *          - "two-sided plane_screen": non-absorbing counter plane
	     *                                      rays intersecting from both sides
	     *                                      are counted
	     *          - "square":                 specular reflection on front surface
	     *                                      defined by direction of surface
	     *                                      normal. total absorption on rear
	     *                                      surface.
	     *          - "triangle":               specular reflection on front surface
	     *                                      defined by direction of surface
	     *                                      normal. total absorption on rear
	     *                                      surface.
	     */
	    if (config_setting_lookup_string(this_t, "name", &S) !=
		CONFIG_TRUE) {
		fprintf(stderr,
			"missing 'name' keyword in 'targets' section %u\n",
			i + 1);
		status = ERR;
	    }
	    if (config_setting_lookup_string(this_t, "type", &type) !=
		CONFIG_TRUE) {
		fprintf(stderr,
			"missing 'type' keyword in 'targets' section %u\n",
			i + 1);
		status = ERR;
	    }

	    /* check target specific settings */

	    if (!strcmp(type, "one-sided plane screen")) {
		/*
		 * one-sided plane screen:
		 *  - array 'point' (point on plane) [x,y,z] / double
		 *  - array 'normal' (normal vector of plane) [x,y,z] / double
		 *  - array 'x' (direction of local x-axis of plane) [x,y,z] / double
		 */
		config_setting_t *point, *normal, *x;

		if ((point =
		     config_setting_get_member(this_t, "point")) == NULL) {
		    fprintf(stderr,
			    "missing 'point' array in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else if (config_setting_is_array(point) == CONFIG_FALSE
			   || config_setting_length(point) != 3) {
		    fprintf(stderr,
			    "setting 'point' in 'targets' section %u is not array with 3 coordinates\n",
			    i + 1);
		    status = ERR;
		}
		/* end keyword 'point' found */
		if ((normal =
		     config_setting_get_member(this_t,
					       "normal")) == NULL) {
		    fprintf(stderr,
			    "missing 'normal' array in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else if (config_setting_is_array(normal) == CONFIG_FALSE
			   || config_setting_length(normal) != 3) {
		    fprintf(stderr,
			    "setting 'normal' in 'targets' section %u is not array with 3 coordinates\n",
			    i + 1);
		    status = ERR;
		}		/* end keyword 'normal' found */
		if ((x = config_setting_get_member(this_t, "x")) == NULL) {
		    fprintf(stderr,
			    "missing 'x' array in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else if (config_setting_is_array(normal) == CONFIG_FALSE
			   || config_setting_length(normal) != 3) {
		    fprintf(stderr,
			    "setting 'x' in 'targets' section %u is not array with 3 coordinates\n",
			    i + 1);
		    status = ERR;
		}		/* end keyword 'x' found */
	    }
	    /* end 'one-sided plane_screen' */
	    else if (!strcmp(type, "two-sided plane screen")) {
		/*
		 * two-sided plane screen:
		 *  - array 'point' (point on plane) [x,y,z] / double
		 *  - array 'normal' (normal vector of plane) [x,y,z] / double
		 *  - array 'x' (direction of local x-axis of plane) [x,y,z] / double
		 */
		config_setting_t *point, *normal, *x;

		if ((point =
		     config_setting_get_member(this_t, "point")) == NULL) {
		    fprintf(stderr,
			    "missing 'point' array in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else if (config_setting_is_array(point) == CONFIG_FALSE
			   || config_setting_length(point) != 3) {
		    fprintf(stderr,
			    "setting 'point' in 'targets' section %u is not array with 3 coordinates\n",
			    i + 1);
		    status = ERR;
		}
		/* end keyword 'point' found */
		if ((normal =
		     config_setting_get_member(this_t,
					       "normal")) == NULL) {
		    fprintf(stderr,
			    "missing 'normal' array in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else if (config_setting_is_array(normal) == CONFIG_FALSE
			   || config_setting_length(normal) != 3) {
		    fprintf(stderr,
			    "setting 'normal' in 'targets' section %u is not array with 3 coordinates\n",
			    i + 1);
		    status = ERR;
		}		/* end keyword 'normal' found */
		if ((x = config_setting_get_member(this_t, "x")) == NULL) {
		    fprintf(stderr,
			    "missing 'x' array in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else if (config_setting_is_array(normal) == CONFIG_FALSE
			   || config_setting_length(normal) != 3) {
		    fprintf(stderr,
			    "setting 'x' in 'targets' section %u is not array with 3 coordinates\n",
			    i + 1);
		    status = ERR;
		}		/* end keyword 'x' found */
	    } /* end 'two-sided plane_screen' */
	    else if (!strcmp(type, "square")) {
		/*
		 * square
		 *  - array 'point' (corner point of square) [x,y,z] / double
		 *  - array 'x' (direction of local x-axis of plane) [x,y,z] / double
		 *  - array 'y' (direction of local y-axis of plane) [x,y,z] / double
		 *
		 * NOTE: normal of plane is defined by 'X' cross 'Y'. only rays
		 *       anti-parallel to normal are reflected. ray impiging
		 *       parallel to square hit its back side and are absorbed
		 */
		config_setting_t *point, *X, *Y;

		if ((point =
		     config_setting_get_member(this_t, "point")) == NULL) {
		    fprintf(stderr,
			    "missing 'point' array in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else if (config_setting_is_array(point) == CONFIG_FALSE
			   || config_setting_length(point) != 3) {
		    fprintf(stderr,
			    "setting 'point' in 'targets' section %u is not array with 3 coordinates\n",
			    i + 1);
		    status = ERR;
		}
		/* end keyword 'point' found */
		if ((X = config_setting_get_member(this_t, "x")) == NULL) {
		    fprintf(stderr,
			    "missing 'x' array in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else if (config_setting_is_array(X) == CONFIG_FALSE
			   || config_setting_length(X) != 3) {
		    fprintf(stderr,
			    "setting 'x' in 'targets' section %u is not array with 3 coordinates\n",
			    i + 1);
		    status = ERR;
		}		/* end keyword 'x' found */
		if ((Y = config_setting_get_member(this_t, "y")) == NULL) {
		    fprintf(stderr,
			    "missing 'y' array in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else if (config_setting_is_array(Y) == CONFIG_FALSE
			   || config_setting_length(Y) != 3) {
		    fprintf(stderr,
			    "setting 'y' in 'targets' section %u is not array with 3 coordinates\n",
			    i + 1);
		    status = ERR;
		}		/* end keyword 'y' found */
	    }			/* end 'square' */
	    else if (!strcmp(type, "triangle")) {
		/*
		 * triangle
		 *  - array 'P1' (corner point of triangle) [x,y,z] / double
		 *  - array 'P2' (corner point of triangle) [x,y,z] / double
		 *  - array 'P3' (corner point of triangle) [x,y,z] / double
		 *
		 * NOTE: normal of triangle is defined by ('P2'-'P1') cross ('P3'-'P1').
		 *       only rays anti-parallel to normal are reflected. ray impiging
		 *       parallel to the triangle hit its back side and are absorbed
		 */
		config_setting_t *P1, *P2, *P3;

		if ((P1 = config_setting_get_member(this_t, "P1")) == NULL) {
		    fprintf(stderr,
			    "missing 'P1' array in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else if (config_setting_is_array(P1) == CONFIG_FALSE
			   || config_setting_length(P1) != 3) {
		    fprintf(stderr,
			    "setting 'P1' in 'targets' section %u is not array with 3 coordinates\n",
			    i + 1);
		    status = ERR;
		}
		/* end keyword 'P1' found */
		if ((P2 = config_setting_get_member(this_t, "P2")) == NULL) {
		    fprintf(stderr,
			    "missing 'P2' array in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else if (config_setting_is_array(P2) == CONFIG_FALSE
			   || config_setting_length(P2) != 3) {
		    fprintf(stderr,
			    "setting 'P2' in 'targets' section %u is not array with 3 coordinates\n",
			    i + 1);
		    status = ERR;
		}
		/* end keyword 'P2' found */
		if ((P3 = config_setting_get_member(this_t, "P3")) == NULL) {
		    fprintf(stderr,
			    "missing 'P3' array in 'targets' section %u\n",
			    i + 1);
		    status = ERR;
		} else if (config_setting_is_array(P3) == CONFIG_FALSE
			   || config_setting_length(P3) != 3) {
		    fprintf(stderr,
			    "setting 'P3' in 'targets' section %u is not array with 3 coordinates\n",
			    i + 1);
		    status = ERR;
		}
		/* end keyword 'P3' found */
	    }			/* end 'triangle' */
	}			/* end 'this_t', check next target */
    }				/* end 'targets' section present */
    return status;
}

void dump_data(FILE * f, double *data, const size_t n_data,
	       const size_t n_items)
{
    /*
     * writes 'n_data' lines with 'n_items' items per line
     * to file 'f' (tab separated).
     */
    size_t i, j;

    for (i = 0; i < n_data; i++) {
	const size_t N = i * n_items;

	for (j = 0; j < n_items; j++) {
	    fprintf(f, "%g", data[N + j]);
	    if (j == n_items - 1)
		fputc('\n', f);
	    else
		fputc('\t', f);
	}
    }
}

void shrink_memory(double **data, size_t * n_data, size_t * n_alloc,
		   const size_t N)
{
    /*
     * shrink memory to minimum (BLOCK_SIZE)
     * 'N' doubles per data set
     */
    double *t = (double *) realloc(*data, N * BLOCK_SIZE * sizeof(double));

    *data = t;
    *n_data = 0;
    *n_alloc = BLOCK_SIZE;
}

void try_increase_memory(double **data, size_t * n_data, size_t * n_alloc,
			 const size_t N, FILE * dump_file, int *dump_flag,
			 const int n_targets)
{
    /*
     * we try to increase the size of the '**data' buffer by
     * BLOCK_SIZE*'N' ('N' items per data set). if
     * memory can not be allocated of if the buffer already
     * has reached the maximum size MAX_BLOCK_SIZE*'N', '**data'
     * is written to the 'dump_file' and the size of the buffer
     * is decreased to BLOCK_SIZE*'N'. this initiates a dump cycle
     * as all targets will dump their data and decrease their
     * buffer too during the next call of 'interception()' in
     * 'run_simulation()'
     */
    if (*n_alloc == MAX_BLOCK_SIZE) {	/* max size reached, initiate dump cycle */
	dump_data(dump_file, *data, *n_data, N);
	shrink_memory(data, n_data, n_alloc, N);

	*dump_flag = n_targets - 1;
    } else {			/* try to increase buffer */
	const size_t n = *n_data + BLOCK_SIZE;
	double *t = (double *) realloc(*data, N * n * sizeof(double));
	if (t) {		/* success, update state */
	    *data = t;
	    *n_alloc = n;
	} else {		/* memory exhausted, initiate dump cycle */
	    dump_data(dump_file, *data, *n_data, N);
	    shrink_memory(data, n_data, n_alloc, N);

	    *dump_flag = n_targets - 1;
	}
    }
}

void cross_product(const double a[3], const double b[3], double result[3])
{
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

void g2l(const double *mat, const double *origin, const double *g,
	 double *l)
/*
 * expresses vector 'vec' (global) in local coordinates
 *     l(x, y, z) = MT (g(x, y, z) - o(x, y, z))
 */
{
    int i;
    double t[3];

    for (i = 0; i < 3; i++)
	t[i] = g[i] - origin[i];

    for (i = 0; i < 3; i++)
	l[i] = cblas_ddot(3, t, 1, &mat[3 * i], 1);
}

void l2g(const double *mat, const double *origin, const double *l,
	 double *g)
/*
 * expresses vector 'vec' (local) in global coordinates
 *     g(x, y, z) = M l(x, y, z) + o(x, y, z)
 */
{
    int i;

    for (i = 0; i < 3; i++)
	g[i] = cblas_ddot(3, &mat[i], 3, l, 1) + origin[i];

}

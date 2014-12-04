/*	targets.c
 *
 * Copyright (C) 2010,2011,2012,2013,2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_cblas.h>

#include "io_utils.h"
#include "targets.h"

target_t *target_alloc(const target_type_t * type,
		       config_setting_t * this_t, const int file_mode)
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

    (T->type->init_state) (T->state, this_t, file_mode);	/* initialize data structures */

    return T;
}

void target_free(target_t * T)
{
    if (T->type->free_state)
	(T->type->free_state) (T->state);
    free(T->state);
    free(T);
}

double *icpt(const target_t * T, ray_t * ray)
{
    return (T->type->get_intercept) (T->state, ray);
}

ray_t *out_ray(const target_t * T, ray_t * ray, double *hit,
	       const gsl_rng * r)
{
    return (T->type->get_out_ray) (T->state, ray, hit, r);
}

void init_PTDT(const target_t * T)
{
    (T->type->init_PTDT) (T->state);
}

void flush_PTDT_outbuf(const target_t * T)
{
    if (T->type->flush_PTDT_outbuf)
	(T->type->flush_PTDT_outbuf) (T->state);
}

void free_PTDT(void *p)
{
    PTDT_t *data = (PTDT_t *) p;

    if (data->buf)
	free(data->buf);
    free(data);
}


int check_targets(config_t * cfg)
{
    int status = NO_ERR;
    const config_setting_t *t = config_lookup(cfg, "targets");

    if (t == NULL) {
	fprintf(stderr, "missing 'targets' keyword\n");
	status += ERR;
    } else {			/* 'targets' section present */
	const int n_targets = config_setting_length(t);
	int i;

	if (n_targets == 0) {
	    fprintf(stderr, "empty 'targets' section\n");
	    status += ERR;
	}

	for (i = 0; i < n_targets; ++i) {
	    config_setting_t *this_t =
		config_setting_get_elem(t, (unsigned int) i);

	    const char *type = NULL;

	    /*
	     * keywords common to all targets
	     * 'name':  identifier / string
	     * 'type':  type of target / string
	     *          - "one-sided plane_screen": non-absorbing counter plane
	     *                                      only rays intersecting anti-parallel
	     *                                      to plane's normal vector are counted
	     *          - "two-sided plane_screen": non-absorbing counter plane
	     *                                      rays intersecting from both sides
	     *                                      are counted
	     *          - "rectangle":              reflection on front surface
	     *                                      defined by direction of surface
	     *                                      normal. total absorption on rear
	     *                                      surface.
	     *          - "triangle":               reflection on front surface
	     *                                      defined by direction of surface
	     *                                      normal. total absorption on rear
	     *                                      surface.
	     *          - "ellipsoid"             : ellipsoid. reflection on
	     *                                      inside surface, total absorption
	     *                                      on outside surface.
	     *          - "disk":                   reflection on front surface
	     *                                      defined by direction of surface
	     *                                      normal. total absorption on rear
	     *                                      surface.
	     *          - "annulus":                reflection on front surface
	     *                                      defined by direction of surface
	     *                                      normal. total absorption on rear
	     *                                      surface.
	     *          - "window":                 round window. reflectivity defined
	     *                                      by index of refraction. total
	     *                                      absorbtion on cylinder walls
	     *                                      (inside and outside). absorption
	     *                                      inside governed by absorption
	     *                                      spectrum.
	     */
	    status += check_string("targets", this_t, "name", i);

	    status +=
		check_return_string("targets", this_t, "type", i, &type);
	    if (!type)
		continue;

	    /* check target specific settings */

	    if (!strcmp(type, "one-sided plane screen")) {
		/*
		 * one-sided plane screen:
		 *  - array 'point' (point on plane) [x,y,z] / double
		 *  - array 'normal' (normal vector of plane) [x,y,z] / double
		 *  - array 'x' (direction of local x-axis of plane) [x,y,z] / double
		 */
		status += check_array("targets", this_t, "point", i);
		status += check_array("targets", this_t, "normal", i);
		status += check_array("targets", this_t, "x", i);

	    } /* end 'one-sided plane_screen' */
	    else if (!strcmp(type, "two-sided plane screen")) {
		/*
		 * two-sided plane screen:
		 *  - array 'point' (point on plane) [x,y,z] / double
		 *  - array 'normal' (normal vector of plane) [x,y,z] / double
		 *  - array 'x' (direction of local x-axis of plane) [x,y,z] / double
		 */
		status += check_array("targets", this_t, "point", i);
		status += check_array("targets", this_t, "normal", i);
		status += check_array("targets", this_t, "x", i);

	    } /* end 'two-sided plane_screen' */
	    else if (!strcmp(type, "rectangle")) {
		/*
		 * rectangle
		 *  - array 'P1' (corner point of rectangle) [x,y,z] / double
		 *  - array 'P2' (corner point of rectangle) [x,y,z] / double
		 *               direction of local x-axis of plane: 'P2'-'P1'
		 *  - array 'P3' (corner point of rectangle) [x,y,z] / double
		 *               direction of local y-axis of plane: 'P2'-'P1'
		 *
		 * NOTE: normal of plane is defined by 'X' cross 'Y'. only rays
		 *       anti-parallel to normal are reflected. ray impinging
		 *       parallel to rectangle hit its back side and are absorbed
		 */
		status += check_array("targets", this_t, "P1", i);
		status += check_array("targets", this_t, "P2", i);
		status += check_array("targets", this_t, "P3", i);
		status +=
		    check_string("targets", this_t, "reflectivity", i);
		status += check_file("targets", this_t, "reflectivity", i);
		status +=
		    check_string("targets", this_t, "reflectivity_model",
				 i);
		status +=
		    check_reflectivity_model("targets", this_t,
					     "reflectivity_model", i);

	    } /* end 'rectangle' */
	    else if (!strcmp(type, "triangle")) {
		/*
		 * triangle
		 *  - array 'P1' (corner point of triangle) [x,y,z] / double
		 *  - array 'P2' (corner point of triangle) [x,y,z] / double
		 *  - array 'P3' (corner point of triangle) [x,y,z] / double
		 *
		 * NOTE: normal of triangle is defined by ('P2'-'P1') cross ('P3'-'P1').
		 *       only rays anti-parallel to normal are reflected. ray impinging
		 *       parallel to the triangle hit its back side and are absorbed
		 */
		status += check_array("targets", this_t, "P1", i);
		status += check_array("targets", this_t, "P2", i);
		status += check_array("targets", this_t, "P3", i);
		status +=
		    check_string("targets", this_t, "reflectivity", i);
		status += check_file("targets", this_t, "reflectivity", i);
		status +=
		    check_string("targets", this_t, "reflectivity_model",
				 i);
		status +=
		    check_reflectivity_model("targets", this_t,
					     "reflectivity_model", i);

	    } /* end 'triangle' */
	    else if (!strcmp(type, "ellipsoid")) {
		/*
		 * ellipsoid:
		 *  - array 'center' (center of ellipsoid) [x,y,z] / double
		 *  - array 'x-axis' (direction of local x-axis) [x,y,z] / double
		 *  - array 'z-axis' (direction of local z-axis) [x,y,z] / double
		 *  - array 'axes' (half axes of ellipsoid) [a,b,c] / double
		 *  - 'z-min', 'z-max' (ellipsoid only valid for 'z-min' <= z <= 'z-max' / doubles
		 */
		status += check_array("targets", this_t, "center", i);
		status += check_array("targets", this_t, "x", i);
		status += check_array("targets", this_t, "z", i);
		status += check_array("targets", this_t, "axes", i);
		status += check_float("targets", this_t, "z_min", i);
		status += check_float("targets", this_t, "z_min", i);
		status +=
		    check_string("targets", this_t, "reflecting_surface",
				 i);
		status +=
		    check_string("targets", this_t, "reflectivity", i);
		status += check_file("targets", this_t, "reflectivity", i);
		status +=
		    check_string("targets", this_t, "reflectivity_model",
				 i);
		status +=
		    check_reflectivity_model("targets", this_t,
					     "reflectivity_model", i);

	    } /* end 'ellipsoid' */
	    else if (!strcmp(type, "disk")) {
		/*
		 * disk:
		 *  - array 'P' (center of disk) [x,y,z] / double
		 *  - array 'N' (direction of local x-axis) [x,y,z] / double
		 *  - 'r' (radius of disk) double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
		 *  - string 'reflectivity' (file name of reflectivity spectrum)
		 */
		status += check_array("targets", this_t, "P", i);
		status += check_array("targets", this_t, "N", i);
		status += check_float("targets", this_t, "r", i);
		status += check_array("targets", this_t, "x", i);
		status +=
		    check_string("targets", this_t, "reflectivity", i);
		status += check_file("targets", this_t, "reflectivity", i);
		status +=
		    check_string("targets", this_t, "reflectivity_model",
				 i);
		status +=
		    check_reflectivity_model("targets", this_t,
					     "reflectivity_model", i);

	    } /* end 'disk' */
	    else if (!strcmp(type, "annulus")) {
		/*
		 * annulus:
		 *  - array 'P' (center of annulus) [x,y,z] / double
		 *  - array 'N' (direction of local x-axis) [x,y,z] / double
		 *  - 'R' (outer radius of annulus) double
		 *  - 'r' (inner radius of annulus) double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
		 *  - string 'reflectivity' (file name of reflectivity spectrum)
		 */
		status += check_array("targets", this_t, "P", i);
		status += check_array("targets", this_t, "N", i);
		status += check_float("targets", this_t, "R", i);
		status += check_float("targets", this_t, "r", i);
		status += check_array("targets", this_t, "x", i);
		status +=
		    check_string("targets", this_t, "reflectivity", i);
		status += check_file("targets", this_t, "reflectivity", i);
		status +=
		    check_string("targets", this_t, "reflectivity_model",
				 i);
		status +=
		    check_reflectivity_model("targets", this_t,
					     "reflectivity_model", i);

	    } /* end 'annulus' */
	    else if (!strcmp(type, "cylinder")) {
		/*
		 * cylinder:
		 * C: center point of one end face
		 *  - array 'C' (center of first end face) [x,y,z] / double
		 *  - array 'a' (cylinder axis / local z-axis) [x,y,z] / double
		 *  - 'r' (radius of cylinder) double
		 *  - 'l' (length of cylinder) double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
		 *  - string 'reflecting_surface' (reflecting surface ["inside"|"outside"])
		 *  - string 'reflectivity' (file name of reflectivity spectrum)
		 *  - string 'reflectivity_model' (name of reflectivity model)
		 */
		status += check_array("targets", this_t, "C", i);
		status += check_array("targets", this_t, "a", i);
		status += check_float("targets", this_t, "r", i);
		status += check_float("targets", this_t, "l", i);
		status += check_array("targets", this_t, "x", i);
		status +=
		    check_string("targets", this_t, "reflecting_surface",
				 i);
		status +=
		    check_string("targets", this_t, "reflectivity", i);
		status += check_file("targets", this_t, "reflectivity", i);
		status +=
		    check_string("targets", this_t, "reflectivity_model",
				 i);
		status +=
		    check_reflectivity_model("targets", this_t,
					     "reflectivity_model", i);

	    } /* end 'cylinder' */
	    else if (!strcmp(type, "window")) {
		/*
		 * window:
		 *  - array 'C' (center of front surface) [x,y,z] / double
		 *  - array 'a' (direction of cylinder axis) [x,y,z] / double
		 *  - 'd' (thickness of window) double
		 *  - 'r' (radius of window) double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
		 *  - string 'absorptivity' (file name of absorption spectrum)
		 *  - string 'idx_refraction' (file name of dispersion curve)
		 */
		status += check_array("targets", this_t, "C", i);
		status += check_array("targets", this_t, "a", i);
		status += check_float("targets", this_t, "d", i);
		status += check_float("targets", this_t, "r", i);
		status += check_array("targets", this_t, "x", i);
		status +=
		    check_string("targets", this_t, "absorptivity", i);
		status += check_file("targets", this_t, "absorptivity", i);
		status +=
		    check_string("targets", this_t, "idx_refraction", i);
		status +=
		    check_file("targets", this_t, "idx_refraction", i);
		status +=
		    check_reflectivity_model("targets", this_t,
					     "reflectivity_model", i);
	    } /* end 'window' */
	    else if (!strcmp(type, "paraboloid")) {
		/*
		 * paraboloid:
		 *  - array 'vertex' (lowest point of paraboloid) [x,y,z] / double
		 *  - 'focal_length' (focal length of paraboloid) double
		 *  - array 'z' (direction of local z-axis) [x,y,z] / double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
		 *  - 'z-min', 'z-max' (paraboloid only valid for 'z-min' <= z <= 'z-max' / doubles
		 *  - string 'reflectivity' (file name of reflection spectrum)
		 *  - string 'reflectivity_model'
		 *  - string 'reflecting_surface'
		 */
		status += check_array("targets", this_t, "vertex", i);
		status +=
		    check_float("targets", this_t, "focal_length", i);
		status += check_array("targets", this_t, "z", i);
		status += check_array("targets", this_t, "x", i);
		status += check_float("targets", this_t, "z_min", i);
		status += check_float("targets", this_t, "z_min", i);
		status +=
		    check_string("targets", this_t, "reflectivity", i);
		status += check_file("targets", this_t, "reflectivity", i);
		status +=
		    check_string("targets", this_t, "reflectivity_model",
				 i);
		status +=
		    check_reflectivity_model("targets", this_t,
					     "reflectivity_model", i);
		status +=
		    check_string("targets", this_t, "reflecting_surface",
				 i);
	    } /* end 'paraboloid' */
	    else if (!strcmp(type, "sphere")) {
		/*
		 * sphere:
		 *  - array 'origin' (center of sphere) [x,y,z] / double
		 *  - 'radius' (radius of sphere) double
		 *  - array 'z' (direction of local z-axis) [x,y,z] / double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
		 *  - 'z-min', 'z-max' (ellipsoid only valid for 'z-min' <= z <= 'z-max' / doubles
		 *  - string 'reflectivity' (file name of reflection spectrum)
		 *  - string 'reflectivity_model'
		 *  - string 'reflecting_surface'
		 */
		status += check_array("targets", this_t, "origin", i);
		status += check_float("targets", this_t, "radius", i);
		status += check_array("targets", this_t, "z", i);
		status += check_array("targets", this_t, "x", i);
		status += check_float("targets", this_t, "z_min", i);
		status += check_float("targets", this_t, "z_min", i);
		status +=
		    check_string("targets", this_t, "reflectivity", i);
		status += check_file("targets", this_t, "reflectivity", i);
		status +=
		    check_string("targets", this_t, "reflectivity_model",
				 i);
		status +=
		    check_reflectivity_model("targets", this_t,
					     "reflectivity_model", i);
		status +=
		    check_string("targets", this_t, "reflecting_surface",
				 i);
	    }			/* end 'sphere' */
	}			/* end 'this_t', check next target */
    }				/* end 'targets' section present */

    return status;
}

void init_spectrum(const char *f_name, gsl_spline ** refl_spectrum)
{
    FILE *refl_data;
    double *lambda;
    double *refl;
    size_t n_lambda;

    refl_data = fopen(f_name, "r");
    read_data(refl_data, &lambda, &refl, &n_lambda);
    fclose(refl_data);

    *refl_spectrum = gsl_spline_alloc(gsl_interp_linear, n_lambda);

    gsl_spline_init(*refl_spectrum, lambda, refl, n_lambda);

    free(lambda);
    free(refl);
}

static void write_target_header(const int fd, const char *name,
				const char *type, const double *origin,
				const double *Mat)
{
    char string[1024];

    snprintf(string, 1023,
	     "# %s (%s)\n"
	     "#\n"
	     "# Transformations:\n"
	     "#    g(x, y, z) = MT l(x, y, z) + origin(x, y, z)\n"
	     "#    l(x, y, z) = M (g(x, y, z) - origin(x, y, z))\n"
	     "# M:\n"
	     "#   \t% g\t% g\t% g\n"
	     "#   \t% g\t% g\t% g\n"
	     "#   \t% g\t% g\t% g\n"
	     "# origin:\n"
	     "#   \t% g\t% g\t% g\n#\n"
	     "#\n"
	     "# x\ty\t[z]\tpower\tlambda\t\t(z component is missing for plane targets!)\n"
	     "#\n",
	     name, type, Mat[0], Mat[1], Mat[2], Mat[3], Mat[4], Mat[5],
	     Mat[6], Mat[7], Mat[8], origin[0], origin[1], origin[2]);

    write(fd, string, strlen(string));
}

int init_output(const int file_mode, const char *target_type,
		config_setting_t * this_target, double point[], double M[])
{
    int i;
    int fh = -1;

    if (config_setting_lookup_bool(this_target, "no_output", &i) ==
	CONFIG_FALSE || i == 0) {
	const char *name;
	char f_name[256];

	config_setting_lookup_string(this_target, "name", &name);
	snprintf(f_name, 256, "%s.dat", name);
	fh = open(f_name, O_CREAT | O_WRONLY | file_mode,
		  S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);

	/* write header to dump file */
	if (file_mode == O_TRUNC)
	    write_target_header(fh, name, target_type, point, M);

    }

    return fh;
}

void init_refl_model(const config_setting_t * s, char *model,
		     void **refl_model_params)
{
    const char *S;

    config_setting_lookup_string(s, "reflectivity_model", &S);
    if (!strcmp(S, "specular"))
	*model = (char) SPECULAR;
    else if (!strcmp(S, "lambertian"))
	*model = (char) LAMBERTIAN;
    else if (!strcmp(S, "microfacet_gaussian")) {
	double *number;

	*model = (char) MICROFACET_GAUSSIAN;
	number = (double *) malloc(sizeof(double));

	config_setting_lookup_float(s, "microfacet_gaussian_sigma",
				    number);
	*number /= (180.0 * M_PI);	/* degree to radian */
	*refl_model_params = number;
    }

}

char init_reflecting_surface(config_setting_t * this_target)
{
    const char *S;

    config_setting_lookup_string(this_target, "reflecting_surface", &S);
    if (!strcmp(S, "inside"))
	return INSIDE;
    else
	return OUTSIDE;
}


double *init_M(config_setting_t * this_target, const char *x,
	       const char *z)
{
/*
 * generate transform matrix M to convert
 * between local and global coordinates
 * l2g:   g(x, y, z) = M l(x, y, z) + o(x, y, z))
 * g2l:   l(x, y, z) = MT (g(x, y, z) - o(x, y, z))
 *
 * input x,z and calculate y from x cross y = z
 */
    double *M = (double *) malloc(9 * sizeof(double));

    read_vector_normalize(this_target, x, M);
    read_vector_normalize(this_target, z, &M[6]);

    orthonormalize(M, &M[3], &M[6]);

    return M;
}

void per_thread_init(pthread_key_t key, size_t n)
{
    PTDT_t *data = (PTDT_t *) malloc(sizeof(PTDT_t));

    data->buf = (char *) malloc(BUF_SIZE * n);
    data->i = 0;
    data->flag = 0;

    pthread_setspecific(key, data);
}

void per_thread_flush(int fh, pthread_key_t key, pthread_mutex_t * mutex)
{
    PTDT_t *data = pthread_getspecific(key);

    if (data->i != 0)		/* write rest of buffer to file. */
	if (fh != -1) {
	    pthread_mutex_lock(mutex);
	    write(fh, data->buf, sizeof(char) * data->i);
	    pthread_mutex_unlock(mutex);
	}
}

static void free_refl_model(const char model, void *refl_model_params)
/*
 * free's model specific parameters
 */
{

    switch (model) {

    case MICROFACET_GAUSSIAN:
	{
	    free((double *) refl_model_params);
	}

    }

}

void state_free(int fh, double *M, gsl_spline * s, char model, void *p)
{
    if (fh != -1)
	close(fh);

    free(M);

    if (s)
	gsl_spline_free(s);

    if (model)
	free_refl_model(model, p);
}

#define WRITE_UCHAR(VAR,BUF_PTR,BUF_IDX) do { \
	memcpy(&(BUF_PTR[BUF_IDX]), &VAR, sizeof(unsigned char)); \
	BUF_IDX += sizeof(unsigned char); \
} while(0);

#define WRITE_FLOAT(VAR,BUF_PTR,BUF_IDX) do { \
	float f = (float) VAR; \
	memcpy(&(BUF_PTR[BUF_IDX]), &f, sizeof(float)); \
	BUF_IDX += sizeof(float); \
} while(0);

void store_xy(const int fd, ray_t * ray, const double *hit,
	      const double *m, const double *point, PTDT_t * data,
	      pthread_mutex_t * mutex_writefd)
{
    double hit_local[3];

    /* transform to local coordinates */
    g2l(m, point, hit, hit_local);

    if (data->i == BUF_SIZE * (4 * sizeof(float) + sizeof(unsigned char))) {
	/*
	 * buffer full, write to disk
	 */
	pthread_mutex_lock(mutex_writefd);
	write(fd, data->buf, sizeof(char) * data->i);
	pthread_mutex_unlock(mutex_writefd);

	data->i = 0;
    }

    WRITE_FLOAT(hit_local[0], data->buf, data->i);
    WRITE_FLOAT(hit_local[1], data->buf, data->i);
    WRITE_FLOAT(ray->power, data->buf, data->i);
    WRITE_FLOAT(ray->lambda, data->buf, data->i);
    WRITE_UCHAR(ray->n_refl, data->buf, data->i);
}

void store_xyz(const int fd, ray_t * ray, const double *hit,
	       const double *m, const double *point, PTDT_t * data,
	       pthread_mutex_t * mutex_writefd)
{
    double hit_local[3];

    /* transform to local coordinates */
    g2l(m, point, hit, hit_local);

    if (data->i == BUF_SIZE * (5 * sizeof(float) + sizeof(unsigned char))) {
	/*
	 * buffer full, write to disk
	 */
	pthread_mutex_lock(mutex_writefd);
	write(fd, data->buf, sizeof(char) * data->i);
	pthread_mutex_unlock(mutex_writefd);

	data->i = 0;
    }

    WRITE_FLOAT(hit_local[0], data->buf, data->i);
    WRITE_FLOAT(hit_local[1], data->buf, data->i);
    WRITE_FLOAT(hit_local[2], data->buf, data->i);
    WRITE_FLOAT(ray->power, data->buf, data->i);
    WRITE_FLOAT(ray->lambda, data->buf, data->i);
    WRITE_UCHAR(ray->n_refl, data->buf, data->i);
}

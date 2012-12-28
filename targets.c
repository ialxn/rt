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
		    check_string("targets", this_t, "reflecting_surface", i);
		status +=
		    check_string("targets", this_t, "reflectivity", i);
		status += check_file("targets", this_t, "reflectivity", i);
		status += check_string("targets", this_t, "reflectivity_model", i);
		status +=
		    check_reflectivity_model("targets", this_t,
					     "reflectivity_model", i);

	    }			/* end 'cylinder' */
	}			/* end 'this_t', check next target */
    }				/* end 'targets' section present */

    return status;
}

void init_refl_spectrum(const char *f_name, gsl_spline ** refl_spectrum)
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
    const char *name;
    int i;
    int fh = -1;

    config_setting_lookup_string(this_target, "name", &name);
    if (config_setting_lookup_bool(this_target, "no_output", &i) ==
	CONFIG_FALSE || i == 0) {
	char f_name[256];

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

void per_thread_init(pthread_key_t key, size_t n)
{
    PTDT_t *data = (PTDT_t *) malloc(sizeof(PTDT_t));

    data->buf = (float *) malloc(BUF_SIZE * n * sizeof(float));
    data->i = 0;
    data->flag = 0;

    pthread_setspecific(key, data);
}

void per_thread_flush(int fh, pthread_key_t key, pthread_mutex_t mutex)
{
    PTDT_t *data = pthread_getspecific(key);

    if (data->i != 0)		/* write rest of buffer to file. */
	if (fh != -1) {
	    pthread_mutex_lock(&mutex);
	    write(fh, data->buf, sizeof(float) * data->i);
	    pthread_mutex_unlock(&mutex);
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

void state_free(int fh, gsl_spline * s, char model, void *p)
{
    if (fh != -1)
	close(fh);

    if (s)
	gsl_spline_free(s);

    if (model)
	free_refl_model(model, p);
}

static double *validate_intercepts(const ray_t * ray, const double t1,
				   const double t2)
{
/*
 * smallest positive is possible:
 *   target is in front of rays origin (t>0)
 *   AND
 *   require minimum path length to detect
 *   self intersection because rays initially
 *   emitted by source do not set flag
 */
    double t;
    double *intercept = NULL;

    t = GSL_MIN(t1, t2);
    if (t > GSL_SQRT_DBL_EPSILON) {

	intercept = (double *) malloc(3 * sizeof(double));
	a_plus_cb(intercept, ray->orig, t, ray->dir);

    } else {

	t = GSL_MAX(t1, t2);
	if (t > GSL_SQRT_DBL_EPSILON) {

	    intercept = (double *) malloc(3 * sizeof(double));
	    a_plus_cb(intercept, ray->orig, t, ray->dir);

	}
    }

    return intercept;
}

double *intercept_sphere(const ray_t * ray, const double *center,
			 const double radius)
{
/*
 * from http://wiki.cgsociety.org/index.php/Ray_Sphere_Intersection
 *
 * In vector notation, the equations are as follows:
 *
 * Equation for a sphere (C: center, P: point on sphere, r: radius)
 *	|| P - C || = r^2
 *
 * Equation for a line (O: origin, P: pointson line, D: direction unit
 * vector, t: distance from O)
 *
 *	P = O + t D
 *
 * Solution (d_12) of resulting quadratic equation:
 *
 *	A t^2 + B t + C = 0
 *
 *	with
 *
 *	A = D dot D := 1 (unit vector)
 *	B = 2 (O-C) dot D
 *	C = (O-C) dot (O-C) - r^2
 *
 * if discriminatnt "B^2 - 4AC" (i.e. B^2 - 4C) is negative:	miss
 * 						       else:	2 solns
 *								(may coincide)
 *
 * to aviod poor numerical precision if B is close to sqrt(discriminant) use
 *
 *	t1 = q / A (i.e. q)
 *      t2 = C / q
 *
 * with
 *
 *	q = (-B + sqrt(discriminant)) / 2	B < 0
 *	q = (-B - sqrt(discriminant)) / 2	else
 *
 */
    double OminusC[3];
    double B, C;
    double discriminant;

    diff(OminusC, ray->orig, center);

    B = 2.0 * cblas_ddot(3, OminusC, 1, ray->dir, 1);
    C = cblas_ddot(3, OminusC, 1, OminusC, 1) - radius * radius;

    discriminant = B * B - 4.0 * C;

    if (discriminant < 0.0)	/* does not intercept */
	return NULL;
    else if (fabs(discriminant) < GSL_DBL_EPSILON)	/* 1 intercept */
	return validate_intercepts(ray, -B * 0.5, 0.0);
    else {			/* 2 intercepts */
	double q;

	if (B < 0)
	    q = (-B + sqrt(discriminant)) * 0.5;
	else
	    q = (-B - sqrt(discriminant)) * 0.5;

	/*
	 * validate intercepts (smallest positive)
	 */
	return validate_intercepts(ray, q, C / q);
    }
}

double *intercept_plane(const ray_t * ray, const double *plane_normal,
			const double *plane_point, int *hits_front)
{
/*
 * calculate point of interception d
 *
 * d = {(\mathbf{p_0}-\mathbf{l_0})\cdot\mathbf{n} \over \mathbf{l}\cdot\mathbf{n}}
 *
 * with
 *       p_0: point on the plane
 *         n: normal vector of the plane (|n|=1)
 *       l_0: origin of the line
 *         l: unit vector in direction of the line
 *
 * If the line starts outside the plane and is parallel to the plane, there is no intersection.
 * In this case, the above denominator will be zero and the numerator will be non-zero. If the
 * line starts inside the plane and is parallel to the plane, the line intersects the plane
 * everywhere. In this case, both the numerator and denominator above will be zero. In all other
 * cases, the line intersects the plane once and d represents the intersection as the distance
 * along the line from \mathbf{l_0}.
 *
 * returns dynamically allocated intercept
 *     or
 * NULL if:
 *         - ray is parallel to plane (or ray->orig lies within plane. this should never
 *           occur as planar targets set the flag 'last_was_hit'
 *         - plane is not in front (propagation direction) of ray
 *
 * 'hits_front' is 1 if ray hits front of plane (ray is antiparallel to normal vector)
 * and 0 otherwise
 */

    double t1, t3;
    double t2[3];
    double d;
    double *intercept;

    t1 = cblas_ddot(3, ray->dir, 1, plane_normal, 1);	/* l dot n */

    if (fabs(t1) < GSL_SQRT_DBL_EPSILON)	/* line is parallel to target, no hit possible */
	return NULL;

    if (t1 < 0.0)
	*hits_front = 1;
    else
	*hits_front = 0;

    diff(t2, plane_point, ray->orig);	/* p_0 - l_0 */
    t3 = cblas_ddot(3, t2, 1, plane_normal, 1);	/* (p_0 - l_0) dot N */

    if (fabs(t3) < GSL_SQRT_DBL_EPSILON)	/* line does start in target, conservative */
	return NULL;

    /*
     * 'ray' intercepts target plane
     */
    d = t3 / t1;
    if (d < 0.0)		/* target is not in front */
	return NULL;

    intercept = (double *) malloc(3 * sizeof(double));
    a_plus_cb(intercept, ray->orig, d, ray->dir);

    return intercept;

}

extern void store_xy(const int fd, ray_t * ray, const double *hit,
		     const double *m, const double *point, PTDT_t * data,
		     pthread_mutex_t * mutex_writefd)
{
    double hit_local[3];

    /* transform to local coordinates */
    g2l(m, point, hit, hit_local);

    if (data->i == BUF_SIZE * 4) {	/* buffer full, write to disk */

	pthread_mutex_lock(mutex_writefd);
	write(fd, data->buf, sizeof(float) * data->i);
	pthread_mutex_unlock(mutex_writefd);

	data->i = 0;
    }

    data->buf[data->i++] = (float) hit_local[0];
    data->buf[data->i++] = (float) hit_local[1];
    data->buf[data->i++] = (float) ray->power;
    data->buf[data->i++] = (float) ray->lambda;
}

extern void store_xyz(const int fd, ray_t * ray, const double *hit,
		      const double *m, const double *point, PTDT_t * data,
		      pthread_mutex_t * mutex_writefd)
{
    double hit_local[3];

    /* transform to local coordinates */
    g2l(m, point, hit, hit_local);

    if (data->i == BUF_SIZE * 5) {	/* buffer full, write to disk */

	pthread_mutex_lock(mutex_writefd);
	write(fd, data->buf, sizeof(float) * data->i);
	pthread_mutex_unlock(mutex_writefd);

	data->i = 0;
    }

    data->buf[data->i++] = (float) hit_local[0];
    data->buf[data->i++] = (float) hit_local[1];
    data->buf[data->i++] = (float) hit_local[2];
    data->buf[data->i++] = (float) ray->power;
    data->buf[data->i++] = (float) ray->lambda;
}

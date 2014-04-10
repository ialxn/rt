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
#include "reflect.h"
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

void init_PTDT(const target_t * T)
{
    (T->type->init_PTDT) (T->state);
}

void flush_PTDT_outbuf(const target_t * T)
{
    (T->type->flush_PTDT_outbuf) (T->state);
}

void free_PTDT(void *p)
{
    PTDT_t *data = (PTDT_t *) p;

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

	    }			/* end 'annulus' */
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
	*refl_model_params = number;
    }

}

void free_refl_model(const char model, void *refl_model_params)
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
		     const double *m, const double *point, PTDT_t * data)
{
    double hit_local[3];

    /* transform to local coordinates */
    g2l(m, point, hit, hit_local);

    if (data->i == BUF_SIZE * 4) {	/* buffer full, write to disk */
	write(fd, data->buf, sizeof(float) * data->i);
	data->i = 0;
    }

    data->buf[data->i++] = (float) hit_local[0];
    data->buf[data->i++] = (float) hit_local[1];
    data->buf[data->i++] = (float) ray->power;
    data->buf[data->i++] = (float) ray->lambda;
}

extern void store_xyz(const int fd, ray_t * ray, const double *hit,
		      const double *m, const double *point, PTDT_t * data)
{
    double hit_local[3];

    /* transform to local coordinates */
    g2l(m, point, hit, hit_local);

    if (data->i == BUF_SIZE * 5) {	/* buffer full, write to disk */
	write(fd, data->buf, sizeof(float) * data->i);
	data->i = 0;
    }

    data->buf[data->i++] = (float) hit_local[0];
    data->buf[data->i++] = (float) hit_local[1];
    data->buf[data->i++] = (float) hit_local[2];
    data->buf[data->i++] = (float) ray->power;
    data->buf[data->i++] = (float) ray->lambda;
}

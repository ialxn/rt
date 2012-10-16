/*	targets.c
 *
 * Copyright (C) 2010, 2011 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_cblas.h>

#include "io_util.h"
#include "targets.h"

target_t *target_alloc(const target_type_t * type,
		       config_setting_t * this_t, config_t * cfg,
		       const char *file_mode)
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

    (T->type->init_state) (T->state, this_t, cfg, file_mode);	/* initialize data structures */

    return T;
}

void target_free(target_t * T)
{
    (T->type->free_state) (T->state);
    free(T->state);
    free(T);
}

double *interception(const target_t * T, ray_t * in_ray)
{
    return (T->type->get_intercept) (T->state, in_ray);
}

ray_t *out_ray(const target_t * T, ray_t * in_ray, double *hit,
	       const gsl_rng * r)
{
    return (T->type->get_out_ray) (T->state, in_ray, hit, r);
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
	     *          - "rectangle":              specular reflection on front surface
	     *                                      defined by direction of surface
	     *                                      normal. total absorption on rear
	     *                                      surface.
	     *          - "triangle":               specular reflection on front surface
	     *                                      defined by direction of surface
	     *                                      normal. total absorption on rear
	     *                                      surface.
	     *          - "ellipsoid"             : ellipsoid. specular reflection on
	     *                                      inside surface, total absorption
	     *                                      on outside surface.
	     *          - "disk":                   specular reflection on front surface
	     *                                      defined by direction of surface
	     *                                      normal. total absorption on rear
	     *                                      surface.
	     *          - "annulus":                specular reflection on front surface
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

	    }			/* end 'annulus' */
	}			/* end 'this_t', check next target */
    }				/* end 'targets' section present */

    return status;
}

void init_refl_spectrum(const char *f_name, gsl_spline ** spline,
			gsl_interp_accel ** acc)
{
    FILE *spectrum;
    double *lambda;
    double *refl;
    size_t n_lambda;

    spectrum = fopen(f_name, "r");
    read_data(spectrum, &lambda, &refl, &n_lambda);
    fclose(spectrum);

    /* cspline will be used to interpolate */
    *acc = gsl_interp_accel_alloc();
    *spline = gsl_spline_alloc(gsl_interp_cspline, n_lambda);

    gsl_spline_init(*spline, lambda, refl, n_lambda);

    free(lambda);
    free(refl);
}

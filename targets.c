/*	targets.c
 *
 * Copyright (C) 2010,2011,2012,2013,2014,2015,2016,2017 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <errno.h>
#include <string.h>
#include <stdlib.h>

#include "io_utils.h"
#include "likely.h"
#include "targets.h"


/*
 * public functions to access / manipulate targets
 */
target_t *target_alloc(const target_type_t * type,
		       config_setting_t * this_t, const int file_mode,
		       const int keep_closed, const double P_factor)
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

    /* initialize data structures */
    if ((T->type->init_state) (T->state, this_t, file_mode,
			       keep_closed, P_factor) == ERR) {
	target_free(T);
	return NULL;
    }

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



/*
 * utility functions
 */
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
	     *          - "annulus"
	     *          - "cone"
	     *          - "cylinder"
	     *          - "disk"
	     *          - "ellipsoid"
	     *          - "paraboloid"
	     *          - "one-sided plane_screen"
	     *          - "two-sided plane_screen"
	     *          - "rectangle"
	     *          - "sphere"
	     *          - "triangle"
	     *          - "window"
	     */
	    status += check_string("targets", this_t, "name", i);

	    status +=
		check_return_string("targets", this_t, "type", i, &type);
	    if (!type)
		continue;

	    /* check target specific settings */

	    if (!strcmp(type, "annulus")) {
		/*
		 * annulus:
		 *  - array 'P' (center of annulus) [x,y,z] / double
		 *  - array 'N' (direction of local z-axis) [x,y,z] / double
		 *  - 'R' (outer radius of annulus) double
		 *  - 'r' (inner radius of annulus) double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
		 *  - string 'reflectivity' (file name of reflectivity spectrum)
		 *  - string 'reflectivity_model' (name of reflectivity model)
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
	    else if (!strcmp(type, "cone")) {
		/*
		 * sphere:
		 *  - array 'origin' (base point of cone) [x,y,z] / double
		 *  - 'R' (radius of base disk of cone) double
		 *  - 'r' (radius of top disk of cone, can be zero) double
		 *  - 'h' (height of cone) double
		 *  - array 'axis' (cone axis, local z-axis) [x,y,z] / double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
		 *  - string 'reflectivity' (file name of reflection spectrum)
		 *  - string 'reflectivity_model' (name of reflectivity model)
		 *  - string 'reflecting_surface' ["inside"|"outside"]
		 */
		status += check_array("targets", this_t, "origin", i);
		status += check_float("targets", this_t, "R", i);
		status += check_float("targets", this_t, "r", i);
		status += check_float("targets", this_t, "h", i);
		status += check_array("targets", this_t, "axis", i);
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
		status +=
		    check_string("targets", this_t, "reflecting_surface",
				 i);
	    } /* end 'cone' */
	    else if (!strcmp(type, "cpc")) {
		/*
		 * cpc:
		 * origin: center point of exit aperture
		 *  - array 'origin' (center of exit aperture) [x,y,z] / double
		 *  - array 'axis' (cpc axis / local z-axis) [x,y,z] / double
		 *  - 'acceptance_angle' phi / double
		 *  - 'truncation_angle' (smaller than phi) / double
		 *  - 'exit_radius' / double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
		 *  - string 'reflectivity' (file name of reflectivity spectrum)
		 *  - string 'reflectivity_model' (name of reflectivity model)
		 */
		status += check_array("targets", this_t, "origin", i);
		status += check_array("targets", this_t, "axis", i);
		status +=
		    check_float("targets", this_t, "acceptance_angle", i);
		status +=
		    check_float("targets", this_t, "truncation_angle", i);
		status += check_float("targets", this_t, "exit_radius", i);
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
	    } /* end 'cpc' */
	    else if (!strcmp(type, "cylinder")) {
		/*
		 * cylinder:
		 * C: center point of first end face
		 *  - array 'C' (center of first end face) [x,y,z] / double
		 *  - array 'a' (cylinder axis / local z-axis) [x,y,z] / double
		 *  - 'r' (radius of cylinder) double
		 *  - 'l' (length of cylinder) double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
		 *  - string 'reflecting_surface' ["inside"|"outside"]
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
	    else if (!strcmp(type, "disk")) {
		/*
		 * disk:
		 *  - array 'P' (center of disk) [x,y,z] / double
		 *  - array 'N' (direction of local z-axis) [x,y,z] / double
		 *  - 'r' (radius of disk) double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
		 *  - string 'reflectivity' (file name of reflectivity spectrum)
		 *  - string 'reflectivity_model' (name of reflectivity model)
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
	    else if (!strcmp(type, "ellipsoid")) {
		/*
		 * ellipsoid:
		 *  - array 'center' (center of ellipsoid) [x,y,z] / double
		 *  - array 'x-axis' (direction of local x-axis) [x,y,z] / double
		 *  - array 'z-axis' (direction of local z-axis) [x,y,z] / double
		 *  - array 'axes' (half axes of ellipsoid) [a,b,c] / double
		 *  - 'z-min', 'z-max' (ellipsoid only valid for
		 *                      'z-min' <= z <= 'z-max') / doubles
		 *  - string 'reflecting_surface' ["inside"|"outside"]
		 *  - string 'reflectivity' (file name of reflectivity spectrum)
		 *  - string 'reflectivity_model' (name of reflectivity model)
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
	    else if (!strcmp(type, "paraboloid")) {
		/*
		 * paraboloid:
		 *  - array 'vertex' (lowest point of paraboloid) [x,y,z] / double
		 *  - 'focal_length' (focal length of paraboloid) double
		 *  - array 'z' (direction of local z-axis) [x,y,z] / double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
		 *  - 'z-min', 'z-max' (paraboloid only valid for
		 *                      'z-min' <= z <= 'z-max') / doubles
		 *  - string 'reflecting_surface' ["inside"|"outside"]
		 *  - string 'reflectivity' (file name of reflectivity spectrum)
		 *  - string 'reflectivity_model' (name of reflectivity model)
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
	    else if (!strcmp(type, "one-sided plane screen")) {
		/*
		 * one-sided plane screen:
		 *  - array 'point' (point on plane) [x,y,z] / double
		 *  - array 'normal' (normal vector of plane) [x,y,z] / double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
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
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
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
		 *  - string 'reflectivity' (file name of reflectivity spectrum)
		 *  - string 'reflectivity_model' (name of reflectivity model)
		 *
		 * NOTE: normal of triangle is defined by 'X' cross 'Y'.
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
	    else if (!strcmp(type, "sphere")) {
		/*
		 * sphere:
		 *  - array 'origin' (center of sphere) [x,y,z] / double
		 *  - 'radius' (radius of sphere) double
		 *  - array 'z' (direction of local z-axis) [x,y,z] / double
		 *  - array 'x' (direction of local x-axis) [x,y,z] / double
		 *  - 'z-min', 'z-max' (sphere only valid for
		 *                      'z-min' <= z <= 'z-max' / doubles
		 *  - string 'reflecting_surface' ["inside"|"outside"]
		 *  - string 'reflectivity' (file name of reflectivity spectrum)
		 *  - string 'reflectivity_model' (name of reflectivity model)
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
	    } /* end 'sphere' */
	    else if (!strcmp(type, "triangle")) {
		/*
		 * triangle
		 *  - array 'P1' (corner point of triangle) [x,y,z] / double
		 *  - array 'P2' (corner point of triangle) [x,y,z] / double
		 *  - array 'P3' (corner point of triangle) [x,y,z] / double
		 *  - string 'reflectivity' (file name of reflectivity spectrum)
		 *  - string 'reflectivity_model' (name of reflectivity model)
		 *
		 * NOTE: normal of triangle is defined by
		 *              ('P2'-'P1') cross ('P3'-'P1').
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
		 *  - string 'reflectivity_model' (name of reflectivity model)
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
	    }			/* end 'window' */
	}			/* end 'this_t', check next target */
    }				/* end 'targets' section present */

    return status;
}

static void write_target_header(const int fd, const char *name,
				const char *type, const double P_factor,
				const double *origin, const double *Mat)
{
    char string[1024];

    snprintf(string, 1023,
	     "# %s (%s)\n"
	     "#\n"
	     "# One absorbed ray represents %f W\n"
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
	     name, type, P_factor, Mat[0], Mat[1], Mat[2], Mat[3], Mat[4],
	     Mat[5], Mat[6], Mat[7], Mat[8], origin[0], origin[1],
	     origin[2]);

    write(fd, string, strlen(string));
}

int init_output(const char *target_type, config_setting_t * this_target,
		const int file_mode, const double P_factor,
		union fh_t *output, int *flags, double point[], double M[])
{
    int i;

    if (config_setting_lookup_bool(this_target, "no_output", &i) ==
	CONFIG_FALSE || i == 0) {
	const char *name;
	char f_name[256];

	*flags |= OUTPUT_REQUIRED;

	config_setting_lookup_string(this_target, "name", &name);
	snprintf(f_name, 256, "%s.dat", name);
	if ((output->fh =
	     open(f_name, O_CREAT | O_WRONLY | file_mode,
		  S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)) < 0) {
	    fprintf(stderr,
		    "Error %s (errno=%d) opening output for target %s\n",
		    strerror(errno), errno, name);
	    return ERR;
	} else {		/* write header to dump file (only if new file) */
	    if (file_mode == O_TRUNC)
		write_target_header(output->fh, name, target_type,
				    P_factor, point, M);
	}

	if (*flags & KEEP_CLOSED) {
	    close(output->fh);
	    output->fname = strdup(f_name);
	}

    }

    return NO_ERR;
}

void init_spectrum(config_setting_t * this_target, const char *kw,
		   gsl_spline ** spectrum)
{
    FILE *data;
    double *lambda;
    double *values;
    size_t n_lambda;
    const char *f_name;

    config_setting_lookup_string(this_target, kw, &f_name);
    data = fopen(f_name, "r");
    read_data(data, &lambda, &values, &n_lambda);
    fclose(data);

    *spectrum = gsl_spline_alloc(gsl_interp_linear, n_lambda);
    gsl_spline_init(*spectrum, lambda, values, n_lambda);

    free(lambda);
    free(values);
}

void init_refl_model(const config_setting_t * s,
		     refl_func_pointer_t * refl_func,
		     void **refl_func_pars)
{
    const char *S;

    config_setting_lookup_string(s, "reflectivity_model", &S);

    if (!strcmp(S, "specular")) {
	*refl_func = (refl_func_pointer_t) reflect_specular;
	refl_func_pars = NULL;
    } else if (!strcmp(S, "lambertian")) {
	*refl_func = (refl_func_pointer_t) reflect_lambertian;
	refl_func_pars = NULL;
    } else if (!strcmp(S, "microfacet_gaussian")) {
	double *number;

	*refl_func = (refl_func_pointer_t) reflect_microfacet_gaussian;
	number = (double *) malloc(sizeof(double));
	config_setting_lookup_float(s, "microfacet_gaussian_sigma",
				    number);
	*number /= (180.0 * M_PI);	/* degree to radian */
	*refl_func_pars = number;
    }
}

int init_reflecting_surface(config_setting_t * this_target)
{
    const char *S;

    config_setting_lookup_string(this_target, "reflecting_surface", &S);
    if (!strcmp(S, "inside"))
	return 0;
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

static void write_buffer(union fh_t output, const int flags,
			 PTDT_t * data, pthread_mutex_t * mutex)
{
    if (flags & KEEP_CLOSED) {
	int fd;

	pthread_mutex_lock(mutex);
	fd = open(output.fname, O_APPEND | O_WRONLY);
	write(fd, data->buf, sizeof(char) * data->i);
	close(fd);
	pthread_mutex_unlock(mutex);
    } else {
	pthread_mutex_lock(mutex);
	write(output.fh, data->buf, sizeof(char) * data->i);
	pthread_mutex_unlock(mutex);
    }
    data->i = 0;

}

void per_thread_flush(union fh_t output, const int flags,
		      pthread_key_t key, pthread_mutex_t * mutex)
{
    PTDT_t *data = pthread_getspecific(key);

    if (data->i && (flags & OUTPUT_REQUIRED))	/* write rest of buffer to file. */
	write_buffer(output, flags, data, mutex);
}

void state_free(union fh_t output, int flags, double *M,
		gsl_spline * s, refl_func_pointer_t refl_func, void *p)
{
    if (flags & OUTPUT_REQUIRED) {
	if (flags & KEEP_CLOSED)
	    free(output.fname);
	else
	    close(output.fh);
    }

    free(M);

    if (s)
	gsl_spline_free(s);

    /*
     * free's model specific parameters
     */
    if (refl_func == reflect_microfacet_gaussian)
	free((double *) p);

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

void store_xy(union fh_t output, const int flags, ray_t * ray,
	      const double *hit, const double *m, const double *point,
	      PTDT_t * data, pthread_mutex_t * mutex_writefd)
{
    double hit_local[3];

    /* transform to local coordinates */
    g2l(m, point, hit, hit_local);

    if (unlikely		/* buffer full is unlikely */
	(data->i ==
	 BUF_SIZE * (3 * sizeof(float) + sizeof(unsigned char))))
	write_buffer(output, flags, data, mutex_writefd);

    WRITE_FLOAT(hit_local[0], data->buf, data->i);
    WRITE_FLOAT(hit_local[1], data->buf, data->i);
    WRITE_FLOAT(ray->lambda, data->buf, data->i);
    WRITE_UCHAR(ray->n_refl, data->buf, data->i);
}

void store_xyz(union fh_t output, const int flags, ray_t * ray,
	       const double *hit, const double *m, const double *point,
	       PTDT_t * data, pthread_mutex_t * mutex_writefd)
{
    double hit_local[3];

    /* transform to local coordinates */
    g2l(m, point, hit, hit_local);

    if (unlikely		/* buffer full is unlikely */
	(data->i ==
	 BUF_SIZE * (4 * sizeof(float) + sizeof(unsigned char))))
	write_buffer(output, flags, data, mutex_writefd);

    WRITE_FLOAT(hit_local[0], data->buf, data->i);
    WRITE_FLOAT(hit_local[1], data->buf, data->i);
    WRITE_FLOAT(hit_local[2], data->buf, data->i);
    WRITE_FLOAT(ray->lambda, data->buf, data->i);
    WRITE_UCHAR(ray->n_refl, data->buf, data->i);
}

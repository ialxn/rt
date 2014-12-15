/*	intercept.c
 *
 * Copyright (C) 2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_poly.h>

#include "intercept.h"


static int soln_in_range(const double s, const double min,
			 const double max, const double r0[3],
			 const double dir[3])
{
/*
 * solution is 'r0' + 's' x 'dir'
 * we are only interested in z-component
 */
    double z;

    if (s < GSL_SQRT_DBL_EPSILON)
	return 0;		/* negative solution */

    z = r0[2] + s * dir[2];
    if (z < min || z > max)
	return 0;
    else
	return 1;
}

static int find_first_soln_restricted(const int n_solns,
				      const double x_small,
				      const double x_large,
				      const double z_min,
				      const double z_max,
				      const double *l_orig,
				      const double *l_dir,
				      double *l_intercept)
{
/*
 * calculate first intercept (in local system) that fullfills the restriction
 *      z_min <= z_component_of_solution <= z_max
 * if none ist found return 0 (and 'l_intercept' is not modified) and 1
 * with 'l_intercept' set otherwise
 */
    if (n_solns == 0)
	return 0;

    if (x_small < GSL_SQRT_DBL_EPSILON && x_large < GSL_SQRT_DBL_EPSILON)
	return 0;		/* none valid */

    if (soln_in_range(x_small, z_min, z_max, l_orig, l_dir))
	a_plus_cb(l_intercept, l_orig, x_small, l_dir);	/* smaller valid */
    else if (soln_in_range(x_large, z_min, z_max, l_orig, l_dir))
	a_plus_cb(l_intercept, l_orig, x_large, l_dir);	/* larger valid */
    else
	return 0;		/* none valid */

    return 1;
}

static double *find_first_soln(const int n_solns, const double x_small,
			       const double x_large, const ray_t * ray)
{
    /*
     * first intercept corresponds to smallest positive solution.
     * require minimum travels corresponding to GSL_SQRT_DBL_EPSILON
     * to detect self intersections.
     */
    double *intercept;
    double x;

    if (n_solns == 0)
	return NULL;

    if (x_small > GSL_SQRT_DBL_EPSILON)
	x = x_small;
    else if (x_large > GSL_SQRT_DBL_EPSILON)
	x = x_large;
    else			/* both solutions are invalid */
	return NULL;

    intercept = (double *) malloc(3 * sizeof(double));
    a_plus_cb(intercept, ray->orig, x, ray->dir);

    return intercept;
}

static double *test_cyl_intercept(const double x, const double *orig,
				  const double *dir, const double *c,
				  const double *a, const double l)
{
/*
 * check if intercept is between two faces
 * based on file LinePosition3d.m from same source.
 *
 * projection of vector 'intercept' - 'c' onto 'a' must fullfill 0 < x < 'd'
 * if we exclude faces
 */
    double *intercept = NULL;

    if (x > GSL_DBL_EPSILON) {	/* cylider in front */
	double t1, t3[3];

	intercept = (double *) malloc(3 * sizeof(double));
	a_plus_cb(intercept, orig, x, dir);
	diff(t3, intercept, c);
	t1 = cblas_ddot(3, t3, 1, a, 1);

	if (t1 < 0.0 || t1 > l) {	/* out of bounds */
	    free(intercept);
	    intercept = NULL;
	}
    }
    return intercept;
}



/*
 * surface normals
 */
void cyl_surf_normal(double *const intercept, const double *C,
		     const double *a, const double r, double *const normal)
{
/*
 * calculate surface normal 'normal' at 'icpt' on cylinder wall
 */
    double IC[3];
    double t;

    diff(IC, intercept, C);
    t = cblas_ddot(3, IC, 1, a, 1);
    a_plus_cb(normal, IC, -t, a);

    cblas_dscal(3, 1.0 / r, normal, 1);
}

void ell_surf_normal(const double *point, const double *axes,
		     double *const normal)
{
/*
 * this normal is for ellisoid cetered at (0,0,0) i.e.
 *
 *      x^2/a^2 + y^2/b^2 + z^2/c^2 - 1 = 0
 *
 * and points outwards! (factor 2.0 left out)
 */
    int i;
    double norm;

    for (i = 0, norm = 0.0; i < 3; i++) {
	normal[i] = point[i] / axes[i];
	norm += normal[i] * normal[i];
    }
    norm = sqrt(norm);
    cblas_dscal(3, 1.0 / norm, normal, 1);
}

void par_surf_normal(const double *point, const double foc2,
		     double *const normal)
{
/*
 * this normal {x/2a, y/2a, -1} is for paraboloid centered at (0,0,0) i.e.
 * x^2/4a + y^2/4a -z = 0 and points outwards!
 * foc2: 1/2a
 * point, normal: local coordinate system
 */
    normal[0] = foc2 * point[0];
    normal[1] = foc2 * point[1];
    normal[2] = -1.0;

    normalize(normal);
}

void sph_surf_normal(const double *point, double *normal)
{
/*
 * surface normal of sphere at 'point' is 'point' in local system
 * with origin at cetnter of sphere
 */
    memcpy(normal, point, 3 * sizeof(double));
    normalize(normal);
}



/*
 * intercepts
 */
double *intercept_cylinder(const ray_t * ray, const double *c,
			   const double *a, const double r, const double l,
			   int *hits_outside)
{
/*
 * returns first (closes to origin of ray) intercept (x,y,z) between
 * ray and cylinder wall.
 * cylinder is defined by:
 *     - center point 'c' of first face
 *     - normalized vector 'a' pointing in direction of second face.
 *     - radius 'r' of cylinder
 *     - length 'l' of the cylinder.
 *
 * based on:
 * http://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d/content/geom3d/geom3d/
 * file: intersectLineCylinder.m
 */
    double t1, t3[3];
    double e[3], f[3];
    double A, B, C;
    double x_small, x_large;
    int n_solns;
    double *intercept;

    /*
     * Starting point of the line:              l0 = line(1:3)';
     *                                               = ray->origin
     * Direction vector of the line:            dl = line(4:6)';
     *                                               = ray->dir
     * Starting position of the cylinder:       c0 = cylinder(1:3)';
     *                                               = c
     * Direction vector of the cylinder:        dc = cylinder(4:6)' - c0;
     *                                               = a*l
     * Radius of the cylinder:                   r = cylinder(7);
     *
     * Resolution of a quadratic equation to find the increment
     * Substitution of parameters
     * e = dl - (dot(dl,dc)/dot(dc,dc))*dc;
     * f = (l0-c0) - (dot(l0-c0,dc)/dot(dc,dc))*dc;
     * Note: 'a' is normalized. 'd' * 'a' -> dc
     */

    t1 = cblas_ddot(3, ray->dir, 1, a, 1);	/* dot(dl,dc)/dot(dc,dc) */
    a_plus_cb(e, ray->dir, -t1 * l, a);

    diff(t3, ray->orig, c);	/* l0-c0 */
    t1 = cblas_ddot(3, t3, 1, a, 1);	/* dot((l0-c0),dc)/dot(dc,dc) */
    a_plus_cb(f, t3, -t1 * l, a);

    /*
     * Coefficients of 2-nd order equation
     * A = dot(e, e);
     * B = 2*dot(e,f);
     * C = dot(f,f) - r^2;
     */
    A = cblas_ddot(3, e, 1, e, 1);
    B = 2.0 * cblas_ddot(3, e, 1, f, 1);
    C = cblas_ddot(3, f, 1, f, 1) - r * r;

    n_solns = gsl_poly_solve_quadratic(A, B, C, &x_small, &x_large);

    if (n_solns == 0)
	return NULL;

    if ((intercept =
	 test_cyl_intercept(x_small, ray->orig, ray->dir, c, a,
			    l)) == NULL)
	/*
	 * smaller solution not valid, try larger instead
	 */
	intercept =
	    test_cyl_intercept(x_large, ray->orig, ray->dir, c, a, l);

    /*
     * if valid intercept found, test if outside surface is hit
     * where 'ray->dir' dot "surface_normal" is negative
     */
    if (intercept) {
	cyl_surf_normal(intercept, c, a, r, t3);
	if (cblas_ddot(3, ray->dir, 1, t3, 1) < 0.0)
	    *hits_outside = 1;
	else
	    *hits_outside = 0;
    }
    return intercept;
}

double *intercept_ellipsoid(const ray_t * ray, const double *M,
			    const double center[3], const double axes[3],
			    const double z_min, const double z_max,
			    int *hits_outside)
{


    int i;
    double r_O[3], r_N[3];	/* origin, direction of ray in local system */
    double N_ellipsoid[3];
    double A = 0.0, B = 0.0, C = -1.0;
    double x_small, x_large;
    int n_solns;
    double l_intercept[3];
    double *intercept;


    /*
     * calculate point of interception D
     */
    /*
     * transform 'ray' from global to local system
     * origin 'ray': rotate / translate by origin of local system
     * dir 'ray': rotate only
     */
    g2l(M, center, ray->orig, r_O);
    g2l_rot(M, ray->dir, r_N);

    /*
     * solve quadratic equation
     */
    for (i = 0; i < 3; i++) {
	A += r_N[i] * r_N[i] / axes[i];
	B += 2.0 * r_O[i] * r_N[i] / axes[i];
	C += r_O[i] * r_O[i] / axes[i];
    }

    n_solns = gsl_poly_solve_quadratic(A, B, C, &x_small, &x_large);
    if (!find_first_soln_restricted
	(n_solns, x_small, x_large, z_min, z_max, r_O, r_N, l_intercept))
	return NULL;

    /*
     * valid intercept found, test if outside surface is hit
     * where 'ray->dir' dot "surface_normal" is negative
     */
    ell_surf_normal(l_intercept, axes, N_ellipsoid);
    if (cblas_ddot(3, r_N, 1, N_ellipsoid, 1) < 0.0)
	*hits_outside = 1;
    else
	*hits_outside = 0;

    /* convert to global coordinates, origin is 'state->center' */
    intercept = (double *) malloc(3 * sizeof(double));
    l2g(M, center, l_intercept, intercept);

    return intercept;
}

double *intercept_paraboloid(const ray_t * ray, const double *M,
			     const double vertex[3], const double foc2,
			     const double foc4, const double z_min,
			     const double z_max, int *hits_outside)
{
    double r_O[3], r_N[3];
    double l_intercept[3];
    double *intercept;
    double N_paraboloid[3];
    double A, B, C;
    double x_small, x_large;
    int n_solns;

    /*
     * transform 'ray' from global to local system
     */
    g2l(M, vertex, ray->orig, r_O);
    g2l_rot(M, ray->dir, r_N);
    /*
     * calculate coefficients of quadratic equation
     */
    A = foc4 * (r_N[0] * r_N[0] + r_N[1] * r_N[1]);
    B = foc2 * (r_N[0] * r_O[0] + r_N[1] * r_O[1]) - r_N[2];
    C = foc4 * (r_O[0] * r_O[0] + r_O[1] * r_O[1]) - r_O[2];

    n_solns = gsl_poly_solve_quadratic(A, B, C, &x_small, &x_large);
    if (!find_first_soln_restricted
	(n_solns, x_small, x_large, z_min, z_max, r_O, r_N, l_intercept))
	return NULL;

    /*
     * valid intercept found, test if outside surface is hit
     * where 'ray->dir' dot "surface_normal" is negative
     */
    par_surf_normal(l_intercept, foc2, N_paraboloid);
    if (cblas_ddot(3, r_N, 1, N_paraboloid, 1) < 0.0)
	*hits_outside = 1;
    else
	*hits_outside = 0;

    /* convert to global coordinates, origin is 'state->vertex' */
    intercept = (double *) malloc(3 * sizeof(double));
    l2g(M, vertex, l_intercept, intercept);

    return intercept;
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

double *intercept_sphere(const ray_t * ray, const double *M,
			 const double *center, const double radius,
			 const double z_min, const double z_max,
			 int *hits_outside)
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
 */
    double *intercept;
    double B, C;
    double x_small, x_large;
    int n_solns;

    if (M == NULL) {		/* called from 'virtual_target_solid_sphere */
	double OminusC[3];

	diff(OminusC, ray->orig, center);

	B = 2.0 * cblas_ddot(3, OminusC, 1, ray->dir, 1);
	C = cblas_ddot(3, OminusC, 1, OminusC, 1) - radius * radius;

	n_solns = gsl_poly_solve_quadratic(1.0, B, C, &x_small, &x_large);
	intercept = find_first_soln(n_solns, x_small, x_large, ray);

    } else {
	double r_O[3], r_N[3];	/* origin, direction of ray in local system */
	double l_intercept[3];
	double N_sphere[3];

	g2l(M, center, ray->orig, r_O);
	g2l_rot(M, ray->dir, r_N);

	B = 2.0 * cblas_ddot(3, r_O, 1, r_N, 1);
	C = cblas_ddot(3, r_O, 1, r_O, 1) - radius * radius;

	n_solns = gsl_poly_solve_quadratic(1.0, B, C, &x_small, &x_large);
	if (!find_first_soln_restricted
	    (n_solns, x_small, x_large, z_min, z_max, r_O, r_N,
	     l_intercept))
	    return NULL;

	sph_surf_normal(l_intercept, N_sphere);
	if (cblas_ddot(3, r_N, 1, N_sphere, 1) < 0)
	    *hits_outside = 1;
	else
	    *hits_outside = 0;

	intercept = (double *) malloc(3 * sizeof(double));
	l2g(M, center, l_intercept, intercept);

    }

    return intercept;
}

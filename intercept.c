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
#include <gsl/gsl_cblas.h>

#include "intercept.h"

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

static double *test_cyl_intercept(const double x, const double *orig,
				  const double *dir, const double *c,
				  const double *a, const double l)
{
/*
 * check if intercept is between two faces
 * based on file LinePosition3d.m from same source.
 *
 * projection of vector 'icpt' - 'c' onto 'a' must fullfill 0 < x < 'd'
 * if we exclude faces
 */
    double *icpt = NULL;

    if (x > GSL_DBL_EPSILON) {	/* cylider in front */
	double t1, t3[3];

	icpt = (double *) malloc(3 * sizeof(double));
	a_plus_cb(icpt, orig, x, dir);
	diff(t3, icpt, c);
	t1 = cblas_ddot(3, t3, 1, a, 1);

	if (t1 < 0.0 || t1 > l) {	/* out of bounds */
	    free(icpt);
	    icpt = NULL;
	}
    }
    return icpt;
}



/*
 * surface normals
 */
void cyl_surf_normal(double *const icpt, const double *C,
		     const double *a, const double r, double *const normal)
{
/*
 * calculate surface normal 'normal' at 'icpt' on cylinder wall
 */
    double IC[3];
    double t;

    diff(IC, icpt, C);
    t = cblas_ddot(3, IC, 1, a, 1);
    a_plus_cb(normal, IC, -t, a);

    cblas_dscal(3, 1.0 / r, normal, 1);
}


/*
 * intercepts
 */
extern double *intercept_cylinder(const ray_t * ray, const double *c,
				  const double *a, const double r,
				  const double l, int *hits_outside)
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
    double delta;
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

    /*
     * compute discriminant
     */
    delta = B * B - 4.0 * A * C;

    /*
     * check existence of solution(s)
     */
    if (delta < 0.0)		/* no interception */
	return NULL;
    else {
	/*
	 * 2 (maybe degenerate) solutions (intercepts with infinite cylinder).
	 * identify the ones intercepting the cylinder between
	 * its end faces (0, 1, or 2 solutions). projection of vector
	 * 'icpt' - 'c' onto 'a' must fullfill 0 < x < 'd'
	 * pick the one closest to origin of ray and assign it to icpt
	 */
	const double s = sqrt(delta);
	const double AA = 1.0 / (2.0 * A);
	const double x1 = (-B + s) * AA;
	const double x2 = (-B - s) * AA;
	double xmin, xmax;

	/*
	 * determine the correct solution corresponding
	 * to the smalles positive 'x' within the length
	 * of the cylinder.
	 */
	if (x1 > x2) {
	    xmax = x1;
	    xmin = x2;
	} else {
	    xmax = x2;
	    xmin = x1;
	}

	if ((intercept =
	     test_cyl_intercept(xmin, ray->orig, ray->dir, c, a,
				l)) == NULL)
	    /*
	     * xmin not valid, try xmax instead
	     */
	    intercept =
		test_cyl_intercept(xmax, ray->orig, ray->dir, c, a, l);

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


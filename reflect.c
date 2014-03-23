/*	reflect.c
 *
 * Copyright (C) 2014 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <string.h>
#include <gsl/gsl_cblas.h>

#include "reflect.h"


static void reflect_specular(ray_t * r, const double N[3],
			     const double P[3])
/*
 * return specularly reflected ray 'r'. surface normal of reflecting surface
 * is 'N'. incoming ray has been determined before to intersect at 'P'. thus
 * origin or refelcted ray will be 'P'.
 */
{
    const double t = cblas_ddot(3, N, 1, r->dir, 1);	/* 'N' dot 'r' */

    cblas_daxpy(3, -2.0 * t, N, 1, r->dir, 1);	/* 'r' - 2 * 'N' dot 'r' * 'N' */
    memcpy(r->orig, P, 3 * sizeof(double));	/* update origin */
}


void reflect(ray_t * r, const double N[3], const double P[3],
	     const char model)
{

    switch (model) {

    case SPECULAR:
	reflect_specular(r, N, P);
	break;
    }

}

/*	sources.c
 *
 * Copyright (C) 2010,2011,2012,2013 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <string.h>

#include "io_utils.h"
#include "sources.h"

source_t *source_alloc(const source_type_t * T, config_setting_t * this_s)
{
    source_t *S;

    /* generic part */
    if ((S = (source_t *) malloc(sizeof(source_t))) == NULL) {
	fprintf(stderr, "failed to allocate space for 'source_t'\n");
	return NULL;
    }

    if ((S->state = malloc(T->size)) == NULL) {
	free(S);
	fprintf(stderr,
		"failed to allocate space for generic state of source\n");
	return NULL;
    };

    S->type = T;

    (S->type->init_state) (S->state, this_s);	/* initialize data structures */

    return S;
}

void source_free(source_t * S)
{
    (S->type->free_state) (S->state);
    free(S->state);
    free(S);
}

ray_t *emit_ray(const source_t * S, const gsl_rng * r)
{
    return (S->type->emit_ray) (S->state, r);
}

const char *get_source_type(const source_t * S)
{
    return S->type->type;
}

const char *get_source_name(const source_t * S)
{
    return (S->type->get_source_name) (S->state);
}

double get_ppr(const source_t * S)
{
    return (S->type->get_source_ppr) (S->state);
}

void init_rays_remain(const source_t * S)
{
    (S->type->init_rays_remain) (S->state);
}


int check_sources(config_t * cfg)
{
    int status = NO_ERR;
    const config_setting_t *s = config_lookup(cfg, "sources");

    if (s == NULL) {
	fprintf(stderr, "missing 'sources' keyword\n");
	status += ERR;
    } else {			/* 'sources' section present */
	int i;
	const int n_sources = config_setting_length(s);

	if (n_sources == 0) {
	    fprintf(stderr, "empty 'sources' section\n");
	    status += ERR;
	}

	for (i = 0; i < n_sources; ++i) {
	    const char *type = NULL;
	    config_setting_t *this_s =
		config_setting_get_elem(s, (unsigned int) i);

	    /*
	     * keywords common to all sources
	     * 'name':  identifier / string
	     * 'type':  type of source / string
	     *          - "uniform point_source": uniform point source
	     * 'power': power [W] of source / double
	     * 'n_rays': number of rays used for this source / int
	     * 'spectrum': name of file containing spectrum of source
	     */
	    status += check_string("sources", this_s, "name", i);
	    status += check_float("sources", this_s, "power", i);
	    status += check_int("sources", this_s, "n_rays", i);
	    status += check_string("sources", this_s, "spectrum", i);
	    status += check_file("sources", this_s, "spectrum", i);

	    status +=
		check_return_string("sources", this_s, "type", i, &type);
	    if (!type)
		continue;

	    /* check source specific settings */
	    if (!strcmp(type, "uniform point source")) {
		/*
		 * uniform point source:
		 *  - array 'origin' [x,y,z] / double
		 */
		status += check_array("sources", this_s, "origin", i);

	    } /* end 'uniform point source' */
	    else if (!strcmp(type, "sphere")) {
		/*
		 * sphere:
		 *  - array 'origin' [x,y,z] / double
		 *  - 'radius' / double
		 */
		status += check_array("sources", this_s, "origin", i);
		status += check_float("sources", this_s, "radius", i);

	    } /* end 'sphere' */
	    else if (!strcmp(type, "spot source")) {
		/*
		 * spot source:
		 *  - array 'origin' [x,y,z] / double
		 *  - array 'direction' [x,y,z] / double
		 *  - 'theta' / double
		 */
		status += check_array("sources", this_s, "origin", i);
		status += check_array("sources", this_s, "direction", i);
		status += check_float("sources", this_s, "theta", i);

	    }			/* end 'uniform point source' */
	}			/* end 'this_s', check next source */
    }

    return status;
}

static double *calc_CDF(const double *I, const double *lambda,
			const size_t n_lambda)
/*
 * calculate normalized cumulative spectrum (CDF).
 * this will be the cumulative distribution function
 * needed to obtain a random wavelength.
 *
 * careful: we have to make sure that the interval
 *          of the random varialble (0..1) is
 *          mapped to the interval (lambda_min
 *          to lambda_max). thus
 *          0 <= lambda_min < l <= lambda_max-lambda_min
 *          'lambda' has been shifted by 'lambda_min'
 *          by caller (this will be reversed
 *          in 'get_new_ray()' by adding 'lambda_min')
 *          make the CDF fall into the range 0..1
 *          after normalization and shift.
 */
{
    size_t i;
    double *CDF = (double *) malloc(n_lambda * sizeof(double));

    /* calculate CDF. include offset by 'I[0]' */
    for (i = 1, CDF[0] = 0.0; i < n_lambda; i++)
	CDF[i] = CDF[i - 1] + (I[i] - I[0]) * (lambda[i] - lambda[i - 1]);

    /* normalize 'CDF' */
    for (i = 0; i < n_lambda; i++)
	CDF[i] /= CDF[n_lambda - 1];

    return CDF;
}

void init_spectrum(const char *f_name, gsl_spline ** spline,
		   double *lambda_min)
{
    FILE *spectrum;
    double *lambda;
    double *I;
    double *CDF;
    size_t n_lambda;
    size_t i;

    spectrum = fopen(f_name, "r");
    read_data(spectrum, &lambda, &I, &n_lambda);
    fclose(spectrum);

    /* shift 'lambda' by lambda_min */
    *lambda_min = lambda[0];
    for (i = 0; i < n_lambda; i++)
	lambda[i] -= *lambda_min;

    CDF = calc_CDF(I, lambda, n_lambda);

    /* cspline will be used to interpolate */
    *spline = gsl_spline_alloc(gsl_interp_cspline, n_lambda);

    /* 'CDF' -> x and 'lambda' -> y */
    gsl_spline_init(*spline, CDF, lambda, n_lambda);

    free(lambda);
    free(I);
    free(CDF);
}

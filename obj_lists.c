/*	obj_list.c
 *
 * Copyright (C) 2010,2011,2012,2013,2014,2015,2016 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <stdlib.h>
#include <string.h>

#include "obj_lists.h"


source_list_t *init_sources(config_t * cfg, int *n_sources,
			    const double P_factor)
{
    int i;
    int n_del = 0;
    source_list_t *s_list;
    const config_setting_t *sources = config_lookup(cfg, "sources");

    s_list = (source_list_t *) malloc(sizeof(source_list_t));
    INIT_LIST_HEAD(&s_list->list);

    *n_sources = config_setting_length(sources);
    for (i = 0; i < *n_sources; ++i) {	/* iterate through all sources */
	source_list_t *new_entry;
	source_t *new_source;
	config_setting_t *this_source =
	    config_setting_get_elem(sources, (unsigned int) i);
	const char *name;	/* name (identifier) of 'this_source' */
	const char *type;	/* type (identifier) of 'this_source' */

	config_setting_lookup_string(this_source, "name", &name);

	/*
	 * allocate the correct source type.
	 * because the 'type' listed in the configuration is a string
	 * and not a C-type we need the 'if() else if() else()' ladder
	 * going over all possibilities and using the correct (first)
	 * argument in the call to 'source_alloc()'. 
	 */
	config_setting_lookup_string(this_source, "type", &type);
	if (!strcmp(type, "arc"))
	    new_source = source_alloc(source_arc, this_source, P_factor);
	else if (!strcmp(type, "solid cone"))
	    new_source =
		source_alloc(source_solid_cone, this_source, P_factor);
	else if (!strcmp(type, "solid cylinder"))
	    new_source =
		source_alloc(source_solid_cylinder, this_source, P_factor);
	else if (!strcmp(type, "solid sphere"))
	    new_source =
		source_alloc(source_solid_sphere, this_source, P_factor);
	else if (!strcmp(type, "sphere"))
	    new_source =
		source_alloc(source_sphere, this_source, P_factor);
	else if (!strcmp(type, "spot source"))
	    new_source = source_alloc(source_spot, this_source, P_factor);
	else if (!strcmp(type, "uniform point source"))
	    new_source =
		source_alloc(source_uniform_point_source, this_source,
			     P_factor);
	else {
	    fprintf(stderr,
		    "    unknown source type (%s) found. Ignoring source %s\n",
		    type, name);
	    ++n_del;
	    continue;
	}

	new_entry = (source_list_t *) malloc(sizeof(source_list_t));
	new_entry->s = new_source;
	list_add_tail(&new_entry->list, &s_list->list);

    }
    (*n_sources) -= n_del;
    return s_list;
}


static void add_targets(target_list_t * t_list, config_t * cfg,
			int *n_targets, const int file_mode,
			const int keep_closed, const double P_factor)
{

    int i;
    int n_del = 0;
    const config_setting_t *targets = config_lookup(cfg, "targets");

    *n_targets = config_setting_length(targets);

    for (i = 0; i < *n_targets; ++i) {	/* iterate through all targets */
	target_list_t *new_entry;
	target_t *new_target;
	config_setting_t *this_target =
	    config_setting_get_elem(targets, (unsigned int) i);
	const char *name;	/* name (identifier) of 'this_target' */
	const char *type;	/* type (identifier) of 'this_target' */

	config_setting_lookup_string(this_target, "name", &name);

	/*
	 * allocate the correct target type.
	 * because the 'type' listed in the configuration is a string
	 * and not a C-type we need the 'if() else if() else()' ladder
	 * going over all possibilities and using the correct (first)
	 * argument in the call to 'target_alloc()'. 
	 */
	config_setting_lookup_string(this_target, "type", &type);
	if (!strcmp(type, "annulus"))
	    new_target =
		target_alloc(target_annulus, this_target, file_mode,
			     keep_closed, P_factor);
	else if (!strcmp(type, "cone"))
	    new_target =
		target_alloc(target_cone, this_target, file_mode,
			     keep_closed, P_factor);
	else if (!strcmp(type, "cylinder"))
	    new_target =
		target_alloc(target_cylinder, this_target, file_mode,
			     keep_closed, P_factor);
	else if (!strcmp(type, "disk"))
	    new_target =
		target_alloc(target_disk, this_target, file_mode,
			     keep_closed, P_factor);
	else if (!strcmp(type, "ellipsoid"))
	    new_target =
		target_alloc(target_ellipsoid, this_target, file_mode,
			     keep_closed, P_factor);
	else if (!strcmp(type, "paraboloid"))
	    new_target =
		target_alloc(target_paraboloid, this_target, file_mode,
			     keep_closed, P_factor);
	else if (!strcmp(type, "one-sided plane screen"))
	    new_target =
		target_alloc(target_plane_screen_one_sided, this_target,
			     file_mode, keep_closed, P_factor);
	else if (!strcmp(type, "two-sided plane screen"))
	    new_target =
		target_alloc(target_plane_screen_two_sided, this_target,
			     file_mode, keep_closed, P_factor);
	else if (!strcmp(type, "rectangle"))
	    new_target =
		target_alloc(target_rectangle, this_target, file_mode,
			     keep_closed, P_factor);
	else if (!strcmp(type, "sphere"))
	    new_target =
		target_alloc(target_sphere, this_target, file_mode,
			     keep_closed, P_factor);
	else if (!strcmp(type, "triangle"))
	    new_target =
		target_alloc(target_triangle, this_target, file_mode,
			     keep_closed, P_factor);
	else if (!strcmp(type, "window"))
	    new_target =
		target_alloc(target_window, this_target, file_mode,
			     keep_closed, P_factor);
	else {
	    fprintf(stderr,
		    "    unknown target type (%s) found. Ignoring target %s\n",
		    type, name);
	    ++n_del;
	    continue;
	}

	new_entry = (target_list_t *) malloc(sizeof(target_list_t));
	new_entry->t = new_target;
	list_add_tail(&new_entry->list, &t_list->list);
    }
    (*n_targets) -= n_del;
}


static void add_virtual_targets(target_list_t * t_list, config_t * cfg)
{
    /*
     * add the virtual targets that correspond to solid sources here
     * curently implemented:
     *          - solid cylinder
     *          - solid sphere
     */
    const config_setting_t *sources = config_lookup(cfg, "sources");
    int n_sources = config_setting_length(sources);
    int i;

    for (i = 0; i < n_sources; ++i) {	/* iterate through all sources */
	target_list_t *new_entry;
	target_t *new_virtual_target;
	const char *type;
	config_setting_t *this_source =
	    config_setting_get_elem(sources, (unsigned int) i);

	config_setting_lookup_string(this_source, "type", &type);
	if (!strcmp(type, "solid cylinder"))
	    new_virtual_target =
		target_alloc(virtual_target_solid_cylinder, this_source, 0,
			     0, 0);
	else if (!strcmp(type, "solid cone"))
	    new_virtual_target =
		target_alloc(virtual_target_solid_cone, this_source, 0, 0,
			     0);
	else if (!strcmp(type, "solid sphere"))
	    new_virtual_target =
		target_alloc(virtual_target_solid_sphere, this_source, 0,
			     0, 0);
	else			/* transparent source found, no virtual target added */
	    continue;

	new_entry = (target_list_t *) malloc(sizeof(target_list_t));
	new_entry->t = new_virtual_target;
	list_add_tail(&new_entry->list, &t_list->list);
    }
}


target_list_t *init_targets(config_t * cfg, int *n_targets,
			    const int file_mode, const int keep_closed,
			    const double P_factor)
{

    target_list_t *t_list =
	(target_list_t *) malloc(sizeof(target_list_t));

    INIT_LIST_HEAD(&t_list->list);

    add_targets(t_list, cfg, n_targets, file_mode, keep_closed, P_factor);
    add_virtual_targets(t_list, cfg);

    return t_list;
}


void source_list_free(source_list_t * s)
{
    struct list_head *pos, *pos_t;

    list_for_each_safe(pos, pos_t, &s->list) {
	source_list_t *this = list_entry(pos, source_list_t, list);

	source_free(this->s);
	list_del(pos);
	free(this);
    }
}


void target_list_free(target_list_t * t)
{
    struct list_head *pos, *pos_t;

    list_for_each_safe(pos, pos_t, &t->list) {
	target_list_t *this = list_entry(pos, target_list_t, list);

	target_free(this->t);
	list_del(pos);
	free(this);
    }
}

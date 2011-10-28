/*	obj_list.c
 *
 * Copyright (C) 2010, 2011 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <stdlib.h>
#include <string.h>

#include "obj_lists.h"


source_list_t *init_sources(config_t * cfg, int *n_sources)
{
    int i;
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
	if (!strcmp(type, "uniform point source"))
	    new_source =
		source_alloc(source_uniform_point_source, this_source,
			     NULL);
	else if (!strcmp(type, "sphere"))
	    new_source = source_alloc(source_sphere, this_source, NULL);
	else if (!strcmp(type, "spot source"))
	    new_source = source_alloc(source_spot, this_source, NULL);
	else {
	    fprintf(stderr,
		    "Unknown source type (%s) found. Ignoring source %s\n",
		    type, name);
	    continue;
	}

	new_entry = (source_list_t *) malloc(sizeof(source_list_t));
	new_entry->s = new_source;
	list_add_tail(&new_entry->list, &s_list->list);

    }
    return s_list;
}

target_list_t *init_targets(config_t * cfg, int *n_targets,
			    const char *file_mode)
{
    int i;
    target_list_t *t_list;
    const config_setting_t *targets = config_lookup(cfg, "targets");

    *n_targets = config_setting_length(targets);

    t_list = (target_list_t *) malloc(sizeof(target_list_t));
    INIT_LIST_HEAD(&t_list->list);

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
	if (!strcmp(type, "one-sided plane screen"))
	    new_target =
		target_alloc(target_plane_screen_one_sided, this_target,
			     NULL, file_mode);
	else if (!strcmp(type, "two-sided plane screen"))
	    new_target =
		target_alloc(target_plane_screen_two_sided, this_target,
			     NULL, file_mode);
	else if (!strcmp(type, "rectangle"))
	    new_target =
		target_alloc(target_rectangle, this_target, NULL,
			     file_mode);
	else if (!strcmp(type, "triangle"))
	    new_target =
		target_alloc(target_triangle, this_target, NULL,
			     file_mode);
	else if (!strcmp(type, "ellipsoid"))
	    new_target =
		target_alloc(target_ellipsoid, this_target, NULL,
			     file_mode);
	else if (!strcmp(type, "annulus"))
	    new_target =
		target_alloc(target_annulus, this_target, NULL, file_mode);
	else if (!strcmp(type, "disk"))
	    new_target =
		target_alloc(target_disk, this_target, NULL, file_mode);
	else {
	    fprintf(stderr,
		    "Unknown target type (%s) found. Ignoring target %s\n",
		    type, name);
	    continue;
	}

	new_entry = (target_list_t *) malloc(sizeof(target_list_t));
	new_entry->t = new_target;
	list_add_tail(&new_entry->list, &t_list->list);

	if (file_mode[0] == 'w') {	/* write header to new dump_file */
	    char string[256];
	    double *mat;
	    int idx;

	    snprintf(string, 256, "# %s (%s)\n",
		     get_target_name(new_target),
		     get_target_type(new_target));
	    dump_string(new_target, string);

	    snprintf(string, 256, "# transformation matrix M:\n");
	    dump_string(new_target, string);
	    snprintf(string, 256,
		     "#    g(x, y, z) = M l(x, y, z) + o(x, y, z))\n");
	    dump_string(new_target, string);
	    snprintf(string, 256,
		     "#    l(x, y, z) = MT (g(x, y, z) - o(x, y, z))\n");
	    dump_string(new_target, string);

	    mat = M(new_target);
	    for (idx = 0; idx < 3; idx++) {
		const int offset = 3 * idx;

		snprintf(string, 256, "#   %g\t%g\t%g\n", mat[offset],
			 mat[offset + 1], mat[offset + 2]);
		dump_string(new_target, string);
	    }
	    snprintf(string, 256,
		     "#\n# x\ty\t[z]\tpower\tlambda\t\t(z component is missing for plane targets!)\n#\n");
	    dump_string(new_target, string);
	}

    }
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

/*	obj_list.c
 *
 * Copyright (C) 2010 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <stdlib.h>
#include <string.h>

#include "obj_lists.h"

source_list_t *init_sources(config_t *cfg)
{
    int i;
    source_list_t *s_list;
    const config_setting_t *sources = config_lookup(cfg, "sources");
    const int n_sources = config_setting_length(sources);

    s_list = (source_list_t *) malloc(sizeof(source_list_t));
    INIT_LIST_HEAD(&s_list->list);

    for (i = 0; i < n_sources; ++i) {	/* iterate through all sources */
	source_list_t *new_entry;
	source_t *new_source;
	config_setting_t *this_source;
	const char *name;	/* name (identifier) of 'this_source' */
	const char *type;	/* type (identifier) of 'this_source' */

	this_source = config_setting_get_elem(sources, (unsigned int) i);
	config_setting_lookup_string(this_source, "name", &name);

	/*
	 * allocate the correct source type.
	 * because the 'type' listed in the configuration is a string
	 * and not a C-type we need the 'if() else if() else()' ladder
	 * going over all possibilities and using the correct (first)
	 * argument in the call to 'source_alloc()'. 
	 */
	config_setting_lookup_string(this_source, "type", &type);
	if (strstr(type, "uniform point source"))
	    new_source =
		source_alloc(source_uniform_point_source, cfg, name);
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

target_list_t *init_targets(config_t *cfg, int *n_targets)
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
	config_setting_t *this_target;
	const char *name;	/* name (identifier) of 'this_target' */
	const char *type;	/* type (identifier) of 'this_target' */

	this_target = config_setting_get_elem(targets, (unsigned int) i);
	config_setting_lookup_string(this_target, "name", &name);

	/*
	 * allocate the correct target type.
	 * because the 'type' listed in the configuration is a string
	 * and not a C-type we need the 'if() else if() else()' ladder
	 * going over all possibilities and using the correct (first)
	 * argument in the call to 'target_alloc()'. 
	 */
	config_setting_lookup_string(this_target, "type", &type);
	if (strstr(type, "one-sided plane screen"))
	    new_target =
		target_alloc(target_plane_screen_one_sided, cfg, name);
	else if (strstr(type, "two-sided plane screen"))
	    new_target =
		target_alloc(target_plane_screen_two_sided, cfg, name);
	else {
	    fprintf(stderr,
		    "Unknown target type (%s) found. Ignoring target %s\n",
		    type, name);
	    continue;
	}

	new_entry = (target_list_t *) malloc(sizeof(target_list_t));
	new_entry->t = new_target;
	list_add_tail(&new_entry->list, &t_list->list);

    }
    return t_list;
}

void source_list_free(source_list_t *s)
{
    struct list_head *pos, *pos_t;

    list_for_each_safe(pos, pos_t, &s->list) {
	source_list_t *this = list_entry(pos, source_list_t, list);

	source_free(this->s);
	list_del(pos);
	free(this);
    }
}


void target_list_free(target_list_t *t)
{
    struct list_head *pos, *pos_t;

    list_for_each_safe(pos, pos_t, &t->list) {
	target_list_t *this = list_entry(pos, target_list_t, list);

	target_free(this->t);
	list_del(pos);
	free(this);
    }
}

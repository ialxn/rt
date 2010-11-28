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

source_list_t *init_sources(config_t * cfg)
{
    source_list_t *s_list;
    const config_setting_t *s = config_lookup(cfg, "sources");
    const int count = config_setting_length(s);

    int i;

    s_list = (source_list_t *) malloc(sizeof(source_list_t));
    INIT_LIST_HEAD(&(*s_list).list);

    for (i = 0; i < count; ++i) {
	source_list_t *new_elem;
	source_t *new_s;
	config_setting_t *this_s;
	const char *name;

	this_s = config_setting_get_elem(s, (unsigned int) i);
	config_setting_lookup_string(this_s, "name", &(name));

	new_s = source_alloc(source_uniform_point_source, cfg, name);

	new_elem = (source_list_t *) malloc(sizeof(source_list_t));
	new_elem->s = new_s;
	list_add_tail(&(new_elem->list), &(*s_list).list);

    }

    return (s_list);
}

target_list_t *init_targets(config_t * cfg)
{
    target_list_t *t_list;
    const config_setting_t *t = config_lookup(cfg, "targets");
    const int count = config_setting_length(t);

    int i;

    t_list = (target_list_t *) malloc(sizeof(target_list_t));
    INIT_LIST_HEAD(&(*t_list).list);

    for (i = 0; i < count; ++i) {
	target_list_t *new_elem;
	target_t *new_t;
	config_setting_t *this_t;
	const char *name;

	new_t = (target_t *) malloc(sizeof(target_t));

	this_t = config_setting_get_elem(t, (unsigned int) i);
	config_setting_lookup_string(this_t, "name", &(name));

	new_t->name = strdup(name);
	new_elem = (target_list_t *) malloc(sizeof(target_list_t));
	new_elem->t = new_t;
	list_add_tail(&(new_elem->list), &(*t_list).list);

    }

    return (t_list);
}

void source_list_free(source_list_t * s)
{
    struct list_head *pos, *pos_t;

    list_for_each_safe(pos, pos_t, &(*s).list) {
	source_list_t *this = list_entry(pos, source_list_t, list);

	source_free(this->s);
	list_del(pos);

	free(this);

    }
}


void target_list_free(target_list_t * t)
{
    struct list_head *pos, *pos_t;

    list_for_each_safe(pos, pos_t, &(*t).list) {
	target_list_t *this = list_entry(pos, target_list_t, list);

	free(this->t->name);
	free(this->t);
	list_del(pos);

	free(this);

    }
}

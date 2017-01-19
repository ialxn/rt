/*	likely.h
 *
 * Copyright (C) 2015,2016,2017 Ivo Alxneit
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
#ifndef __LIKELY_H__
#define __LIKELY_H__


#ifdef __GNUC__
#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)
#else
#define likely(x)       !!(x)
#define unlikely(x)     !!(x)
#endif


#endif				/* __LIKELY_H__ */

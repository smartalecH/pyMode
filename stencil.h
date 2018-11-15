
/*
    wgms3d - a full-vectorial finite-difference mode solver.

    Copyright (C) 2005-2012  Michael Krause <m.krause@tu-harburg.de>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _STENCIL_H
#define _STENCIL_H

/*  This file defines the discretization stencil. Most of the stuff in
 *  "#if 0" here is experimental. In particular, support for stencils
 *  larger than nine points is not fully functional right now. */

/* NDIRS >= NDO ! */

struct direction {
    double dr;
    double dz;
};

#define DIRS5 \
{ \
    { 0, n }, \
    { e, n }, \
    { e, 0 }, \
    { 0, -s }, \
    { -w, 0 } \
};

#define DIRS8 \
{ \
    { 0, n }, \
    { e, n }, \
    { e, 0 }, \
    { e, -s }, \
    { 0, -s }, \
    { -w, -s }, \
    { -w, 0 }, \
    { -w, n } \
};

#define DIRS24 \
{ \
    { 0, n }, \
    { e, n }, \
    { e, 0 }, \
    { e, -s }, \
    { 0, -s }, \
    { -w, -s }, \
    { -w, 0 }, \
    { -w, n }, \
    { 0, nn }, \
    { e, nn }, \
    { ee, nn }, \
    { ee, n }, \
    { ee, 0 }, \
    { ee, -s }, \
    { ee, -ss }, \
    { e, -ss }, \
    { 0, -ss }, \
    { -w, -ss }, \
    { -ww, -ss }, \
    { -ww, -s }, \
    { -ww, 0 }, \
    { -ww, n }, \
    { -ww, nn }, \
    { -w, nn } \
};

#define DIRS16 \
{ \
    { 0, n }, \
    { e, n }, \
    { e, 0 }, \
    { e, -s }, \
    { 0, -s }, \
    { -w, -s }, \
    { -w, 0 }, \
    { -w, n }, \
    { 0, nn }, \
    { ee, nn }, \
    { ee, 0 }, \
    { ee, -ss }, \
    { 0, -ss }, \
    { -ww, -ss }, \
    { -ww, 0 }, \
    { -ww, nn } \
};

#define DIRS16b \
{ \
    { 0, n }, \
    { e, n }, \
    { e, 0 }, \
    { e, -s }, \
    { 0, -s }, \
    { -w, -s }, \
    { -w, 0 }, \
    { -w, n }, \
    { 0, nn }, \
    { ee, nn }, \
    { ee, 0 }, \
    { ee, -ss }, \
    { 0, -ss }, \
    { -ww, -ss }, \
    { -ww, 0 }, \
    { -ww, nn },				\
    { e, nn },				\
    { ee, -s },				\
};

/* The point listed in dirs[i] corresponds to stencil point number dirs_to_stencil_map[i]. */

#define DIRS5MAP { 0, 1, 2, 3, 5, 7 }

#define DIRS8MAP { 0, 1, 2, 3, 4, 5, 6, 7, 8 }

#define DIRS24MAP { 0, 1, 2, 3, 4, 5, 6, 7, 8, \
	    9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 }

#define DIRS16MAP { 0, 1, 2, 3, 4, 5, 6, 7, 8, \
	    9, 11, 13, 15, 17, 19, 21, 23 }

#define DIRS16bMAP { 0, 1, 2, 3, 4, 5, 6, 7, 8, \
	    9, 11, 13, 15, 17, 19, 21, 23, 10, 14 }


/* ------------------------------------------------------------------------------------------- */

#if 0 /* second-order, five-point stencil (P + N + E + S + W + NE) */
#define NUM_GHOST_POINTS 1
#define ENABLE_STANDARD_DIFFOPS 0
#define NDO 5
#define NSP 9
#define DIRS DIRS5
#define NDIRS 5
static const int dirs_to_stencil_map[NDIRS+1] = DIRS5MAP;
#endif

/* ------------------------------------------------------------------------------------------- */

#if 1 /* working second-order, nine-point stencil version */
#define NUM_GHOST_POINTS 1
#define ENABLE_STANDARD_DIFFOPS 1
#define NDO 5
#define NSP 9
#define DIRS DIRS8
#define NDIRS 8
static const int dirs_to_stencil_map[NDIRS+1] = DIRS8MAP;
#endif

/* ------------------------------------------------------------------------------------------- */

#if 0 /* ??? second-order, twenty-five-point stencil */
#define NUM_GHOST_POINTS 2
#define ENABLE_STANDARD_DIFFOPS 0
#define NDO 5
#define NSP 25
#define DIRS DIRS24
#define NDIRS 24
static const int dirs_to_stencil_map[NDIRS+1] = DIRS24MAP;
#endif

/* ------------------------------------------------------------------------------------------- */

#if 0 /* "third-order" (ansatz only), twenty-five-point stencil version */
#define NUM_GHOST_POINTS 2
#define ENABLE_STANDARD_DIFFOPS 0
#define NDO 9
#define NSP 25
#define DIRS DIRS24
#define NDIRS 24
static const int dirs_to_stencil_map[NDIRS+1] = DIRS24MAP;
#endif

/* ------------------------------------------------------------------------------------------- */

#if 0 /* fourth-order version, 25-point stencil, BUGGY. */
#define NUM_GHOST_POINTS 2
#define ENABLE_STANDARD_DIFFOPS 0
#define NDO 14
#define NSP 25
#define DIRS DIRS24
#define NDIRS 24
static const int dirs_to_stencil_map[NDIRS+1] = DIRS24MAP;
#endif

/* ------------------------------------------------------------------------------------------- */

#if 0
#define NUM_GHOST_POINTS 2
#define ENABLE_STANDARD_DIFFOPS 0
#define NDO 14
#define NSP 25
#define DIRS DIRS16
#define NDIRS 16
static const int dirs_to_stencil_map[NDIRS+1] = DIRS16MAP;
#endif

#if 0
#define NUM_GHOST_POINTS 2
#define ENABLE_STANDARD_DIFFOPS 0
#define NDO 14
#define NSP 25
#define DIRS DIRS16b
#define NDIRS 18
static const int dirs_to_stencil_map[NDIRS+1] = DIRS16bMAP;
#endif

#if 0
#define NUM_GHOST_POINTS 2
#define ENABLE_STANDARD_DIFFOPS 0
#define NDO 5
#define NSP 25
#define DIRS DIRS16
#define NDIRS 16
static const int dirs_to_stencil_map[NDIRS+1] = DIRS16MAP;
#endif

#if 0
#define NUM_GHOST_POINTS 2
#define ENABLE_STANDARD_DIFFOPS 0
#define NDO 5
#define NSP 25
#define DIRS DIRS16b
#define NDIRS 18
static const int dirs_to_stencil_map[NDIRS+1] = DIRS16bMAP;
#endif

/* For a given field component, the value of the field, e.g. H^r, is
 * considered as well as NDO derivatives D_i (1 <= i <= NDO).
 *
 * D_0 = id
 *
 * D_1 = \partial / \partial r
 * D_2 = \partial / \partial z
 * D_3 = \partial / \partial rr
 * D_4 = \partial / \partial rz
 * D_5 = \partial / \partial zz
 *
 * D_6 = \partial / \partial rrr
 * D_7 = \partial / \partial rrz
 * D_8 = \partial / \partial rzz
 * D_9 = \partial / \partial zzz
 *
 * D_10 = \partial / \partial rrrr
 * D_11 = \partial / \partial rrrz
 * D_12 = \partial / \partial rrzz
 * D_13 = \partial / \partial rzzz
 * D_14 = \partial / \partial zzzz
 *
 */

/*

  The vector of derivatives (vector f in the paper) has the following
  12 components:

       H^r              0
       H^r_r            1
       H^r_z            2
       H^r_rr           3
       H^r_rz           4
       H^r_zz           5

       H^z              6
       H^z_r            7
       H^z_z            8
       H^z_rr           9
       H^z_rz          10
       H^z_zz          11

*/

#define NF (NDO+1)
#define NV (2*NF)

#endif

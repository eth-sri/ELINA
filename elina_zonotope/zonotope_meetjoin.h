/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
 *  This software is distributed under GNU Lesser General Public License Version 3.0.
 *  For more information, see the ELINA project website at:
 *  http://elina.ethz.ch
 *
 *  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
 *  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
 *  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
 *  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 *  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY     
 *  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
 *  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
 *  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
 *  CONTRACT, TORT OR OTHERWISE).
 *
 */


#ifndef _ZONOTOPE_MEETJOIN_H_
#define _ZONOTOPE_MEETJOIN_H_

#include "zonotope.h"
#include "zonotope_internal.h"


#ifdef __cplusplus
extern "C" {
#endif

/* Meet and Join */
/*****************/
/* 1.Meet */
zonotope_t* zonotope_meet(elina_manager_t* man, bool destructive, zonotope_t* z1, zonotope_t* z2);
zonotope_t* zonotope_meet_array(elina_manager_t* man, zonotope_t** tab, size_t size);
zonotope_t* zonotope_meet_lincons_array(elina_manager_t* man,
		bool destructive,
		zonotope_t* z,
		elina_lincons0_array_t* array);
zonotope_t* zonotope_meet_tcons_array(elina_manager_t* man,
		bool destructive,
		zonotope_t* z,
		elina_tcons0_array_t* array);

/* 2.Join */
zonotope_t* zonotope_join(elina_manager_t* man, bool destructive, zonotope_t* z1, zonotope_t* z2);
zonotope_t* zonotope_join_array(elina_manager_t* man, zonotope_t** tab, size_t size);

#ifdef __cplusplus
}
#endif

#endif

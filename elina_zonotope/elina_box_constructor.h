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

#ifndef _ELINA_BOX_CONSTRUCTOR_H_
#define _ELINA_BOX_CONSTRUCTOR_H_

#include "elina_box_internal.h"
#include "elina_box_representation.h"
#include "elina_linearize.h"
#include "elina_box_meetjoin.h"

#ifdef __cplusplus
extern "C" {
#endif

elina_box_t* elina_box_bottom(elina_manager_t* man, size_t intdim, size_t realdim);
elina_box_t* elina_box_top(elina_manager_t* man, size_t intdim, size_t realdim);
elina_box_t* elina_box_of_box(elina_manager_t* man,
		  size_t intdim, size_t realdim,
		  elina_interval_t** tinterval);
elina_dimension_t elina_box_dimension(elina_manager_t* man, elina_box_t* a);
bool elina_box_is_bottom(elina_manager_t* man, elina_box_t* a);
bool elina_box_is_top(elina_manager_t* man, elina_box_t* a);
bool elina_box_is_leq(elina_manager_t* man, elina_box_t* a, elina_box_t* b);
bool elina_box_is_eq(elina_manager_t* man, elina_box_t* a, elina_box_t* b);
elina_interval_t* elina_box_bound_dimension(elina_manager_t* man,
				   elina_box_t* a, elina_dim_t dim);
elina_interval_t* elina_box_bound_linexpr(elina_manager_t* man,
				 elina_box_t* a, elina_linexpr0_t* expr);

elina_lincons0_array_t elina_box_to_lincons_array(elina_manager_t* man, elina_box_t* a);

#ifdef __cplusplus
}
#endif

#endif

/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
 *  This software is distributed under GNU Lesser General Public License
 * Version 3.0. For more information, see the ELINA project website at:
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

#include <string.h>
#include <stdio.h>

#include "elina_box_internal.h"


void elina_box_internal_init(elina_box_internal_t* intern)
{
  intern->bound_linexpr_internal_itv = elina_interval_alloc();
  intern->bound_linexpr_internal_itv2 = elina_interval_alloc();
  intern->bound_linexpr_itv= elina_interval_alloc();
  intern->meet_lincons_internal_itv = elina_interval_alloc();
  intern->meet_lincons_internal_itv2 = elina_interval_alloc();
  intern->meet_lincons_internal_itv3 = elina_interval_alloc();
  intern->meet_lincons_internal_bound = elina_scalar_alloc();
}
void elina_box_internal_clear(elina_box_internal_t* intern)
{
  elina_interval_free(intern->bound_linexpr_internal_itv);
  elina_interval_free(intern->bound_linexpr_internal_itv2);
  elina_interval_free(intern->bound_linexpr_itv);
  elina_interval_free(intern->meet_lincons_internal_itv);
  elina_interval_free(intern->meet_lincons_internal_itv2);
  elina_interval_free(intern->meet_lincons_internal_itv3);
  elina_scalar_free(intern->meet_lincons_internal_bound);
}

elina_box_internal_t* elina_box_internal_alloc(void)
{
  elina_box_internal_t* intern = malloc(sizeof(elina_box_internal_t));
  elina_box_internal_init(intern);
  return intern;
}
void elina_box_internal_free(elina_box_internal_t* intern)
{
  elina_box_internal_clear(intern);
  free(intern);
}

elina_manager_t* elina_box_manager_alloc(void)
{
  size_t i;
  elina_box_internal_t* itv;
  elina_manager_t* man;
  void** funptr;

  if (!elina_fpu_init()) {
    fprintf(stderr,"elina_box_manager_alloc cannot change the FPU rounding mode\n");
  }

  itv = elina_box_internal_alloc();
  man = elina_manager_alloc("elina box", "1.0 with DOUBLE ",
			 itv, (void (*)(void*))elina_box_internal_free);
  funptr = man->funptr;

  funptr[ELINA_FUNID_COPY] = &elina_box_copy;
  funptr[ELINA_FUNID_FREE] = &elina_box_free;
  //funptr[ELINA_FUNID_ASIZE] = &elina_box_size;
  //funptr[ELINA_FUNID_MINIMIZE] = &elina_box_minimize;
  //funptr[ELINA_FUNID_CANONICALIZE] = &elina_box_canonicalize;
  //funptr[ELINA_FUNID_HASH] = &elina_box_hash;
  //funptr[ELINA_FUNID_APPROXIMATE] = &elina_box_approximate;
  funptr[ELINA_FUNID_FPRINT] = &elina_box_fprint;
  //funptr[ELINA_FUNID_FPRINTDIFF] = &elina_box_fprintdiff;
  //funptr[ELINA_FUNID_FDUMP] = &elina_box_fdump;
  //funptr[ELINA_FUNID_SERIALIZE_RAW] = &elina_box_serialize_raw;
  //funptr[ELINA_FUNID_DESERIALIZE_RAW] = &elina_box_deserialize_raw;
  funptr[ELINA_FUNID_BOTTOM] = &elina_box_bottom;
  funptr[ELINA_FUNID_TOP] = &elina_box_top;
  funptr[ELINA_FUNID_OF_BOX] = &elina_box_of_box;
  funptr[ELINA_FUNID_DIMENSION] = &elina_box_dimension;
  funptr[ELINA_FUNID_IS_BOTTOM] = &elina_box_is_bottom;
  funptr[ELINA_FUNID_IS_TOP] = &elina_box_is_top;
  funptr[ELINA_FUNID_IS_LEQ] = &elina_box_is_leq;
  funptr[ELINA_FUNID_IS_EQ] = &elina_box_is_eq;
  //funptr[ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED] = &elina_box_is_dimension_unconstrained;
  //funptr[ELINA_FUNID_SAT_INTERVAL] = &elina_box_sat_interval;
  //funptr[ELINA_FUNID_SAT_LINCONS] = &elina_box_sat_lincons;
  //funptr[ELINA_FUNID_SAT_TCONS] = &elina_box_sat_tcons;
  funptr[ELINA_FUNID_BOUND_DIMENSION] = &elina_box_bound_dimension;
  funptr[ELINA_FUNID_BOUND_LINEXPR] = &elina_box_bound_linexpr;
  // funptr[ELINA_FUNID_BOUND_TEXPR] = &elina_box_bound_texpr;
  // funptr[ELINA_FUNID_TO_BOX] = &elina_box_to_box;
  // funptr[ELINA_FUNID_TO_LINCONS_ARRAY] = &elina_box_to_lincons_array;
  // funptr[ELINA_FUNID_TO_TCONS_ARRAY] = &elina_box_to_tcons_array;
  // funptr[ELINA_FUNID_TO_GENERATOR_ARRAY] = &elina_box_to_generator_array;
  // funptr[ELINA_FUNID_MEET] = &elina_box_meet;
  // funptr[ELINA_FUNID_MEET_ARRAY] = &elina_box_meet_array;
  funptr[ELINA_FUNID_MEET_LINCONS_ARRAY] = &elina_box_meet_lincons_array;
  //funptr[ELINA_FUNID_MEET_TCONS_ARRAY] = &elina_box_meet_tcons_array;
  funptr[ELINA_FUNID_JOIN] = &elina_box_join;
  // funptr[ELINA_FUNID_JOIN_ARRAY] = &elina_box_join_array;
  // funptr[ELINA_FUNID_ADD_RAY_ARRAY] = &elina_box_add_ray_array;
  // funptr[ELINA_FUNID_ASSIGN_LINEXPR_ARRAY] = &elina_box_assign_linexpr_array;
  // funptr[ELINA_FUNID_SUBSTITUTE_LINEXPR_ARRAY] =
  // &elina_box_substitute_linexpr_array; funptr[ELINA_FUNID_ASSIGN_TEXPR_ARRAY]
  // = &elina_box_assign_texpr_array; funptr[ELINA_FUNID_SUBSTITUTE_TEXPR_ARRAY]
  // = &elina_box_substitute_texpr_array;
  funptr[ELINA_FUNID_ADD_DIMENSIONS] = &elina_box_add_dimensions;
  funptr[ELINA_FUNID_REMOVE_DIMENSIONS] = &elina_box_remove_dimensions;
  // funptr[ELINA_FUNID_PERMUTE_DIMENSIONS] = &elina_box_permute_dimensions;
  // funptr[ELINA_FUNID_FORGET_ARRAY] = &elina_box_forget_array;
  // funptr[ELINA_FUNID_EXPAND] = &elina_box_expand;
  // funptr[ELINA_FUNID_FOLD] = &elina_box_fold;
  // funptr[ELINA_FUNID_WIDENING] = &elina_box_widening;
  // funptr[ELINA_FUNID_CLOSURE] = &elina_box_closure;

  for (i=0; i<ELINA_EXC_SIZE; i++){
    elina_manager_set_abort_if_exception(man, i, false);
  }
  return man;
}

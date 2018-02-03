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

#include <stdlib.h>
#include <unistd.h>

#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>


#include "zonotope_internal.h"
//#include "zonotope_itv_utils.h"

#include "zonotope_representation.h"
#include "zonotope_constructor.h"
#include "zonotope_meetjoin.h"
#include "zonotope_assign.h"
#include "zonotope_resize.h"
#include "zonotope_otherops.h"

double zonotope_copy_time=0;
double zonotope_is_equal_time=0;
double zonotope_is_lequal_time=0;
double zonotope_permute_dimension_time=0;
double zonotope_add_dimension_time=0;
double zonotope_remove_dimension_time=0;
double zonotope_top_time=0;
double zonotope_bottom_time=0;
double zonotope_join_time=0;
double zonotope_free_time=0;
double zonotope_forget_array_time=0;
double zonotope_meet_lincons_time=0;
double zonotope_to_box_time=0;
double zonotope_of_box_time=0;
double zonotope_is_top_time=0;
double zonotope_is_bottom_time=0;
double zonotope_assign_linexpr_time=0;

elina_manager_t* zonotope_manager_alloc(void)
{
	//CALL();
	elina_manager_t* man;
	void** funptr;

	
	elina_manager_t* manNS = elina_box_manager_alloc();
	zonotope_internal_t *zt = zonotope_internal_alloc(manNS);

	man = elina_manager_alloc("Zonotope",/* Library name */
			"1.0", /* version */
			zt, /* internal structure */
			(void (*)(void*))zonotope_internal_free /* free function for internal */
			);

	funptr = man->funptr;

	/* Internal representation */
	/***************************/
	/* 1.Memory */
	funptr[ELINA_FUNID_COPY] = &zonotope_copy;
	funptr[ELINA_FUNID_FREE] = &zonotope_free;
	/*funptr[ELINA_FUNID_SIZE] = &zonotope_size;*/
	//funptr[ELINA_FUNID_ASIZE] = &zonotope_size;
	/* 2.Control of internal representation */
	//funptr[ELINA_FUNID_MINIMIZE] = &zonotope_minimize;
	//funptr[ELINA_FUNID_CANONICALIZE] = &zonotope_canonicalize;
	//funptr[ELINA_FUNID_HASH] = &zonotope_hash;
	//funptr[ELINA_FUNID_APPROXIMATE] = &zonotope_approximate;
	/* 3.Printing */
	funptr[ELINA_FUNID_FPRINT] = &zonotope_fprint;
	//funptr[ELINA_FUNID_FPRINTDIFF] = &zonotope_fprintdiff;
	//funptr[ELINA_FUNID_FDUMP] = &zonotope_fdump;
	/* 4.Serialisation */
	//funptr[ELINA_FUNID_SERIALIZE_RAW] = &zonotope_serialize_raw;
	//funptr[ELINA_FUNID_DESERIALIZE_RAW] = &zonotope_deserialize_raw;

	/* Constructors */
	/****************/
	/* 1.Basic constructors */
	funptr[ELINA_FUNID_BOTTOM] = &zonotope_bottom;
	funptr[ELINA_FUNID_TOP] = &zonotope_top;
	funptr[ELINA_FUNID_OF_BOX] = &zonotope_of_box;
	/* 2.Accessors */
	funptr[ELINA_FUNID_DIMENSION] = &zonotope_dimension;
	/* 3.Tests */
	funptr[ELINA_FUNID_IS_BOTTOM] = &zonotope_is_bottom;
	funptr[ELINA_FUNID_IS_TOP] = &zonotope_is_top;
	//funptr[ELINA_FUNID_IS_LEQ] = &zonotope_is_leq;
	funptr[ELINA_FUNID_IS_EQ] = &zonotope_is_eq;
	//funptr[ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED] = &zonotope_is_dimension_unconstrained;
	//funptr[ELINA_FUNID_SAT_TCONS] = &zonotope_sat_tcons; /*  */
	//funptr[ELINA_FUNID_SAT_INTERVAL] = &zonotope_sat_interval;
	//funptr[ELINA_FUNID_SAT_LINCONS] = &zonotope_sat_lincons;
	/* 4.Extraction of properties */
        // funptr[ELINA_FUNID_BOUND_TEXPR] = &zonotope_bound_texpr; /*  */
        // funptr[ELINA_FUNID_BOUND_DIMENSION] = &zonotope_bound_dimension;
        // funptr[ELINA_FUNID_BOUND_LINEXPR] = &zonotope_bound_linexpr;
        funptr[ELINA_FUNID_TO_BOX] = &zonotope_to_box;
	//funptr[ELINA_FUNID_TO_TCONS_ARRAY] = &zonotope_to_tcons_array; /*  */
	funptr[ELINA_FUNID_TO_LINCONS_ARRAY] = &zonotope_to_lincons_array;
	//funptr[ELINA_FUNID_TO_GENERATOR_ARRAY] = &zonotope_to_generator_array;

	/* Meet and Join */
	/*****************/
	/* 1.Meet */
	//funptr[ELINA_FUNID_MEET] = &zonotope_meet; /* */
	//funptr[ELINA_FUNID_MEET_ARRAY] = &zonotope_meet_array; /*  */
	funptr[ELINA_FUNID_MEET_LINCONS_ARRAY] = &zonotope_meet_lincons_array; /*  */
	//funptr[ELINA_FUNID_MEET_TCONS_ARRAY] = &zonotope_meet_tcons_array; /*  */
	/* 2.Join */
	funptr[ELINA_FUNID_JOIN] = &zonotope_join;
	//funptr[ELINA_FUNID_JOIN_ARRAY] = &zonotope_join_array;

	//funptr[ELINA_FUNID_ADD_RAY_ARRAY] = &zonotope_add_ray_array;

	/* Assign and Substitute */
	/*************************/
	funptr[ELINA_FUNID_ASSIGN_LINEXPR_ARRAY] = &zonotope_assign_linexpr_array;
	//funptr[ELINA_FUNID_SUBSTITUTE_LINEXPR_ARRAY] = &zonotope_substitute_linexpr_array;
	//funptr[ELINA_FUNID_ASSIGN_TEXPR_ARRAY] = &zonotope_assign_texpr_array;
	//funptr[ELINA_FUNID_SUBSTITUTE_TEXPR_ARRAY] = &zonotope_substitute_texpr_array;

	/* Resize dimensions */
	/*********************/
	funptr[ELINA_FUNID_ADD_DIMENSIONS] = &zonotope_add_dimensions;
	funptr[ELINA_FUNID_REMOVE_DIMENSIONS] = &zonotope_remove_dimensions;
	funptr[ELINA_FUNID_PERMUTE_DIMENSIONS] = &zonotope_permute_dimensions;

	/* Other functions */
	/*******************/
	funptr[ELINA_FUNID_FORGET_ARRAY] = &zonotope_forget_array;
	//funptr[ELINA_FUNID_EXPAND] = &zonotope_expand;
	//funptr[ELINA_FUNID_FOLD] = &zonotope_fold;
	//funptr[ELINA_FUNID_WIDENING] = &zonotope_widening;
	//funptr[ELINA_FUNID_CLOSURE] = &zonotope_closure;

	/* ?! */
	/******/
	/*
	   ELINA_FUNID_CHANGE_ENVIRONMENT
	   ELINA_FUNID_RENAME_ARRAY
	   ELINA_FUNID_SIZE2
	 */
	man->option.abort_if_exception[ELINA_EXC_INVALID_ARGUMENT] = false;
	return man;
}

/* back pointer to our internal structure from the manager */
zonotope_internal_t* zonotope_init_from_manager(elina_manager_t* man, elina_funid_t funid)
{
  // printf("funid: %d %d\n",funid,funid==ELINA_FUNID_UNKNOWN);
  // fflush(stdout);
  zonotope_internal_t *pr = (zonotope_internal_t *)man->internal;
  pr->funid = funid;
  if (!(pr->man))
    pr->man = man;
  return pr;
}

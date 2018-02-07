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

#ifndef _ELINA_BOX_H_
#define _ELINA_BOX_H_

#include <stdio.h>
#if defined (HAS_APRON)
#include "apron_wrapper.h"
#else
#include "elina_scalar.h"
#include "elina_interval.h"
#include "elina_dimension.h"
#include "elina_manager.h"
#include "elina_linexpr0.h"
#include "elina_lincons0.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* OCaml interface */
elina_manager_t* elina_box_manager_alloc(void);

typedef struct elina_box_t elina_box_t;

/* ********************************************************************** */
/* I. General management */
/* ********************************************************************** */

/* ============================================================ */
/* 1 Memory */
/* ============================================================ */

elina_box_t* elina_box_copy(elina_manager_t* man, elina_box_t* a);
  /* Return a copy of a box value, on
     which destructive update does not affect the initial value. */

void elina_box_free(elina_manager_t* man, elina_box_t* a);
  /* Free all the memory used by the box value */


/* ============================================================ */
/* 2 Printing */
/* ============================================================ */

void elina_box_fprint(FILE* stream,
	       elina_manager_t* man,
	       elina_box_t* a,
	       char** name_of_dim);
  /* Print the box value in a pretty way, using function
     name_of_dim to name dimensions */

/* ============================================================ */
/* 3 Basic constructors */
/* ============================================================ */

/* We assume that dimensions [0..intdim-1] correspond to integer variables, and
   dimensions [intdim..intdim+realdim-1] to real variables */

elina_box_t* elina_box_bottom(elina_manager_t* man, size_t intdim, size_t realdim);
  /* Create a bottom (empty) value */

elina_box_t* elina_box_top(elina_manager_t* man, size_t intdim, size_t realdim);
  /* Create a top (universe) value */

elina_box_t* elina_box_of_box(elina_manager_t* man,
		  size_t intdim, size_t realdim,
		  elina_interval_t** tinterval);
  /* Abstract an hypercube defined by the array of intervals
     of size intdim+realdim */

/* ============================================================ */
/* 4 Accessors */
/* ============================================================ */

elina_dimension_t elina_box_dimension(elina_manager_t* man, elina_box_t* a);
  /* Return the dimensionality of the box values */

/* ============================================================ */
/* 5 Tests */
/* ============================================================ */

bool elina_box_is_bottom(elina_manager_t* man, elina_box_t* a);
bool elina_box_is_top(elina_manager_t* man, elina_box_t* a);

bool elina_box_is_leq(elina_manager_t* man, elina_box_t* a1, elina_box_t* a2);
  /* inclusion check */

bool elina_box_is_eq(elina_manager_t* man, elina_box_t* a1, elina_box_t* a2);
  /* equality check */


/* ============================================================ */
/* 6 Extraction of properties */
/* ============================================================ */

elina_interval_t* elina_box_bound_dimension(elina_manager_t* man,
				elina_box_t* a, elina_dim_t dim);
  /* Returns the interval taken by a linear expression
     over the box value */

elina_interval_t* elina_box_bound_linexpr(elina_manager_t* man,
				 elina_box_t* a, elina_linexpr0_t* expr);
  /* Returns the interval taken by a linear expression
     over the box value */
elina_lincons0_array_t elina_box_to_lincons_array(elina_manager_t* man, elina_box_t* a);

/* ============================================================ */
/* 7 Meet and Join */
/* ============================================================ */

elina_box_t* elina_box_meet(elina_manager_t* man, bool destructive, elina_box_t* a1, elina_box_t* a2);
    
elina_box_t* elina_box_join(elina_manager_t* man, bool destructive, elina_box_t* a1, elina_box_t* a2);
  /* Join of 2 box values */


elina_box_t* elina_box_meet_lincons_array(elina_manager_t* man,
			      bool destructive,
			      elina_box_t* a,
			      elina_lincons0_array_t* array);
  /* Meet of an box value with a set of constraints */


/* ============================================================ */
/* 8 Change and permutation of dimensions */
/* ============================================================ */

elina_box_t* elina_box_add_dimensions(elina_manager_t* man,
			  bool destructive,
			  elina_box_t* a,
			  elina_dimchange_t* dimchange,
			  bool project);
elina_box_t* elina_box_remove_dimensions(elina_manager_t* man,
			     bool destructive,
			     elina_box_t* a,
			     elina_dimchange_t* dimchange);

/* ============================================================ */
/* 9 Assignment */
/* ============================================================ */

elina_box_t* elina_box_assign_linexpr_array(elina_manager_t* man,
                                    bool destructive,
                                    elina_box_t* a,
                                    elina_dim_t* tdim,
                                    elina_linexpr0_t** texpr,
                                    size_t size,
                                    elina_box_t* dest);

/* ============================================================ */
/* 10 Projection */
/* ============================================================ */

elina_box_t* elina_box_forget_array(elina_manager_t* man,
                            bool destructive,
                            elina_box_t* a, elina_dim_t* tdim, size_t size,
                            bool project);

/* ============================================================ */
/* 11 Widening */
/* ============================================================ */
elina_box_t* elina_box_widening(elina_manager_t* man,
                                elina_box_t* a1, elina_box_t* a2);

#ifdef __cplusplus
}
#endif

#endif

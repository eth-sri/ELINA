/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
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


/* ********************************************************************** */
/* opt_pk.h: Interface of the ELINA linear relation library  */
/* ********************************************************************** */


#ifndef __OPT_PK_H__
#define __OPT_PK_H__


#include "elina_int.h"
#include "elina_rat.h"
#include "comp_list.h"


#ifdef __cplusplus
extern "C" {
#endif

#if defined (HAS_APRON)
#include "apron_wrapper.h"
//#include "num.h"
//#include "numint.h"
//#include "numrat.h"
//#include "bound.h"

//#include "itv.h"
//#include "itv_linexpr.h"
//#include "itv_linearize.h"

#else
#include "elina_coeff.h"
#include "elina_dimension.h"
#include "elina_linexpr0.h"
#include "elina_texpr0.h"
#include "elina_lincons0.h"
#include "elina_tcons0.h"
#include "elina_manager.h"

#endif

#include "elina_generic.h"
#include "elina_linearize_texpr.h"


/* The invariant of the representation of a polyhedron is the following: if the
   polyhedron is empty, then C==F==satC==satF==0. Otherwise, we have
   (C || F) && (satC || satF || !(C && F)).
   This means that a non-empty polyhedron has a minimal representation minimal
   if and only if C && F if and only if satC || satF. */

typedef enum opt_pk_status_t {
  opt_pk_status_conseps=0x1,
  opt_pk_status_consgauss=0x2,
  opt_pk_status_gengauss=0x4,
  opt_pk_status_minimaleps=0x8
} opt_pk_status_t;


struct opt_pk_t {
  /* private data: do not use directly ! */
  struct opt_matrix_t* C;
  struct opt_matrix_t* F;
  struct opt_satmat_t* satC;
  struct opt_satmat_t* satF;
  unsigned short int intdim;
  unsigned short int realdim;
  size_t nbeq;
  size_t nbline;
  bool is_minimized;
  opt_pk_status_t status;
};

typedef struct opt_pk_t opt_pk_t;
typedef struct opt_pk_internal_t opt_pk_internal_t;

typedef struct opt_pk_array_t {
	opt_pk_t ** poly;
	array_comp_list_t * acl;
	unsigned short int maxcols;
	bool is_bottom;
}opt_pk_array_t;




/* ============================================================ */
/* A. Constructor for ELINA manager (to be freed with elina_manager_free). */
/* ============================================================ */

elina_manager_t* opt_pk_manager_alloc(bool strict);
  /* Allocate a ELINA manager for convex polyhedra.

     If the Boolean parameter is true, abstract values generated with the
     manager can have strict constraints (like x>0). Otherwise they are defined
     using only loose constraints. Managers and abstract values in strict or
     loose mode are incompatible.
  */

/* ============================================================ */
/* B. Options */
/* ============================================================ */

opt_pk_internal_t* opt_pk_manager_get_internal(elina_manager_t* man);

/* For setting options when one has a elina_manager_t object, one can use the
   ELINA function elina_manager_get_internal with a cast. */

void opt_pk_set_approximate_max_coeff_size(opt_pk_internal_t* opk, size_t size);
//void opt_pk_print(elina_manager_t* man, opt_pk_t* po, char** name_of_dim);

/* ============================================================ */
/* D. Conversions */
/* ============================================================ */


//elina_abstract0_t* opt_pk_to_abstract0(elina_manager_t* man, opt_pk_t* poly);
  

/* ============================================================ */
/* D. Constructor and destructor for internal manager */
/* ============================================================ */

/* Allocates pk and initializes it with a default size */
struct opt_pk_internal_t* opt_pk_internal_alloc(bool strict);
/* Clear and free pk */
void opt_pk_internal_free(opt_pk_internal_t* pk);

/* ********************************************************************** */
/* I. General management */
/* ********************************************************************** */


/* ============================================================ */
/* I.1 Memory */
/* ============================================================ */

opt_pk_array_t* opt_pk_copy(elina_manager_t* man, opt_pk_array_t* o);
  /* Return a copy of an abstract value, on
     which destructive update does not affect the initial value. */

void opt_pk_free(elina_manager_t* man, opt_pk_array_t* o);
  /* Free all the memory used by the abstract value */

size_t opt_pk_size(elina_manager_t* man, opt_pk_array_t* o);
  /* Return the abstract size of a polyhedron, which is the number of
     coefficients of its current representation, possibly redundant. */


/* ============================================================ */
/* I.2 Control of internal representation */
/* ============================================================ */

void opt_pk_minimize(elina_manager_t* man, opt_pk_array_t* o);
  /* Minimize the size of the representation of a.
     This may result in a later recomputation of internal information.
  */

void opt_pk_array_canonicalize(elina_manager_t* man, opt_pk_array_t* o);
  /* Put the polyhedron with minimized constraints and frames.  If in addition
     the integer man->option->canonicalize.algorithm is strictly positive,
     normalize equalities and lines, and also strict constraints */


/* ============================================================ */
/* I.3 Printing */
/* ============================================================ */

void opt_pk_array_fprint(FILE* stream,
	       elina_manager_t* man,
	       opt_pk_array_t* o,
	       char** name_of_dim);
  /* Print the abstract value in a pretty way, using function
     name_of_dim to name dimensions */

elina_membuf_t opt_pk_array_serialize_raw(elina_manager_t* man, opt_pk_array_t* a);

opt_pk_array_t* opt_pk_array_deserialize_raw(elina_manager_t* man, void* p, size_t* size);

/* ********************************************************************** */
/* II. Constructor, accessors, tests and property extraction */
/* ********************************************************************** */

/* ============================================================ */
/* II.1 Basic constructors */
/* ============================================================ */

/* We assume that dimensions [0..intdim-1] correspond to integer variables, and
   dimensions [intdim..intdim+realdim-1] to real variables */

opt_pk_array_t* opt_pk_bottom(elina_manager_t* man, size_t intdim, size_t realdim);
  /* Create a bottom (empty) value */

opt_pk_array_t* opt_pk_top(elina_manager_t* man, size_t intdim, size_t realdim);
  /* Create a top (universe) value */


opt_pk_array_t* opt_pk_of_box(elina_manager_t* man,
		size_t intdim, size_t realdim,
		elina_interval_t** tinterval);
  /* Abstract an hypercube defined by the array of intervals
     of size intdim+realdim */

/* ============================================================ */
/* II.2 Accessors */
/* ============================================================ */

elina_dimension_t opt_pk_dimension(elina_manager_t* man, opt_pk_array_t* o);
/* Return the total number of dimensions of the abstract values */

/* ============================================================ */
/* II.3 Tests */
/* ============================================================ */

bool opt_pk_is_bottom(elina_manager_t* man, opt_pk_array_t* o);
  /* Emptiness test
     algorithm >= 0: strict behaviour, compute canonical form if necessary
     algorithm < 0: lazy behaviour, always cheap
  */
bool opt_pk_is_top(elina_manager_t* man, opt_pk_array_t* o);
  /* Universe test
     algorithm >= 0: strict behaviour, compute canonical form if necessary
     algorithm < 0: lazy behaviour, always cheap
  */

bool opt_pk_is_leq(elina_manager_t* man, opt_pk_array_t* o1, opt_pk_array_t* o2);
  /* Inclusion test:
     Is always strict
     algorithm > 0: (nearly always) compute canonical forms
     algorithm <= 0: compute dual representations only if necessary
  */

bool opt_pk_is_eq(elina_manager_t* man, opt_pk_array_t* o1, opt_pk_array_t* o2);
  /* Equality test:
     Is always strict
     Use algorithm field of is_leq.
  */

bool opt_pk_sat_lincons(elina_manager_t* man, opt_pk_array_t* oa, elina_lincons0_t* lincons);
  /* Satisfiability of a linear constraint
     Is always strict
     algorithm > 0: (nearly always) compute canonical form
     algorithm <= 0: compute dual representation only if necessary
  */

bool opt_pk_sat_tcons(elina_manager_t* man, opt_pk_array_t* oa, elina_tcons0_t* cons);
  /* Satisfiability of a tree expression constraint. */

  /* Inclusion of a dimension in an interval
     Is always strict
     algorithm > 0: (nearly always) compute canonical form
     algorithm <= 0: compute dual representation only if necessary
  */

bool opt_pk_is_dimension_unconstrained(elina_manager_t* man, opt_pk_array_t* oa,
				   elina_dim_t dim);
  /* Is a dimension unconstrained ?
     Is always strict
     algorithm > 0: compute canonical form
     algorithm <= 0: compute dual representation only if necessary
  */

/* ============================================================ */
/* II.4 Extraction of properties */
/* ============================================================ */

elina_interval_t* opt_pk_bound_linexpr(elina_manager_t* man,
				opt_pk_array_t* o, elina_linexpr0_t* expr);
  /* Returns the interval taken by a linear expression
     over the abstract value.

     algorithm > 0: compute canonical form
     algorithm <= 0: compute dual representation only if necessary
  */

//elina_interval_t* opt_pk_bound_texpr(elina_manager_t* man,
//			      opt_pk_t* o, elina_texpr0_t* expr);
  /* Returns the interval taken by a tree expression
     over the abstract value. */

elina_interval_t* opt_pk_bound_dimension(elina_manager_t* man,
				  opt_pk_array_t* o, elina_dim_t dim);
  /* Returns the interval taken by the dimension
     over the abstract value

     algorithm > 0: compute canonical form
     algorithm <= 0: compute dual representation only if necessary
  */

elina_lincons0_array_t opt_pk_to_lincons_array(elina_manager_t* man, opt_pk_array_t* o);
  /* Converts an abstract value to a polyhedra
     (conjunction of linear constraints).

     Always consider canonical form */

//elina_tcons0_array_t opt_pk_to_tcons_array(elina_manager_t* man, opt_pk_t* a);
  /* Converts an abstract value to a
     conjunction of tree expressions constraints. */

elina_interval_t** opt_pk_to_box(elina_manager_t* man, opt_pk_array_t* o);
  /* Converts an abstract value to an interval/hypercube.
     The size of the resulting array is pk_dimension(man,a).  This
     function can be reimplemented by using pk_bound_linexpr

     algorithm >= 0: compute canonical form
     algorithm < 0: compute dual representation only if necessary
  */



/* ********************************************************************** */
/* III. Operations */
/* ********************************************************************** */

/* ============================================================ */
/* III.1 Meet and Join */
/* ============================================================ */

opt_pk_array_t* opt_pk_meet(elina_manager_t* man, bool destructive, opt_pk_array_t* o1, opt_pk_array_t* o2);

opt_pk_array_t* opt_pk_join(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, opt_pk_array_t* ob);
  /* Meet and Join of 2 abstract values */

//opt_pk_t* opt_pk_meet_array(elina_manager_t* man, opt_pk_t** tab, size_t size);

//opt_pk_t* opt_pk_join_array(elina_manager_t* man, opt_pk_t** tab, size_t size);
  /* Meet and Join of a non empty array of abstract values. */

opt_pk_array_t* opt_pk_meet_lincons_array(elina_manager_t* man,
			    bool destructive, opt_pk_array_t* o,
			    elina_lincons0_array_t* array);
opt_pk_array_t* opt_pk_meet_tcons_array(elina_manager_t* man,
			  bool destructive, opt_pk_array_t* oa,
			  elina_tcons0_array_t* array);
  /* Meet of an abstract value with a set of constraints. */


  /* Generalized time elapse operator */

/* ============================================================ */
/* III.2 Assignement and Substitutions */
/* ============================================================ */
opt_pk_array_t* opt_pk_assign_linexpr_array(elina_manager_t* man,
			      bool destructive, opt_pk_array_t* oa,
			      elina_dim_t* tdim,
			      elina_linexpr0_t** texpr,
			      size_t size,
			      opt_pk_array_t* dest);

opt_pk_array_t* opt_pk_substitute_linexpr_array(elina_manager_t* man,
				  bool destructive, opt_pk_array_t* oa,
				  elina_dim_t* tdim,
				  elina_linexpr0_t** texpr,
				  size_t size,
				  opt_pk_array_t* dest);

opt_pk_array_t* opt_pk_assign_texpr_array(elina_manager_t* man,
			    bool destructive, opt_pk_array_t* oa,
			    elina_dim_t* tdim,
			    elina_texpr0_t** texpr,
			    size_t size,
			    opt_pk_array_t* dest);

  /* Parallel Assignement and Substitution of several dimensions by interval
     expressons. */

/* ============================================================ */
/* III.3 Projections */
/* ============================================================ */

opt_pk_array_t* opt_pk_forget_array(elina_manager_t* man,
		      bool destructive, opt_pk_array_t* o,
		      elina_dim_t* tdim, size_t size,
		      bool project);

/* ============================================================ */
/* III.4 Change and permutation of dimensions */
/* ============================================================ */

opt_pk_array_t* opt_pk_add_dimensions(elina_manager_t* man,
			bool destructive, opt_pk_array_t* o,
			elina_dimchange_t* dimchange,
			bool project);

opt_pk_array_t* opt_pk_remove_dimensions(elina_manager_t* man,
			   bool destructive, opt_pk_array_t* o,
			   elina_dimchange_t* dimchange);

opt_pk_array_t* opt_pk_permute_dimensions(elina_manager_t* man,
			    bool destructive,
			    opt_pk_array_t* o,
			    elina_dimperm_t* permutation);

void remove_block_and_factor(opt_pk_array_t *op, comp_list_t *cl);

/* ============================================================ */
/* III.5 Expansion and folding of dimensions */
/* ============================================================ */

opt_pk_array_t* opt_pk_expand(elina_manager_t* man,
		bool destructive, opt_pk_array_t* o,
		elina_dim_t dim,
		size_t n);
  /* Expand the dimension dim into itself + n additional dimensions.
     It results in (n+1) unrelated dimensions having same
     relations with other dimensions. The (n+1) dimensions are put as follows:

     - original dimension dim

     - if the dimension is integer, the n additional dimensions are put at the
       end of integer dimensions; if it is real, at the end of the real
       dimensions.
  */

opt_pk_array_t* opt_pk_fold(elina_manager_t* man,
	      bool destructive, opt_pk_array_t* o,
	      elina_dim_t* tdim,
	      size_t size);
  /* Fold the dimensions in the array tdim of size n>=1 and put the result
     in the first dimension in the array. The other dimensions of the array
     are then removed (using pk_permute_remove_dimensions). */

/* ============================================================ */
/* III.6 Widening */
/* ============================================================ */

/* Widening */

opt_pk_array_t* opt_pk_widening(elina_manager_t* man, opt_pk_array_t* o1, opt_pk_array_t* o2);

/* ============================================================ */
/* III.7 Closure operation */
/* ============================================================ */

/* Returns the topological closure of a possibly opened abstract value */




#ifdef __cplusplus
}
#endif

#endif

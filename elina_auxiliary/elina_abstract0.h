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

/* ************************************************************************* */
/* elina_abstract0.h: generic operations on numerical abstract values */
/* ************************************************************************* */


#ifndef _ELINA_ABSTRACT0_H_
#define _ELINA_ABSTRACT0_H_

typedef struct elina_abstract0_t elina_abstract0_t;

#include "elina_manager.h"
#include "elina_texpr0.h"
#include "elina_tcons0.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Generic abstract value at level 0 */
struct elina_abstract0_t {
  void* value;       /* Abstract value of the underlying library */
  elina_manager_t* man; /* Used to identify the effective type of value */
};

/* ********************************************************************** */
/* I. General management */
/* ********************************************************************** */

/* ============================================================ */
/* I.1 Memory */
/* ============================================================ */

elina_abstract0_t* elina_abstract0_copy(elina_manager_t* man, elina_abstract0_t* a);
  /* Return a copy of an abstract value, on
     which destructive update does not affect the initial value. */

void elina_abstract0_free(elina_manager_t* man, elina_abstract0_t* a);
  /* Free all the memory used by the abstract value */

size_t elina_abstract0_size(elina_manager_t* man, elina_abstract0_t* a);
  /* Return the abstract size of an abstract value (see elina_manager_t) */


/* ============================================================ */
/* I.2 Control of internal representation */
/* ============================================================ */


void elina_abstract0_minimize(elina_manager_t* man, elina_abstract0_t* a);
  /* Minimize the size of the representation of a.
     This may result in a later recomputation of internal information.
  */

void elina_abstract0_canonicalize(elina_manager_t* man, elina_abstract0_t* a);
  /* Put the abstract value in canonical form. (not yet clear definition) */

int elina_abstract0_hash(elina_manager_t* man, elina_abstract0_t* a);
  /* Return an hash value for the abstract value.  Two abstract values in
     canonical from (according to @code{elina_abstract1_canonicalize}) and
     considered as equal by the function elina_abstract0_is_eq should be given the
     same hash value (this implies more or less a canonical form).
  */


void elina_abstract0_approximate(elina_manager_t* man, elina_abstract0_t* a, int algorithm);
  /* Perform some transformation on the abstract value, guided by the
     field algorithm.

     The transformation may lose information.  The argument "algorithm"
     overrides the field algorithm of the structure of type elina_funopt_t
     associated to elina_abstract0_approximate (commodity feature). */

/* ============================================================ */
/* I.3 Printing */
/* ============================================================ */


void elina_abstract0_fprint(FILE* stream,
			 elina_manager_t* man, elina_abstract0_t* a, char** name_of_dim);
  /* Print the abstract value in a pretty way, using function name_of_dim to
     name dimensions. If name_of_dim==NULL, use x0, x1, ... */

void elina_abstract0_fprintdiff(FILE* stream,
			     elina_manager_t* man,
			     elina_abstract0_t* a1, elina_abstract0_t* a2,
			     char** name_of_dim);
  /* Print the difference between a1 (old value) and a2 (new value),
     using function name_of_dim to name dimensions.
     The meaning of difference is library dependent. */


void elina_abstract0_fdump(FILE* stream, elina_manager_t* man, elina_abstract0_t* a);
  /* Dump the internal representation of an abstract value,
     for debugging purposes */

/* ============================================================ */
/* I.4 Serialization */
/* ============================================================ */


elina_membuf_t elina_abstract0_serialize_raw(elina_manager_t* man, elina_abstract0_t* a);
  /* Allocate a memory buffer (with malloc), output the abstract value in raw
     binary format to it and return a pointer on the memory buffer and the size
     of bytes written.  It is the user responsibility to free the memory
     afterwards (with free). */


elina_abstract0_t* elina_abstract0_deserialize_raw(elina_manager_t* man, void* ptr, size_t* size);
  /* Return the abstract value read in raw binary format from the input stream
     and store in size the number of bytes read */

/* ********************************************************************** */
/* II. Constructor, accessors, tests and property extraction */
/* ********************************************************************** */

/* ============================================================ */
/* II.1 Basic constructors */
/* ============================================================ */

/* We assume that dimensions [0..intdim-1] correspond to integer variables, and
   dimensions [intdim..intdim+realdim-1] to real variables */

elina_abstract0_t* elina_abstract0_bottom(elina_manager_t* man, size_t intdim, size_t realdim);
  /* Create a bottom (empty) value */


elina_abstract0_t* elina_abstract0_top(elina_manager_t* man, size_t intdim, size_t realdim);
  /* Create a top (universe) value */

elina_abstract0_t* elina_abstract0_of_box(elina_manager_t* man,
				    size_t intdim, size_t realdim,
				    elina_interval_t** tinterval);
  /* Abstract an hypercube defined by the array of intervals
     of size intdim+realdim */

/* ============================================================ */
/* II.2 Accessors */
/* ============================================================ */

elina_dimension_t elina_abstract0_dimension(elina_manager_t* man, elina_abstract0_t* a);
  /* Return the dimensionality of the abstract values */

/* ============================================================ */
/* II.3 Tests */
/* ============================================================ */

/* In abstract tests,

   - true means that the predicate is certainly true.

   - false means by default don't know (an exception has occurred, or the exact
     computation was considered too expensive to be performed).

     However, if the flag exact in the manager is true, then false means really
     that the predicate is false.
*/

bool elina_abstract0_is_bottom(elina_manager_t* man, elina_abstract0_t* a);
bool elina_abstract0_is_top(elina_manager_t* man, elina_abstract0_t* a);


bool elina_abstract0_is_leq(elina_manager_t* man, elina_abstract0_t* a1, elina_abstract0_t* a2);
  /* inclusion check */
bool elina_abstract0_is_eq(elina_manager_t* man, elina_abstract0_t* a1, elina_abstract0_t* a2);
  /* equality check */

bool elina_abstract0_sat_lincons(elina_manager_t* man, elina_abstract0_t* a, elina_lincons0_t* lincons);
bool elina_abstract0_sat_tcons(elina_manager_t* man, elina_abstract0_t* a, elina_tcons0_t* tcons);
  /* does the abstract value satisfy the constraint ? */
bool elina_abstract0_sat_interval(elina_manager_t* man, elina_abstract0_t* a,
			      elina_dim_t dim, elina_interval_t* interval);
  /* is the dimension included in the interval in the abstract value ? */

bool elina_abstract0_is_dimension_unconstrained(elina_manager_t* man,
					     elina_abstract0_t* a, elina_dim_t dim);
  /* is the dimension unconstrained in the abstract value ?  If it is the case,
     we have forget(man,a,dim) == a */

/* ============================================================ */
/* II.4 Extraction of properties */
/* ============================================================ */

elina_interval_t* elina_abstract0_bound_linexpr(elina_manager_t* man,
					  elina_abstract0_t* a, elina_linexpr0_t* expr);
elina_interval_t* elina_abstract0_bound_texpr(elina_manager_t* man,
					elina_abstract0_t* a, elina_texpr0_t* expr);
  /* Returns the interval taken by the expression
     over the abstract value */


elina_interval_t* elina_abstract0_bound_dimension(elina_manager_t* man,
					    elina_abstract0_t* a, elina_dim_t dim);
  /* Returns the interval taken by the dimension
     over the abstract value */


elina_lincons0_array_t elina_abstract0_to_lincons_array(elina_manager_t* man, elina_abstract0_t* a);
  /* Converts an abstract value to a polyhedra
     (conjunction of linear constraints).

     The constraints are normally guaranteed to be really linear (without intervals) */

elina_tcons0_array_t elina_abstract0_to_tcons_array(elina_manager_t* man, elina_abstract0_t* a);
  /* Converts an abstract value to conjunction of tree expressions constraints.

     The constraints are normally guaranteed to be scalar (without intervals) */

elina_interval_t** elina_abstract0_to_box(elina_manager_t* man, elina_abstract0_t* a);
  /* Converts an abstract value to an interval/hypercube.
     The size of the resulting array is elina_abstract0_dimension(man,a).  This
     function can be reimplemented by using elina_abstract0_bound_linexpr */




/* ********************************************************************** */
/* III. Operations */
/* ********************************************************************** */

/* ============================================================ */
/* III.1 Meet and Join */
/* ============================================================ */


elina_abstract0_t* elina_abstract0_meet(elina_manager_t* man,
				  bool destructive, elina_abstract0_t* a1, elina_abstract0_t* a2);

elina_abstract0_t* elina_abstract0_join(elina_manager_t* man,
				  bool destructive, elina_abstract0_t* a1, elina_abstract0_t* a2);
  /* Meet and Join of 2 abstract values */


elina_abstract0_t* elina_abstract0_meet_array(elina_manager_t* man,
					elina_abstract0_t** tab, size_t size);

elina_abstract0_t* elina_abstract0_join_array(elina_manager_t* man,
					elina_abstract0_t** tab, size_t size);
  /* Meet and Join of an array of abstract values.
     Raises an [[exc_invalid_argument]] exception if [[size==0]]
     (no way to define the dimensionality of the result in such a case */


elina_abstract0_t*
elina_abstract0_meet_lincons_array(elina_manager_t* man,
				bool destructive, elina_abstract0_t* a, elina_lincons0_array_t* array);
elina_abstract0_t*
elina_abstract0_meet_tcons_array(elina_manager_t* man,
				bool destructive, elina_abstract0_t* a, elina_tcons0_array_t* array);
  /* Meet of an abstract value with a set of constraints */


/* ============================================================ */
/* III.2 Assignment and Substitutions */
/* ============================================================ */

elina_abstract0_t*
elina_abstract0_assign_linexpr_array(elina_manager_t* man,
				  bool destructive,
				  elina_abstract0_t* org,
				  elina_dim_t* tdim, elina_linexpr0_t** texpr, size_t size,
				  elina_abstract0_t* dest);
elina_abstract0_t*
elina_abstract0_assign_texpr_array(elina_manager_t* man,
				bool destructive,
				elina_abstract0_t* org,
				elina_dim_t* tdim, elina_texpr0_t** texpr, size_t size,
				elina_abstract0_t* dest);
elina_abstract0_t*
elina_abstract0_substitute_linexpr_array(elina_manager_t* man,
				      bool destructive,
				      elina_abstract0_t* org,
				      elina_dim_t* tdim, elina_linexpr0_t** texpr, size_t size,
				      elina_abstract0_t* dest);
elina_abstract0_t*
elina_abstract0_substitute_texpr_array(elina_manager_t* man,
				    bool destructive,
				    elina_abstract0_t* org,
				    elina_dim_t* tdim, elina_texpr0_t** texpr, size_t size,
				    elina_abstract0_t* dest);

  /* Parallel Assignment and Substitution of several dimensions by
     linear/tree expressions in abstract value org.

     dest is an optional argument. If not NULL, semantically speaking,
     the result of the transformation is intersected with dest. This is
     useful for precise backward transformations in lattices like intervals or
     octagons. */

/* ============================================================ */
/* III.3 Projections */
/* ============================================================ */

elina_abstract0_t*
elina_abstract0_forget_array(elina_manager_t* man,
			  bool destructive,
			  elina_abstract0_t* a, elina_dim_t* tdim, size_t size,
			  bool project);

/* ============================================================ */
/* III.4 Change and permutation of dimensions */
/* ============================================================ */

elina_abstract0_t*
elina_abstract0_add_dimensions(elina_manager_t* man,
			    bool destructive,
			    elina_abstract0_t* a,elina_dimchange_t* dimchange,
			    bool project);
elina_abstract0_t*
elina_abstract0_remove_dimensions(elina_manager_t* man,
			       bool destructive,
			       elina_abstract0_t* a, elina_dimchange_t* dimchange);
  /* Size of the permutation is supposed to be equal to
     the dimension of the abstract value */
elina_abstract0_t*
elina_abstract0_permute_dimensions(elina_manager_t* man,
				bool destructive,
				elina_abstract0_t* a, elina_dimperm_t* perm);

/* ============================================================ */
/* III.5 Expansion and folding of dimensions */
/* ============================================================ */


elina_abstract0_t*
elina_abstract0_expand(elina_manager_t* man,
		    bool destructive,
		    elina_abstract0_t* a, elina_dim_t dim, size_t n);
  /* Expand the dimension dim into itself + n additional dimensions.
     It results in (n+1) unrelated dimensions having same
     relations with other dimensions. The (n+1) dimensions are put as follows:

     - original dimension dim

     - if the dimension is integer, the n additional dimensions are put at the
       end of integer dimensions; if it is real, at the end of the real
       dimensions.
  */


elina_abstract0_t*
elina_abstract0_fold(elina_manager_t* man,
		  bool destructive,
		  elina_abstract0_t* a, elina_dim_t* tdim, size_t size);
  /* Fold the dimensions in the array tdim of size n>=1 and put the result
     in the first dimension in the array. The other dimensions of the array
     are then removed (using elina_abstract0_permute_remove_dimensions). */

/* ============================================================ */
/* III.6 Widening */
/* ============================================================ */

elina_abstract0_t* elina_abstract0_widening(elina_manager_t* man,
				      elina_abstract0_t* a1, elina_abstract0_t* a2);

/* ============================================================ */
/* III.7 Closure operation */
/* ============================================================ */

/* Returns the topological closure of a possibly opened abstract value */
elina_abstract0_t* elina_abstract0_closure(elina_manager_t* man, bool destructive, elina_abstract0_t* a);

/* ********************************************************************** */
/* IV. Functions offered by the ELINA interface */
/* ********************************************************************** */

/* These functions do not correspond to functions in the underlying library. */

static inline
elina_manager_t* elina_abstract0_manager(elina_abstract0_t* a)
  { return a->man; }
  /* Return a reference to the manager contained in the abstract value.
     The reference should not be freed */

elina_abstract0_t* elina_abstract0_of_lincons_array(elina_manager_t* man,
					      size_t intdim, size_t realdim,
					      elina_lincons0_array_t* array);
elina_abstract0_t* elina_abstract0_of_tcons_array(elina_manager_t* man,
					    size_t intdim, size_t realdim,
					    elina_tcons0_array_t* array);
  /* Abstract a conjunction of tree constraints */

elina_abstract0_t* elina_abstract0_assign_linexpr(elina_manager_t* man,
					    bool destructive,
					    elina_abstract0_t* org,
					    elina_dim_t dim, elina_linexpr0_t* expr,
					    elina_abstract0_t* dest);
elina_abstract0_t* elina_abstract0_assign_texpr(elina_manager_t* man,
					  bool destructive,
					  elina_abstract0_t* org,
					  elina_dim_t dim, elina_texpr0_t* expr,
					  elina_abstract0_t* dest);
elina_abstract0_t* elina_abstract0_substitute_linexpr(elina_manager_t* man,
						bool destructive,
						elina_abstract0_t* org,
						elina_dim_t dim, elina_linexpr0_t* expr,
						elina_abstract0_t* dest);
elina_abstract0_t* elina_abstract0_substitute_texpr(elina_manager_t* man,
					      bool destructive,
					      elina_abstract0_t* org,
					      elina_dim_t dim, elina_texpr0_t* expr,
					      elina_abstract0_t* dest);
  /* Assignment and Substitution of a single dimension by an expression in
     abstract value org.

     dest is an optional argument. If not NULL, semantically speaking,
     the result of the transformation is intersected with dest. This is
     useful for precise backward transformations in lattices like intervals or
     octagons.
 */

/* Applying an elina_dimchange2_t transformation (dimension adding followed by
   dimension removal/projection).  If project is true, the newly added
   dimension are projected on their 0-hyperplane. */
elina_abstract0_t* elina_abstract0_apply_dimchange2(elina_manager_t* man,
					      bool destructive,
					      elina_abstract0_t* a, elina_dimchange2_t* dimchange2,
					      bool project);

/* Widening with threshold */
elina_abstract0_t*
elina_abstract0_widening_threshold(elina_manager_t* man,
				elina_abstract0_t* a1, elina_abstract0_t* a2,
				elina_lincons0_array_t* array);

/* ********************************************************************** */
/* ********************************************************************** */
/* Internal functions */
/* ********************************************************************** */
/* ********************************************************************** */

/* ********************************************************************** */
/* 0. Utility and checking functions */
/* ********************************************************************** */

/* Constructor for elina_abstract0_t */

static inline
elina_abstract0_t* elina_abstract0_cons(elina_manager_t* man, void* value)
{
  elina_abstract0_t* res = (elina_abstract0_t*)malloc(sizeof(elina_abstract0_t));
  res->value = value;
  res->man = elina_manager_copy(man);
  return res;
}
static inline
void _elina_abstract0_free(elina_abstract0_t* a)
{
  void (*ptr)(elina_manager_t*,elina_abstract0_t*) = (void (*) (elina_manager_t*,elina_abstract0_t*))(a->man->funptr[ELINA_FUNID_FREE]);
  ptr(a->man,(elina_abstract0_t*)a->value);
  elina_manager_free(a->man);
  free(a);
}
static inline
elina_abstract0_t* elina_abstract0_cons2(elina_manager_t* man, bool destructive, elina_abstract0_t* oldabs, void* newvalue)
{
  if (destructive){
    if (oldabs->man != man){
      elina_manager_free(oldabs->man);
      oldabs->man = elina_manager_copy(man);
    }
    oldabs->value = newvalue;
    return oldabs;
  }
  else {
    return elina_abstract0_cons(man,newvalue);
  }
}

/* ====================================================================== */
/* 0.1 Checking typing w.r.t. manager */
/* ====================================================================== */

/*
  These functions return true if everything is OK, otherwise they raise an
  exception in the manager and return false.
*/

/* One abstract value */

void elina_abstract0_checkman1_raise(elina_funid_t funid, elina_manager_t* man, elina_abstract0_t* a);

static inline
bool elina_abstract0_checkman1(elina_funid_t funid, elina_manager_t* man, elina_abstract0_t* a)
{
  if (man->library != a->man->library){
    elina_abstract0_checkman1_raise(funid,man,a);
    return false;
  }
  else
    return true;
}

/* Two abstract values */
bool elina_abstract0_checkman2(elina_funid_t funid,
			    elina_manager_t* man, elina_abstract0_t* a1, elina_abstract0_t* a2);

/* Array of abstract values */
bool elina_abstract0_checkman_array(elina_funid_t funid,
				 elina_manager_t* man, elina_abstract0_t** tab, size_t size);

/* ====================================================================== */
/* 0.2 Checking compatibility of arguments: abstract values */
/* ====================================================================== */

/* Getting dimensions without checks */
static inline
elina_dimension_t _elina_abstract0_dimension(elina_abstract0_t* a)
{
  elina_dimension_t (*ptr)(elina_manager_t*,...) = (elina_dimension_t (*) (elina_manager_t*,...))(a->man->funptr[ELINA_FUNID_DIMENSION]);
  return ptr(a->man,a->value);
}

/* Check that the 2 abstract values have the same dimensionality */
bool elina_abstract0_check_abstract2(elina_funid_t funid, elina_manager_t* man,
				  elina_abstract0_t* a1, elina_abstract0_t* a2);

/* Check that the array of abstract values have the same dimensionality.*/
bool elina_abstract0_check_abstract_array(elina_funid_t funid, elina_manager_t* man,
				       elina_abstract0_t** tab, size_t size);

/* ====================================================================== */
/* 0.3 Checking compatibility of arguments: dimensions */
/* ====================================================================== */

/* Check that the dimension makes sense in the given dimensionality */
void elina_abstract0_check_dim_raise(elina_funid_t funid, elina_manager_t* man,
				  elina_dimension_t dimension, elina_dim_t dim,
				  const char* prefix);
static inline
bool elina_abstract0_check_dim(elina_funid_t funid, elina_manager_t* man,
			    elina_dimension_t dimension, elina_dim_t dim)
{
  if (dim>=dimension.intdim+dimension.realdim){
    elina_abstract0_check_dim_raise(funid,man,dimension,dim,
				 "incompatible dimension for the abstract value");
    return false;
  } else {
    return true;
  }
}

/* Check that the array of dimensions make sense in the given dimensionality */
bool elina_abstract0_check_dim_array(elina_funid_t funid, elina_manager_t* man,
				  elina_dimension_t dimension, elina_dim_t* tdim, size_t size);

/* ====================================================================== */
/* 0.4 Checking compatibility of arguments: expressions */
/* ====================================================================== */

void elina_abstract0_check_expr_raise(elina_funid_t funid, elina_manager_t* man,
				   elina_dimension_t dimension,
				   elina_dim_t dim,
				   char* prefix);

/* Check that the linear expression makes sense in the given dimensionality */
elina_dim_t elina_abstract0_check_linexpr_check(elina_dimension_t dimension,
					  elina_linexpr0_t* expr);
bool elina_abstract0_check_linexpr(elina_funid_t funid, elina_manager_t* man,
				elina_dimension_t dimension,
				elina_linexpr0_t* expr);

/* Check that the tree expression makes sense in the given dimensionality */
elina_dim_t elina_abstract0_check_texpr_check(elina_dimension_t dimension,
					elina_texpr0_t* expr);
bool elina_abstract0_check_texpr(elina_funid_t funid, elina_manager_t* man,
			      elina_dimension_t dimension,
			      elina_texpr0_t* expr);

/* ====================================================================== */
/* 0.5 Checking compatibility of arguments: array of expressions/constraints/generators */
/* ====================================================================== */

/* Check that array of linear expressions makes sense in the given dimensionality */
bool elina_abstract0_check_linexpr_array(elina_funid_t funid, elina_manager_t* man,
				      elina_dimension_t dimension,
				      elina_linexpr0_t** texpr, size_t size);

/* Check that array of linear constraint makes sense in the given dimensionality */
bool elina_abstract0_check_lincons_array(elina_funid_t funid, elina_manager_t* man,
				      elina_dimension_t dimension,
				      elina_lincons0_array_t* array);


/* Check that array of tree expressions makes sense in the given dimensionality */
bool elina_abstract0_check_texpr_array(elina_funid_t funid, elina_manager_t* man,
				    elina_dimension_t dimension,
				    elina_texpr0_t** texpr, size_t size);

/* Check that array of tree constraint makes sense in the given dimensionality */
bool elina_abstract0_check_tcons_array(elina_funid_t funid, elina_manager_t* man,
				    elina_dimension_t dimension,
				    elina_tcons0_array_t* array);


elina_abstract0_t*
elina_abstract0_meetjoin(elina_funid_t funid,
		      /* either meet or join */
		      elina_manager_t* man, bool destructive,
		      elina_abstract0_t* a1, elina_abstract0_t* a2);
elina_abstract0_t*
elina_abstract0_asssub_linexpr(elina_funid_t funid,
			    /* either assign or substitute */
			    elina_manager_t* man,
			    bool destructive,
			    elina_abstract0_t* a,
			    elina_dim_t dim, elina_linexpr0_t* expr,
			    elina_abstract0_t* dest);
elina_abstract0_t*
elina_abstract0_asssub_linexpr_array(elina_funid_t funid,
				  /* either assign or substitute */
				  elina_manager_t* man,
				  bool destructive,
				  elina_abstract0_t* a,
				  elina_dim_t* tdim, elina_linexpr0_t** texpr, size_t size,
				  elina_abstract0_t* dest);
elina_abstract0_t*
elina_abstract0_asssub_texpr(elina_funid_t funid,
			    /* either assign or substitute */
			  elina_manager_t* man,
			  bool destructive,
			  elina_abstract0_t* a,
			  elina_dim_t dim, elina_texpr0_t* expr,
			  elina_abstract0_t* dest);
elina_abstract0_t*
elina_abstract0_asssub_texpr_array(elina_funid_t funid,
				/* either assign or substitute */
				elina_manager_t* man,
				bool destructive,
				elina_abstract0_t* a,
				elina_dim_t* tdim, elina_texpr0_t** texpr, size_t size,
				elina_abstract0_t* dest);
#ifdef __cplusplus
}
#endif

#endif

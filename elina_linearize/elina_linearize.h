/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
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


#ifndef _ELINA_LINEARIZE_H_
#define _ELINA_LINEARIZE_H_

#include "elina_manager.h"
#include "elina_rat.h"
#include "elina_abstract0.h"
#include "elina_linexpr0_arith.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ********************************************************************** */
/* Auxiliary functions for add, multiplications etc. */
/* ********************************************************************** */

bool elina_interval_sqrt(elina_interval_t *dst, elina_interval_t *src, elina_scalar_discr_t discr);

void linearize_elina_lincons0_array(elina_lincons0_array_t* array, bool meet, elina_scalar_discr_t discr);

bool quasilinearize_elina_linexpr0(elina_linexpr0_t* expr, elina_interval_t** env, bool for_meet_inequality, elina_scalar_discr_t discr);

bool quasilinearize_elina_lincons0(elina_lincons0_t* cons, elina_interval_t** env, bool meet, elina_scalar_discr_t discr);

bool quasilinearize_elina_lincons0_array(elina_lincons0_array_t* array, elina_interval_t ** env, bool meet, elina_scalar_discr_t discr);

void elina_lincons0_reduce_integer(elina_lincons0_t* cons, size_t intdim, elina_scalar_discr_t discr);

char elina_lincons0_array_reduce_integer(elina_lincons0_array_t* array, size_t intdim, elina_scalar_discr_t discr);

void elina_lincons0_set_bool(elina_lincons0_t* cons, bool value, elina_scalar_discr_t discr);

char eval_elina_cstlincons0(elina_lincons0_t* cons);

void elina_lincons0_array_reinit(elina_lincons0_array_t* array, size_t size);

/* ********************************************************************** */
/* I. Evaluation of interval linear expressions */
/* ********************************************************************** */

bool elina_interval_eval_elina_linexpr0(elina_interval_t * itv, elina_linexpr0_t* expr, elina_interval_t** env, elina_scalar_discr_t discr);

elina_interval_t* eval_elina_linexpr0(elina_manager_t* man,
		 elina_abstract0_t* abs,
		 elina_linexpr0_t* expr,
		 elina_scalar_discr_t discr,
		 bool* pexact);

/* These functions are dedicated to implementors of domains. They offer generic
   default implementations for some of the operations required by the ELINA
   API, when there is no more specific and efficient implementation for the
   domain being implemented.

   To use them, the function allocating manager, which is specific to the domain,
   should put the corresponding pointers in the virtual table to them.

   They manipulated "unboxed" abstract values, which are native to the
   underlying library: they are not yet boxed with the manager in the type
   elina_abstract0_t.
*/

/* The following functions use the given abstract value for transforming
   interval linear expressions (resp. constraints, arrays of expressions,
   arrays od constraints) in quasilinear corresponding objects.

   They use to_box and dimension (and is_bottom if NDEBUG is undefined) generic
   functions.

   - discr allows to choose the type of scalars used for computations and for
     the result.

   - pexact is a pointer to a Boolean, which is set to true if all the
     conversions and computations were exact.

   For elina_quasilinearize_XXX functions, if the argument does not need any modification,
   then it is returned itself.

   Calling elina_linearize_linexpr0_array is more efficient than calling N times
   elina_linearize_linexpr0 because the conversion of abstract value to bounding
   boxes is done only once, as well as other internal allocations.
*/

/* ********************************************************************** */
/* II. Quasilinearization of interval linear expressions */
/* ********************************************************************** */

elina_linexpr0_t*
elina_quasilinearize_linexpr0(elina_manager_t* man,
			   void* abs,
			   elina_linexpr0_t* linexpr0,
			   bool* pexact,
			   elina_scalar_discr_t discr);

elina_lincons0_t
elina_quasilinearize_lincons0(elina_manager_t* man,
			   void* abs,
			   elina_lincons0_t* lincons0,
			   bool* pexact,
			   elina_scalar_discr_t discr,
			   bool meet);

elina_linexpr0_t**
elina_quasilinearize_linexpr0_array(elina_manager_t* man,
				 void* abs,
				 elina_linexpr0_t** texpr, size_t size,
				 bool* pexact,
				 elina_scalar_discr_t discr);

elina_lincons0_array_t
elina_quasilinearize_lincons0_array(elina_manager_t* man,
				 void* abs,
				 elina_lincons0_array_t* array,
				 bool* pexact,
				 elina_scalar_discr_t discr,
				 bool linearize,
				 bool meet);

/* ********************************************************************** */
/* III. Evaluation of tree expressions */
/* ********************************************************************** */

elina_interval_t* elina_eval_texpr0(elina_manager_t* man,
			      elina_abstract0_t* abs,
			      elina_texpr0_t* expr,
			      elina_scalar_discr_t discr,
			      bool* pexact);

/* ********************************************************************** */
/* IV. Interval linearization of linear tree expressions */
/* ********************************************************************** */

/* Linearize a tree expression that is (syntaxically) interval linear with
   exact arithmetic.

   Compared to elina_intlinearize_texpr0() function below, this functions does
   not require a bounding box for dimensions.

   If the precondition is violated, returns NULL.
*/


elina_linexpr0_t* elina_intlinearize_texpr0_intlinear(elina_manager_t* man,
						elina_texpr0_t* expr,
						elina_scalar_discr_t discr);

/* ********************************************************************** */
/* V. Interval linearization of tree expressions */
/* ********************************************************************** */

elina_linexpr0_t* elina_intlinearize_texpr0(elina_manager_t* man,
				      elina_abstract0_t* abs,
				      elina_texpr0_t* expr,
				      bool* pexact,
				      elina_scalar_discr_t discr,
				      bool quasilinearize);

elina_linexpr0_t** elina_intlinearize_texpr0_array(elina_manager_t* man,
					     elina_abstract0_t* abs,
					     elina_texpr0_t** texpr, size_t size,
					     bool* pexact,
					     elina_scalar_discr_t discr,
					     bool quasilinearize);

elina_lincons0_t elina_intlinearize_tcons0(elina_manager_t* man,
				     elina_abstract0_t* abs,
				     elina_tcons0_t* cons,
				     bool* pexact,
				     elina_scalar_discr_t discr,
				     bool quasilinearize, bool meet);

elina_lincons0_array_t elina_intlinearize_tcons0_array(elina_manager_t* man,
						 elina_abstract0_t* abs,
						 elina_tcons0_array_t* array,
						 bool* pexact,
						 elina_scalar_discr_t discr,
						 elina_linexpr_type_t type, bool meet,
						 bool boxize, size_t kmax, bool intervalonly);



#ifdef __cplusplus
}
#endif

#endif

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



#ifndef _ELINA_GENERIC_H_
#define _ELINA_GENERIC_H_

#if defined (HAS_APRON)
#include "apron_wrapper.h"
#else

#include "elina_manager.h"
#include "elina_scalar.h"
#include "elina_abstract0.h"
#endif

#include "elina_linearize.h"

#ifdef __cplusplus
extern "C" {
#endif

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

/* ********************************************************************** */
/* I. Constructors, accessors, tests and property extraction */
/* ********************************************************************** */

bool elina_generic_sat_tcons(elina_manager_t* man, void* abs, elina_tcons0_t* cons,
			  elina_scalar_discr_t discr, bool quasilinearize);
  /* This function implements a generic sat_tcons operation using
     elina_linearize_texpr0 and sat_lincons operations. */

elina_interval_t* elina_generic_bound_texpr(elina_manager_t* man, void* abs, elina_texpr0_t* expr,
				      elina_scalar_discr_t discr, bool quasilinearize);
  /* This function implements a generic bound_texpr operation using to_box and 
     elina_eval_texpr0 operations. */
  
elina_tcons0_array_t elina_generic_to_tcons_array(elina_manager_t* man,
					    void* abs);
  /* This function implements a generic to_tcons_array operation using
     to_lincons_array operation. */

/* ********************************************************************** */
/* II. Operations */
/* ********************************************************************** */

/* ============================================================ */
/*  Meet/Join on arrays of abstract values */
/* ============================================================ */

void* elina_generic_meetjoin_array(bool meet,
				elina_manager_t* man,
				void** tab, size_t size);
  /* This function implements a generic meet/join_array operation using copy and
     meet/join operations. */

static inline
void* elina_generic_meet_array(elina_manager_t* man,
			    void** tab, size_t size);
  /* This function implements a generic meet_array operation using copy and
     meet operations. */
static inline
void* elina_generic_join_array(elina_manager_t* man,
			    void** tab, size_t size);
  /* This function implements a generic join_array operation using copy and
     meet operations. */

/* ============================================================ */
/*  Meet with array of constraints */
/* ============================================================ */

void* elina_generic_meet_quasilinearize_lincons_array(elina_manager_t* man,
						   bool destructive, void* abs, elina_lincons0_array_t* array,
						   elina_scalar_discr_t discr, bool linearize,
						   void* (*meet_lincons_array)(elina_manager_t*, 
									       bool, void*,elina_lincons0_array_t*));

void*
elina_generic_meet_intlinearize_tcons_array(elina_manager_t* man,
					 bool destructive, void* abs, elina_tcons0_array_t* array,
					 elina_scalar_discr_t discr, elina_linexpr_type_t linearize,
					 void* (*meet_lincons_array)(elina_manager_t*,
								     bool, void*,
								     elina_lincons0_array_t*));

/* ============================================================ */
/*  Assignments/Substitutions */
/* ============================================================ */

void* elina_generic_asssub_linexpr_array(bool assign,
				      elina_manager_t* man,
				      bool destructive, void* abs, elina_dim_t* tdim, elina_linexpr0_t** texpr, size_t size,
				      void* dest);
void* elina_generic_asssub_texpr_array(bool assign,
				    elina_manager_t* man,
				    bool destructive, void* abs, elina_dim_t* tdim, elina_texpr0_t** texpr, size_t size,
				    void* dest);
  /*
    These functions implement generic parallel assignment/substitution
    operations by:
    1. introducing primed dimensions
    2. transforming linear expressions into equality constraints relating the
    assigned primed dimension and the linear expression
    If dest!=NULL
      3. introducing primed dimensions in dest
      4. exchanging primed and unprimed dimensions in dest
      5. intersecting the abstract value with the modified dest
    6. intersecting the obtained abstract value with the constraints
    7. exchanging primed and unprimed dimensions
    8. removing the introduced (primed) dimensions
    
   It relies on: is_bottom, copy, dimension, add_dimensions,
   permute_dimensions, remove_dimensions, meet_lincons_array/meet_tcons_array, meet and free
   abstract operations.
   
   Meaning of parameters:
   - assign selects the operation: true means assignment, false substitution
   - The other parameters have the meaning they have for parallel
     assignment/substitution
*/
static inline
void* elina_generic_assign_linexpr_array(elina_manager_t* man,
				      bool destructive, void* abs, elina_dim_t* tdim, elina_linexpr0_t** texpr, size_t size,
				      void* dest);
static inline
void* elina_generic_assign_texpr_array(elina_manager_t* man,
				    bool destructive, void* abs, elina_dim_t* tdim, elina_texpr0_t** texpr, size_t size,
				    void* dest);
  /*
     These functions implement generic parallel assignment operations by
     relying on is_bottom, copy, dimension, add_dimensions, permute_dimensions,
     remove_dimensions, meet_lincons_array or meet_tcons_array abstract
     operations.
  */
static inline
void* elina_generic_substitute_linexpr_array(elina_manager_t* man,
					  bool destructive, void* abs, elina_dim_t* tdim, elina_linexpr0_t** texpr, size_t size,
					  void* dest);
static inline
void* elina_generic_substitute_texpr_array(elina_manager_t* man,
					bool destructive, void* abs, elina_dim_t* tdim, elina_texpr0_t** texpr, size_t size,
					void* dest);
  /*
     These functions implement generic parallel assignment operations by
     relying on is_bottom, copy, dimension, add_dimensions, permute_dimensions,
     remove_dimensions, meet_lincons_array or meet_tcons_array abstract
     operations.
  */

/* ********************************************************************** */
/* III. Inline functions definitions */
/* ********************************************************************** */
static inline 
void* elina_generic_meet_array(elina_manager_t* man,
			    void** tab, size_t size)
{ return elina_generic_meetjoin_array(true,man,tab,size); }

static inline 
void* elina_generic_join_array(elina_manager_t* man,
			    void** tab, size_t size)
{ return elina_generic_meetjoin_array(false,man,tab,size); }

static inline
void* elina_generic_assign_linexpr_array(elina_manager_t* man,
				      bool destructive, void* abs, elina_dim_t* tdim, elina_linexpr0_t** texpr, size_t size,
				      void* dest)
{
  return elina_generic_asssub_linexpr_array(true,
					 man, destructive, abs, tdim, texpr, size,
					 dest);
}
static inline
void* elina_generic_substitute_linexpr_array(elina_manager_t* man,
					  bool destructive, void* abs, elina_dim_t* tdim, elina_linexpr0_t** texpr, size_t size,
					  void* dest)
{
  return elina_generic_asssub_linexpr_array(false,
					 man, destructive, abs, tdim, texpr, size,
					 dest);
}

static inline
void* elina_generic_assign_texpr_array(elina_manager_t* man,
				    bool destructive, void* abs, elina_dim_t* tdim, elina_texpr0_t** texpr, size_t size,
				    void* dest)
{
  return elina_generic_asssub_texpr_array(true,
				       man, destructive, abs, tdim, texpr, size,
				       dest);
}
static inline
void* elina_generic_substitute_texpr_array(elina_manager_t* man,
					bool destructive,
					void* abs,
					elina_dim_t* tdim,
					elina_texpr0_t** texpr,
					size_t size,
					void* dest)
{
  return elina_generic_asssub_texpr_array(false,
				       man, destructive, abs, tdim, texpr, size,
				       dest);
}

#ifdef __cplusplus
}
#endif

#endif

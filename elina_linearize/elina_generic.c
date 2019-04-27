/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
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
/* elina_generic: generic functions for library implementors */
/* ************************************************************************* */

#include "elina_generic.h"

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
/* I. Constructors */
/* ********************************************************************** */

bool elina_generic_sat_tcons(elina_manager_t* man, void* abs, elina_tcons0_t* cons,
			     elina_scalar_discr_t discr,
			     bool quasilinearize)
{
  bool (*is_bottom)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_IS_BOTTOM];
  bool (*sat_lincons)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_SAT_LINCONS];
  bool exact;
  elina_lincons0_t lincons0;
  bool res;
  elina_abstract0_t a0;

  if (is_bottom(man,abs)){
  man->result.flag_exact = man->result.flag_best = true;
    return true;
  }

  a0.value = abs;
  a0.man = man;
  lincons0 = elina_intlinearize_tcons0(man,&a0,cons,&exact,discr,quasilinearize,false);
  res = elina_lincons0_is_sat(&lincons0) || sat_lincons(man,abs,&lincons0);
  elina_lincons0_clear(&lincons0);
  if (!exact){
    man->result.flag_exact = man->result.flag_best = false;
  }
  return res;
}

/*
   This function implements a generic bound_texpr operation using linearisation
   and bound_linexpr
*/
elina_interval_t* elina_generic_bound_texpr(elina_manager_t* man, void* abs, elina_texpr0_t* expr,
				      elina_scalar_discr_t discr, bool quasilinearize)
{
  bool (*is_bottom)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_IS_BOTTOM];
  elina_interval_t* (*bound_linexpr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_BOUND_LINEXPR];
  bool exact;
  elina_linexpr0_t* linexpr0;
  elina_interval_t* res;
  elina_abstract0_t a0;

  if (is_bottom(man,abs)){
    res = elina_interval_alloc();
    elina_interval_set_bottom(res);
    return res;
  }

  a0.value = abs;
  a0.man = man;
  linexpr0 = elina_intlinearize_texpr0(man,&a0,expr,&exact,discr,quasilinearize);
  res = bound_linexpr(man,abs,linexpr0);
  elina_linexpr0_free(linexpr0);
  if (!exact){
    man->result.flag_exact = man->result.flag_best = false;
  }
  return res;
}

elina_tcons0_array_t elina_generic_to_tcons_array(elina_manager_t* man,
					    void* abs)
{
  elina_lincons0_array_t (*to_lincons_array)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_TO_LINCONS_ARRAY];
  elina_lincons0_array_t array = to_lincons_array(man,abs);
  elina_tcons0_array_t res = elina_tcons0_array_make(array.size);
  size_t i;
  for (i=0; i<array.size; i++){
    res.p[i] = elina_tcons0_from_lincons0(&array.p[i]);
  }
  elina_lincons0_array_clear(&array);
  return res;
}


/* ********************************************************************** */
/* II. Operations */
/* ********************************************************************** */

/* ============================================================ */
/*  Meet/Join on arrays of abstract values */
/* ============================================================ */

/*
   This function implements a generic meet/join_array operation using copy and meet
   operations.
*/
void* elina_generic_meetjoin_array(bool meet,
				elina_manager_t* man,
				void** tab, size_t size)
{
  void* (*copy)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_COPY];
  void* (*meetjoin)(elina_manager_t*,...) = man->funptr[meet ? ELINA_FUNID_MEET : ELINA_FUNID_JOIN];
  size_t i;
  void* res;
  bool exact,best;
  if (size==1){
    return copy(man,tab[0]);
  }
  else {
    res = meetjoin(man,false,tab[0],tab[1]);
    exact = man->result.flag_exact;
    best =  man->result.flag_best;
    for (i=2; i<size; i++){
      res = meetjoin(man,true,res,tab[i]);
      exact = exact && man->result.flag_exact;
      best =  best && man->result.flag_best;
    }
    man->result.flag_exact = exact;
    man->result.flag_best = best;
    return res;
  }
}

/* ============================================================ */
/*  Meet with array of constraints */
/* ============================================================ */

void*
elina_generic_meet_quasilinearize_lincons_array(elina_manager_t* man,
					     bool destructive,
					     void* abs,
					     elina_lincons0_array_t* array,
					     elina_scalar_discr_t discr,
					     bool linearize,
					     void* (*meet_lincons_array)(elina_manager_t*,
									 bool, void*,
									 elina_lincons0_array_t*))
{
  bool exact;
  elina_lincons0_array_t array2;
  void* res;
  void* (*copy)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_COPY];
  bool (*is_bottom)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_IS_BOTTOM];

  man->result.flag_exact = man->result.flag_best = true;

  if (is_bottom(man,abs) || array->size==0){
    res = destructive ? abs : copy(man, abs);
  }
  else {
    array2 = elina_quasilinearize_lincons0_array(man,abs,array,&exact,
					      discr,linearize,true);
    res = meet_lincons_array(man,destructive,abs,&array2);
    if (!exact){
      man->result.flag_exact = man->result.flag_best = false;
    }
    if (array2.p!=array->p){
      elina_lincons0_array_clear(&array2);
    }
  }
  return res;
}

void*
elina_generic_meet_intlinearize_tcons_array(elina_manager_t* man,
					 bool destructive,
					 void* abs,
					 elina_tcons0_array_t* array,
					 elina_scalar_discr_t discr,
					 elina_linexpr_type_t linearize,
					 void* (*meet_lincons_array)(elina_manager_t*,
								     bool, void*,
								     elina_lincons0_array_t*))
{
  bool exact;
  elina_lincons0_array_t array2;
  void* res;
  void* (*copy)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_COPY];
  bool (*is_bottom)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_IS_BOTTOM];
  elina_abstract0_t a0;

  man->result.flag_exact = man->result.flag_best = true;

  if (is_bottom(man,abs) || array->size==0){
    res = destructive ? abs : copy(man,abs);
  }
  else {
    a0.value = abs;
    a0.man = man;
    array2 = elina_intlinearize_tcons0_array(man,&a0,array,&exact,discr,linearize,true,true,2,false);
    res = meet_lincons_array(man,destructive,abs,&array2);
    if (!exact){
      man->result.flag_exact = man->result.flag_best = false;
    }
	
    elina_lincons0_array_clear(&array2);
  }
  return res;
}

/* ============================================================ */
/*  Assignments/Substitutions */
/* ============================================================ */

/*
   This function implements generic parallel assignment/substitution
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
   permute_dimensions, remove_dimensions, meet_lincons_array, meet and free
   abstract operations.

   Meaning of parameters:
   - assign selects the operation: true means assignment, false substitution
   - The other parameters have the meaning they have for parallel
     assignment/substitution
*/

void* elina_generic_asssub_linexpr_array(bool assign,
				      elina_manager_t* man,
				      bool destructive, void* abs, elina_dim_t* tdim, elina_linexpr0_t** texpr, size_t size,
				      void* dest)
{
  bool (*is_bottom)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_IS_BOTTOM];
  void* (*copy)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_COPY];
  void* (*add_dimensions)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_ADD_DIMENSIONS];
  void* (*permute_dimensions)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_PERMUTE_DIMENSIONS];
  void* (*remove_dimensions)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_REMOVE_DIMENSIONS];
  void* (*meet_lincons_array)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_MEET_LINCONS_ARRAY];
  void* (*meet)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_MEET];
  void (*elina_free)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_FREE];
  elina_dimension_t (*dimension)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_DIMENSION];
  size_t i;
  elina_dimension_t d, dsup;
  elina_dimchange_t dimchange;
  elina_dimperm_t permutation;
  elina_lincons0_array_t array;
  void* abs2;
  bool exact,best;

  if (is_bottom(man,abs)){
    man->result.flag_exact = man->result.flag_best = true;
    return destructive ? abs : copy(man,abs);
  }
  /* 1. Compute the number of integer and real dimensions assigned */
  d = dimension(man,abs);
  dsup.intdim = 0;
  dsup.realdim = 0;
  for (i=0; i<size; i++){
    if (tdim[i]<d.intdim)
      dsup.intdim++;
    else
      dsup.realdim++;
  }
  /* 2. Build dimchange (for addition of primed dimensions) */
  elina_dimchange_init(&dimchange,dsup.intdim,dsup.realdim);
  for (i=0;i<dsup.intdim;i++)
    dimchange.dim[i]=d.intdim;
  for (i=dsup.intdim;i<dsup.intdim+dsup.realdim;i++)
    dimchange.dim[i]=d.intdim+d.realdim;

  /* 3. Build permutation (exchanging primed and unprimed dimensions) */
  elina_dimperm_init(&permutation,d.intdim+d.realdim+dsup.intdim+dsup.realdim);
  elina_dimperm_set_id(&permutation);
  {
    int index_int = 0;
    int index_real = 0;
    for (i=0; i<size; i++){
      elina_dim_t dim = tdim[i];
      elina_dim_t dimp;
      if (dim<d.intdim){
	dimp = d.intdim+index_int;
	index_int++;
      } else {
	dim += dsup.intdim;
	dimp = d.intdim+dsup.intdim+d.realdim+index_real;
	index_real++;
      }
      permutation.dim[dim] = dimp;
      permutation.dim[dimp] = dim;
    }
  }
  /* 4. Add primed dimensions to abstract value */
  abs2 = add_dimensions(man,destructive,abs,&dimchange,false);
  exact = man->result.flag_exact;
  best =  man->result.flag_best;
  /* From now, work by side-effect on abs2 */

  /* 5. Build constraints system
     An assignment x'_i := a_ij x_j + b_i becomes
     an equality constraint -x'_i + a_ij x_j + b_i = 0
     Primed and unprimed dimensiosn permuted if dest!=NULL
  */
  array = elina_lincons0_array_make(size);
  for (i=0; i<size; i++){
    elina_dim_t dim = tdim[i];
    if (dim>=d.intdim) dim += dsup.intdim;
    elina_dim_t dimp = permutation.dim[dim];
    elina_linexpr0_t* expr = elina_linexpr0_add_dimensions(texpr[i],&dimchange);
    elina_linexpr0_set_coeff_scalar_double(expr,dimp,-1.0);
    elina_lincons0_t cons = elina_lincons0_make(ELINA_CONS_EQ,expr,NULL);
    array.p[i] = cons;
  }
  /* 6. Permute unprimed and primed dimensions if !assign */
  if (!assign){
    abs2 = permute_dimensions(man,true,abs2,&permutation);
    exact = exact && man->result.flag_exact;
    best =  best && man->result.flag_best;
  }
  /* 7. If dest!=NULL, perform intersection */
  if (dest!=NULL){
    void* dest2 = add_dimensions(man,false,dest,&dimchange,false);
    exact = exact && man->result.flag_exact;
    best =  best && man->result.flag_best;

    if (assign){
      dest2 = permute_dimensions(man,true,dest2,&permutation);
      exact = exact && man->result.flag_exact;
      best =  best && man->result.flag_best;
    }
    abs2 = meet(man,true,abs2,dest2);
    exact = exact && man->result.flag_exact;
    best =  best && man->result.flag_best;

    elina_free(man,dest2);
  }
  /* 8. Perform meet of abs2 with constraints */
  abs2 = meet_lincons_array(man,true,abs2,&array);
  exact = exact && man->result.flag_exact;
  best =  best && man->result.flag_best;

  /* 9. Permute unprimed and primed dimensions if assign */
  if (assign){
    abs2 = permute_dimensions(man,true,abs2,&permutation);
    exact = exact && man->result.flag_exact;
    best =  best && man->result.flag_best;
  }
  /* 10. Remove extra dimensions */
  elina_dimchange_add_invert(&dimchange);
  abs2 = remove_dimensions(man,true,abs2,&dimchange);
  exact = exact && man->result.flag_exact;
  best =  best && man->result.flag_best;

  /* 11. Free allocated objects */
  elina_dimperm_clear(&permutation);
  elina_dimchange_clear(&dimchange);
  elina_lincons0_array_clear(&array);

  man->result.flag_exact = exact;
  man->result.flag_best = best;
  return abs2;
}

void* elina_generic_asssub_texpr_array(bool assign,
				    elina_manager_t* man,
				    bool destructive, void* abs, elina_dim_t* tdim, elina_texpr0_t** texpr, size_t size,
				    void* dest)
{
  bool (*is_bottom)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_IS_BOTTOM];
  void* (*copy)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_COPY];
  void* (*add_dimensions)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_ADD_DIMENSIONS];
  void* (*permute_dimensions)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_PERMUTE_DIMENSIONS];
  void* (*remove_dimensions)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_REMOVE_DIMENSIONS];
  void* (*meet_tcons_array)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_MEET_TCONS_ARRAY];
  void* (*meet)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_MEET];
  void (*elina_free)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_FREE];
  elina_dimension_t (*dimension)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_DIMENSION];
  size_t i;
  elina_dimension_t d, dsup;
  elina_dimchange_t dimchange;
  elina_dimperm_t permutation;
  elina_tcons0_array_t array;
  void* abs2;
  bool exact,best;

  if (is_bottom(man,abs)){
    man->result.flag_exact = man->result.flag_best = true;
    return destructive ? abs : copy(man,abs);
  }
  /* 1. Compute the number of integer and real dimensions assigned */
  d = dimension(man,abs);
  dsup.intdim = 0;
  dsup.realdim = 0;
  for (i=0; i<size; i++){
    if (tdim[i]<d.intdim)
      dsup.intdim++;
    else
      dsup.realdim++;
  }
  /* 2. Build dimchange (for addition of primed dimensions) */
  elina_dimchange_init(&dimchange,dsup.intdim,dsup.realdim);
  for (i=0;i<dsup.intdim;i++)
    dimchange.dim[i]=d.intdim;
  for (i=dsup.intdim;i<dsup.intdim+dsup.realdim;i++)
    dimchange.dim[i]=d.intdim+d.realdim;

  /* 3. Build permutation (exchanging primed and unprimed dimensions) */
  elina_dimperm_init(&permutation,d.intdim+d.realdim+dsup.intdim+dsup.realdim);
  elina_dimperm_set_id(&permutation);
  {
    int index_int = 0;
    int index_real = 0;
    for (i=0; i<size; i++){
      elina_dim_t dim = tdim[i];
      elina_dim_t dimp;
      if (dim<d.intdim){
	dimp = d.intdim+index_int;
	index_int++;
      } else {
	dim += dsup.intdim;
	dimp = d.intdim+dsup.intdim+d.realdim+index_real;
	index_real++;
      }
      permutation.dim[dim] = dimp;
      permutation.dim[dimp] = dim;
    }
  }

  /* 4. Add primed dimensions to abstract value */
  abs2 = add_dimensions(man,destructive,abs,&dimchange,false);
  exact = man->result.flag_exact;
  best =  man->result.flag_best;
  /* From now, work by side-effect on abs2 */

  /* 5. Build constraints system
     An assignment x'_i := a_ij x_j + b_i becomes
     an equality constraint -x'_i + a_ij x_j + b_i = 0
  */
  array = elina_tcons0_array_make(size);
  for (i=0; i<size; i++){
    elina_dim_t dim = tdim[i];
    if (dim>=d.intdim) dim += dsup.intdim;
    elina_dim_t dimp = permutation.dim[dim];
    elina_texpr0_t* expr = elina_texpr0_add_dimensions(texpr[i],&dimchange);
    expr = elina_texpr0_binop(ELINA_TEXPR_SUB, expr,elina_texpr0_dim(dimp),
			   ELINA_RTYPE_REAL, ELINA_RDIR_RND);
    elina_tcons0_t cons = elina_tcons0_make(ELINA_CONS_EQ,expr,NULL);
    array.p[i] = cons;
  }
  /* 6. Permute unprimed and primed dimensions if !assign */
  if (!assign){
    abs2 = permute_dimensions(man,true,abs2,&permutation);
    exact = exact && man->result.flag_exact;
    best =  best && man->result.flag_best;
  }
  /* 7. If dest!=NULL, perform intersection */
  if (dest!=NULL){
    void* dest2 = add_dimensions(man,false,dest,&dimchange,false);
    exact = exact && man->result.flag_exact;
    best =  best && man->result.flag_best;

    if (assign){
      dest2 = permute_dimensions(man,true,dest2,&permutation);
      exact = exact && man->result.flag_exact;
      best =  best && man->result.flag_best;
    }
    abs2 = meet(man,true,abs2,dest2);
    exact = exact && man->result.flag_exact;
    best =  best && man->result.flag_best;

    elina_free(man,dest2);
  }
  /* 8. Perform meet of abs2 with constraints */
  abs2 = meet_tcons_array(man,true,abs2,&array);
  exact = exact && man->result.flag_exact;
  best =  best && man->result.flag_best;

  /* 9. Permute unprimed and primed dimensions if assign */
  if (assign){
    abs2 = permute_dimensions(man,true,abs2,&permutation);
    exact = exact && man->result.flag_exact;
    best =  best && man->result.flag_best;
  }
  /* 10. Remove extra dimensions */
  elina_dimchange_add_invert(&dimchange);
  abs2 = remove_dimensions(man,true,abs2,&dimchange);
  exact = exact && man->result.flag_exact;
  best =  best && man->result.flag_best;

  /* 11. Free allocated objects */
  elina_dimperm_clear(&permutation);
  elina_dimchange_clear(&dimchange);
  elina_tcons0_array_clear(&array);
  man->result.flag_exact = exact;
  man->result.flag_best = best;
  return abs2;
}

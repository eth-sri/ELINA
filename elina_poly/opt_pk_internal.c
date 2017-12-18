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

/* ************************************************************************* */
/* internal.c: internal manager */
/* ************************************************************************* */

#include "opt_pk_config.h"
#include "opt_pk_internal.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"

/* ********************************************************************** */
/* I. Constructor and destructor for internal */
/* ********************************************************************** */


#if defined(TIMING)
	double copy_time = 0;
	double is_lequal_time = 0;
	double permute_dimension_time = 0;
	double add_dimension_time = 0;
	double remove_dimension_time = 0;
	double top_time = 0;
	double bottom_time = 0;
	double meet_time = 0;
	double join_time = 0;
	double widening_time = 0;
	double free_time = 0;
	double forget_array_time = 0;
	double meet_lincons_time = 0;
	double poly_to_box_time = 0;
	double is_top_time = 0;
	double is_bottom_time = 0;
	double expand_time = 0;
	double fold_time = 0;
	double sat_lincons_time = 0;
	double assign_linexpr_time = 0;
	double substitute_linexpr_time = 0;
	double bound_dimension_time = 0;
	double opt_conversion_time = 0;
	long int join_count = 0;
	double poly_is_unconstrained_time = 0;
#endif

/* Initialize opk with size maxdims */
void opt_pk_internal_init(opt_pk_internal_t* opk, size_t maxdims)
{
  size_t i;

  opk->exn = ELINA_EXC_NONE;

  opk->maxdims = maxdims;
  opk->maxcols = maxdims+3;
  
  opk->vector_numintp = opt_vector_alloc(opk->maxcols);
  //opk->vector_mpqp = malloc( (maxdims+3)*sizeof(mpq_t));
  //for (i=0; i<maxdims+3; i++)
    //mpq_init(pk->vector_mpqp[i]);

  opk->vector_tmp = opt_vector_alloc(5);

  opk->matrix_dimp = malloc(opk->maxdims*sizeof(elina_dim_t));
  opk->matrix_acc = 0;
  opk->matrix_prod = 0;

  /* opk->cherni_bitstringp = bitstring_alloc(bitindex_size(pk->maxrows));*/
  opk->cherni_intp = (int*)malloc(opk->maxcols * sizeof(int));
  opk->cherni_prod = 0;

  opk->poly_numrat = (elina_rat_t *)calloc(1,sizeof(elina_rat_t));
  opk->poly_numintp = opt_vector_alloc(opk->maxcols);
  opk->poly_numintp2 = opt_vector_alloc(opk->maxcols);
  opk->poly_dimp = malloc(opk->maxdims*sizeof(elina_dim_t));
  opk->poly_dimp2 = malloc(opk->maxdims*sizeof(elina_dim_t));
  opk->poly_fold_dimp = malloc(opk->maxdims*sizeof(elina_dim_t));
  opk->poly_prod = 0;
}

/* Allocates opk and initializes it with a default size */
opt_pk_internal_t* opt_pk_internal_alloc(bool strict)
{
  opt_pk_internal_t* opk = (opt_pk_internal_t*)malloc(sizeof(opt_pk_internal_t));

  opk->strict = strict;
  opk->dec = strict ? 3 : 2;
  opk->max_coeff_size = 0;
  opk->approximate_max_coeff_size = 2;

  opt_pk_internal_init(opk,10);

  return opk;
}

/* Clear opk */
void opt_pk_internal_clear(opt_pk_internal_t* opk)
{
  size_t i;
  if (opk->vector_numintp) opt_vector_free(opk->vector_numintp,opk->maxcols);
  opk->vector_numintp = NULL;
  if (opk->vector_tmp) opt_vector_free(opk->vector_tmp,5);
  opk->vector_tmp = NULL;
  if (opk->matrix_dimp) free(opk->matrix_dimp);
  opk->matrix_dimp = NULL;
   opk->matrix_acc = 0;
   opk->matrix_prod = 0;

  
  if (opk->cherni_intp) free(opk->cherni_intp);
  opk->cherni_intp = NULL;
   
  opk->cherni_prod = 0;

  free(opk->poly_numrat);
  if (opk->poly_numintp) opt_vector_free(opk->poly_numintp, opk->maxcols);
  opk->poly_numintp = NULL; 

  if (opk->poly_numintp2) opt_vector_free(opk->poly_numintp2, opk->maxcols);
  opk->poly_numintp2 = NULL;

  if (opk->poly_dimp) free(opk->poly_dimp);
  opk->poly_dimp = NULL;
  if (opk->poly_dimp2) free(opk->poly_dimp2);
  opk->poly_dimp2 = NULL;
  if (opk->poly_fold_dimp) free(opk->poly_fold_dimp);
  opk->poly_fold_dimp = NULL;
  

  opk->poly_prod=0;
  opk->maxdims = 0;
  opk->maxrows = 0;
  opk->maxcols = 0;
}

/* Clear and free pk */
void opt_pk_internal_free(opt_pk_internal_t* opk)
{
  opt_pk_internal_clear(opk);
  free(opk);
}



/* Reallocate opk only if a bigger dimension is required */
void opt_pk_internal_realloc_lazy(opt_pk_internal_t* opk, size_t maxdims)
{
  if (maxdims > opk->maxdims){
     //printf("GG: %d %d\n",opk->maxcols,maxdims);
    //fflush(stdout);
    opt_pk_internal_clear(opk);
    //printf("KK\n");
    //fflush(stdout);
    opt_pk_internal_init(opk,maxdims);
    //printf("KK\n");
    //fflush(stdout);
  }
}

/* ********************************************************************** */
/* II. Options */
/* ********************************************************************** */


void opt_pk_set_approximate_max_coeff_size(opt_pk_internal_t* opk, size_t size){  
  opk->approximate_max_coeff_size = size;
}


/* ********************************************************************** */
/* III. Initialization from manager */
/* ********************************************************************** */

elina_manager_t* opt_pk_manager_alloc(bool strict)
{
  size_t i;
  opt_pk_internal_t* opk;
  elina_manager_t* man;
  void** funptr;

  opk = opt_pk_internal_alloc(strict);
  opt_pk_set_approximate_max_coeff_size(opk, 1);
  man = elina_manager_alloc(
      strict ? "optpoly, strict mode" : "optpoly, loose mode",
      "1.0 with Long Long Int", opk, (void (*)(void *))opt_pk_internal_free);
  funptr = man->funptr;
  
  funptr[ELINA_FUNID_COPY] = &opt_pk_copy;
  funptr[ELINA_FUNID_FREE] = &opt_pk_free;
  funptr[ELINA_FUNID_ASIZE] = &opt_pk_size;
  funptr[ELINA_FUNID_MINIMIZE] = &opt_pk_minimize;
  funptr[ELINA_FUNID_CANONICALIZE] = &opt_pk_array_canonicalize;
  //funptr[ELINA_FUNID_HASH] = &opt_pk_hash;
  //funptr[ELINA_FUNID_APPROXIMATE] = &opt_pk_approximate;
  funptr[ELINA_FUNID_FPRINT] = &opt_pk_array_fprint;
  funptr[ELINA_FUNID_BOTTOM] = &opt_pk_bottom;
  funptr[ELINA_FUNID_TOP] = &opt_pk_top;
  funptr[ELINA_FUNID_OF_BOX] = &opt_pk_of_box;
  funptr[ELINA_FUNID_DIMENSION] = &opt_pk_dimension;
  funptr[ELINA_FUNID_IS_BOTTOM] = &opt_pk_is_bottom;
  funptr[ELINA_FUNID_IS_TOP] = &opt_pk_is_top;
  funptr[ELINA_FUNID_IS_LEQ] = &opt_pk_is_leq;
  funptr[ELINA_FUNID_IS_EQ] = &opt_pk_is_eq;
  funptr[ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED] = &opt_pk_is_dimension_unconstrained;
  //funptr[ELINA_FUNID_SAT_INTERVAL] = &opt_pk_sat_interval;
  funptr[ELINA_FUNID_SAT_LINCONS] = &opt_pk_sat_lincons;
  funptr[ELINA_FUNID_SAT_TCONS] = &opt_pk_sat_tcons;
  funptr[ELINA_FUNID_BOUND_DIMENSION] = &opt_pk_bound_dimension;
  funptr[ELINA_FUNID_BOUND_LINEXPR] = &opt_pk_bound_linexpr;
  //funptr[ELINA_FUNID_BOUND_TEXPR] = &opt_pk_bound_texpr;
  funptr[ELINA_FUNID_TO_BOX] = &opt_pk_to_box;
  funptr[ELINA_FUNID_TO_LINCONS_ARRAY] = &opt_pk_to_lincons_array;
  funptr[ELINA_FUNID_MEET] = &opt_pk_meet;
  //funptr[ELINA_FUNID_MEET_ARRAY] = &opt_pk_meet_array;
  funptr[ELINA_FUNID_MEET_LINCONS_ARRAY] = &opt_pk_meet_lincons_array;
  funptr[ELINA_FUNID_MEET_TCONS_ARRAY] = &opt_pk_meet_tcons_array;
  funptr[ELINA_FUNID_JOIN] = &opt_pk_join;
  //funptr[ELINA_FUNID_JOIN_ARRAY] = &opt_pk_join_array;
  funptr[ELINA_FUNID_ASSIGN_LINEXPR_ARRAY] = &opt_pk_assign_linexpr_array;
  funptr[ELINA_FUNID_SUBSTITUTE_LINEXPR_ARRAY] = &opt_pk_substitute_linexpr_array;
  funptr[ELINA_FUNID_ASSIGN_TEXPR_ARRAY] = &opt_pk_assign_texpr_array;
  funptr[ELINA_FUNID_ADD_DIMENSIONS] = &opt_pk_add_dimensions;
  funptr[ELINA_FUNID_REMOVE_DIMENSIONS] = &opt_pk_remove_dimensions;
  funptr[ELINA_FUNID_PERMUTE_DIMENSIONS] = &opt_pk_permute_dimensions;
  funptr[ELINA_FUNID_FORGET_ARRAY] = &opt_pk_forget_array;
  funptr[ELINA_FUNID_EXPAND] = &opt_pk_expand;
  funptr[ELINA_FUNID_FOLD] = &opt_pk_fold;
  funptr[ELINA_FUNID_WIDENING] = &opt_pk_widening;

  elina_manager_set_abort_if_exception(man, ELINA_EXC_TIMEOUT, false);
  elina_manager_set_abort_if_exception(man, ELINA_EXC_OUT_OF_SPACE, false);
  elina_manager_set_abort_if_exception(man, ELINA_EXC_OVERFLOW, false);
  
  return man;
}


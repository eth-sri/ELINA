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


/* ************************************************************************* */
/* opt_pk_internal.h: internal manager */
/* ************************************************************************* */


#ifndef _OPT_PK_INTERNAL_H_
#define _OPT_PK_INTERNAL_H_

#include "opt_pk_config.h"
#include "rdtsc.h"
#include "opt_pk.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(TIMING)
	#define start_timing()				\
		tsc_counter start, end;				\
		double cycles;					\
		CPUID();					\
		RDTSC(start)						

  	#define record_timing(counter)		\
		RDTSC(end);				\
  		CPUID();				\
  		cycles = (double)(COUNTER_DIFF(end, start));	\
  		counter += cycles	
	extern double copy_time;
	extern double is_equal_time;
	extern double is_lequal_time;
	extern double permute_dimension_time;
	extern double add_dimension_time;
	extern double remove_dimension_time;
	extern double top_time;
	extern double bottom_time;
	extern double meet_time;
	extern double join_time;
	extern double widening_time;
	extern double free_time;
	extern double forget_array_time;
	extern double meet_lincons_time;
	extern double poly_to_box_time;
	extern double is_top_time;
	extern double is_bottom_time;
	extern double expand_time;
	extern double fold_time;
	extern double sat_lincons_time;
	extern double assign_linexpr_time;
	extern double substitute_linexpr_time;
	extern double bound_dimension_time;
	extern double opt_conversion_time;
	extern long int join_count;
	extern double poly_is_unconstrained_time;
#endif
/* ********************************************************************** */
/* I. Types */
/* ********************************************************************** */

/* These variables are used by various functions.  The prefix XXX_
   indicates that the variable is used by the module XXX. */

struct opt_pk_internal_t {
  enum elina_exc_t exn;

  bool strict;
  size_t dec;

  size_t maxdims;
  size_t maxcols;
  size_t maxrows;

  elina_funid_t funid;
  elina_funopt_t* funopt;

  size_t max_coeff_size; /* Used for overflow exception in vector_combine */
  size_t approximate_max_coeff_size;

  opt_numint_t * vector_numintp; /* of size maxcols */

  //mpq_t* vector_mpqp; /* of size maxdims+3 */

  opt_numint_t * vector_tmp;    /* of size 5 */
  
  elina_dim_t* matrix_dimp;                /* of size maxdims */
  opt_numint_t matrix_acc;
  opt_numint_t matrix_prod;

  /* bitstring_t* cherni_bitstringp; */ /* of size maxrows */
  int* cherni_intp;                /* of size maxcols */
  opt_numint_t cherni_prod;             

  
  //bound_t poly_bound;
  elina_rat_t * poly_numrat;
  opt_numint_t * poly_numintp;            /* of size maxcols */
  opt_numint_t * poly_numintp2;           /* of size maxcols */
  /* bitstring_t* poly_bitstringp; */    /* of size maxrows */
  elina_dim_t* poly_dimp;                /* of size maxdims */
  elina_dim_t* poly_dimp2;               /* of size maxdims */
  elina_dim_t* poly_fold_dimp;               /* of size maxdims */
  //struct matrix_t* poly_matspecial; 
  opt_numint_t poly_prod; 
};

/* ********************************************************************** */
/* A. Constructor and destructor for internal */
/* ********************************************************************** */

opt_pk_internal_t* opt_pk_internal_alloc(bool strict);
  /* Allocates pk and initializes it with a default size */
void opt_pk_internal_free(opt_pk_internal_t* pk);
  /* Clear and free pk */
void opt_pk_internal_realloc_lazy(opt_pk_internal_t* pk, size_t maxdims);
  /* Reallocate pk only if a bigger dimension is required */

static inline opt_pk_internal_t* opt_pk_init_from_manager(elina_manager_t* man, elina_funid_t funid);
  /* Initializes some fields of pk from manager */

/* ********************************************************************** */
/* Definition of inline functions */
/* ********************************************************************** */
static inline opt_pk_internal_t* opt_pk_init_from_manager(elina_manager_t* man, elina_funid_t funid)
{
	//printf("id: %d\n",funid);
	//fflush(stdout);
  opt_pk_internal_t* opk = (opt_pk_internal_t*)man->internal;
  opk->funid = funid;
  opk->funopt = &man->option.funopt[funid];
  man->result.flag_exact = man->result.flag_best = false;
  return opk;
}


#ifdef __cplusplus
}
#endif

#endif

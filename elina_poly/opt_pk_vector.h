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


/* ********************************************************************** */
/* opt_pk_vector.h: operations on vectors */
/* ********************************************************************** */


/* This header file defines operations on vectors. A vector is just an array of
   elements of type pkint_t, as a consequence functions need to be given also
   their size. */

#ifndef __OPT_PK_VECTOR_H__
#define __OPT_PK_VECTOR_H__

#include "opt_pk_config.h"
#include "opt_pk_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Basic Operations */
opt_numint_t* opt_vector_alloc(unsigned short int size);
void opt_vector_realloc(opt_numint_t** ov, unsigned short int size, unsigned short int nsize);
void opt_vector_free(opt_numint_t* q, unsigned short int size);
void opt_vector_clear(opt_numint_t * ov, unsigned short int size);
void opt_vector_copy(opt_numint_t * dst, opt_numint_t* src, unsigned short int size);
void opt_vector_print(opt_numint_t * ov, unsigned short int size);

/* Normalization */
bool opt_vector_normalize(opt_pk_internal_t* opk,
		      opt_numint_t* ov, unsigned short int size);
bool opt_vector_normalize_expr(opt_pk_internal_t* opk,
			   opt_numint_t* ov, unsigned short int size);
bool opt_vector_normalize_constraint(opt_pk_internal_t* opk,
				 opt_numint_t* ov, 
				unsigned short int intdim, unsigned short int realdim);
bool opt_vector_normalize_constraint_int(opt_pk_internal_t* opk,
				     opt_numint_t* ov, 
				     unsigned short int intdim, unsigned short int realdim);

opt_numint_t * opt_vector_scalar_product(opt_pk_internal_t *opk, opt_numint_t * src, opt_numint_t s, unsigned short int size);

void opt_vector_sum(opt_pk_internal_t * opk, opt_numint_t * ov1, opt_numint_t *ov2, unsigned short int size);


/* Comparison */
int opt_vector_compare_coeff(opt_pk_internal_t* opk,
		   opt_numint_t* ov1, opt_numint_t* ov2, unsigned short int size);

int opt_vector_compare(opt_pk_internal_t* opk,
		   opt_numint_t * ov1, opt_numint_t * ov2, unsigned short int size);

/* Combination and Algebraic Operations */
void opt_vector_combine(opt_pk_internal_t* opk,
		    opt_numint_t* ov1, opt_numint_t* ov2,
		    opt_numint_t* ov3, size_t k, unsigned short intsize, bool add);

opt_numint_t * opt_vector_neg(opt_pk_internal_t *opk, opt_numint_t *ov, unsigned short int size);

void opt_vector_mul_scalar(opt_numint_t *dst, opt_numint_t *src, opt_numint_t prod, unsigned short int size);

void opt_vector_add(opt_numint_t *dst, opt_numint_t *op1, opt_numint_t *op2, unsigned short int size);

/*opt_numint_t opt_vector_product(opt_pk_internal_t* opk,
		    opt_numint_t* q1, opt_numint_t* q2, size_t size);*/

opt_numint_t opt_vector_product_strict(opt_pk_internal_t* opk,
			   opt_numint_t* r1, opt_numint_t* r2, unsigned short int size);

//opt_numint_t opt_vector_product_strict_comp(opt_pk_internal_t* opk,
//			   opt_numint_t* q1, opt_numint_t* q2, unsigned short int * ind_map1,
//			   unsigned short int *ind_map2, unsigned short int size1, 
//			   unsigned short int size2, unsigned short int size);

void opt_vector_permute_dimensions(opt_pk_internal_t* opk,
			       opt_numint_t* nov, opt_numint_t* ov, unsigned short int size,
			       elina_dim_t* permut);

void opt_vector_add_dimensions(opt_pk_internal_t* opk,
			   opt_numint_t* nov, 
			   opt_numint_t* ov, unsigned short int size,
			   elina_dimchange_t* dimchange);

void opt_vector_remove_dimensions(opt_pk_internal_t* pk,
			      opt_numint_t* nov, 
			      opt_numint_t* ov, unsigned short int size,
			      elina_dimchange_t* dimchange);

/* Predicates that can be useful for users */
bool opt_vector_is_null(opt_pk_internal_t* opk,
		    opt_numint_t* ov, unsigned short int size);

bool opt_vector_is_null_expr(opt_pk_internal_t *opk,
			     opt_numint_t* ov, unsigned short int size);

bool opt_vector_is_positivity_constraint(opt_pk_internal_t* opk,
				     opt_numint_t* ov, unsigned short int size);
bool opt_vector_is_dummy_constraint(opt_pk_internal_t* opk,
				  opt_numint_t* ov, unsigned short int size);

bool opt_vector_is_integer(opt_pk_internal_t* opk,
		       opt_numint_t* ov,
		       size_t intdim, size_t realdim);

/* Functions meant to be internal */

opt_numint_t* _opt_vector_alloc_int(unsigned short int size);

void opt_vector_simplify(opt_pk_internal_t *opk, opt_numint_t * ov, unsigned short int size);

opt_numint_t opt_vector_gcd(opt_pk_internal_t* opk,
		opt_numint_t* ov, unsigned short int size);

unsigned short int * build_index_vector(opt_numint_t *ov, unsigned short int size);

opt_numint_t opt_vector_product(opt_pk_internal_t* opk,
		    opt_numint_t* q1, opt_numint_t* q2, size_t size);

opt_numint_t opt_vector_product_with_index(opt_pk_internal_t* opk,
		    opt_numint_t* q1, opt_numint_t* q2, unsigned short int *nz);

opt_numint_t * opt_map_vector(opt_numint_t *q, unsigned short int *map, 
			      unsigned short int comp_size, unsigned short int size);

bool is_vertex(opt_numint_t * v);
bool is_ray(opt_numint_t * v);
bool is_line(opt_numint_t * v);

bool opt_vector_equal(opt_pk_internal_t* opk,
		   opt_numint_t* ov1, opt_numint_t* ov2,
		   unsigned short int size, unsigned short int* ind);


/* Bounding the value of an (interval) linear expression (itv_linexpr) in a
   generator vector
*/
void opt_vector_bound_elina_linexpr0(opt_pk_internal_t* opk,
			      elina_rat_t *inf, elina_rat_t *sup,
			      elina_linexpr0_t* linexpr,
			      opt_numint_t* vec, size_t size);

#ifdef __cplusplus
}
#endif

#endif

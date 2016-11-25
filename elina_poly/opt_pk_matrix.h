/*
	Copyright 2016 Software Reliability Lab, ETH Zurich

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

		http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
*/



/* ********************************************************************** */
/* opt_pk_matrix.h: operations on matrices */
/* ********************************************************************** */


/*
A matrix is represented in the following manner: the coefficients are stored
in an private array of pkint_t p_init of size
_maxrows*nbcolumns. To access to elements, one use an array of
pointers p, the $i^{\mbox{\scriptsize nth}}$ element of which points
to the $i^{\mbox{\scriptsize nth}}$ row of the matrix. This array is
initialized by the constructor. The advantage of this representation is to be
able to exchange easily rows of the matrix by exchanging the pointers,
without having to allocate at construction _maxrows arrays for each
rows.  nbrows indicates that only the first nbrows rows are
used.
*/

#ifndef __OPT_PK_MATRIX_H__
#define __OPT_PK_MATRIX_H__

#include "opt_pk_internal.h"
#include "opt_pk_satmat.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct opt_matrix_t {
  /* public part */
  opt_numint_t ** p;     /* array of pointers to rows */
  size_t nbrows;      /* number of effective rows */
  unsigned short int nbcolumns;   /* size of rows */

  /* private part */
  size_t  _maxrows;   /* number of rows allocated */
  bool _sorted;
} opt_matrix_t;

/* Normal functions */

/* information about private part */
static inline size_t opt_matrix_get_maxrows(opt_matrix_t* mat)
{ return mat->_maxrows; }

static inline bool opt_matrix_is_sorted(opt_matrix_t* mat)
{ return mat->_sorted; }

/* Basic Operations */
opt_matrix_t* opt_matrix_alloc(size_t nbrows, unsigned short int nbcols, bool s);
void      opt_matrix_resize_rows(opt_matrix_t* mat, size_t nbrows);
void      opt_matrix_resize_rows_lazy(opt_matrix_t* mat, size_t nbrows);
//void      matrix_minimize(matrix_t* mat);
void      opt_matrix_free(opt_matrix_t* mat);
void      opt_matrix_clear(opt_matrix_t* mat);
void      opt_matrix_print(opt_matrix_t* mat);
void      opt_matrix_fprint(FILE* stream, opt_matrix_t* mat);
opt_matrix_t* opt_matrix_copy(opt_matrix_t* mat);
bool      opt_matrix_equal(opt_matrix_t* oma, opt_matrix_t* omb);

/* Operations on rows */
void opt_matrix_normalize_row(opt_pk_internal_t* opk,
			  opt_matrix_t* oc, size_t l);
void opt_matrix_combine_rows(opt_pk_internal_t* opk,
			 opt_matrix_t* oc, size_t l1, size_t l2, size_t l3, size_t k,bool add);
int opt_matrix_compare_rows(opt_pk_internal_t* opk,
			opt_matrix_t* oc, size_t l1, size_t l2);
void opt_matrix_exch_rows(opt_matrix_t* om, size_t l1, size_t l2);
void opt_matrix_move_rows(opt_matrix_t* om, size_t destrow, size_t orgrow, size_t size);
void opt_matrix_union(opt_matrix_t *op1, opt_matrix_t *op2);

/* Normalization of rows */
bool opt_matrix_normalize_constraint(opt_pk_internal_t* opk,
				   opt_matrix_t* oc, 
				   size_t intdim, size_t realdim);
//bool matrix_normalize_constraint_int(pk_internal_t* pk,
//				       matrix_t* mat, 
//				       size_t intdim, size_t realdim);

/* Sorting & Merging */
void opt_matrix_sort_rows(opt_pk_internal_t* opk,
		      opt_matrix_t* om);

void opt_matrix_sort_rows_from(opt_pk_internal_t* opk,
		      opt_matrix_t* om, size_t start, size_t num);

void opt_matrix_sort_rows_with_sat(opt_pk_internal_t* opk,
			       opt_matrix_t* mat, opt_satmat_t* sat);

opt_matrix_t* opt_matrix_append(opt_matrix_t* oma, opt_matrix_t* omb);
void opt_matrix_append_with(opt_matrix_t* oma, opt_matrix_t* omb);
void opt_matrix_revappend_with(opt_matrix_t* oma, opt_matrix_t* omb);

opt_matrix_t* opt_matrix_merge_sort(opt_pk_internal_t* opk,
			    opt_matrix_t* oma, opt_matrix_t* omb);
void opt_matrix_merge_sort_with(opt_pk_internal_t* opk,
			    opt_matrix_t* mat, opt_matrix_t* cmat);


void opt_matrix_project_var(opt_pk_internal_t *opk, 
			    opt_matrix_t * oc, 
			    size_t start, unsigned short int nbcolumns, 
			    elina_dim_t *tdim, size_t dimsup);

opt_matrix_t* opt_matrix_add_dimensions(opt_pk_internal_t* opk,
				bool destructive,
				opt_matrix_t* oc,
				elina_dimchange_t* dimchange, bool project);

opt_matrix_t* opt_matrix_remove_dimensions(opt_pk_internal_t* pk,
				   bool destructive,
				   opt_matrix_t* oc,
				   elina_dimchange_t* dimchange);

opt_matrix_t* opt_matrix_permute_dimensions(opt_pk_internal_t* opk,
				    bool destructive,
				    opt_matrix_t* oc,
				    elina_dim_t* permutation);

void opt_matrix_rearrange(opt_matrix_t *oc, size_t nbeq);

size_t opt_generator_rearrange(opt_matrix_t *F, opt_satmat_t * satF);

size_t opt_matrix_gauss_elimination(opt_pk_internal_t *opk, opt_matrix_t *oc, size_t nbeq);

void gauss_backsubstitute(opt_pk_internal_t* opk, opt_matrix_t* oc, size_t rank);

void opt_matrix_backsubstitute(opt_pk_internal_t *opk, opt_matrix_t *oc, size_t rank);


bool opt_matrix_is_bottom(opt_pk_internal_t *opk, opt_matrix_t * oc);


itv_t* opt_matrix_to_box(opt_pk_internal_t* opk,
		     opt_matrix_t* oc);



/* Predicates that can be useful for users */
//static inline 
//bool matrix_is_row_dummy_constraint(pk_internal_t* pk,
//				    matrix_t* mat,
//				    size_t l)
//{
//  return vector_is_dummy_constraint(pk,
//				    mat->p[l],mat->nbcolumns);
//}


/* Functions meant to be internal */
opt_matrix_t* _opt_matrix_alloc_int(size_t nr, unsigned short int nc, bool s);


opt_matrix_t* opt_matrix_assign_variable(opt_pk_internal_t* opk,
				 bool destructive, opt_matrix_t* mat,
				 elina_dim_t dim, opt_numint_t* tab);

opt_matrix_t* opt_matrix_assign_variables(opt_pk_internal_t* opk,
				  opt_matrix_t* mat,
				  elina_dim_t* tdim,
				  opt_numint_t** tvec,
				  size_t size);

opt_matrix_t* opt_matrix_substitute_variable(opt_pk_internal_t* opk,
				     bool destructive,
				     opt_matrix_t* mat,
				     elina_dim_t dim, opt_numint_t* tab);

void opt_generator_init(opt_pk_internal_t *opk, opt_matrix_t * mat, unsigned short int comp_size, size_t start);

/*******************************
	Remove common generators
********************************/
void remove_common_gen(opt_pk_internal_t *opk, opt_matrix_t * F, size_t start);

void opt_generator_bound_dimension(opt_pk_internal_t* opk,
			    itv_t itv,
			    elina_dim_t dim,
			    opt_matrix_t* of);

itv_t* opt_generator_to_box(opt_pk_internal_t* opk,
		     opt_matrix_t* of);

#ifdef __cplusplus
}
#endif

#endif

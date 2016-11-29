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
/* opt_pk_user.h: conversions with interface datatypes */
/* ********************************************************************** */

/* This header file defines operations on vectors. A vector is just an array of
   elements of type opt_pkint_t, as a consequence functions need to be given also
   their size. */

#ifndef __OPT_PK_USER_H__
#define __OPT_PK_USER_H__

#include "opt_pk_config.h"
#include "opt_pk_internal.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif


static void print_pk_array(opt_pk_array_t * op){
	unsigned short int nbcolumns = op->maxcols;
	array_comp_list_t * acl = op->acl;	
	unsigned short int num_comp = acl->size;
	opt_pk_t ** poly_array = op->poly;
	print_array_comp_list(acl,nbcolumns);
	unsigned short int k;
	for(k = 0; k < num_comp; k++){
		opt_pk_t * poly = poly_array[k];
		opt_matrix_t * mat = poly->C;
		opt_matrix_fprint(stdout,mat);
		fprintf(stdout,"\n");
	}
	fflush(stdout);
}

bool opt_vector_set_dim_bound(opt_pk_internal_t* opk,
			  opt_numint_t* vec,
			  elina_dim_t dim,
			  numrat_t numrat,
			  int mode,
			  size_t intdim, size_t realdim,
			  bool integer);


/* Fills the vector with the quasi-linear expression (itv_linexpr) */
void opt_vector_set_itv_linexpr(opt_pk_internal_t* opk,
			    opt_numint_t* ov,
			    itv_linexpr_t* expr,
			    size_t dim,
			    int mode);

/* Fills the vector(s) with the linear constraint cons */
void opt_vector_set_itv_lincons(opt_pk_internal_t* opk,
			    opt_numint_t* ov,
			    itv_lincons_t* cons,
			    size_t intdim, size_t realdim,
			    bool integer);

/* Fills the vector(s) with the linear constraint cons for testing
   satisfiability. Returns false if unsatisfiable
 */
bool opt_vector_set_itv_lincons_sat(opt_pk_internal_t* opk,
				opt_numint_t* ov,
				itv_lincons_t* cons,
				size_t intdim, size_t realdim,
				bool integer);



bool opt_matrix_set_itv_lincons_array(opt_pk_internal_t* opk,
				  opt_matrix_t** mat,
				  itv_lincons_array_t* array,
				  size_t intdim, size_t realdim,
				  bool integer);



elina_lincons0_t opt_lincons0_of_vector(opt_pk_internal_t* opk,
				 opt_numint_t* ov, unsigned short int * ca,
				 unsigned short int v_size, unsigned short int size);

#ifdef __cplusplus
}
#endif

#endif

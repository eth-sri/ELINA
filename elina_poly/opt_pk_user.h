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

elina_rat_t * elina_scalar_set_rat(elina_scalar_t * scalar);

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
			  elina_scalar_t * numrat,
			  int mode,
			  size_t intdim, size_t realdim,
			  bool integer);


/* Fills the vector with the quasi-linear expression (elina_linexpr) */
void opt_vector_set_elina_linexpr0(opt_pk_internal_t* opk,
			    opt_numint_t* ov,
			    elina_linexpr0_t* expr,
			    size_t dim,
			    int mode);

/* Fills the vector(s) with the linear constraint cons */
void opt_vector_set_elina_lincons0(opt_pk_internal_t* opk,
			    opt_numint_t* ov,
			    elina_lincons0_t* cons,
			    size_t intdim, size_t realdim,
			    bool integer);

/* Fills the vector(s) with the linear constraint cons for testing
   satisfiability. Returns false if unsatisfiable
 */
bool opt_vector_set_elina_lincons0_sat(opt_pk_internal_t* opk,
				opt_numint_t* ov,
				elina_lincons0_t* cons,
				size_t intdim, size_t realdim,
				bool integer);



bool opt_matrix_set_elina_lincons0_array(opt_pk_internal_t* opk,
				  opt_matrix_t** mat,
				  elina_lincons0_array_t* array,
				  size_t intdim, size_t realdim,
				  bool integer);



elina_lincons0_t opt_lincons0_of_vector(opt_pk_internal_t* opk,
				 opt_numint_t* ov, unsigned short int * ca,
				 unsigned short int v_size, unsigned short int size);



#ifdef __cplusplus
}
#endif

#endif

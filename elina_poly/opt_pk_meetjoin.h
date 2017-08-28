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


#ifndef _OPT_PK_MEETJOIN_H_
#define _OPT_PK_MEETJOIN_H_

#include "opt_pk_config.h"
#include "opt_pk.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ********************************************************************** */
/* I. Meet/Join */
/* ********************************************************************** */

bool opt_poly_meet_matrix(bool meet, bool lazy,
		      elina_manager_t* man,
		      opt_pk_t* op, 
		      opt_pk_t* o, opt_matrix_t* mat);

bool opt_poly_meet_elina_lincons_array(bool lazy,
				 elina_manager_t* man,
				 opt_pk_t* op, opt_pk_t* oa,
				 elina_lincons0_array_t* array);

void opt_poly_meet(bool meet,
	       bool lazy,
	       elina_manager_t* man,
	       opt_pk_array_t* op, opt_pk_array_t* oa, opt_pk_array_t* ob);

comp_list_t * lincons0_to_comp_list(opt_pk_internal_t * opk, elina_lincons0_t * cons);

comp_list_t * linexpr0_to_comp_list(opt_pk_internal_t * opk, elina_linexpr0_t * expr);

elina_linexpr0_t * copy_linexpr0_with_comp_list(opt_pk_internal_t *opk, elina_linexpr0_t * src, 
					     unsigned short int * ca, unsigned short int comp_size);


void copy_lincons0_with_comp_list(opt_pk_internal_t *opk, elina_lincons0_t * dst, 
				  elina_lincons0_t * src, unsigned short int * ca, unsigned short int comp_size);

bool is_linexpr_zero(elina_linexpr0_t * expr);

int elina_coeff_sgn(elina_coeff_t * coeff);

#ifdef __cplusplus
}
#endif

#endif

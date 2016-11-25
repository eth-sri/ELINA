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

bool opt_poly_meet_itv_lincons_array(bool lazy,
				 elina_manager_t* man,
				 opt_pk_t* op, opt_pk_t* oa,
				 itv_lincons_array_t* array);

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

#ifdef __cplusplus
}
#endif

#endif

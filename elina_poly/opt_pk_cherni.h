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

#ifndef __OPT_PK_CHERNI_H_
#define __OPT_PK_CHERNI_H_

#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_satmat.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"

#ifdef __cplusplus
extern "C"{
#endif

/* ********************************************************************** */
/* Conversion algorithm */
/* ********************************************************************** */

size_t opt_cherni_conversion(opt_pk_internal_t* opk,
			 opt_matrix_t* con, size_t start,
			 opt_matrix_t* ray, opt_satmat_t* osc, size_t nbline);


int  opt_cherni_simplify(opt_pk_internal_t* opk,
		     opt_matrix_t* con, opt_matrix_t* ray,
		     opt_satmat_t* satf, const size_t nbline);

void opt_cherni_minimize(opt_pk_internal_t* opk,
		     bool con_to_ray,
		     opt_pk_t* op);


void opt_cherni_add_and_minimize(opt_pk_internal_t* opk, 
			     bool con_to_ray,
			     opt_pk_t* op,
			     size_t start);


#ifdef __cplusplus
}
#endif

#endif



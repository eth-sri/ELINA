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
/* opt_pk_constructor.h: constructors and accessors */
/* ********************************************************************** */



#ifndef _OPT_PK_CONSTRUCTOR_H_
#define _OPT_PK_CONSTRUCTOR_H_

#include "opt_pk_config.h"
#include "opt_pk_vector.h"

#include "opt_pk_matrix.h"
#include "opt_pk.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Fill the first (opk->dec-1) rows of the matrix with the constraints of the
   universe polyhedron */
void opt_matrix_fill_constraint_top(opt_pk_internal_t* opk, opt_matrix_t* oc, size_t start);


void opt_poly_set_bottom(opt_pk_internal_t* opk, opt_pk_array_t* op);
void opt_poly_set_top(opt_pk_internal_t* opk, opt_pk_array_t* op);

#ifdef __cplusplus
}
#endif

#endif

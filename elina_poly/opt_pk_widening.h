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

#ifndef __OPT_PK_WIDENING_H
#define __OPT_PK_WIDENING_H

#include "opt_pk_config.h"
#include "opt_pk_internal.h"
#include "opt_pk.h"

#ifdef __cplusplus
extern "C" {
#endif

bool is_vectors_equal_comp_list(opt_pk_internal_t *opk, opt_numint_t * v1, 
				opt_numint_t * v2, unsigned short int * ca1, 
				unsigned short int * ca2, unsigned short int comp_size1, 
				unsigned short int comp_size2);

void vector_copy_comp_list(opt_pk_internal_t *opk, opt_numint_t * dst, opt_numint_t * src, 
			   unsigned short int * ind_map, unsigned short int comp_size);

#ifdef __cplusplus
}
#endif

#endif

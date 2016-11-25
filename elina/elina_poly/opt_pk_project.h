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


#ifndef __OPT_PK_PROJECT_H__
#define __OPT_PK_PROJECT_H__

#include "opt_pk_config.h"
#include "opt_pk_internal.h"
#include "opt_pk_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

void opt_poly_projectforget_array(bool project,
			      elina_manager_t* man,	
			      opt_pk_t* op, opt_pk_t* oa, 
			      elina_dim_t* tdim, size_t size, bool destructive);

opt_matrix_t * extreme_projection(opt_pk_internal_t *opk, 
				  opt_matrix_t *oc, 
				  elina_dim_t * tdim, size_t size, size_t proj);


#ifdef __cplusplus
}
#endif

#endif

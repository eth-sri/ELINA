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
/* opt_pk_assign.h: Assignements and Substitutions */
/* ********************************************************************** */



#ifndef _OPT_PK_ASSIGN_H_
#define _OPT_PK_ASSIGN_H_

#include "opt_pk_config.h"
#include "opt_pk.h"

#ifdef __cplusplus
extern "C" {
#endif

opt_pk_array_t* opt_poly_asssub_linexpr_det(bool assign, elina_manager_t* man,
			      bool destructive,
			      opt_pk_array_t* oa,
			      elina_dim_t dim, elina_linexpr0_t* linexpr);

opt_pk_array_t* opt_poly_asssub_linexpr_array_det(elina_manager_t* man,
				    bool destructive,
				    opt_pk_array_t* pa,
				    elina_dim_t* tdim, elina_linexpr0_t** texpr, 
				    size_t size);

#ifdef __cplusplus
}
#endif

#endif

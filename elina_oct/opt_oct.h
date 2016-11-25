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


#ifndef __OPT_OCT_H_INCLUDED__
#define __OPT_OCT_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#include "ap_generic.h"
#include "ap_coeff.h"
#include "ap_dimension.h"
#include "ap_expr0.h"
#include "ap_manager.h"



ap_manager_t* opt_oct_manager_alloc(void);

/* Enlarge each bound by epsilon times the maximum finite bound in 
     the octagon */

ap_abstract0_t* 
ap_abstract0_opt_oct_add_epsilon(ap_manager_t* man, 
			     ap_abstract0_t* a, 
			     ap_scalar_t* epsilon);
  
/* Enlarge each bound from a1 by epsilon times the maximum finite bound in 
     a2. Only those bounds in a1 that are not stable in a2 are enlared. */

ap_abstract0_t* 
ap_abstract0_opt_oct_add_epsilon_bin(ap_manager_t* man, 
				 ap_abstract0_t* a1, 
				 ap_abstract0_t* a2, 
				 ap_scalar_t* epsilon);
  


/* Widening with threshold.
     array is assumed to contain nb thresholds, sorted in increasing order. */
ap_abstract0_t* 
ap_abstract0_opt_oct_widening_thresholds(ap_manager_t* man,
				     ap_abstract0_t* a1,
				     ap_abstract0_t* a2,
				     ap_scalar_t** arr,
				     size_t nb );


/* Standard narrowing: refine only +oo constraint */
ap_abstract0_t* 
ap_abstract0_opt_oct_narrowing( ap_manager_t* man,
			    ap_abstract0_t* a1,
			    ap_abstract0_t* a2 );
  
#ifdef __cplusplus
 }
#endif

#endif

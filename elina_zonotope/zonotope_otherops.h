/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
 *  This software is distributed under GNU Lesser General Public License
 * Version 3.0. For more information, see the ELINA project website at:
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

#ifndef _ZONOTOPE_OTHEROPS_H_
#define _ZONOTOPE_OTHEROPS_H_

//#include "elina_var.h"
//#include "elina_abstract1.h"
#include "zonotope_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Other functions */
/*******************/

int elina_manager_zonotope_get_nsym(elina_manager_t* man);
 
//void elina_abstract1_aff_build(elina_manager_t* man, elina_abstract1_t * abstract, elina_var_t var, unsigned int index, elina_interval_t *itv, bool isunion);
//void elina_abstract1_ns_meet_lincons_array(elina_manager_t* man, elina_abstract1_t* abstract1, elina_lincons0_array_t* lincons);
//void elina_abstract1_ns_meet_box_array(elina_manager_t* man, elina_abstract1_t* abstract1, elina_interval_t** box, size_t size);

zonotope_t* zonotope_forget_array(elina_manager_t* man,
		bool destructive, zonotope_t* a,
		elina_dim_t* tdim, size_t size,
		bool project);

//zonotope_t zonotope_expand(elina_manager_t* man,
//		bool destructive, zonotope_t* a,
//		elina_dim_t var,
//		elina_dim_t* tvar, size_t size);

//zonotope_t zonotope_fold(elina_manager_t* man,
//		bool destructive, zonotope_t* a,
//		elina_dim_t* tvar, size_t size);

//zonotope_t* zonotope_widening(elina_manager_t* man,
//		  zonotope_t* a1, 
//		  zonotope_t* a2);

//zonotope_t zonotope_closure(elina_manager_t* man, bool destructive, zonotope_t* a);

/* retourne [0,+oo] */
/* for internal use */
/* pour test */
//elina_interval_t* zonotope_create_pos(elina_manager_t* man);

#ifdef __cplusplus
}
#endif

#endif


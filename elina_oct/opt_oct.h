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



#ifndef __OPT_OCT_H_INCLUDED__
#define __OPT_OCT_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#if defined (HAS_APRON)
#include "ap_generic.h"
#include "ap_coeff.h"
#include "ap_dimension.h"
#include "ap_expr0.h"
#include "ap_manager.h"

#else
#include "elina_generic.h"
#include "elina_coeff.h"
#include "elina_dimension.h"
#include "elina_manager.h"
#endif


elina_manager_t* opt_oct_manager_alloc(void);

/* Enlarge each bound by epsilon times the maximum finite bound in 
     the octagon */

elina_abstract0_t* 
elina_abstract0_opt_oct_add_epsilon(elina_manager_t* man, 
			     elina_abstract0_t* a, 
			     elina_scalar_t* epsilon);
  
/* Enlarge each bound from a1 by epsilon times the maximum finite bound in 
     a2. Only those bounds in a1 that are not stable in a2 are enlared. */

elina_abstract0_t* 
elina_abstract0_opt_oct_add_epsilon_bin(elina_manager_t* man, 
				 elina_abstract0_t* a1, 
				 elina_abstract0_t* a2, 
				 elina_scalar_t* epsilon);
  


/* Widening with threshold.
     array is assumed to contain nb thresholds, sorted in increasing order. */
elina_abstract0_t* 
elina_abstract0_opt_oct_widening_thresholds(elina_manager_t* man,
				     elina_abstract0_t* a1,
				     elina_abstract0_t* a2,
				     elina_scalar_t** arr,
				     size_t nb );


/* Standard narrowing: refine only +oo constraint */
elina_abstract0_t* 
elina_abstract0_opt_oct_narrowing(elina_manager_t* man,
			    elina_abstract0_t* a1,
			    elina_abstract0_t* a2 );
  
#ifdef __cplusplus
 }
#endif

#endif

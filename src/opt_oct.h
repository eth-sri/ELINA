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

ap_abstract0_t* 
ap_abstract0_opt_oct_add_epsilon(ap_manager_t* man, 
			     ap_abstract0_t* a, 
			     ap_scalar_t* epsilon);
  /* Enlarge each bound by epsilon times the maximum finite bound in 
     the octagon */

ap_abstract0_t* 
ap_abstract0_opt_oct_add_epsilon_bin(ap_manager_t* man, 
				 ap_abstract0_t* a1, 
				 ap_abstract0_t* a2, 
				 ap_scalar_t* epsilon);
  /* Enlarge each bound from a1 by epsilon times the maximum finite bound in 
     a2. Only those bounds in a1 that are not stable in a2 are enlared. */

#ifdef __cplusplus
 }
#endif

#endif

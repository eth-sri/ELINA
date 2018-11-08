/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
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

#ifndef _ELINA_BOX_MEETJOIN_H_
#define _ELINA_BOX_MEETJOIN_H_

#include "elina_box_internal.h"
#include "elina_linearize_texpr.h"

#ifdef __cplusplus
extern "C" {
#endif

elina_box_t* elina_box_meet(elina_manager_t* man, bool destructive, elina_box_t* a1, elina_box_t* a2);
elina_box_t* elina_box_join(elina_manager_t* man, bool destructive, elina_box_t* a1, elina_box_t* a2);
elina_box_t* elina_box_meet_lincons_array(elina_manager_t* man,
			      bool destructive,
			      elina_box_t* a,
			      elina_lincons0_array_t* array);
    
elina_box_t* elina_box_widening(elina_manager_t* man,
                        elina_box_t* a1, elina_box_t* a2);

bool elina_double_interval_eval_elina_linexpr0(double * itv_inf, double * itv_sup, elina_linexpr0_t* expr, double* env_inf, double * env_sup, elina_scalar_discr_t discr);

bool elina_double_boxize_lincons0_array(double * res_inf, double * res_sup, bool* tchange,
				      elina_lincons0_array_t* array,
				      double* env_inf, double * env_sup, size_t intdim,
				      size_t kmax,
				      bool intervalonly, elina_scalar_discr_t discr);

void elina_double_interval_mul(double *a_inf, double *a_sup, double b_inf, double b_sup, double c_inf, double c_sup);

void elina_double_interval_div(double *a_inf, double *a_sup, double b_inf, double b_sup, double c_inf, double c_sup);

#ifdef __cplusplus
}
#endif

#endif

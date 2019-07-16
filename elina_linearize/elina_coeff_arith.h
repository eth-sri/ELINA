/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
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


#ifndef _ELINA_COEFF_ARITH_H_
#define _ELINA_COEFF_ARITH_H_

#include "elina_interval_arith.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "elina_interval_arith.h"
#if defined (HAS_APRON)
#include "apron_wrapper.h"
#else
#include "elina_coeff.h"
#endif
/***********************************************
	Arithmetic on elina_coeff
************************************************/
void elina_coeff_mul_scalar(elina_coeff_t * dst, elina_coeff_t * src, elina_scalar_t * mul, elina_scalar_discr_t discr);

void elina_coeff_add(elina_coeff_t * op, elina_coeff_t * op1, elina_coeff_t * op2, elina_scalar_discr_t discr);

void elina_coeff_sub_num(elina_coeff_t * dst, elina_coeff_t * src, elina_scalar_t * sub, elina_scalar_discr_t discr);

void elina_coeff_mul_interval(elina_coeff_t * dst, elina_coeff_t *src, elina_interval_t * interval, elina_scalar_discr_t discr);

void elina_interval_set_elina_coeff(elina_interval_t *interval, elina_coeff_t *coeff);

#ifdef __cplusplus
}
#endif

#endif

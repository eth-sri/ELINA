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

#ifndef _ELINA_INTERVAL_ARITH_H_
#define _ELINA_INTERVAL_ARITH_H_

#include "elina_scalar_arith.h"

#if defined (HAS_APRON)
#include "apron_wrapper.h"
#else
#include "elina_interval.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ************************************************************************* */
/* Interval Arithmetic */
/* ************************************************************************* */

void elina_interval_add(elina_interval_t * op, elina_interval_t *op1, elina_interval_t *op2, elina_scalar_discr_t discr);

void elina_interval_mul_scalar(elina_interval_t * dst, elina_interval_t * src, elina_scalar_t * mul, elina_scalar_discr_t discr);

void elina_interval_mul(elina_interval_t *op, elina_interval_t *op1, elina_interval_t *op2, elina_scalar_discr_t discr);

void elina_interval_sub(elina_interval_t* a, elina_interval_t *b, elina_interval_t* c, elina_scalar_discr_t discr);

bool elina_interval_sqrt(elina_interval_t *dst, elina_interval_t *src, elina_scalar_discr_t discr);

void elina_interval_pow(elina_interval_t *a, elina_interval_t *b, elina_interval_t *n, elina_scalar_discr_t discr);

void elina_interval_abs(elina_interval_t *a, elina_interval_t *b,  elina_scalar_discr_t discr); 

void elina_interval_trunc(elina_interval_t *a, elina_interval_t *b, elina_scalar_discr_t discr);

void elina_interval_ceil(elina_interval_t *a, elina_interval_t *b, elina_scalar_discr_t discr);

void elina_interval_floor(elina_interval_t *a, elina_interval_t *b, elina_scalar_discr_t discr);

void elina_interval_to_int(elina_interval_t *a, elina_interval_t *b, elina_scalar_discr_t discr);

void elina_interval_to_float(elina_interval_t *a, elina_interval_t *b, elina_scalar_discr_t discr);

void elina_interval_to_double(elina_interval_t *a, elina_interval_t *b, elina_scalar_discr_t discr);

void elina_interval_div(elina_interval_t *a, elina_interval_t *b, elina_interval_t *c, elina_scalar_discr_t discr);

void elina_interval_mod(elina_interval_t *a, elina_interval_t *b, elina_interval_t *c, bool is_int, elina_scalar_discr_t discr);

bool elina_interval_is_int(elina_interval_t *a, elina_scalar_discr_t discr);

void elina_interval_magnitude(elina_scalar_t *a, elina_interval_t *b);

void elina_interval_range_rel(elina_scalar_t *a, elina_interval_t *b, elina_scalar_discr_t discr);

void elina_interval_enlarge_bound(elina_interval_t *a, elina_interval_t *b, elina_scalar_t *c,elina_scalar_discr_t discr);

bool elina_interval_canonicalize(elina_interval_t *a, bool integer, elina_scalar_discr_t discr);

void elina_interval_set_to_int(elina_interval_t *a, elina_int_t i, elina_int_t j, elina_scalar_discr_t discr);

void elina_interval_mul_2exp(elina_interval_t *a, elina_interval_t *b, int i, elina_scalar_discr_t discr);

void elina_interval_convert(elina_interval_t * a, elina_scalar_discr_t discr);

#ifdef __cplusplus
}
#endif

#endif

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


#ifndef _ELINA_SCALAR_ARITH_H_
#define _ELINA_SCALAR_ARITH_H_


#include <float.h>
#include "gmp.h"
#include "mpfr.h"

#if defined (HAS_APRON)
#include "apron_wrapper.h"
#else
#include "elina_scalar.h"
#endif

#include "elina_rat.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ************************************************************************* */
/* Scalar Arithmetic */
/* ************************************************************************* */

void elina_scalar_add_uint(elina_scalar_t *a, elina_scalar_t * b, unsigned long int c, elina_scalar_discr_t discr);

void elina_scalar_add(elina_scalar_t *op, elina_scalar_t *op1, elina_scalar_t *op2, elina_scalar_discr_t discr);

void elina_scalar_mul(elina_scalar_t *op, elina_scalar_t *op1, elina_scalar_t *op2, elina_scalar_discr_t discr);

void elina_scalar_sub_uint(elina_scalar_t *a, elina_scalar_t * b, unsigned long int c, elina_scalar_discr_t discr);

void elina_scalar_sub(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_t *c, elina_scalar_discr_t discr);

void elina_scalar_div_2(elina_scalar_t * dst, elina_scalar_t *src);

void elina_scalar_div(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_t *c, elina_scalar_discr_t discr);

void elina_scalar_max(elina_scalar_t *op, elina_scalar_t * op1, elina_scalar_t *op2);

void elina_scalar_min(elina_scalar_t *op, elina_scalar_t * op1, elina_scalar_t *op2);

void elina_scalar_sqrt(elina_scalar_t *up, elina_scalar_t *down, elina_scalar_t *b,elina_scalar_discr_t discr);

void int_set_elina_scalar(long int *a, elina_scalar_t * scalar);

void elina_scalar_pow(elina_scalar_t *up, elina_scalar_t *down, elina_scalar_t *b, unsigned long n,elina_scalar_discr_t discr);

void elina_scalar_ceil(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_discr_t discr);

void elina_scalar_floor(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_discr_t discr);

void elina_scalar_trunc(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_discr_t discr);

void elina_scalar_to_double(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_discr_t discr);

void elina_scalar_to_float(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_discr_t discr);

bool elina_scalar_is_integer(elina_scalar_t *a);

void elina_scalar_set_to_int(elina_scalar_t *a, elina_int_t i, elina_scalar_discr_t discr);

void elina_scalar_mul_2exp(elina_scalar_t *a, elina_scalar_t *b, int c , elina_scalar_discr_t discr);

void elina_scalar_convert(elina_scalar_t *a, elina_scalar_discr_t discr);

#ifdef __cplusplus
}
#endif

#endif

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


#ifndef _ELINA_LINEXPR0_ARITH_H_
#define _ELINA_LINEXPR0_ARITH_H_

#if defined (HAS_APRON)
#include "apron_wrapper.h"
#else
#include "elina_linexpr0.h"
#endif

#include "elina_coeff_arith.h"

#ifdef __cplusplus
extern "C" {
#endif

void elina_linexpr0_reinit(elina_linexpr0_t* expr, size_t size);

void elina_linexpr0_init(elina_linexpr0_t* expr, size_t size);

void elina_linexpr0_clear(elina_linexpr0_t* e);

void elina_linexpr0_neg(elina_linexpr0_t* expr);

void elina_linexpr0_scale(elina_linexpr0_t* expr, elina_interval_t *interval, elina_scalar_discr_t discr);

void elina_linexpr0_add(elina_linexpr0_t **res, elina_linexpr0_t **exprA, elina_linexpr0_t **exprB, elina_scalar_discr_t discr);

void elina_linexpr0_sub(elina_linexpr0_t **res, elina_linexpr0_t **exprA, elina_linexpr0_t **exprB, elina_scalar_discr_t discr);

void elina_linexpr0_div(elina_linexpr0_t* expr, elina_interval_t *interval, elina_scalar_discr_t discr);

#ifdef __cplusplus
}
#endif

#endif

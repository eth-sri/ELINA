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

/* ************************************************************************* */
/* elina_coeff.h: coefficients, that are either scalars or intervals */
/* ************************************************************************* */


#ifndef _ELINA_COEFF_H_
#define _ELINA_COEFF_H_

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "elina_config.h"
#include "elina_scalar.h"
#include "elina_interval.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum elina_coeff_discr_t {
  ELINA_COEFF_SCALAR,
  ELINA_COEFF_INTERVAL
} elina_coeff_discr_t;
  /* Discriminant for coefficients */

typedef struct elina_coeff_t {
  elina_coeff_discr_t discr; /* discriminant for coefficient */
  union {
    elina_scalar_t* scalar;       /* cst (normal linear expression) */
    elina_interval_t* interval;   /* interval (quasi-linear expression) */
  } val;
} elina_coeff_t;

/* ====================================================================== */
/* Basics */
/* ====================================================================== */

elina_coeff_t* elina_coeff_alloc(elina_coeff_discr_t elina_coeff_discr);
  /* Initialization, specifying the type of the coefficient */
void elina_coeff_reinit(elina_coeff_t* coeff, elina_coeff_discr_t elina_coeff_discr, elina_scalar_discr_t elina_scalar_discr);
  /* Changing the type of scalar(s) and the type of the coefficient */
void elina_coeff_free(elina_coeff_t* a);
  /* Free a coefficient */
void elina_coeff_fprint(FILE* stream, elina_coeff_t* a);
static inline void elina_coeff_print(elina_coeff_t *a) {
  elina_coeff_fprint(stdout, a);
}
  /* Printing */

void elina_coeff_reduce(elina_coeff_t* coeff);
  /* If the coefficient is an interval [a;a], convert it to a scalar */
static inline void elina_coeff_swap(elina_coeff_t *a, elina_coeff_t *b) {
  elina_coeff_t t = *a;
  *a = *b;
  *b = t;
}
  /* Exchange */

/* ====================================================================== */
/* Assignments */
/* ====================================================================== */

void elina_coeff_set(elina_coeff_t* a, elina_coeff_t* b);
  /* Assignment */
void elina_coeff_set_scalar(elina_coeff_t* coeff, elina_scalar_t* scalar);
void elina_coeff_set_scalar_mpq(elina_coeff_t* coeff, mpq_t mpq);
void elina_coeff_set_scalar_int(elina_coeff_t* coeff, long int num);
void elina_coeff_set_scalar_frac(elina_coeff_t* coeff, long int num, unsigned long int den);
void elina_coeff_set_scalar_double(elina_coeff_t* coeff, double num);
void elina_coeff_set_scalar_mpfr(elina_coeff_t* coeff, mpfr_t mpfr);
  /* Assign a coefficient of type SCALAR, with resp.
     - a coeff
     - a rational of type mpq_t, converted to type MPQ
     - an integer, converted to type MPQ
     - a rational, converted to type MPQ
     - a double, converted to type DOUBLE
     - a MPFR, converted to type MPFR
  */
void elina_coeff_set_interval(elina_coeff_t* coeff, elina_interval_t* itv);
void elina_coeff_set_interval_scalar(elina_coeff_t* coeff, elina_scalar_t* inf, elina_scalar_t* sup);
void elina_coeff_set_interval_mpq(elina_coeff_t* coeff, mpq_t inf, mpq_t sup);
void elina_coeff_set_interval_int(elina_coeff_t* coeff, long int inf, long int sup);
void elina_coeff_set_interval_frac(elina_coeff_t* coeff,
                                  long int numinf, unsigned long int deninf,
                                  long int numsup, unsigned long int densup);
void elina_coeff_set_interval_double(elina_coeff_t* coeff, double inf, double sup);
void elina_coeff_set_interval_top(elina_coeff_t* coeff);
void elina_coeff_set_interval_mpfr(elina_coeff_t* coeff, mpfr_t inf, mpfr_t sup);
  /* Assign a coefficient of type INTERVAL, with resp.
     - an interval of coeff
     - an interval of rationals of type MPQ
     - an interval of integers, converted to type MPQ
     - an interval of rationals, converted to type MPQ
     - an interval of double, converted to type DOUBLE
     - an interval of MPFR, converted to type MPFR
     - a top interval (type not precised).
  */

/* ====================================================================== */
/* Combined allocation and assignment */
/* ====================================================================== */
elina_coeff_t* elina_coeff_alloc_set(elina_coeff_t* coeff);
elina_coeff_t* elina_coeff_alloc_set_scalar(elina_scalar_t* scalar);
elina_coeff_t* elina_coeff_alloc_set_interval(elina_interval_t* interval);

/* ====================================================================== */
/* Tests */
/* ====================================================================== */

int elina_coeff_cmp(elina_coeff_t* coeff1, elina_coeff_t* coeff2);
  /* Non Total Comparison:
     - If the 2 coefficients are both scalars, corresp. to elina_scalar_cmp
     - If the 2 coefficients are both intervals, corresp. to elina_interval_cmp
     - otherwise, -3 if the first is a scalar, 3 otherwise
  */
bool elina_coeff_equal(elina_coeff_t* coeff1, elina_coeff_t* coeff2);
  /* Equality */

bool elina_coeff_zero(elina_coeff_t* coeff);
  /* Return true iff coeff is a zero scalar or an interval with zero bounds */
bool elina_coeff_equal_int(elina_coeff_t* coeff, int i);
  /* Return true iff coeff is a scalar equals to i or an interval with bounds equal to i */

/* ====================================================================== */
/* Other operations */
/* ====================================================================== */
void elina_coeff_neg(elina_coeff_t* a, elina_coeff_t* b);
  /* Negation */

long elina_coeff_hash(elina_coeff_t* coeff);
  /* Hash code */

/* ====================================================================== */
/* FOR INTERNAL USE ONLY */
/* ====================================================================== */
void elina_coeff_init(elina_coeff_t* coeff, elina_coeff_discr_t elina_coeff_discr);
void elina_coeff_init_set(elina_coeff_t* coeff, elina_coeff_t* coeff2);
void elina_coeff_clear(elina_coeff_t* coeff);

#ifdef __cplusplus
}
#endif

#endif

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

/* ************************************************************************* */
/* elina_interval.h: intervals */
/* ************************************************************************* */


#ifndef _ELINA_INTERVAL_H_
#define _ELINA_INTERVAL_H_

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "elina_config.h"
#include "elina_scalar.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct elina_interval_t {
  elina_scalar_t* inf;
  elina_scalar_t* sup;
} elina_interval_t;

/* ====================================================================== */
/* Basics */
/* ====================================================================== */
elina_interval_t* elina_interval_alloc(void);
  /* Initialization, using DOUBLE as default type for scalars */
void elina_interval_reinit(elina_interval_t* interval, elina_scalar_discr_t elina_scalar_discr);
  /* Change the type of scalars */
void elina_interval_free(elina_interval_t* interval);
  /* Free an interval */
void elina_interval_fprint(FILE* stream, elina_interval_t* a);

void elina_interval_print(elina_interval_t* a);

  /* Printing */
void elina_interval_swap(elina_interval_t* a, elina_interval_t* b);

  /* Exchange */

/* ====================================================================== */
/* Assignments */
/* ====================================================================== */

void elina_interval_set(elina_interval_t* interval, elina_interval_t* interval2);
  /* Assignment */
void elina_interval_set_scalar(elina_interval_t* interval, elina_scalar_t* inf, elina_scalar_t* sup);
void elina_interval_set_mpq(elina_interval_t* interval, mpq_t inf, mpq_t sup);
void elina_interval_set_int(elina_interval_t* interval, long int inf, long int sup);
void elina_interval_set_frac(elina_interval_t* interval, long int numinf, unsigned long int deninf, long int numsup, unsigned long int densup);
void elina_interval_set_double(elina_interval_t* interval, double inf, double sup);
void elina_interval_set_mpfr(elina_interval_t* interval, mpfr_t inf, mpfr_t sup);
  /* Assignment from resp.
     - two scalars
     - two rationals of type MPQ
     - two integers, giving [inf,dup]
     - two rationals, giving [numinf/deninf,numsup/densup]
     - two double values
     - two MPFR floating-point numbers
  */
void elina_interval_set_top(elina_interval_t* interval);
  /* Assignment to universe interval [-oo,oo],
     does not change the type of scalars */
void elina_interval_set_bottom(elina_interval_t* interval);
  /* Assignment to empty interval [1,-1],
     does not change the type of scalars */

/* ====================================================================== */
/* Combined allocation and assignments */
/* ====================================================================== */

elina_interval_t* elina_interval_alloc_set(elina_interval_t* interval);
  /* Assignment */

/* ====================================================================== */
/* Tests */
/* ====================================================================== */

bool elina_interval_is_top(elina_interval_t* interval);
  /* Is it the universe interval ? */
bool elina_interval_is_bottom(elina_interval_t* interval);
  /* Is it an empty interval ? */
bool elina_interval_is_leq(elina_interval_t* i1, elina_interval_t* i2);
  /* Inclusion test */
int elina_interval_cmp(elina_interval_t* i1, elina_interval_t* i2);
  /* Comparison:
     0: equality
     -1: i1 included in i2
     +1: i2 included in i1
     -2: i1.inf less than i2.inf
     +2: i1.inf greater than i2.inf
  */
bool elina_interval_equal(elina_interval_t* i1, elina_interval_t* i2);
bool elina_interval_equal_int(elina_interval_t* i, int b);
  /* Equality */

/* ====================================================================== */
/* Other operations */
/* ====================================================================== */

void elina_interval_neg(elina_interval_t* a, elina_interval_t* b);
  /* Negation */
long elina_interval_hash(elina_interval_t* itv);
  /* Hash code */

/* ====================================================================== */
/* Array of intervals */
/* ====================================================================== */

elina_interval_t** elina_interval_array_alloc(size_t size);
  /* Allocating an array of intervals, initialized with [0,0] values */
void elina_interval_array_free(elina_interval_t** array, size_t size);
  /* Clearing and freeing an array of intervals */

#ifdef __cplusplus
}
#endif

#endif

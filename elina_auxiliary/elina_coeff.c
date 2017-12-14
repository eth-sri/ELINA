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
/* elina_coeff.c: coefficients, that are either scalars or intervals */
/* ************************************************************************* */

#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include "elina_coeff.h"

/* ====================================================================== */
/* Basics */
/* ====================================================================== */

/* FOR INTERNAL USE ONLY */
void elina_coeff_init(elina_coeff_t* coeff, elina_coeff_discr_t coeff_discr)
{
  coeff->discr = coeff_discr;
  switch (coeff_discr){
  case ELINA_COEFF_SCALAR:
    coeff->val.scalar = elina_scalar_alloc();
    break;
  case ELINA_COEFF_INTERVAL:
    coeff->val.interval = elina_interval_alloc();
    break;
  }
}
/* FOR INTERNAL USE ONLY */
void elina_coeff_init_set(elina_coeff_t* coeff, elina_coeff_t* coeff2)
{
  coeff->discr = coeff2->discr;
  switch (coeff2->discr){
  case ELINA_COEFF_SCALAR:
    coeff->val.scalar = elina_scalar_alloc_set(coeff2->val.scalar);
    break;
  case ELINA_COEFF_INTERVAL:
    coeff->val.interval = elina_interval_alloc_set(coeff2->val.interval);
    break;
  }
}
/* FOR INTERNAL USE ONLY */
void elina_coeff_clear(elina_coeff_t* coeff)
{
  switch (coeff->discr){
  case ELINA_COEFF_SCALAR:
    elina_scalar_free(coeff->val.scalar);
    break;
  case ELINA_COEFF_INTERVAL:
    elina_interval_free(coeff->val.interval);
    break;
  }
}

elina_coeff_t* elina_coeff_alloc(elina_coeff_discr_t coeff_discr)
{
  elina_coeff_t* coeff = malloc(sizeof(elina_coeff_t));
  coeff->discr = coeff_discr;
  switch (coeff_discr){
  case ELINA_COEFF_SCALAR:
    coeff->val.scalar = elina_scalar_alloc();
    break;
  case ELINA_COEFF_INTERVAL:
    coeff->val.interval = elina_interval_alloc();
    break;
  }
  return coeff;
}
void elina_coeff_reinit(elina_coeff_t* coeff, elina_coeff_discr_t coeff_discr, elina_scalar_discr_t elina_scalar_discr)
{
  switch (coeff->discr){
  case ELINA_COEFF_SCALAR:
    if (coeff_discr != coeff->discr){
      elina_scalar_free(coeff->val.scalar);
      coeff->val.interval = elina_interval_alloc();
      elina_interval_reinit(coeff->val.interval,elina_scalar_discr);
    }
    else {
      elina_scalar_reinit(coeff->val.scalar,elina_scalar_discr);
    }
    break;
  case ELINA_COEFF_INTERVAL:
    if (coeff_discr != coeff->discr){
      elina_interval_free(coeff->val.interval);
      coeff->val.scalar = elina_scalar_alloc();
      elina_scalar_reinit(coeff->val.scalar,elina_scalar_discr);
    }
    else {
      elina_interval_reinit(coeff->val.interval,elina_scalar_discr);
    }
    break;
  }
  coeff->discr = coeff_discr;
}

void elina_coeff_free(elina_coeff_t* coeff)
{
  switch (coeff->discr){
  case ELINA_COEFF_SCALAR:
    elina_scalar_free(coeff->val.scalar);
    break;
  case ELINA_COEFF_INTERVAL:
    elina_interval_free(coeff->val.interval);
    break;
  }
  free(coeff);
}

void elina_coeff_fprint(FILE* stream, elina_coeff_t* a)
{
  switch(a->discr){
  case ELINA_COEFF_SCALAR:
    elina_scalar_fprint(stream,a->val.scalar);
    break;
  case ELINA_COEFF_INTERVAL:
    elina_interval_fprint(stream,a->val.interval);
    break;
  }
}

void elina_coeff_reduce(elina_coeff_t* coeff)
{
  if (coeff->discr==ELINA_COEFF_INTERVAL){
    if (elina_scalar_equal(coeff->val.interval->inf,coeff->val.interval->sup)){
      /* We cheat with the good rules */
      elina_scalar_t* scalar = coeff->val.interval->inf;
      elina_scalar_free(coeff->val.interval->sup);
      free(coeff->val.interval);
      coeff->val.scalar = scalar;
      coeff->discr=ELINA_COEFF_SCALAR;
    }
  }
}

/* ====================================================================== */
/* Combined allocation and assignment */
/* ====================================================================== */
elina_coeff_t* elina_coeff_alloc_set_scalar(elina_scalar_t* scalar)
{
  elina_coeff_t* coeff = malloc(sizeof(elina_coeff_t));
  coeff->discr = ELINA_COEFF_SCALAR;
  coeff->val.scalar = elina_scalar_alloc_set(scalar);
  return coeff;
}
elina_coeff_t* elina_coeff_alloc_set_interval(elina_interval_t* interval)
{
  elina_coeff_t* coeff = malloc(sizeof(elina_coeff_t));
  coeff->discr = ELINA_COEFF_INTERVAL;
  coeff->val.interval = elina_interval_alloc_set(interval);
  return coeff;
}
elina_coeff_t* elina_coeff_alloc_set(elina_coeff_t* coeff)
{
  switch (coeff->discr){
  case ELINA_COEFF_SCALAR:
    return elina_coeff_alloc_set_scalar(coeff->val.scalar);
  case ELINA_COEFF_INTERVAL:
    return elina_coeff_alloc_set_interval(coeff->val.interval);
    break;
  default:
    return NULL;
  }
}

/* ====================================================================== */
/* Assignments */
/* ====================================================================== */

void elina_coeff_set(elina_coeff_t* a, elina_coeff_t* b)
{

  if(a->discr != b->discr){
    elina_coeff_reinit(a,b->discr,ELINA_SCALAR_DOUBLE);
  }
  switch(b->discr){
  case ELINA_COEFF_SCALAR:
    elina_scalar_set(a->val.scalar,b->val.scalar);
    break;
  case ELINA_COEFF_INTERVAL:
    elina_interval_set(a->val.interval,b->val.interval);
    break;
  }
}
void elina_coeff_set_scalar(elina_coeff_t* coeff, elina_scalar_t* scalar)
{
  elina_coeff_reinit(coeff,ELINA_COEFF_SCALAR,scalar->discr); 
  elina_scalar_set(coeff->val.scalar,scalar); 
}
void elina_coeff_set_scalar_mpq(elina_coeff_t* coeff, mpq_t mpq)
{
  elina_coeff_reinit(coeff,ELINA_COEFF_SCALAR,ELINA_SCALAR_MPQ); 
  elina_scalar_set_mpq(coeff->val.scalar,mpq); 
}
void elina_coeff_set_scalar_mpfr(elina_coeff_t* coeff, mpfr_t mpfr)
{
  elina_coeff_reinit(coeff,ELINA_COEFF_SCALAR,ELINA_SCALAR_MPFR); 
  elina_scalar_set_mpfr(coeff->val.scalar,mpfr); 
}
void elina_coeff_set_scalar_int(elina_coeff_t* coeff, long int num)
{
  elina_coeff_reinit(coeff,ELINA_COEFF_SCALAR,ELINA_SCALAR_MPQ); 
  elina_scalar_set_int(coeff->val.scalar,num); 
}
void elina_coeff_set_scalar_frac(elina_coeff_t* coeff, long int num, unsigned long int den)
{
  elina_coeff_reinit(coeff,ELINA_COEFF_SCALAR,ELINA_SCALAR_MPQ); 
  elina_scalar_set_frac(coeff->val.scalar,num,den); 
}
void elina_coeff_set_scalar_double(elina_coeff_t* coeff, double num)
{ 
  elina_coeff_reinit(coeff,ELINA_COEFF_SCALAR,ELINA_SCALAR_DOUBLE); 
  elina_scalar_set_double(coeff->val.scalar,num); 
}
void elina_coeff_set_interval(elina_coeff_t* coeff, elina_interval_t* itv)
{
  elina_coeff_reinit(coeff,ELINA_COEFF_INTERVAL,ELINA_SCALAR_DOUBLE); 
  elina_interval_set(coeff->val.interval,itv); 
}
void elina_coeff_set_interval_scalar(elina_coeff_t* coeff, elina_scalar_t* inf, elina_scalar_t* sup)
{
  elina_coeff_reinit(coeff,ELINA_COEFF_INTERVAL,ELINA_SCALAR_DOUBLE);
  elina_interval_set_scalar(coeff->val.interval,inf,sup);
}
void elina_coeff_set_interval_mpq(elina_coeff_t* coeff, mpq_t inf, mpq_t sup)
{
  elina_coeff_reinit(coeff,ELINA_COEFF_INTERVAL,ELINA_SCALAR_MPQ); 
  elina_interval_set_mpq(coeff->val.interval,inf,sup);
}
void elina_coeff_set_interval_mpfr(elina_coeff_t* coeff, mpfr_t inf, mpfr_t sup)
{
  elina_coeff_reinit(coeff,ELINA_COEFF_INTERVAL,ELINA_SCALAR_MPFR); 
  elina_interval_set_mpfr(coeff->val.interval,inf,sup);
}
void elina_coeff_set_interval_int(elina_coeff_t* coeff, long int inf, long int sup)
{
  elina_coeff_reinit(coeff,ELINA_COEFF_INTERVAL,ELINA_SCALAR_MPQ); 
  elina_interval_set_int(coeff->val.interval,inf,sup); 
}
void elina_coeff_set_interval_frac(elina_coeff_t* coeff,
				long int numinf, unsigned long int deninf, 
				long int numsup, unsigned long int densup)
{
  elina_coeff_reinit(coeff,ELINA_COEFF_INTERVAL,ELINA_SCALAR_MPQ); 
  elina_interval_set_frac(coeff->val.interval,numinf,deninf,numsup,densup); 
}
void elina_coeff_set_interval_double(elina_coeff_t* coeff, double inf, double sup)
{
  elina_coeff_reinit(coeff,ELINA_COEFF_INTERVAL,ELINA_SCALAR_DOUBLE); 
  elina_interval_set_double(coeff->val.interval,inf,sup); 
}
void elina_coeff_set_interval_top(elina_coeff_t* coeff)
{
  if (coeff->discr == ELINA_COEFF_SCALAR)
    elina_coeff_reinit(coeff,ELINA_COEFF_INTERVAL,ELINA_SCALAR_DOUBLE);
  elina_interval_set_top(coeff->val.interval);
}

/* ====================================================================== */
/* Tests */
/* ====================================================================== */

int elina_coeff_cmp(elina_coeff_t* coeff1, elina_coeff_t* coeff2)
{
  if (coeff1->discr==coeff2->discr){
    switch (coeff1->discr){
    case ELINA_COEFF_SCALAR:
      return elina_scalar_cmp(coeff1->val.scalar,coeff2->val.scalar);
    case ELINA_COEFF_INTERVAL:
      return elina_interval_cmp(coeff1->val.interval,coeff2->val.interval);
    default:
      abort();
      return 0;
    }
  } 
  else {
    return (coeff1->discr==ELINA_COEFF_SCALAR) ? -3 : 3;
  }
}
bool elina_coeff_equal(elina_coeff_t* coeff1, elina_coeff_t* coeff2)
{
  if (coeff1->discr==coeff2->discr){
    switch (coeff1->discr){
    case ELINA_COEFF_SCALAR:
      return elina_scalar_equal(coeff1->val.scalar,coeff2->val.scalar);
    case ELINA_COEFF_INTERVAL:
      return elina_interval_equal(coeff1->val.interval,coeff2->val.interval);
    default:
      abort();
    }
  }
  else
    return false;
}
bool elina_coeff_zero(elina_coeff_t* coeff)
{
  switch (coeff->discr){
  case ELINA_COEFF_SCALAR:
    return elina_scalar_sgn(coeff->val.scalar)==0;
  case ELINA_COEFF_INTERVAL:
    return (elina_scalar_sgn(coeff->val.interval->inf)==0) && (elina_scalar_sgn(coeff->val.interval->sup)==0);
  default:
    abort();
  }
}
bool elina_coeff_equal_int(elina_coeff_t* coeff, int i)
{
  switch (coeff->discr){
  case ELINA_COEFF_SCALAR:
    return elina_scalar_equal_int(coeff->val.scalar,i)==0;
  case ELINA_COEFF_INTERVAL:
    return (elina_scalar_equal_int(coeff->val.interval->inf,i)==0) && (elina_scalar_equal_int(coeff->val.interval->sup,i)==0);
  default:
    abort();
  }
}

/* ====================================================================== */
/* Other operations */
/* ====================================================================== */

void elina_coeff_neg(elina_coeff_t* a, elina_coeff_t* b)
{
  elina_coeff_set(a,b);
  switch(b->discr){
  case ELINA_COEFF_SCALAR:
    elina_scalar_neg(a->val.scalar,b->val.scalar);
    break;
  case ELINA_COEFF_INTERVAL:
    elina_interval_neg(a->val.interval,b->val.interval);
    break;
  default:
    abort();
  }
}

/* Hash */

long elina_coeff_hash(elina_coeff_t* coeff)
{
  switch (coeff->discr){
  case ELINA_COEFF_SCALAR:
    return elina_scalar_hash(coeff->val.scalar);
  case ELINA_COEFF_INTERVAL:
    return elina_interval_hash(coeff->val.interval);
  default:
    abort();
    return 0;
  }
}

void elina_coeff_print(elina_coeff_t* a)
{ elina_coeff_fprint(stdout,a); }


void elina_coeff_swap(elina_coeff_t* a, elina_coeff_t* b)
{ elina_coeff_t t = *a; *a = *b; *b = t; }


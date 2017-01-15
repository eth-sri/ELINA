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
/* elina_interval.c: intervals */
/* ************************************************************************* */


#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include "elina_interval.h"

/* ====================================================================== */
/* Basics */
/* ====================================================================== */

elina_interval_t* elina_interval_alloc()
{
  elina_interval_t* itv = malloc(sizeof(elina_interval_t));
  itv->inf = elina_scalar_alloc();
  itv->sup = elina_scalar_alloc();
  return itv;
}
void elina_interval_reinit(elina_interval_t* itv,  elina_scalar_discr_t d)
{
  elina_scalar_reinit(itv->inf,d);
  elina_scalar_reinit(itv->sup,d);
}
void elina_interval_free(elina_interval_t* itv)
{
  elina_scalar_free(itv->inf);
  elina_scalar_free(itv->sup);
  free(itv);
}
void elina_interval_fprint(FILE* stream, elina_interval_t* a)
{
  fprintf(stream,"[");
  elina_scalar_fprint(stream,a->inf);
  fprintf(stream,",");
  elina_scalar_fprint(stream,a->sup);
  fprintf(stream,"]");
}

/* ====================================================================== */
/* Assignments */
/* ====================================================================== */

void elina_interval_set(elina_interval_t* interval, elina_interval_t* interval2)
{
  elina_scalar_set(interval->inf,interval2->inf);
  elina_scalar_set(interval->sup,interval2->sup); 
}
void elina_interval_set_scalar(elina_interval_t* interval, elina_scalar_t* inf, elina_scalar_t* sup)
{
  elina_scalar_set(interval->inf,inf);
  elina_scalar_set(interval->sup,sup);
}
void elina_interval_set_mpq(elina_interval_t* interval, mpq_t inf, mpq_t sup)
{
  elina_scalar_set_mpq(interval->inf,inf);
  elina_scalar_set_mpq(interval->sup,sup);
}
void elina_interval_set_int(elina_interval_t* interval, long int inf, long int sup)
{
  elina_scalar_set_int(interval->inf,inf);
  elina_scalar_set_int(interval->sup,sup);
}
void elina_interval_set_frac(elina_interval_t* interval, long int numinf, unsigned long int deninf, long int numsup, unsigned long int densup)
{
  elina_scalar_set_frac(interval->inf,numinf,deninf);
  elina_scalar_set_frac(interval->sup,numsup,densup);
}
void elina_interval_set_double(elina_interval_t* interval, double inf, double sup)
{
  elina_scalar_set_double(interval->inf,inf);
  elina_scalar_set_double(interval->sup,sup);
}
void elina_interval_set_mpfr(elina_interval_t* interval, mpfr_t inf, mpfr_t sup)
{
  elina_scalar_set_mpfr(interval->inf,inf);
  elina_scalar_set_mpfr(interval->sup,sup);
}
void elina_interval_set_top(elina_interval_t* interval)
{
  elina_scalar_set_infty(interval->inf,-1);
  elina_scalar_set_infty(interval->sup,+1);
}
void elina_interval_set_bottom(elina_interval_t* interval)
{
  switch (interval->inf->discr) {
  case ELINA_SCALAR_DOUBLE: interval->inf->val.dbl = 1.; break;
  case ELINA_SCALAR_MPQ:    mpq_set_si(interval->inf->val.mpq,1,1); break;
  case ELINA_SCALAR_MPFR:   mpfr_set_si(interval->inf->val.mpfr,1,GMP_RNDU); break;
  default:               abort();
  }
  switch (interval->sup->discr) {
  case ELINA_SCALAR_DOUBLE: interval->sup->val.dbl = -1.; break;
  case ELINA_SCALAR_MPQ:    mpq_set_si(interval->sup->val.mpq,-1,1); break;
  case ELINA_SCALAR_MPFR:   mpfr_set_si(interval->sup->val.mpfr,-1,GMP_RNDD); break;
  default:               abort();
  }
}

/* ====================================================================== */
/* Combined allocation and assignments */
/* ====================================================================== */

elina_interval_t* elina_interval_alloc_set(elina_interval_t* interval)
{
  elina_interval_t* itv = malloc(sizeof(elina_interval_t));
  itv->inf = elina_scalar_alloc_set(interval->inf);
  itv->sup = elina_scalar_alloc_set(interval->sup);
  return itv;
}

/* ====================================================================== */
/* Tests */
/* ====================================================================== */

bool elina_interval_is_top(elina_interval_t* interval)
{
  return elina_scalar_infty(interval->inf)<0 && elina_scalar_infty(interval->sup)>0;
}
bool elina_interval_is_bottom(elina_interval_t* interval)
{
  return elina_scalar_cmp(interval->inf,interval->sup)>0;
}
bool elina_interval_is_leq(elina_interval_t* itv1, elina_interval_t* itv2)
{
  int sinf = elina_scalar_cmp(itv1->inf,itv2->inf);
  int ssup = elina_scalar_cmp(itv1->sup,itv2->sup);
  return (sinf>=0 && ssup<=0) || elina_interval_is_bottom(itv1);
}
int elina_interval_cmp(elina_interval_t* itv1, elina_interval_t* itv2)
{
  int sinf = elina_scalar_cmp(itv1->inf,itv2->inf);
  int ssup = elina_scalar_cmp(itv1->sup,itv2->sup);

  if (sinf==0 && ssup==0) return 0;
  else if (sinf>=0 && ssup<=0) return -1;
  else if (sinf<=0 && ssup>=0) return 1;
  else {
    bool b1 = elina_interval_is_bottom(itv1);
    bool b2 = elina_interval_is_bottom(itv2);
    if (b1 && b2) return 0;
    else if (b1) return -1;
    else if (b2) return 1;
    else return sinf > 0 ? 2 : -2;
  }
}
bool elina_interval_equal(elina_interval_t* itv1, elina_interval_t* itv2)
{  
  bool inf = elina_scalar_equal(itv1->inf,itv2->inf);
  bool sup = elina_scalar_equal(itv1->sup,itv2->sup);
  if (inf && sup) return true;
  else return elina_interval_is_bottom(itv1) && elina_interval_is_bottom(itv2);
}
bool elina_interval_equal_int(elina_interval_t* itv, int b)
{  
  return elina_scalar_equal_int(itv->inf,b) && elina_scalar_equal_int(itv->sup,b);
}

/* ====================================================================== */
/* Other operations */
/* ====================================================================== */

void elina_interval_neg(elina_interval_t* a, elina_interval_t* b)
{
  if (a==b){
    elina_scalar_swap(a->inf,a->sup);
    elina_scalar_neg(a->inf,a->inf);
    elina_scalar_neg(a->sup,a->sup);
  } else {
    elina_scalar_neg(a->inf,b->sup);
    elina_scalar_neg(a->sup,b->inf);
  }
}

long elina_interval_hash(elina_interval_t* itv)
{
  if (elina_interval_is_bottom(itv)) return 0;
  else return 5*elina_scalar_hash(itv->inf) + 7*elina_scalar_hash(itv->sup);
}


/* ====================================================================== */
/* Array of intervals */
/* ====================================================================== */

elina_interval_t** elina_interval_array_alloc(size_t size)
{
  size_t i;

  elina_interval_t** array = malloc(size*sizeof(elina_interval_t*));
  for (i=0;i<size;i++){
    array[i] = elina_interval_alloc();
  }
  return array;
}

void elina_interval_array_free(elina_interval_t** array, size_t size)
{
  size_t i;
  for (i=0; i<size; i++){
    elina_interval_free(array[i]);
  }
  free(array);
}

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

#include "elina_interval_arith.h"

void elina_interval_add(elina_interval_t * op, elina_interval_t *op1, elina_interval_t *op2, elina_scalar_discr_t discr){
	elina_scalar_add(op->inf,op1->inf,op2->inf,discr);
	elina_scalar_add(op->sup,op1->sup,op2->sup,discr);
}


void elina_interval_mul_scalar(elina_interval_t * dst, elina_interval_t * src, elina_scalar_t * mul, elina_scalar_discr_t discr ){
	elina_scalar_t *dsup, *dinf;
 	elina_scalar_t *ssup, *sinf;
	dsup = dst->sup;
	dinf = dst->inf;
	ssup = src->sup;
	sinf = src->inf;
	elina_scalar_mul(dsup,ssup,mul,discr);
	elina_scalar_mul(dinf,sinf,mul,discr);
	if(elina_scalar_cmp(dsup,dinf)<0){
		elina_scalar_swap(dsup,dinf);
	}
}

void elina_interval_mul(elina_interval_t *op, elina_interval_t *op1, elina_interval_t *op2, elina_scalar_discr_t discr){
	elina_interval_t *tmp1, *tmp2;
	tmp1 = elina_interval_alloc();
	tmp2 = elina_interval_alloc();
	elina_interval_mul_scalar(tmp1,op1,op2->inf,discr);
	elina_interval_mul_scalar(tmp2,op1,op2->sup,discr);
	elina_scalar_min(op->inf,tmp1->inf,tmp2->inf);
	elina_scalar_max(op->sup,tmp1->sup,tmp2->sup);
	elina_interval_free(tmp1);
	elina_interval_free(tmp2);
}

/* ********************************************************************** */
/* Square root */
/* ********************************************************************** */


bool elina_interval_sqrt(elina_interval_t *dst, elina_interval_t *src, elina_scalar_discr_t discr)
{
  bool exact = true;
  if (elina_interval_is_bottom(src) || elina_scalar_sgn(src->sup)<0) {
    /* empty result */
    elina_interval_set_bottom(dst);
    return true;
  }
  /* lower bound */
  elina_scalar_t * tmp = elina_scalar_alloc();
  if (elina_scalar_sgn(src->inf)<=0) {
    elina_scalar_set_to_int(dst->inf,0,discr);
  }
  else {
    
    elina_scalar_sqrt(tmp,dst->inf,src->inf,discr);
    exact = exact && elina_scalar_equal(tmp,dst->inf);
  }
  /* upper bound */
  elina_scalar_sqrt(dst->sup,tmp,src->sup,discr);
  exact = exact && elina_scalar_equal(tmp,dst->sup);
  elina_scalar_free(tmp);
  return exact;
}


/* ********************************************************************** */
/* Subtraction */
/* ********************************************************************** */

void elina_interval_sub(elina_interval_t* a, elina_interval_t *b, elina_interval_t* c, elina_scalar_discr_t discr)
{
  elina_scalar_t * tmp1 = elina_scalar_alloc();
  elina_scalar_t * tmp2 = elina_scalar_alloc();
  elina_scalar_neg(tmp1,c->inf);
  elina_scalar_neg(tmp2,c->sup);
  elina_scalar_add(a->inf,b->inf,tmp2,discr);
  elina_scalar_add(a->sup,b->sup,tmp1,discr);
  elina_scalar_free(tmp1);
  elina_scalar_free(tmp2);
}


void elina_interval_abs(elina_interval_t *a, elina_interval_t *b, elina_scalar_discr_t discr)
{
  if (elina_scalar_sgn(b->inf)>=0){
    /* positive interval */
    elina_interval_set(a,b);
  }
  else if (elina_scalar_sgn(b->sup)<=0){
    /* negative interval */
    elina_interval_neg(a,b);
  }
  else {
    elina_scalar_max(a->sup, b->inf, b->sup);
    elina_scalar_set_to_int(a->inf,0,discr);
  }
}


/* ********************************************************************** */
/* Power */
/* ********************************************************************** */

void elina_interval_pow(elina_interval_t *a, elina_interval_t *b, elina_interval_t *n,elina_scalar_discr_t discr)
{
  long x;
  if (elina_interval_is_bottom(b) || elina_interval_is_bottom(n)) {
    elina_interval_set_bottom(a);
    return;
  }
  /* ensures that the exponent is a singleton */
 
  if (elina_scalar_infty(n->sup) || !elina_scalar_equal(n->inf, n->sup)) {
    elina_interval_set_top(a);
    return;
  }
  /* ensures that the exponent is a positive integer, stores it in x */
  int_set_elina_scalar(&x, n->sup);
  elina_scalar_t * scalar = elina_scalar_alloc();
  elina_scalar_set_to_int(scalar, x,discr);
  if (elina_scalar_cmp(scalar, n->sup) || x < 0) {
    elina_interval_set_top(a);
    return;
  }
  elina_interval_t * interval = elina_interval_alloc();
  if (x & 1) elina_interval_set(interval, b);
  else elina_interval_abs(interval, b, discr);
  elina_scalar_pow(a->sup, scalar, interval->sup, x,discr);
  elina_scalar_pow(scalar, a->inf, interval->inf, x,discr);
  elina_scalar_free(scalar);
  elina_interval_free(interval);
}


/* ********************************************************************** */
/* Truncate */
/* ********************************************************************** */
void elina_interval_trunc(elina_interval_t *a, elina_interval_t *b, elina_scalar_discr_t discr){
   elina_scalar_trunc(a->sup,b->sup, discr); 
   elina_scalar_trunc(a->inf,b->inf,discr); 
}

/* ********************************************************************** */
/* Ceil */
/* ********************************************************************** */
void elina_interval_ceil(elina_interval_t *a, elina_interval_t *b, elina_scalar_discr_t discr){
 elina_scalar_ceil(a->sup,b->sup, discr); 
 elina_scalar_ceil(a->inf,b->inf, discr); 
}

/* ********************************************************************** */
/* Floor */
/* ********************************************************************** */
void elina_interval_floor(elina_interval_t *a, elina_interval_t *b, elina_scalar_discr_t discr){
 elina_scalar_floor(a->sup,b->sup, discr); 
 elina_scalar_floor(a->inf,b->inf, discr); 
}


/* ********************************************************************** */
/* convert to integer */
/* ********************************************************************** */
void elina_interval_to_int(elina_interval_t *a, elina_interval_t *b, elina_scalar_discr_t discr){ 
	elina_scalar_ceil(a->sup,b->sup,discr); 
	elina_scalar_floor(a->inf,b->inf,discr); 
}


/* ********************************************************************** */
/* convert to single precision float */
/* ********************************************************************** */
void elina_interval_to_float(elina_interval_t *a, elina_interval_t *b, elina_scalar_discr_t discr){ 
	elina_scalar_to_float(a->sup,b->sup,discr); 
	elina_scalar_to_float(a->inf,b->inf,discr); 
}

/* ********************************************************************** */
/* convert to double precision float */
/* ********************************************************************** */
void elina_interval_to_double(elina_interval_t *a, elina_interval_t *b, elina_scalar_discr_t discr){
	
 	elina_scalar_to_double(a->sup,b->sup,discr); 
 	elina_scalar_to_double(a->inf,b->inf,discr); 
	
}


/* ====================================================================== */
/* Division */
/* ====================================================================== */

/* Assume that both intervals are positive */
static void elina_interval_divpp(elina_interval_t *a, elina_interval_t *b, elina_interval_t *c, elina_scalar_discr_t discr){
  elina_scalar_t * tmp = elina_scalar_alloc();
  elina_scalar_set(tmp,c->inf);
  elina_scalar_div(a->inf,b->inf,c->sup, discr);
  elina_scalar_div(a->sup,b->sup,tmp, discr);
  elina_scalar_free(tmp);
}

/* Assume that both intervals are negative */
static void elina_interval_divnn(elina_interval_t *a, elina_interval_t *b, elina_interval_t *c, elina_scalar_discr_t discr){
  elina_scalar_t *scalar = elina_scalar_alloc();
  elina_scalar_set(scalar,b->inf);
  elina_scalar_div(a->inf,b->sup,c->inf, discr);
  elina_scalar_div(a->sup,scalar,c->sup, discr);
  elina_scalar_free(scalar);
}

/* Assume that b is positive and c negative */
static void elina_interval_divpn(elina_interval_t *a, elina_interval_t *b, elina_interval_t *c, elina_scalar_discr_t discr){
  elina_scalar_t *scalar = elina_scalar_alloc();
  elina_scalar_set(scalar,b->sup);
  elina_scalar_div(scalar,scalar,c->sup, discr);
  elina_scalar_div(a->sup,b->inf,c->inf, discr);
  elina_scalar_set(a->inf,scalar);
  elina_scalar_free(scalar);
}

/* Assume that b is negative and c positive */
static void elina_interval_divnp(elina_interval_t *a, elina_interval_t *b, elina_interval_t *c, elina_scalar_discr_t discr){
  elina_scalar_div(a->inf,b->inf,c->inf, discr);
  elina_scalar_div(a->sup,b->sup,c->sup, discr);
}

/* Assume that interval c is positive */
static
void elina_interval_divp(elina_interval_t *a, elina_interval_t *b, elina_interval_t *c, elina_scalar_discr_t discr){
  if (elina_scalar_sgn(b->inf)>=0){
    /* b is positive */
    elina_interval_divpp(a,b,c,discr);
  }
  else if (elina_scalar_sgn(b->sup)<=0){
    /* b is negative */
    elina_interval_divnp(a,b,c,discr);
  }
  else {
    /* 0 is in the middle of b: one divides b by c->inf */
    elina_scalar_div(a->inf,b->inf,c->inf, discr);
    elina_scalar_div(a->sup,b->sup,c->inf, discr);
  }
}
/* Assume that interval c is negative */
static void elina_interval_divn(elina_interval_t *a, elina_interval_t *b, elina_interval_t *c, elina_scalar_discr_t discr){
  if (elina_scalar_sgn(b->inf)>=0){
    /* b is positive */
    elina_interval_divpn(a,b,c,discr);
  }
  else if (elina_scalar_sgn(b->sup)<=0){
    /* b is negative */
    elina_interval_divnn(a,b,c,discr);
  }
  else {
    /* 0 is in the middle of b: one cross-divide b by c->sup */
    if (a!=b) {
      elina_scalar_div(a->inf,b->sup,c->sup, discr);
      elina_scalar_div(a->sup,b->inf,c->sup, discr);
    }
    else {
	elina_scalar_t * tmp = elina_scalar_alloc();
        elina_scalar_div(tmp,b->sup,c->sup, discr);
        elina_scalar_div(a->sup,b->inf,c->sup, discr);
        elina_scalar_set(a->inf,tmp);
	elina_scalar_free(tmp);
    }
  }
}

void elina_interval_div(elina_interval_t *a, elina_interval_t *b, elina_interval_t *c, elina_scalar_discr_t discr)
{
  if (elina_scalar_sgn(c->inf)>0){
    /* c is positive */
    elina_interval_divp(a,b,c,discr);
  }
  else if (elina_scalar_sgn(c->sup)<0){
    /* c is negative */
    elina_interval_divn(a,b,c,discr);
  }
  else if (elina_scalar_sgn(b->inf)==0 && elina_scalar_sgn(b->sup)==0){
    /* b is [0,0] */
    elina_interval_set(a,b);
  }
  else {
    elina_interval_set_top(a);
  }
}


/* ********************************************************************** */
/* Modulo */
/* ********************************************************************** */
void elina_interval_mod(elina_interval_t *a, elina_interval_t *b, elina_interval_t *c, bool is_int, elina_scalar_discr_t discr)
{
  elina_interval_t * interval = elina_interval_alloc();
  elina_interval_t * interval2 = elina_interval_alloc();
  /* b-|c|*trunc(b/|c|) */ 
  elina_interval_abs(interval, c, discr);
  if (!elina_scalar_sgn(interval->inf)) elina_interval_set_top(a);
  else {
    elina_interval_div(interval2, b, interval,discr);
    elina_interval_trunc(interval2, interval2,discr);
    elina_interval_mul(interval2, interval2, interval,discr);
    if (is_int) elina_scalar_sub_uint(interval->sup,interval->sup,1,discr);
    if (elina_scalar_sgn(b->sup)<0) {
      /* [-max|c|,0] */
      elina_scalar_neg(interval->inf, interval->sup);
      elina_scalar_set_to_int(interval->sup, 0,discr);
    }
    else if (elina_scalar_sgn(b->inf)>0){
      /* [-max|c|,max|c|] */
      elina_scalar_neg(interval->inf, interval->sup);
    }
    else{
      /* [0,max|c|] */
      elina_scalar_set_to_int(interval->inf, 0,discr);
    }
    elina_interval_sub(a, b, interval2,discr);
    elina_scalar_max(a->inf, a->inf, interval->inf);
    elina_scalar_min(a->sup, a->sup, interval->sup);
  }
  elina_interval_free(interval);
  elina_interval_free(interval2);
}

/* ********************************************************************** */
/* Is Integer */
/* ********************************************************************** */
bool elina_interval_is_int(elina_interval_t *a, elina_scalar_discr_t discr)
{
  elina_scalar_t *tmp = elina_scalar_alloc();
  elina_scalar_trunc(tmp,a->sup,discr);
  if (elina_scalar_cmp(tmp,a->sup)){
      elina_scalar_free(tmp);
      return false;
  }
  elina_scalar_trunc(tmp,a->inf,discr);
  bool res = !elina_scalar_cmp(tmp,a->inf);
  elina_scalar_free(tmp);
  return res;
}


void elina_interval_magnitude(elina_scalar_t *a, elina_interval_t *b)
{
  if (elina_scalar_sgn(b->inf)>=0) elina_scalar_set(a,b->sup);
  else if (elina_scalar_sgn(b->sup)<=0) elina_scalar_neg(a,b->inf);
  else
    elina_scalar_max(a, b->inf, b->sup);
}

void elina_interval_range_rel(elina_scalar_t *a, elina_interval_t *b, elina_scalar_discr_t discr)
{
  elina_scalar_reinit(a,discr);
  elina_scalar_add(a,b->sup,b->inf,discr);
  if (!elina_scalar_infty(a)) {
     elina_scalar_t * tmp = elina_scalar_alloc();
     elina_interval_magnitude(tmp,b);
     elina_scalar_div_2(tmp,tmp);
     elina_scalar_div(a,a,tmp,discr);
     elina_scalar_free(tmp);
  }
}

void elina_interval_enlarge_bound(elina_interval_t *a, elina_interval_t *b, elina_scalar_t *c, elina_scalar_discr_t discr)
{
  elina_scalar_t * tmp = elina_scalar_alloc();
  //elina_scalar_neg(tmp,c);
  elina_scalar_neg(tmp,c);
  //elina_scalar_add(a->inf,b->inf,tmp,discr);
  elina_scalar_add(a->inf,b->inf,tmp,discr);
  //elina_scalar_neg(tmp,tmp);
  elina_scalar_add(a->sup,b->sup,c,discr);
  elina_scalar_free(tmp);
}

bool elina_interval_canonicalize(elina_interval_t *a, bool integer, elina_scalar_discr_t discr)
{
  bool exc;

  if (integer){
    elina_scalar_floor(a->inf,a->inf,discr);
    elina_scalar_floor(a->sup,a->sup,discr);
  }
  if (elina_scalar_infty(a->inf) || elina_scalar_infty(a->sup)) return false;

  /* Check that it is not bottom */
  exc = false;
  if (elina_scalar_cmp(a->sup,a->inf) < 0)
    exc = true;
  return exc;
}

void elina_interval_set_to_int(elina_interval_t *a, elina_int_t i, elina_int_t j, elina_scalar_discr_t discr){
	elina_scalar_set_to_int(a->inf,i,discr);
	elina_scalar_set_to_int(a->sup,j,discr);
}

void elina_interval_mul_2exp(elina_interval_t *a, elina_interval_t *b, int i, elina_scalar_discr_t discr){
	elina_scalar_mul_2exp(a->inf,b->inf,i,discr);
	elina_scalar_mul_2exp(a->sup,b->sup,i,discr);
}

void elina_interval_convert(elina_interval_t * a, elina_scalar_discr_t discr){
	elina_scalar_convert(a->inf,discr);
	elina_scalar_convert(a->sup,discr);
}

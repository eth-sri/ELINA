/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
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


/* ********************************************************************** */
/* Rational numbers for ELINA */
/* ********************************************************************** */

#ifndef _ELINA_RAT_H_
#define _ELINA_RAT_H_

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include "gmp.h"
#include "mpfr.h"

#if defined (HAS_APRON)
#include "apron_wrapper.h"
#else
#include "elina_scalar.h"
#endif

#include "elina_int.h"


#ifdef __cplusplus
extern "C" {
#endif

typedef struct elina_rat_t {
  elina_int_t n; /* numerator */
  elina_int_t d; /* denominator, >=0 */
} elina_rat_t;

static inline int elina_rat_sgn(elina_rat_t *r){
	return elina_int_sgn(r->n);
}

static inline bool elina_rat_infty(elina_rat_t *r){
	return r->d==0;
}

/*  int -> elina_rat */
static inline bool elina_rat_set_int(elina_rat_t *a, long long int num){
	a->n = num;
	a->d = 1;
	return true;
}

/* mpq -> elina_rat */
static inline bool elina_rat_set_mpq(elina_rat_t* a, mpq_t b)
{ 
  elina_int_set_mpz(&a->n,mpq_numref(b));
  elina_int_set_mpz(&a->d,mpq_denref(b));
  return true;
}

/* double -> elina_rat */
static inline bool elina_rat_set_double_tmp(elina_rat_t* a, double k, mpq_t mpq)
{
  if (!isfinite(k)) {elina_rat_set_int(a,0); return false; }
  mpq_set_d(mpq,k);
  elina_rat_set_mpq(a,mpq);
  return true;
}
static inline bool elina_rat_set_double(elina_rat_t* a, double k)
{
  mpq_t mpq;
  mpq_init(mpq);
  elina_rat_set_double_tmp(a,k,mpq);
  mpq_clear(mpq);
  return true;
}

/* mpfr -> elina_rat */
static inline bool elina_rat_set_mpfr(elina_rat_t* a, mpfr_t b)
{
  mpq_t q;
  mp_exp_t e;
  if (!mpfr_number_p(b)) { elina_rat_set_int(a,0); return false; }
  mpq_init(q);
  /* XXX might fail if scaled exponent not representable in mp_exp_t */
  e = mpfr_get_z_exp(mpq_numref(q),b);
  mpz_set_si(mpq_denref(q),1);
  if (e>0) mpq_mul_2exp(q,q,e);
  if (e<0) mpq_div_2exp(q,q,-e);
  elina_rat_set_mpq(a,q);
  mpq_clear(q);
  return true;
}


/* ====================================================================== */
/* Rational operations */
/* ====================================================================== */

static inline void elina_rat_canonicalize(elina_rat_t * r)
{
  if (r->d){
    elina_int_t pgcd;
    pgcd = elina_int_gcd(r->n,r->d);
    if (pgcd==0 || (pgcd==-1 && (r->d==ELINA_INT_MIN || r->n==ELINA_INT_MIN))) {
      fprintf(stderr,"overflow in elina_rat_canonicalize\n");
      return; 
    }
    r->n /= pgcd;
    r->d /= pgcd;
  }
  else {
    r->n = 1;
  }
}


static inline bool elina_rat_set_int2(elina_rat_t *a, elina_int_t i, elina_int_t j)
{ 
  a->n = i;
  a->d = j;
  elina_rat_canonicalize(a);
  return true;
}

static inline int elina_rat_cmp(elina_rat_t *a, elina_rat_t *b)
{ 
  elina_int_t aa,bb;
  if(a->d==0){
	int sgn_a = elina_int_sgn(a->n);
	if(b->d==0){
		int sgn_b = elina_int_sgn(b->n);
		return elina_int_cmp(sgn_a,sgn_b);
	}
	else{
		if(sgn_a<0){
			return -1;
		}
		else{
			return 1;
		}
	}
  }
  else{
	if(b->d==0){
		int sgn_b = elina_int_sgn(b->n);
		if(sgn_b<0){
			return 1;
		}
		else{
			return -1;
		}
	}
	else{
		elina_int_t d = elina_int_lcm(a->d,b->d);
  		aa = a->n * (d / a->d);
  		bb = b->n * (d / b->d);  
  		return elina_int_cmp(aa,bb);
	}
  }
  
}

static inline bool elina_rat_set_elina_scalar(elina_rat_t *a, elina_scalar_t* b)
{
  assert (elina_scalar_infty(b)==0);
  switch (b->discr){
  case ELINA_SCALAR_MPQ:
    return elina_rat_set_mpq(a,b->val.mpq);
  case ELINA_SCALAR_DOUBLE:
    return elina_rat_set_double(a,b->val.dbl);
  case ELINA_SCALAR_MPFR:
    return elina_rat_set_mpfr(a,b->val.mpfr);
  default: abort(); 
  }
}

static inline bool mpq_set_elina_rat(mpq_t a, elina_rat_t *b)
{
  mpz_set_elina_int(mpq_numref(a), b->n);
  mpz_set_elina_int(mpq_denref(a), b->d);
  return true;
}

static inline bool double_set_elina_rat_tmp(double* a, elina_rat_t* b, 
					 mpq_t mpq, mpfr_t mpfr)
{
  mpq_set_elina_rat(mpq,b);
  int res = mpfr_set_q(mpfr,mpq,GMP_RNDU);
  *a = mpfr_get_d(mpfr,GMP_RNDU); /* should be exact */
  return (res==0);
}

static inline bool double_set_elina_rat(double* a, elina_rat_t * b)
{
  mpq_t mpq;
  mpfr_t mpfr;
  mpq_init(mpq);
  mpfr_init2(mpfr,53);
  bool res = double_set_elina_rat_tmp(a,b,mpq,mpfr);
  mpq_clear(mpq);
  mpfr_clear(mpfr);
  return res;
}
static inline bool mpfr_set_elina_rat(mpfr_t a, elina_rat_t *b)
{ 
  int r = mpfr_set_si(a,b->n,GMP_RNDU);
  return !mpfr_div_si(a,a,b->d,GMP_RNDU) && !r;
}

static inline bool elina_scalar_set_elina_rat(elina_scalar_t* a, elina_rat_t *b)
{
  elina_scalar_reinit(a,ELINA_SCALAR_MPQ);
  return mpq_set_elina_rat(a->val.mpq,b);
}

static inline int elina_scalar_cmp_elina_rat(elina_scalar_t * a, elina_rat_t * b){
	elina_scalar_t * tmp = elina_scalar_alloc();
	bool res;
	switch(a->discr){
		case ELINA_SCALAR_MPQ:
			elina_scalar_init(tmp,ELINA_SCALAR_MPQ);
			mpq_set_elina_rat(tmp->val.mpq,b);
			break;
		case ELINA_SCALAR_DOUBLE:
			double_set_elina_rat(&tmp->val.dbl,b);
			break;
		case ELINA_SCALAR_MPFR:
			elina_scalar_init(tmp,ELINA_SCALAR_MPFR);
			mpfr_set_elina_rat(tmp->val.mpfr,b);
			break;
	}
	res = elina_scalar_cmp(a,tmp);
	elina_scalar_free(tmp);
	
	return res;
}

static inline void elina_rat_div(elina_rat_t* a, elina_rat_t* b, elina_rat_t* c)
{
  elina_int_t d;
  d = b->d * c->n;
  if (d<0) {
    a->n = - b->n * c->d;
    a->d = - d;
  }
  else {
    a->n = b->n * c->d;
    a->d = d;
  }
  elina_rat_canonicalize(a);   
}

static inline void elina_rat_inv(elina_rat_t *a, elina_rat_t *b)
{ 
  if (a!=b){
    a->n = b->d;
    a->d = b->n;
  }
  else{
	elina_int_t tmp = a->n;
	a->n = a->d;
	a->d = tmp;
  }
  if (a->d<0){
    a->n = -a->n;
    a->d = -a->d;
  }
}

static inline void elina_rat_sub_uint(elina_rat_t *a, elina_rat_t *b, unsigned long int c)
{ 
  a->n = b->n - (elina_int_t)c * b->d;
  a->d = b->d;
  elina_rat_canonicalize(a);
}


static inline void elina_rat_set_infty(elina_rat_t * r, char sgn){
	if(sgn<0){
		r->n = -1;
		r->d = 0;
	}
	else{
		r->n = 1;
		r->d = 0;
	}
}

static elina_rat_t * elina_scalar_set_rat(elina_scalar_t * scalar){
	elina_rat_t * res = (elina_rat_t *)calloc(1,sizeof(elina_rat_t));
	switch (scalar->discr){
 	 	case ELINA_SCALAR_MPQ:
		    if (!mpz_sgn(mpq_denref(scalar->val.mpq))){
			int sgn = mpz_sgn(mpq_numref(scalar->val.mpq));
			if(sgn>0){
				res->n = 1;
			}
			else if(sgn < 0){
				res->n = -1;
			}
			else{
				abort();
				return NULL;
			}
		      	res->d = 0;
		      return res;
		    }
		    else {
		      elina_rat_set_mpq(res,scalar->val.mpq);
		      return res;
		    }
		    break;
		case ELINA_SCALAR_DOUBLE:
		    if (scalar->val.dbl==INFINITY || scalar->val.dbl==-INFINITY) {
		      if (scalar->val.dbl>0){ 
			  res->n = 1;
		      }
		      else {
			  res->n = -1;
		      }
		      res->d = 0;
		      return res;
		    }
		    else {
		        elina_rat_set_double(res,scalar->val.dbl);
			return res;
		    }
		    break;
		case ELINA_SCALAR_MPFR:
		    if (mpfr_inf_p(scalar->val.mpfr)) {
		      if (mpfr_sgn(scalar->val.mpfr)>0) {
			  res->n=1;
		      }
		      else {
			  res->n = -1;
		      }
			res->d = 0;
		      return res;
		    }
		    else {
		        elina_rat_set_mpfr(res,scalar->val.mpfr);
			return res;
		    }
		    break;
		default:
		    abort();
		    return NULL;
	}
}

static inline void elina_rat_neg(elina_rat_t *a, elina_rat_t* b){
	a->n = -b->n; 
	a->d = b->d; 
}

static inline void elina_rat_mul(elina_rat_t* a, elina_rat_t* b, elina_rat_t* c)
{
  a->n = b->n * c->n;
  a->d = b->d * c->d;
  elina_rat_canonicalize(a);   
}

static inline void elina_rat_add(elina_rat_t* a, elina_rat_t* b, elina_rat_t* c)
{ 
  elina_int_t d = elina_int_lcm(b->d,c->d);
  a->n = b->n * (d / b->d) + (d / c->d) * c->n; 
  a->d = d; 
  elina_rat_canonicalize(a); 
}


static inline void elina_rat_min(elina_rat_t *a, elina_rat_t *b, elina_rat_t *c){ 
       if(elina_rat_cmp(b,c)<=0){
	  a->n = b->n;
	  a->d = b->d;
       }
       else{
	  a->n = c->n;
	  a->d = c->d;
       }
	
}
static inline void elina_rat_max(elina_rat_t *a, elina_rat_t *b, elina_rat_t *c){ 
	if(elina_rat_cmp(b,c)>=0){
	   a->n = b->n;
	   a->d = b->d;
	}
	else{
	   a->n = c->n;
	   a->d = c->d;
	}
}

#ifdef __cplusplus
}
#endif

#endif

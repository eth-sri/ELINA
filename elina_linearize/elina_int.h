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
/* Integers for ELINA */
/* ********************************************************************** */

#ifndef _ELINA_INT_H_
#define _ELINA_INT_H_

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "gmp.h"
#include "mpfr.h"

#if defined (HAS_APRON)
#include "apron_wrapper.h"
#else
#include "elina_scalar.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef long long int elina_int_t;
#define ELINA_INT_ZERO 0LL
#define ELINA_INT_ONE 1LL
#define ELINA_INT_MAX LLONG_MAX
#define ELINA_INT_MIN LLONG_MIN


static inline int elina_int_sgn(elina_int_t a)
{ return (a==ELINA_INT_ZERO ? 0 : (a>ELINA_INT_ZERO ? 1 : -1)); }

static inline elina_int_t elina_int_fdiv_q(elina_int_t a, elina_int_t b)
{
  elina_int_t q;
  if (elina_int_sgn(a)*elina_int_sgn(b)<0 && a % b) q = a / b - 1;
  else q = a / b;
  return q;
}

static inline elina_int_t elina_int_cdiv_q(elina_int_t a, elina_int_t b)
{
  elina_int_t q;
  if (elina_int_sgn(a)*elina_int_sgn(b)>0 && a % b) q = a / b + 1;
  else q = a / b;
  return q;
}


/* mpz -> elina_int */
static inline bool elina_int_set_mpz(elina_int_t *a, mpz_t b)
{
   *a = mpz_get_si(b);
  return true;
}

static inline elina_int_t elina_int_abs(elina_int_t a){
	if(a<0){
		return -a;
	}
	return a;
}

static inline elina_int_t elina_gcd_aux2(elina_int_t a, elina_int_t b)
{ /* a is supposed to be greater than b */
  elina_int_t t;
  while (b!=ELINA_INT_ZERO && a!=b) {
    t = b;
    b = a % b;
    a = t;
  }
  return a;
}
static inline elina_int_t elina_gcd_aux(elina_int_t a, elina_int_t b)
{
  a = elina_int_abs(a);
  b = elina_int_abs(b);
  return (a>=b) ? elina_gcd_aux2(a,b) : elina_gcd_aux2(b,a);
}
static inline elina_int_t elina_int_gcd(elina_int_t b,  elina_int_t c)
{ return elina_gcd_aux(b,c); }


static inline int elina_int_cmp(elina_int_t a, elina_int_t b)
{ return (a==b ? 0 : (a>b ? 1 : -1)); }

static inline elina_int_t elina_lcm_aux(elina_int_t a, elina_int_t b)
{
  a = elina_int_abs(a);
  b = elina_int_abs(b);
  return a / elina_gcd_aux(a,b) * b;
}

static inline elina_int_t elina_int_lcm(elina_int_t b, elina_int_t c)
{ return elina_lcm_aux(b,c); }

static inline bool mpz_set_elina_int(mpz_t a, elina_int_t b)
{
  mpz_set_si(a,b);
  return true;
}

static inline bool mpq_set_elina_int(mpq_t a, elina_int_t b)
{
  if (sizeof(elina_int_t)==sizeof(long int)) {
    mpq_set_si(a,b,1);
    return true;
  }
  else {
    mpz_set_ui(mpq_denref(a),1);
    return mpz_set_elina_int(mpq_numref(a),b);
  }
}


static inline bool int_set_mpq_tmp(long int* a, mpq_t b, 
				      mpz_t q, mpz_t r)
{ 
  mpz_cdiv_qr(q,r,mpq_numref(b),mpq_denref(b));
  *a = mpz_get_si(q);
  return (mpz_sgn(r)==0);
}

static inline bool int_set_mpq(long int* a, mpq_t b)
{ 
  mpz_t q,r;
  mpz_init(q); mpz_init(r);
  bool res = int_set_mpq_tmp(a,b,q,r);
  mpz_clear(q); mpz_clear(r);
  return res;
}


static inline bool int_set_double(long int* a, double b)
{
  double c;
  c = ceil(b);
  if (!isfinite(c)) { *a = 0; return false; }
  *a = (long int)c;
  return (b==c);
}

static inline bool int_set_mpfr(long int* a, mpfr_t b)
{
  if (!mpfr_number_p(b)) { *a = 0; return false; }
  *a = mpfr_get_si(b,GMP_RNDU);
  return mpfr_integer_p(b);
}

#ifdef __cplusplus
}
#endif

#endif

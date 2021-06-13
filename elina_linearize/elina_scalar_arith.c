/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
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

#include "elina_scalar_arith.h"

/* ************************************************************************* */
/*  convert scalar */
/* ************************************************************************* */

void elina_scalar_convert(elina_scalar_t *a, elina_scalar_discr_t discr){
	if(a->discr==discr){
		return;
	}
	mpq_t mpq;
        mpfr_t mpfr;
	double d;
	switch(discr){
		case ELINA_SCALAR_MPQ:
			if(a->discr==ELINA_SCALAR_DOUBLE){
				d = a->val.dbl;
				elina_scalar_reinit(a,discr);
				if (d==INFINITY) { 
					 mpz_set_si(mpq_numref(a->val.mpq),1);
    					 mpz_set_ui(mpq_denref(a->val.mpq), 0);
				}
				else if(d==-INFINITY){
					 mpz_set_si(mpq_numref(a->val.mpq),-1);
    					 mpz_set_ui(mpq_denref(a->val.mpq), 0);
				}
				else{
					mpq_set_d(a->val.mpq,d);
				}		
			}
			else{
				mpfr_init(mpfr);
				mpfr_set(mpfr,a->val.mpfr,GMP_RNDU);
				elina_scalar_reinit(a,discr);
				mp_exp_t e;
  				if (!mpfr_number_p(mpfr)) { 
					mpq_set_si(a->val.mpq,1,0);
				}
				else{
  					e = mpfr_get_z_exp(mpq_numref(a->val.mpq),mpfr);
  					mpz_set_si(mpq_denref(a->val.mpq),1);
  					if (e>0) mpq_mul_2exp(a->val.mpq,a->val.mpq,e);
  					if (e<0) mpq_div_2exp(a->val.mpq,a->val.mpq,-e);
				}
				mpfr_clear(mpfr);
			}
			break;	
		case ELINA_SCALAR_DOUBLE:
			if(a->discr==ELINA_SCALAR_MPQ){
				mpq_init(mpq);
				mpq_set(mpq,a->val.mpq);
				elina_scalar_reinit(a,discr);
				mpfr_t mpfr;
  				mpfr_init2(mpfr,DBL_MANT_DIG);
				mpfr_set_q(mpfr,mpq,GMP_RNDU);
				a->val.dbl = mpfr_get_d(mpfr,GMP_RNDU);
  				mpfr_clear(mpfr);
				mpq_clear(mpq);
			}
			else{
				mpfr_init(mpfr);
				mpfr_set(mpfr,a->val.mpfr,GMP_RNDU);
				elina_scalar_reinit(a,discr);
				a->val.dbl = mpfr_get_d(mpfr,GMP_RNDU);
				mpfr_clear(mpfr);
			}
			break;
		case ELINA_SCALAR_MPFR:
			if(a->discr==ELINA_SCALAR_MPQ){
				mpq_init(mpq);
				mpq_set(mpq,a->val.mpq);
				elina_scalar_reinit(a,discr);
				mpfr_set_q(a->val.mpfr,mpq,GMP_RNDU);
				mpq_clear(mpq);
			}
			else{
				d = a->val.dbl;
				elina_scalar_reinit(a,discr);
				mpfr_set_d(a->val.mpfr,d,GMP_RNDU); 
			}
			break;
	}
	
}

static inline void elina_scalar_add_uint_mpq(mpq_t a, mpq_t b, unsigned long int c)
{
  if (a==b) {
    mpq_t tmp;
    mpq_init(tmp);
    mpq_set_ui(tmp,c,1);
    mpq_add(a,b,tmp);
    mpq_clear(tmp);
  }
  else {
    mpq_set_ui(a,c,1);
    mpq_add(a,a,b);
  }
}

static inline void elina_scalar_add_uint_double(double* a, double* b, unsigned long int c){
	 *a = *b + (elina_int_t)c; 
}

static inline void elina_scalar_add_uint_mpfr(mpfr_t a, mpfr_t b, unsigned long int c){
  mpfr_add_ui(a,b,c,GMP_RNDU);
}

void elina_scalar_add_uint(elina_scalar_t *a, elina_scalar_t * b, unsigned long int c, elina_scalar_discr_t discr)
{
  if(a!=b){
  	elina_scalar_reinit(a,discr);
  }
  if (elina_scalar_infty(b)){
	 elina_scalar_set_infty(a,elina_scalar_sgn(b));
  }
  else { 
	switch(discr){
		case ELINA_SCALAR_MPQ:
			elina_scalar_add_uint_mpq(a->val.mpq,b->val.mpq,c);
			break;
		case ELINA_SCALAR_DOUBLE:
			elina_scalar_add_uint_double(&a->val.dbl,&b->val.dbl,c);
			break;
		case ELINA_SCALAR_MPFR:
			elina_scalar_add_uint_mpfr(a->val.mpfr,b->val.mpfr,c);
			break;
	}
  }
}


void elina_scalar_add(elina_scalar_t *op, elina_scalar_t *op1, elina_scalar_t *op2, elina_scalar_discr_t discr){
	elina_scalar_convert(op1,discr);
	elina_scalar_convert(op2,discr);
	if(op!=op1){
	   	elina_scalar_reinit(op,discr);	
	}
	
	int res1 = elina_scalar_infty(op1);
	int res2 = elina_scalar_infty(op2);
	if(res1){
		elina_scalar_set_infty(op,res1);
		return;
	}
	if(res2){
		elina_scalar_set_infty(op,res2);
		return;
	}
	//elina_scalar_reinit(op1,op2->discr);
	switch(discr){
		case ELINA_SCALAR_MPQ:
			mpq_add(op->val.mpq,op1->val.mpq,op2->val.mpq);
			break;
		case ELINA_SCALAR_DOUBLE:
			op->val.dbl = op1->val.dbl+op2->val.dbl;
			break;
		case ELINA_SCALAR_MPFR:
			mpfr_add(op->val.mpfr,op1->val.mpfr,op2->val.mpfr,GMP_RNDU);
			break;
	}
}



/* ************************************************************************* */
/* Multiplication */
/* ************************************************************************* */
void elina_scalar_mul(elina_scalar_t *op, elina_scalar_t *op1, elina_scalar_t *op2, elina_scalar_discr_t discr){
	elina_scalar_convert(op1,discr);
	elina_scalar_convert(op2,discr);
	if(op!=op1){
		elina_scalar_reinit(op,discr);
	}
	int sgn1 = elina_scalar_sgn(op1);
	int sgn2 = elina_scalar_sgn(op2);
	if(!sgn1){
		elina_scalar_set(op,op1);
		return;
	}
	if(!sgn2){
		elina_scalar_set(op,op2);
		return;
	}
	int res1 = elina_scalar_infty(op1);
	int res2 = elina_scalar_infty(op2);
	if(res1 || res2){
		int res = sgn1*sgn2;
		elina_scalar_set_infty(op,res);
		return;
	}
	
	//elina_scalar_reinit(op1,op2->discr);
		switch(discr){
			case ELINA_SCALAR_MPQ:
				mpq_mul(op->val.mpq,op1->val.mpq,op2->val.mpq);
				break;
			case ELINA_SCALAR_DOUBLE:
				op->val.dbl = op1->val.dbl*op2->val.dbl;
				break;
			case ELINA_SCALAR_MPFR:
				mpfr_mul(op->val.mpfr,op1->val.mpfr,op2->val.mpfr,GMP_RNDU);
				break;
		}
}

/* ********************************************************************** */
/* Subtraction */
/* ********************************************************************** */

static inline void elina_scalar_sub_uint_mpq(mpq_t a, mpq_t b, unsigned long int c)
{
  if (a==b) {
    mpq_t tmp;
    mpq_init(tmp);
    mpq_set_ui(tmp,c,1);
    mpq_sub(a,b,tmp);
    mpq_clear(tmp);
  }
  else {
    mpq_set_ui(a,c,1);
    mpq_sub(a,b,a);
  }
}

static inline void elina_scalar_sub_uint_double(double *a, double *b, unsigned long int c){
	 *a = *b - (elina_int_t)c; 
}

static inline void elina_scalar_sub_uint_mpfr(mpfr_t a, mpfr_t b, unsigned long int c){ 
	mpfr_sub_ui(a,b,c,GMP_RNDU); 
}

void elina_scalar_sub_uint(elina_scalar_t *a, elina_scalar_t * b, unsigned long int c, elina_scalar_discr_t discr)
{
  if(a!=b){
  	elina_scalar_reinit(a,discr);
  }
  if (elina_scalar_infty(b)){
	 elina_scalar_set_infty(a,elina_scalar_sgn(b));
  }
  else { 
	switch(discr){
		case ELINA_SCALAR_MPQ:
			elina_scalar_sub_uint_mpq(a->val.mpq,b->val.mpq,c);
			break;
		case ELINA_SCALAR_DOUBLE:
			elina_scalar_sub_uint_double(&a->val.dbl,&b->val.dbl,c);
			break;
		case ELINA_SCALAR_MPFR:
			elina_scalar_sub_uint_mpfr(a->val.mpfr,b->val.mpfr,c);
			break;
	}
  }
}

void elina_scalar_sub(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_t *c, elina_scalar_discr_t discr){
        if(a!=b){
	   	elina_scalar_reinit(a,c->discr);	
	}
	elina_scalar_neg(c,c);
	elina_scalar_add(a,b,c,discr);
	elina_scalar_neg(c,c);
}

/* ********************************************************************** */
/* Division */
/* ********************************************************************** */


void elina_scalar_div_2(elina_scalar_t * dst, elina_scalar_t *src){
	if(dst!=src){
		elina_scalar_reinit(dst,src->discr);
	}
	switch(src->discr){
		case ELINA_SCALAR_MPQ:
			mpq_div_2exp(dst->val.mpq,src->val.mpq,1);
			break;
		case ELINA_SCALAR_DOUBLE:
			dst->val.dbl = src->val.dbl/2;
			break;
		case ELINA_SCALAR_MPFR:
			mpfr_div_2ui(dst->val.mpfr,src->val.mpfr,1,GMP_RNDU);
			break;
	}
}


void elina_scalar_div(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_t *c, elina_scalar_discr_t discr)
{
  if(a!=b && a != c){
  	elina_scalar_reinit(a,discr);
  }
  elina_scalar_convert(b,discr);
  elina_scalar_convert(c,discr);
  int sgnb = elina_scalar_sgn(b);
  int sgnc = elina_scalar_sgn(c);
  if (!sgnb || elina_scalar_infty(c)) {
	elina_scalar_set_to_int(a,0,discr);
	return;
  }
  else if (!sgnc){
	 elina_scalar_set_infty(a,sgnb);
	return;
  }
  else if(elina_scalar_infty(b)){
	int res = sgnb *sgnc;
	elina_scalar_set_infty(a,res);
	return;
  }
  else{
   	switch(discr){
		case ELINA_SCALAR_MPQ:
			mpq_div(a->val.mpq,b->val.mpq,c->val.mpq);
			break;
		case ELINA_SCALAR_DOUBLE:
			a->val.dbl =  b->val.dbl/c->val.dbl;
			break;
		case ELINA_SCALAR_MPFR:
			mpfr_div(a->val.mpfr,b->val.mpfr,c->val.mpfr,GMP_RNDU);
			break;
	}
  }
}

void elina_scalar_max(elina_scalar_t *op, elina_scalar_t * op1, elina_scalar_t *op2){
	bool res = elina_scalar_cmp(op1,op2);	
	if(res>=0){
		elina_scalar_set(op,op1);
        }
	else{
		elina_scalar_set(op,op2);
	}
}

void elina_scalar_min(elina_scalar_t *op, elina_scalar_t * op1, elina_scalar_t * op2){
	bool res = elina_scalar_cmp(op1,op2);	
	if(res>=0){
		elina_scalar_set(op,op2);
        }
	else{
		elina_scalar_set(op,op1);
	}
}


/* ********************************************************************** */
/* Square Root */
/* ********************************************************************** */

static inline void elina_scalar_sqrt_mpq(mpq_t up, mpq_t down, mpq_t b){
	mpz_t tmp; 
  	int perfect;
  	assert(mpq_sgn(b)>=0);
 	mpz_init(tmp);
  	mpz_mul(tmp,mpq_numref(b),mpq_denref(b));
  	perfect = mpz_perfect_square_p(tmp);
  	mpz_sqrt(mpq_numref(down),tmp);
  	if (perfect) mpz_set(mpq_numref(up),mpq_numref(down));
  	else mpz_add_ui(mpq_numref(up),mpq_numref(down),1);
  	mpz_set(mpq_denref(up),mpq_denref(b));
  	mpz_set(mpq_denref(down),mpq_denref(b));
  	mpq_canonicalize(up);
  	mpq_canonicalize(down);
  	mpz_clear(tmp);
}

static inline void elina_scalar_sqrt_double(double *up, double *down, double *b){
  double x;
  assert(*b>=0);
  x = sqrt(*b);
  assert(x*x>=*b); /* assumes round towards +oo! */
  if (x*x==*b) *down = x;
  else *down = nextafter(x,0);
  *up = x;
}


static inline void elina_scalar_sqrt_mpfr(mpfr_t up, mpfr_t down, mpfr_t b)
{ 
  mpfr_sqrt(up,b,GMP_RNDU);
  mpfr_sqrt(down,b,GMP_RNDD);
}

void elina_scalar_sqrt(elina_scalar_t *up, elina_scalar_t *down, elina_scalar_t *b, elina_scalar_discr_t discr)
{
  if(up!=b){
  	elina_scalar_reinit(up,discr);
  }
  if(down!=b){
  	elina_scalar_reinit(down,discr);
  }
  if(b->discr!=discr){
  	elina_scalar_convert(b,discr);
  }
  //elina_scalar_reinit(b,discr);
  if (elina_scalar_infty(b)) {
    elina_scalar_set_infty(up,1);
    elina_scalar_set_infty(down,1);
  }
  else {
    switch(discr){
	case ELINA_SCALAR_MPQ:
		elina_scalar_sqrt_mpq(up->val.mpq,down->val.mpq,b->val.mpq);
		break;
	case ELINA_SCALAR_DOUBLE:
		elina_scalar_sqrt_double(&up->val.dbl,&down->val.dbl,&b->val.dbl);
		break;
	case ELINA_SCALAR_MPFR:
		elina_scalar_sqrt_mpfr(up->val.mpfr,down->val.mpfr,b->val.mpfr);
		break;
    }
  }
  
}



void int_set_elina_scalar(long int *a, elina_scalar_t * scalar){
	switch(scalar->discr){
		case ELINA_SCALAR_MPQ:
			int_set_mpq(a, scalar->val.mpq);
			break;
		case ELINA_SCALAR_DOUBLE:
			int_set_double(a,scalar->val.dbl);
			break;
		case ELINA_SCALAR_MPFR:
			int_set_mpfr(a,scalar->val.mpfr);
			break;
	}
}

/* ********************************************************************** */
/* Power */
/* ********************************************************************** */

static inline int elina_scalar_pow_mpq(mpq_t up, mpq_t down, mpq_t b, unsigned long int n)
{
  mpz_pow_ui(mpq_numref(up), mpq_numref(b), n);
  mpz_pow_ui(mpq_denref(up), mpq_denref(b), n);
  mpq_canonicalize(up);
  mpz_set(mpq_numref(down), mpq_numref(up));
  mpz_set(mpq_denref(down), mpq_denref(up));
 return 0;
}


static inline int elina_scalar_pow_double(double * up, double *down, double *b, unsigned long int n)
{
  /* we cannot rely on pow, we implement a simple fast power */
  double u, d, fu, fd;
  int sign;
  if (*b < 0) { fu = fd = -*b; sign = n&1;}
  else { fu = fd = *b; sign = 0; }
  u = d = 1;
  while (n) {
    if (n & 1) { u *= fu; d = - ((-d) * fd); }
    fu = fu * fu;
    fd = - ((-fd) * fd);
    n >>= 1;
  }
  if (sign) {
    *up = -d;
    *down = -u;
  }
  else {
    *up = u;
    *down = d;
  }
  return !(isfinite(u) && isfinite(d));
}

static inline int  elina_scalar_pow_mpfr(mpfr_t up, mpfr_t down, mpfr_t b, unsigned long int n)
{
  mpfr_pow_ui(up, b, n, GMP_RNDU);
  mpfr_pow_ui(down, b, n, GMP_RNDD);
  return 0;
}

void elina_scalar_pow(elina_scalar_t *up, elina_scalar_t *down, elina_scalar_t *b, unsigned long n, elina_scalar_discr_t discr)
{
  if(up!=b){
  	elina_scalar_reinit(up,discr);
  }
  if(down!=b){
  	elina_scalar_reinit(down,discr);
  }
  //elina_scalar_reinit(b,discr);
  if (elina_scalar_infty(b)) {
    elina_scalar_set_infty(up, 1);
    elina_scalar_set_infty(down, 1);
  }
  else{
    bool res = true;
    switch(discr){
	case ELINA_SCALAR_MPQ:
		res = elina_scalar_pow_mpq(up->val.mpq,down->val.mpq,b->val.mpq,n);
		break;
	case ELINA_SCALAR_DOUBLE:
		res = elina_scalar_pow_double(&up->val.dbl,&down->val.dbl,&b->val.dbl,n);
		break;
	case ELINA_SCALAR_MPFR:
		res = elina_scalar_pow_mpfr(up->val.mpfr,down->val.mpfr,b->val.mpfr,n);
		break;
     }
     if (res) {
       elina_scalar_set_infty(up,1);
       if (n & 1) elina_scalar_set_infty(down,-1);
       else elina_scalar_set_to_int(down, 0,discr);
    }
  }
}

/* ********************************************************************** */
/* Truncate */
/* ********************************************************************** */
static inline void elina_scalar_trunc_mpq(mpq_t a, mpq_t b)
{
  mpz_tdiv_q(mpq_numref(a),mpq_numref(b),mpq_denref(b));
  mpz_set_si(mpq_denref(a),1);
}

static inline void elina_scalar_trunc_double(double * a, double *b){ 
	*a = trunc(*b); 
}

static inline void elina_scalar_trunc_mpfr(mpfr_t a, mpfr_t b){ 
	mpfr_rint_trunc(a,b,GMP_RNDU); 
}

void elina_scalar_trunc(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_discr_t discr){
  if(a!=b){
  	elina_scalar_reinit(a,discr);
  }
  if(b->discr!=discr){
  	elina_scalar_convert(b,discr);
  }
  //elina_scalar_reinit(b,discr);
  if (elina_scalar_infty(b)){
	 elina_scalar_set_infty(a,elina_scalar_sgn(b));
  }
  else { 
	switch(discr){
	case ELINA_SCALAR_MPQ:
		elina_scalar_trunc_mpq(a->val.mpq,b->val.mpq);
		break;
	case ELINA_SCALAR_DOUBLE:
		elina_scalar_trunc_double(&a->val.dbl,&b->val.dbl);
		break;
	case ELINA_SCALAR_MPFR:
		elina_scalar_trunc_mpfr(a->val.mpfr,b->val.mpfr);
		break;
       }
  }
}


/* ********************************************************************** */
/* Ceiling */
/* ********************************************************************** */

static inline void elina_scalar_ceil_mpq(mpq_t a, mpq_t b)
{
  mpz_cdiv_q(mpq_numref(a),mpq_numref(b),mpq_denref(b));
  mpz_set_si(mpq_denref(a),1);
}

static inline void elina_scalar_ceil_double(double * a, double *b){ 
	*a = ceil(*b); 
}

static inline void elina_scalar_ceil_mpfr(mpfr_t a, mpfr_t b){ 
	mpfr_rint_ceil(a,b,GMP_RNDU); 
}


void elina_scalar_ceil(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_discr_t discr){
	if(a!=b){
		elina_scalar_reinit(a,discr);
	}
	if(b->discr!=discr){
  		elina_scalar_convert(b,discr);
  	}
  	//elina_scalar_reinit(b,discr);
	 if (elina_scalar_infty(b)){
	 	elina_scalar_set_infty(a,elina_scalar_sgn(b));
  	}
	else { 
		switch(discr){
		case ELINA_SCALAR_MPQ:
			elina_scalar_ceil_mpq(a->val.mpq,b->val.mpq);
			break;
		case ELINA_SCALAR_DOUBLE:
			elina_scalar_ceil_double(&a->val.dbl,&b->val.dbl);
			break;
		case ELINA_SCALAR_MPFR:
			elina_scalar_ceil_mpfr(a->val.mpfr,b->val.mpfr);
			break;
	       }
	}
}



/* ********************************************************************** */
/* Floor */
/* ********************************************************************** */

static inline void elina_scalar_floor_mpq(mpq_t a, mpq_t b)
{
  mpz_fdiv_q(mpq_numref(a),mpq_numref(b),mpq_denref(b));
  mpz_set_si(mpq_denref(a),1);
}

static inline void elina_scalar_floor_double(double * a, double *b){ 
	*a = floor(*b); 
}

static inline void elina_scalar_floor_mpfr(mpfr_t a, mpfr_t b){ 
	mpfr_rint_floor(a,b,GMP_RNDU); 
}


void elina_scalar_floor(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_discr_t discr){
	if(a!=b){
		elina_scalar_reinit(a,discr);
	}
	if(b->discr!=discr){
  		elina_scalar_convert(b,discr);
  	}
  	//elina_scalar_reinit(b,discr);
	 if (elina_scalar_infty(b)){
	 	elina_scalar_set_infty(a,elina_scalar_sgn(b));
  	}
	else { 
		
		switch(discr){
		case ELINA_SCALAR_MPQ:
			elina_scalar_floor_mpq(a->val.mpq,b->val.mpq);
			break;
		case ELINA_SCALAR_DOUBLE:
			elina_scalar_floor_double(&a->val.dbl,&b->val.dbl);
			break;
		case ELINA_SCALAR_MPFR:
			elina_scalar_floor_mpfr(a->val.mpfr,b->val.mpfr);
			break;
	       }
	}
}



/* ********************************************************************** */
/* To single precision float */
/* ********************************************************************** */
static inline bool elina_scalar_fits_float_mpq(mpq_t a){
	return ((int)mpz_sizeinbase(mpq_numref(a),2)-
	  (int)mpz_sizeinbase(mpq_denref(a),2)<126); 
}

static inline bool elina_scalar_fits_float_double(double a){
	 int e;
  	 frexp(a,&e);
  	 return (e<127);
}

static inline bool elina_scalar_fits_float_mpfr(mpfr_t a){
	return mpfr_number_p(a) && mpfr_get_exp(a)<126;
}



bool elina_scalar_fits_float(elina_scalar_t *b, elina_scalar_discr_t discr){
	switch(discr){
		case ELINA_SCALAR_MPQ:
			return elina_scalar_fits_float_mpq(b->val.mpq);
		case ELINA_SCALAR_DOUBLE:
			return elina_scalar_fits_float_double(b->val.dbl);
		case ELINA_SCALAR_MPFR:
			return elina_scalar_fits_float_mpfr(b->val.mpfr);
	}
	return false;
} 


void elina_scalar_to_float(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_discr_t discr)
{
  if (elina_scalar_infty(b) || !elina_scalar_fits_float(b,discr)){
    elina_scalar_set_infty(a,elina_scalar_sgn(b));
  }
  else {
    double d;
    elina_double_set_scalar(&d,b,GMP_RNDU);
    elina_scalar_set_double(a,(double)((float)d));
    elina_scalar_convert(a,discr);
  }
}


/* ********************************************************************** */
/* To double precision float */
/* ********************************************************************** */

static inline bool elina_scalar_fits_double_mpq(mpq_t a){
	return ((int)mpz_sizeinbase(mpq_numref(a),2)-
	  (int)mpz_sizeinbase(mpq_denref(a),2)<1022); 
}


static inline bool elina_scalar_fits_double_mpfr(mpfr_t a){
	 return mpfr_number_p(a) && mpfr_get_exp(a)<1022;
}



bool elina_scalar_fits_double(elina_scalar_t *b, elina_scalar_discr_t discr){
	switch(discr){
		case ELINA_SCALAR_MPQ:
			return elina_scalar_fits_double_mpq(b->val.mpq);
		case ELINA_SCALAR_DOUBLE:
			return true;
		case ELINA_SCALAR_MPFR:
			return elina_scalar_fits_double_mpfr(b->val.mpfr);
	}
	return false;
} 


void elina_scalar_to_double(elina_scalar_t *a, elina_scalar_t *b, elina_scalar_discr_t discr)
{
  if (elina_scalar_infty(b) || !elina_scalar_fits_double(b,discr)){
    elina_scalar_set_infty(a,elina_scalar_sgn(b));
  }
  else {
    double d;
    elina_double_set_scalar(&d,b,GMP_RNDU);
    elina_scalar_set_double(a,d);
    elina_scalar_convert(a,discr);
  }
}

bool elina_scalar_is_integer(elina_scalar_t *a){
	switch(a->discr){
		case ELINA_SCALAR_MPQ:
			return (mpz_cmp_ui(mpq_denref(a->val.mpq),1)==0); 	
		case ELINA_SCALAR_DOUBLE:
			 return ceil(a->val.dbl) == a->val.dbl;
		case ELINA_SCALAR_MPFR:
			return mpfr_integer_p(a->val.mpfr);
	}
	return false;
}

void elina_scalar_set_to_int(elina_scalar_t *a, elina_int_t i, elina_scalar_discr_t discr){
	elina_scalar_reinit(a,discr);
	switch(discr){
		case ELINA_SCALAR_MPQ:
			mpq_set_si(a->val.mpq,i,1);
			break;
		case ELINA_SCALAR_DOUBLE:
			a->val.dbl = (double)i;
			break;
		case ELINA_SCALAR_MPFR:
			mpfr_set_si(a->val.mpfr,i,GMP_RNDU);
			break;
	}
}



/* ********************************************************************** */
/* Multiply by power of two */
/* ********************************************************************** */

static inline void elina_scalar_mul_2exp_mpq(mpq_t a, mpq_t b, int c)
{
  if (c>=0) mpq_mul_2exp(a,b,c);
  else mpq_div_2exp(a,b,-c);
}

static inline void elina_scalar_mul_2exp_double(double* a, double* b, int c){
 *a = ldexp(*b,c);
}

static inline void elina_scalar_mul_2exp_mpfr(mpfr_t a, mpfr_t b, int c){
	mpfr_mul_2si(a,b,c,GMP_RNDU);
}


void elina_scalar_mul_2exp(elina_scalar_t *a, elina_scalar_t *b, int c , elina_scalar_discr_t discr){
	if(a!=b){
		elina_scalar_reinit(a,discr);
	}
	switch(discr){
		case ELINA_SCALAR_MPQ:
			elina_scalar_mul_2exp_mpq(a->val.mpq,b->val.mpq,c);
			break;
		case ELINA_SCALAR_DOUBLE:
			elina_scalar_mul_2exp_double(&a->val.dbl,&b->val.dbl,c);
			break;
		case ELINA_SCALAR_MPFR:
			elina_scalar_mul_2exp_mpfr(a->val.mpfr,b->val.mpfr,c);
			break;
	}
}




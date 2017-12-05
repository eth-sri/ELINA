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

#ifndef _OPT_PK_CONFIG_H_
#define _OPT_PK_CONFIG_H_

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#if defined(VECTOR)

 #include <immintrin.h>
 #include "vector_intrin.h"

#else
 #define v_length 1

#endif



#ifdef __cplusplus
#define HAS_BOOL
extern "C" {
#endif

#ifndef HAS_BOOL
#define HAS_BOOL
typedef char bool;
static const bool false = 0;
static const bool true  = 1;
#endif

#define opt_numint_t long long int
#define min(a,b) a < b ? a : b
/* Extension to the num package */
/* size in words */
static inline size_t opt_numint_size(opt_numint_t a){
	return 1;
}


/* Do not change ! */
static const size_t opt_polka_cst = 1;
static const size_t opt_polka_eps = 2;
static inline opt_numint_t opt_numint_gcd2(opt_numint_t a, opt_numint_t b){
	while(b!=0 && a!=b){
		opt_numint_t tmp = b;
		b = a%b;
		a = tmp; 
	}
	return a;
}

static inline opt_numint_t opt_numint_gcd(opt_numint_t a, opt_numint_t b){
	
	if(a < 0){
		a = -a;
	}
	if(b < 0){
		b = -b;
	}
	return (a < b) ? opt_numint_gcd2(b,a) : opt_numint_gcd2(a, b);
}



static inline int opt_numint_sgn(opt_numint_t num){
	return ((num==0) ? 0 : (num < 0) ? -1 : 1);
}

static inline opt_numint_t opt_numint_abs(opt_numint_t num){
	return ((num < 0) ? -num : num);
}


static inline opt_numint_t opt_numint_lcm(opt_numint_t a, opt_numint_t b){
	opt_numint_t gcd = opt_numint_gcd(a,b);
	return (a/(gcd))*b;
}

static inline opt_numint_t opt_numint_fdiv(opt_numint_t num, opt_numint_t den){
	if((opt_numint_sgn(num) *opt_numint_sgn(den) < 0) && (num%den)){
		return num/den - 1;
	}
	else{
		return num/den;
	}
}

#ifdef __cplusplus
}
#endif

#endif


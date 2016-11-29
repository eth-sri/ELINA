/*
	Copyright 2016 Software Reliability Lab, ETH Zurich

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

		http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
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


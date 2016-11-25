
/* This file is part of the APRON Library, released under LGPL license.  Please
   read the COPYING file packaged in the distribution */

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
//#if defined (NUMINT_NATIVE)
//{ return 1; }
//#elif defined(NUMINT_MPZ)
//{ return mpz_size(a); }
//#else
//#error "Here"
//#endif

/* size in bits */
//#if defined (NUMINT_NATIVE)

//#if defined(NUMINT_LONGINT) || defined(NUMINT_LONGLONGINT)
//#if defined(__GNUC__)
//static inline size_t numint_size2(numint_t a)
//{
//  numint_t x;
//  if (a==0) return 0;
//  numint_abs(x,a);
//  return (size_t)(
//		  sizeof(numint_t)*8 - 
//#if defined(NUMINT_LONGINT)
//		  __builtin_clzl((unsigned long int)(*x))
//#else
//		  __builtin_clzll((unsigned long long int)(*x))
//#endif
//		  );
//}
//#else
//static inline size_t numint_size2(numint_t a)
//{ 
//  numint_t x;
//  size_t nbzero;
//  static size_t table[256] = {
//    8,
//    7,
//    6,6,
//    5,5,5,5,
//    4,4,4,4,4,4,4,4,
//    3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
//    2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
//    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
//    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
//    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
//    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
//    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
//    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
//  }    
//  if (a==0) return 0;
//  numint_abs(x,a);
  /* Iterate over leading zero octets */
//#if defined(NUMINT_LONGINT)
//  unsigned long int mask = (0xffUL << (sizeof(numint_t)-1)*8);
//#else
//  unsigned long long int mask = (0xffULL << (sizeof(numint_t)-1)*8);
//#endif
//  nbzero = 0;
//  while ( (*x & *mask)==0 ){
//    nbzero += 8;
//    *mask = *mask >> 256;
//  }
//  *x = (*x & *mask) >> ((sizeof(numint_t)-1)*8 - nbzero);
//  nbzero += table[*x];
//  return (sizeof(numint_t)*8 - nbzero);
//}
//#endif
//#else
//#error "Here"
//#endif

//#elif defined(NUMINT_MPZ)
//static inline size_t numint_size2(numint_t a)
//{ return mpz_sizeinbase(a,2); }
//#else
//#error "Here"
//#endif


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


/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
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


#ifndef __VECTOR_INTRIN_H
#define __VECTOR_INTRIN_H

#ifdef __cplusplus
extern "C" {
#endif

#if defined(SSE)

#define v_length 2
#define v_double_type __m128d 
#define v_int_type __m128i
#define v_load_double _mm_loadu_pd
#define v_store_double _mm_storeu_pd
#define v_set1_double  _mm_set1_pd
#define v_min_double _mm_min_pd
#define v_max_double _mm_max_pd
#define v_add_double _mm_add_pd
#define v_cmp_double _mm_cmp_pd
#define v_set1_int  _mm_set1_epi64x
#define v_double_to_int _mm_castpd_si128
#define v_test_int _mm_testc_si128

#else

#define v_length 4
#define v_double_type __m256d 
#define v_int_type __m256i
#define v_load_double _mm256_loadu_pd
#define v_store_double _mm256_storeu_pd 
#define v_set1_double _mm256_set1_pd
#define v_min_double _mm256_min_pd
#define v_max_double _mm256_max_pd
#define v_add_double _mm256_add_pd
#define v_cmp_double _mm256_cmp_pd
#define v_set1_int _mm256_set1_epi64x
#define v_double_to_int _mm256_castpd_si256
#define v_test_int _mm256_testc_si256

#endif

#ifdef __cplusplus
}
#endif

#endif

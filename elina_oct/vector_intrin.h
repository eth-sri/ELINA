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

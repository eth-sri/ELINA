/*
	Copyright 2015 Department of Computer Science, ETH Zurich

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


#ifndef __OPT_OCT_CLOSURE_DENSE_SCALAR_H_INCLUDED__
#define __OPT_OCT_CLOSURE_DENSE_SCALAR_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#include "opt_oct_hmat.h"

//void print_dense(double *m, int dim);

//double strong_closure_calc_perf_dense(double cycles, int dim);
bool strong_closure_dense_scalar(opt_oct_mat_t *m, double * temp1, double *temp2, int dim, bool is_int);
bool strengthning_int_dense_scalar(opt_oct_mat_t * result, double *temp, int n);
bool floyd_warshall_dense_scalar(opt_oct_mat_t *m, double * temp1, double *temp2, int dim, bool is_int);
bool strengthning_dense_scalar(opt_oct_mat_t * result, double *temp, int n);


#ifdef __cplusplus
}
#endif

#endif 

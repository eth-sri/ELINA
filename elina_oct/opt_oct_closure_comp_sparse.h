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


#ifndef __OPT_OCT_CLOSURE_COMP_SPARSE_H
#define __OPT_OCT_CLOSURE_COMP_SPARSE_H


#include "comp_list.h"
#include "opt_oct_hmat.h"



bool strong_closure_comp_sparse(opt_oct_mat_t *oo, double *temp1, double *temp2, unsigned short int *index1, unsigned short int *index2, int dim, bool is_int);
bool strengthning_int_comp_sparse(opt_oct_mat_t * oo,  unsigned short int * ind1, double *temp, int n);
void strengthening_comp_list(opt_oct_mat_t *oo,comp_list_t * cd, unsigned short int dim);
bool strengthning_comp_sparse(opt_oct_mat_t *oo, unsigned short int * ind1, double *temp, int n);
void compute_index_comp_sparse(double *result, unsigned short int *ca, unsigned short int comp_size, unsigned short int *index1, unsigned short int *index2, unsigned short int k, int dim);

#endif

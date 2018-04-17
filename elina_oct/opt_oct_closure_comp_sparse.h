/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
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

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

#ifndef __OPT_OCT_CLOSURE_DENSE_H_INCLUDED__
#define __OPT_OCT_CLOSURE_DENSE_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#include "opt_oct_hmat.h"
#include <immintrin.h>
#include "vector_intrin.h"

void print_dense(double *m, int dim);

double strong_closure_calc_perf_dense(double cycles, int dim);
bool strong_closure_dense(opt_oct_mat_t *m, double * temp1, double *temp2, int dim, bool is_int);
bool strengthning_int_dense(opt_oct_mat_t * result, double *temp, int n);
bool floyd_warshall_dense(opt_oct_mat_t *m, double * temp1, double *temp2, int dim, bool is_int);
bool strengthning_dense(opt_oct_mat_t * result, double *temp, int n);

#ifdef __cplusplus
}
#endif

#endif 

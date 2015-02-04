#ifndef __OPT_OCT_CLOSURE_DENSE_H_INCLUDED__
#define __OPT_OCT_CLOSURE_DENSE_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#include "opt_oct_hmat.h"

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

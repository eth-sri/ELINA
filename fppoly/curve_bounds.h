#ifndef __CURVE_BOUNDS_H_INCLUDED__
#define __CURVE_BOUNDS_H_INCLUDED__

// This header is also used in a different package and thus should have
// no additional dependencies.

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

void compute_S_curve_bounds(double x_lb, double x_ub, bool is_sigmoid,
                            double* k_lb, double* b_lb, double* k_ub, double* b_ub);
                            
#ifdef __cplusplus
 }
#endif

#endif

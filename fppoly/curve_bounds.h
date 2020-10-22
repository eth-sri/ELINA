#ifndef __CURVE_BOUNDS_H_INCLUDED__
#define __CURVE_BOUNDS_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#include "fppoly.h"

void compute_S_curve_bounds(double x_lb, double x_ub, bool is_sigmoid,
                            double *k_lb, double *b_lb, double *k_ub,
                            double *b_ub);

#ifdef __cplusplus
}
#endif

#endif

#pragma once

void compute_S_curve_bounds(double x_lb, double x_ub, bool is_sigmoid,
                            double *k_lb, double *b_lb, double *k_ub,
                            double *b_ub);

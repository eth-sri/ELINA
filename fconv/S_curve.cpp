#include "asrt.h"
#include <math.h>

void compute_curve_bounds(double x_bound, bool is_sigmoid,
                     double& k_lb, double& b_lb, double& k_ub, double& b_ub) {
    ASRTF(x_bound != 0, "x_bound cannot be zero.");
    if (!is_sigmoid) {
        double y_bound = tanh(x_bound);
        if (x_bound < 0) {
            k_lb = 1 - y_bound * y_bound;
            b_lb = y_bound - x_bound * k_lb;
            k_ub = y_bound / x_bound;
            b_ub = 0;
        } else {
            k_lb = y_bound / x_bound;
            b_lb = 0;
            k_ub = 1 - y_bound * y_bound;
            b_ub = y_bound - x_bound * k_ub;
        }
    } else {
        // Numerically stable sigmoid
        // http://timvieira.github.io/blog/post/2014/02/11/exp-normalize-trick/
        if (x_bound < 0) {
            double y_bound = exp(x_bound);
            y_bound = y_bound / (1 + y_bound);
            k_lb = y_bound * (1 - y_bound);
            b_lb = y_bound - x_bound * k_lb;
            k_ub = (y_bound - 0.5) / x_bound;
            b_ub = 0.5;
        } else {
            double y_bound = 1 / (1 + exp(-x_bound));
            k_lb = (y_bound - 0.5) / x_bound;
            b_lb = 0.5;
            k_ub = y_bound * (1 - y_bound);
            b_ub = y_bound - x_bound * k_ub;
        }
    }
    b_lb -= 1.0E-5;
    b_ub += 1.0E-5;
}

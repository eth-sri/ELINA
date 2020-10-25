#include "asrt.h"
#include <math.h>
#include <cfenv>
#include <cassert>

inline void upward() { fesetround(FE_UPWARD); }

inline void downward() {
    fesetround(FE_DOWNWARD);
}

inline void round(bool up) { fesetround(up ? FE_UPWARD : FE_DOWNWARD); }

inline double sigm(double x, bool up) {
    round(!up);
    double denom = 1 + exp(-x);
    round(up);
    return 1 / denom;
}

inline double my_tanh(double x, bool up) {
  round(up);
  return tanh(x);
}

// y belongs to [0, 1]
// Returns infimum of the tangent line in the point.
double sigm_tang_k_inf(double y) {
  assert(0 <= y && y <= 1 && "y should be within range.");
  downward();
  return y * (1 - y);
}

// y belongs to [-1, 1]
double tanh_tang_k_inf(double y) {
  assert(-1 <= y && y <= 1 && "y should be within range.");
  upward();
  double temp = y * y;
  downward();
  return 1 - temp;
}

double k_of_S_curve_chord_inf(double b_chord, double x_point, double y_point) {
  round(y_point < b_chord);
  double numer = y_point - b_chord;
  downward();
  double k = numer / x_point;
  assert(k >= 0 && "k expected to be non-negative.");
  return k;
}

double b_of_S_curve_tangent(double k_tang, double x_point, double y_point) {
  assert(k_tang >= 0 && "k of tangent is non-negative.");
  round(x_point < 0);
  double temp = k_tang * x_point;
  round(x_point > 0);
  return y_point - temp;
}

void compute_curve_bounds(double x_bound, bool is_sigm,
                          double& k_lb, double& b_lb, double& k_ub, double& b_ub) {
    ASRTF(x_bound != 0, "x_bound cannot be zero.");

    round(true);
    ASRTF(fegetround() == FE_UPWARD, "making sure fe upward is set correctly.");

    double y_inf, y_sup;
    if (is_sigm) {
      y_inf = sigm(x_bound, false);
      y_sup = sigm(x_bound, true);
    } else {
      y_inf = my_tanh(x_bound, false);
      y_sup = my_tanh(x_bound, true);
    }

    double b_chord = is_sigm ? 0.5 : 0.0;

    if (x_bound > 0) {
      b_lb = b_chord;
      k_lb = k_of_S_curve_chord_inf(b_chord, x_bound, y_inf);
      k_ub = is_sigm ? sigm_tang_k_inf(y_sup) : tanh_tang_k_inf(y_sup);
      b_ub = b_of_S_curve_tangent(k_ub, x_bound, y_sup);
    } else {
      b_ub = b_chord;
      k_ub = k_of_S_curve_chord_inf(b_chord, x_bound, y_sup);
      k_lb = is_sigm ? sigm_tang_k_inf(y_inf) : tanh_tang_k_inf(y_inf);
      b_lb = b_of_S_curve_tangent(k_lb, x_bound, y_inf);
    }

    fesetround(FE_TONEAREST);
//
//    b_lb -= 1.0;
//    b_ub += 1.0;
}

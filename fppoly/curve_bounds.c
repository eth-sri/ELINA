#include "curve_bounds.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double PRECISION = 1.0E-7;

// Numerically stable sigmoid.
inline double sigm(double x) {
  if (x < 0) {
    double v = exp(x);
    return v / (1 + v);
  }
  return x == 0 ? 0.5 : 1 / (1 + exp(-x));
}

// Score proportional to area that has to be minimized.
inline double tangent_bound_area_score(double x, double x_eval1, double x_eval2,
                                       bool is_upper, bool is_sigm) {
  double y, k;
  if (is_sigm) {
    y = sigm(x);
    k = (1 - y) * y;
  } else {
    y = tanh(x);
    k = 1 - y * y;
  }
  double b = y - k * x;
  double score = k * x_eval1 + k * x_eval2 + 2 * b;
  if (!is_upper) {
    score = -score;
  }
  return score;
}

/*
 * Any tangent line in [x1, x2] is sound and the area on [x_eval1, x_eval2] is
 * minimized by doing binary search in that region. Although we don't have proof
 * of function unimodality yet, we observe it in practice. In unexpected case it
 * doesn't hold, we return *good enough* tangent line.
 */
void find_best_tangent_line(double *k, double *b, double x1, double x2,
                            double x_eval1, double x_eval2, bool is_upper,
                            bool is_sigm) {
  double x_left, x_right;
  if (x1 < x2) {
    x_left = x1;
    x_right = x2;
  } else {
    x_left = x2;
    x_right = x1;
  }

  double v_left =
      tangent_bound_area_score(x_left, x_eval1, x_eval2, is_upper, is_sigm);
  double v_right =
      tangent_bound_area_score(x_right, x_eval1, x_eval2, is_upper, is_sigm);

  while (x_right - x_left > PRECISION) {
    double x_left_third = x_left + (x_right - x_left) / 3;
    double x_right_third = x_right - (x_right - x_left) / 3;
    double v_left_third = tangent_bound_area_score(x_left_third, x_eval1,
                                                   x_eval2, is_upper, is_sigm);
    double v_right_third = tangent_bound_area_score(x_right_third, x_eval1,
                                                    x_eval2, is_upper, is_sigm);
    bool unimodality = (v_left_third <= v_left + PRECISION ||
                        v_left_third <= v_right + PRECISION) &&
                       (v_right_third <= v_left + PRECISION ||
                        v_right_third <= v_right + PRECISION);
    // Checking unimodality condition in DEBUG mode.
    // Mostly for curiosity as it is not required for correctness, but required
    // for optimality.
    assert(unimodality && "Expected unimodal.");
    if (v_left_third > v_right_third) {
      x_left = x_left_third;
      v_left = v_left_third;
    } else {
      x_right = x_right_third;
      v_right = v_right_third;
    }
  }
  double x;
  if (v_left <= v_right) {
    x = x_left;
  } else {
    x = x_right;
  }
  double y;
  if (is_sigm) {
    y = sigm(x);
    *k = (1 - y) * y;
  } else {
    y = tanh(x);
    *k = 1 - y * y;
  }
  *b = y - *k * x;
}

// Returns 0 iff line goes through x_convex and touches curve in x.
inline double curve_touching_condition(double x, double x_convex,
                                       bool is_sigm) {
  double y, k, y_convex;
  if (is_sigm) {
    y = sigm(x);
    k = (1 - y) * y;
    y_convex = sigm(x_convex);
  } else {
    y = tanh(x);
    k = 1 - y * y;
    y_convex = tanh(x_convex);
  }
  return fabs(k * (x - x_convex) + y_convex - y);
}

double find_x_star(double x_concave, double x_convex, bool is_sigm) {
  double x_left, x_right;
  if (x_concave > 0) {
    x_left = 0;
    x_right = x_concave;
  } else {
    x_left = x_concave;
    x_right = 0;
  }

  double v_left = curve_touching_condition(x_left, x_convex, is_sigm);
  double v_right = curve_touching_condition(x_right, x_convex, is_sigm);

  // Ternary search algorithm.
  while (x_right - x_left > PRECISION / 10) {
    double x_left_third = x_left + (x_right - x_left) / 3;
    double x_right_third = x_right - (x_right - x_left) / 3;
    double v_left_third =
        curve_touching_condition(x_left_third, x_convex, is_sigm);
    double v_right_third =
        curve_touching_condition(x_right_third, x_convex, is_sigm);

    bool unimodality = (v_left_third <= v_left + PRECISION ||
                        v_left_third <= v_right + PRECISION) &&
                       (v_right_third <= v_left + PRECISION ||
                        v_right_third <= v_right + PRECISION);
    assert(unimodality && "Expected unimodality");
    if (!unimodality) {
      abort();
    }
    if (v_left_third > v_right_third) {
      x_left = x_left_third;
      v_left = v_left_third;
    } else {
      x_right = x_right_third;
      v_right = v_right_third;
    }
  }

  double x, v;
  if (v_left < v_right) {
    x = x_left;
    v = v_left;
  } else {
    x = x_right;
    v = v_right;
  }
  assert(fabs(v) <= PRECISION &&
         "Checking that the touching line has been found.");
  if (fabs(v) > PRECISION) {
    abort();
  }

  return x;
}

void get_optimal_curve_bound(double *k, double *b, double x_lb, double x_ub,
                             bool is_upper, bool is_sigm) {
  assert(x_lb < x_ub && "x_lb < x_ub");
  double y_lb, y_ub;
  if (is_sigm) {
    y_lb = sigm(x_lb);
    y_ub = sigm(x_ub);
  } else {
    y_lb = tanh(x_lb);
    y_ub = tanh(x_ub);
  }

  if ((x_ub <= 0 && is_upper) || (0 <= x_lb && !is_upper)) {
    // In this case the bound is simply the chord.
    *k = (y_ub - y_lb) / (x_ub - x_lb);
    *b = y_ub - *k * x_ub;
    return;
  }
  if (x_ub <= 0 || 0 <= x_lb) {
    // Function is either convex or concave, thus any tangent line is a valid
    // bound. We are looking for the one minimizing the area.
    find_best_tangent_line(k, b, x_lb, x_ub, x_lb, x_ub, is_upper, is_sigm);
    return;
  }

  double x_convex, x_concave;
  if (is_upper) {
    x_convex = x_lb;
    x_concave = x_ub;
  } else {
    x_convex = x_ub;
    x_concave = x_lb;
  }

  double k_chord = (y_ub - y_lb) / (x_ub - x_lb);
  double k_concave;
  if (is_sigm) {
    double y = sigm(x_concave);
    k_concave = (1 - y) * y;
  } else {
    double y = tanh(x_concave);
    k_concave = 1 - y * y;
  }
  if (k_chord <= k_concave) {
    // Chord is sound and optimal.
    *k = k_chord;
    *b = y_ub - *k * x_ub;
    return;
  }

  // Finding the x_star. Condition that line goes through convex x edge and
  // touches S curve.
  double x_star = find_x_star(x_concave, x_convex, is_sigm);

  // The sound region that contains optimal line is [x_star, x_concave].
  // Find the best line there with respect to minimizing area on [x_lb, x_ub].
  find_best_tangent_line(k, b, x_star, x_concave, x_lb, x_ub, is_upper,
                         is_sigm);
}

// Values beyond which I just take a simple bound due to potential numerical
// stability issues.
const double TANH_LIM = 5;
const double SIGM_LIM = 10;

void compute_S_curve_bounds(double x_lb, double x_ub, bool is_sigm,
                            double *k_lb, double *b_lb, double *k_ub,
                            double *b_ub) {
  assert(x_lb <= x_ub && "x_lb <= x_ub");
  if (x_lb > x_ub) {
    abort();
  }

  if (x_ub - x_lb <= 1.0E-5) {
    if (is_sigm) {
      *b_lb = sigm(x_lb);
      *b_ub = sigm(x_ub);
    } else {
      *b_lb = tanh(x_lb);
      *b_ub = tanh(x_ub);
    }
    *k_lb = 0;
    *k_ub = 0;
  } else {
    if ((is_sigm && x_lb <= -SIGM_LIM) || (!is_sigm && x_lb <= -TANH_LIM)) {
      *k_lb = 0;
      *b_lb = is_sigm ? sigm(x_lb) : tanh(x_lb);
    } else {
      get_optimal_curve_bound(k_lb, b_lb, x_lb, x_ub, false, is_sigm);
    }
    if ((is_sigm && x_ub >= SIGM_LIM) || (!is_sigm && x_ub >= TANH_LIM)) {
      *k_ub = 0;
      *b_ub = is_sigm ? sigm(x_ub) : tanh(x_ub);
    } else {
      get_optimal_curve_bound(k_ub, b_ub, x_lb, x_ub, true, is_sigm);
    }
  }

  // Adjusting for numerical soundness (with big safe margin).
  *b_lb -= 1.0E-5;
  *b_ub += 1.0E-5;
}

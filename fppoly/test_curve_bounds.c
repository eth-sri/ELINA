#include "curve_bounds.h"
#include <execinfo.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

const double INF_MARGIN = 1;
const double EPS = 1.0E-7;

double d_rand(double lb, double ub) {
  return lb + (rand() / (RAND_MAX / (ub - lb)));
}

double min(double x1, double x2) {
  if (x1 < x2) {
    return x1;
  }
  return x2;
}

double max(double x1, double x2) {
  if (x1 > x2) {
    return x1;
  }
  return x2;
}

double margin_in_point(double x, double k, double b, bool is_upper,
                       bool is_sigm) {
  double y = is_sigm ? 1 / (1 + exp(-x)) : tanh(x);
  double margin = k * x + b - y;
  if (!is_upper) {
    margin = -margin;
  }
  return margin;
}

double margin_in_concave_segm(double k, double b, double x_lb, double x_ub,
                              bool is_upper, bool is_sigm) {
  if ((is_upper && x_ub <= 0) || (!is_upper && x_lb >= 0)) {
    return INF_MARGIN;
  }

  double x_left, x_right;
  if (is_upper) {
    x_left = max(0.0, x_lb);
    x_right = x_ub;
  } else {
    x_left = x_lb;
    x_right = min(0.0, x_ub);
  }

  double v_left = margin_in_point(x_left, k, b, is_upper, is_sigm);
  double v_right = margin_in_point(x_right, k, b, is_upper, is_sigm);

  while (x_right - x_left > EPS) {
    double x_left_third = x_left + (x_right - x_left) / 3;
    double x_right_third = x_right - (x_right - x_left) / 3;

    double v_left_third =
        margin_in_point(x_left_third, k, b, is_upper, is_sigm);
    double v_right_third =
        margin_in_point(x_right_third, k, b, is_upper, is_sigm);

    bool unimodality = (v_left_third <= v_left + EPS / 20 ||
                        v_left_third <= v_right + EPS / 20) &&
                       (v_right_third <= v_left + EPS / 20 ||
                        v_right_third <= v_right + EPS / 20);
    if (!unimodality) {
      printf("Expected unimodal lb %f ub %f %s for %s", x_lb, x_ub,
             is_upper ? "upper" : "lower", is_sigm ? "Sigm" : "Tanh");
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
  return min(v_left, v_right);
}

double margin_in_convex_segm(double k, double b, double x_lb, double x_ub,
                             bool is_upper, bool is_sigm) {
  if ((is_upper && x_lb > 0) || (!is_upper && x_ub < 0)) {
    return INF_MARGIN;
  }

  double x_left, x_right;
  if (is_upper) {
    x_left = x_lb;
    x_right = min(0.0, x_ub);
  } else {
    x_left = max(0.0, x_lb);
    x_right = x_ub;
  }

  double v_left = margin_in_point(x_left, k, b, is_upper, is_sigm);
  double v_right = margin_in_point(x_right, k, b, is_upper, is_sigm);
  return min(v_left, v_right);
}

double margin_in_segm(double x_lb, double x_ub, double k_lb, double b_lb,
                      double k_ub, double b_ub, bool is_sigm) {
  double lb_concave =
      margin_in_concave_segm(k_lb, b_lb, x_lb, x_ub, false, is_sigm);
  double lb_convex =
      margin_in_convex_segm(k_lb, b_lb, x_lb, x_ub, false, is_sigm);
  double ub_concave =
      margin_in_concave_segm(k_ub, b_ub, x_lb, x_ub, true, is_sigm);
  double ub_convex =
      margin_in_convex_segm(k_ub, b_ub, x_lb, x_ub, true, is_sigm);
  double margins[4] = {lb_concave, lb_convex, ub_concave, ub_convex};
  double min_margin = INF_MARGIN;
  for (int i = 0; i < 4; i++) {
    if (margins[i] < min_margin) {
      min_margin = margins[i];
    }
  }
  return min_margin;
}

void test_S_curve_bounds(double LB, double UB, int N_samples, bool is_sigm) {
  printf("Running test lb %f ub %f samples %d for %s\n", LB, UB, N_samples,
         is_sigm ? "Sigm" : "Tanh");

  for (int i = 0; i < N_samples; i++) {
    double x_lb = d_rand(LB, UB);
    double x_ub = d_rand(LB, UB);
    if (x_lb > x_ub) {
      double tmp = x_ub;
      x_ub = x_lb;
      x_lb = tmp;
    }

    double k_lb, b_lb, k_ub, b_ub;
    compute_S_curve_bounds(x_lb, x_ub, is_sigm, &k_lb, &b_lb, &k_ub, &b_ub);
    double min_margin =
        margin_in_segm(x_lb, x_ub, k_lb, b_lb, k_ub, b_ub, is_sigm);

    if (min_margin < EPS) {
      printf("Failed x_lb %f x_ub %f for %s\n", x_lb, x_ub,
             is_sigm ? "Sigm" : "Tanh");
      abort();
    }
  }
  printf("\tpassed\n");
}

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

int main() {
  signal(SIGSEGV, handler);
  srand((unsigned)time(NULL));

  int N_samples = 1000;
#define N_BOUNDS 8
  double bounds[N_BOUNDS] = {0.001, 0.01, 0.1, 0.5, 1, 10, 100, 1000};
  for (int i = 0; i < N_BOUNDS; i++) {
    test_S_curve_bounds(-bounds[i], bounds[i], N_samples, false);
    test_S_curve_bounds(-bounds[i], bounds[i], N_samples, true);
  }

  return 0;
}

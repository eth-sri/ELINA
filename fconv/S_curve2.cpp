#include "S_curve2.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "asrt.h"

const double PRECISION = 1.0E-7;

const double SMALL = 1.0E-4;

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
        k  = 1 - y * y;
    }
    double b = y - k * x;
    double score = k * x_eval1 + k * x_eval2 + 2 * b;
    if (!is_upper) {
        score = -score;
    }
    return score;
}

/*
 * Any tangent line in [x1, x2] is sound and the area on [x_eval1, x_eval2] is minimized
 * by doing binary search in that region.
 * Although we don't have proof of function unimodality yet, we observe it in practice.
 * In unexpected case it doesn't hold, we return *good enough* tangent line.
 */
void find_best_tangent_line(double *k, double *b,
                            double x1, double x2,
                            double x_eval1, double x_eval2,
                            bool is_upper, bool is_sigm) {
    double x_left, x_right;
    if (x1 < x2) {
        x_left = x1;
        x_right = x2;
    } else {
        x_left = x2;
        x_right = x1;
    }

    double v_left = tangent_bound_area_score(x_left, x_eval1, x_eval2, is_upper, is_sigm);
    double v_right = tangent_bound_area_score(x_right, x_eval1, x_eval2, is_upper, is_sigm);

    while (x_right - x_left > PRECISION) {
        double x_left_third = x_left + (x_right - x_left) / 3;
        double x_right_third = x_right - (x_right - x_left) / 3;
        double v_left_third = tangent_bound_area_score(x_left_third, x_eval1, x_eval2, is_upper, is_sigm);
        double v_right_third = tangent_bound_area_score(x_right_third, x_eval1, x_eval2, is_upper, is_sigm);
        bool unimodality =
                (v_left_third <= v_left + PRECISION || v_left_third <= v_right + PRECISION) &&
                (v_right_third <= v_left + PRECISION || v_right_third <= v_right + PRECISION);
        // Checking unimodality condition in DEBUG mode.
        // Mostly for curiosity as it is not required for correctness, but required for optimality.
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
inline double curve_touching_condition(double x, double x_convex, bool is_sigm) {
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
        double v_left_third = curve_touching_condition(x_left_third, x_convex, is_sigm);
        double v_right_third = curve_touching_condition(x_right_third, x_convex, is_sigm);

        bool unimodality =
                (v_left_third <= v_left + PRECISION || v_left_third <= v_right + PRECISION) &&
                (v_right_third <= v_left + PRECISION || v_right_third <= v_right + PRECISION);
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
    assert(fabs(v) <= PRECISION && "Checking that the touching line has been found.");
    if (fabs(v) > PRECISION) {
        abort();
    }

    return x;
}

void get_optimal_curve_bound(double* k, double *b,
                             double x_lb, double x_ub,
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
        *b = y_lb - *k * x_lb;
        return;
    }
    if (x_ub <= 0 || 0 <= x_lb) {
        // Function is either convex or concave, thus any tangent line is a valid bound.
        // We are looking for the one minimizing the area.
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

    // Finding the x_star. Condition that line goes through convex x edge and touches S curve.
    double x_star = find_x_star(x_concave, x_convex, is_sigm);

    // The sound region that contains optimal line is [x_star, x_concave].
    // Find the best line there with respect to minimizing area on [x_lb, x_ub].
    find_best_tangent_line(k, b, x_star, x_concave, x_lb, x_ub, is_upper, is_sigm);
}


// It x is greater - it should give upper bound. If x is smaller - lower bound.
void new_S_curve_tang_bound(double& k, double& b, double x, bool is_sigm, bool upper) {
//    check_round();
    ASRTF(x != 0, "x should equal zero.");

    double y = is_sigm ? sigm(x) : tanh(x);
    if (upper) {
        y += SMALL;
    } else {
        y -= SMALL;
    }
    k = is_sigm ? y * (1 - y) : (1 - y * y);
    k = max(0.0, k - SMALL);
    b = y - k * x;
    if (upper) {
        b += SMALL;
    } else {
        b -= SMALL;
    }
}


void new_S_curve_chord_bound(double& k, double& b, double x_lb, double x_ub, bool is_sigm, bool upper) {
    ASRTF(x_ub - x_lb >= 0.0001, "x_lb and x_ub are too close.");

    double y_lb, y_ub;
    if (is_sigm) {
        y_lb = sigm(x_lb);
        y_ub = sigm(x_ub);
    } else {
        y_lb = tanh(x_lb);
        y_ub = tanh(x_ub);
    }

    if (upper) {
        y_lb += SMALL;
        y_ub += SMALL;
    } else {
        y_lb -= SMALL;
        y_ub -= SMALL;
    }

    k = (y_ub - y_lb) / (x_ub - x_lb);
    if (upper) {
        k += SMALL;
    } else {
        k -= SMALL;
    }

    b = y_lb - k * x_lb;
    if (upper) {
        b += SMALL;
    } else {
        b -= SMALL;
    }
}


// Values beyond which I just take a simple bound due to potential numerical stability issues.
const double TANH_LIM = 3;
const double SIGM_LIM = 5;

void compute_S_curve_bounds(double x_lb, double x_ub, bool is_sigm,
                            double* k_lb, double* b_lb,
                            double* k_ub, double* b_ub) {
    assert(x_lb <= x_ub && "x_lb <= x_ub");
    if (x_lb > x_ub) {
        abort();
    }

    double limit = is_sigm ? SIGM_LIM : TANH_LIM;

    if (x_ub - x_lb <= 1.0E-3 || x_lb >= limit || x_ub <= -limit || (x_lb <= -limit && limit <= x_ub)) {
        if (is_sigm) {
            *b_lb = sigm(x_lb);
            *b_ub = sigm(x_ub);
        } else {
            *b_lb = tanh(x_lb);
            *b_ub = tanh(x_ub);
        }
        *k_lb = 0;
        *k_ub = 0;
    }  else if (x_lb >= 0) {
        new_S_curve_chord_bound(*k_lb, *b_lb, x_lb, x_ub, is_sigm, false);
        new_S_curve_tang_bound(*k_ub, *b_ub, x_ub, is_sigm, true);
    } else if (x_ub <= 0) {
        new_S_curve_tang_bound(*k_lb, *b_lb, x_lb, is_sigm, false);
        new_S_curve_chord_bound(*k_ub, *b_ub, x_lb, x_ub, is_sigm, true);
    } else {
        if (x_lb <= -limit) {
            *k_lb = 0;
            *b_lb = is_sigm ? sigm(x_lb) : tanh(x_lb);
        } else if (abs(x_lb) >= abs(x_ub)) {
            new_S_curve_tang_bound(*k_lb, *b_lb, x_lb, is_sigm, false);
        } else {
            double k1, b1;
            double k2, b2;
            new_S_curve_chord_bound(k1, b1, x_lb, x_ub, is_sigm, false);
            new_S_curve_tang_bound(k2, b2, x_lb, is_sigm, false);

            if (k1 <= k2) {
                *k_lb = k1;
                *b_lb = b1;
            } else {
                *k_lb = k2;
                *b_lb = b2;
            }
        }

        if (x_ub >= limit) {
            *k_ub = 0;
            *b_ub = is_sigm ? sigm(x_ub) : tanh(x_ub);
        } else if (abs(x_ub) >= abs(x_lb)) {
            new_S_curve_tang_bound(*k_ub, *b_ub, x_ub, is_sigm, true);
        }
        else {
            double k1, b1;
            double k2, b2;
            new_S_curve_chord_bound(k1, b1, x_lb, x_ub, is_sigm, true);
            new_S_curve_tang_bound(k2, b2, x_ub, is_sigm, true);

            if (k1 <= k2) {
                *k_ub = k1;
                *b_ub = b1;
            } else {
                *k_ub = k2;
                *b_ub = b2;
            }
        }
    }


//    }
//        {
//        if (x_lb <= -limit) {
//            *k_lb = 0;
//            *b_lb = is_sigm ? sigm(x_lb) : tanh(x_lb);
//        } else {
//            get_optimal_curve_bound(k_lb, b_lb, x_lb, x_ub, false, is_sigm);
//        }
//        if (x_ub >= limit) {
//            *k_ub = 0;
//            *b_ub = is_sigm ? sigm(x_ub) : tanh(x_ub);
//        } else {
//            get_optimal_curve_bound(k_ub, b_ub, x_lb, x_ub, true, is_sigm);
//        }
//    }

    // Adjusting for numerical soundness (with big safe margin).
    *b_lb -= SMALL;
    *b_ub += SMALL;

    double y_lb, y_ub;
    if (is_sigm) {
        y_lb = sigm(x_lb);
        y_ub = sigm(x_ub);
    } else {
        y_lb = tanh(x_lb);
        y_ub = tanh(x_ub);
    }

    double lb1 = *k_lb * x_lb + *b_lb;
    double ub1 = *k_ub * x_lb + *b_ub;
    double lb2 = *k_lb * x_ub + *b_lb;
    double ub2 = *k_ub * x_ub + *b_ub;

    double rel = lb1 - y_lb;
    rel = max(rel, y_lb - ub1);
    rel = max(rel, lb2 - y_ub);
    rel = max(rel, y_ub - ub2);
    if (rel > 0) {
        cout << "rel " << rel << " x_lb " << x_lb << " x_ub " << x_ub << endl;
        *b_lb -= 2 * rel;
        *b_ub += 2 * rel;
    }
}

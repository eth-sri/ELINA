#include "asrt.h"
#include <math.h>
#include <cfenv>
#include <cassert>
#include "fconv.h"
#include <stdio.h>
#include <iostream>
#include <limits>
#include <iomanip>

constexpr int DIGITS = numeric_limits<double>::digits10;

void check_round() {
    fesetround(FE_UPWARD);
    ASRTF(fegetround() == FE_UPWARD, "Making sure that setting floating rounding mode works correctly.");
    fesetround(FE_DOWNWARD);
    ASRTF(fegetround() == FE_DOWNWARD, "Making sure that setting floating rounding mode works correctly.");
}

inline void round(bool up) {
    fesetround(up ? FE_UPWARD : FE_DOWNWARD);
}

inline void downward() {
    fesetround(FE_DOWNWARD);
}

inline void upward() {
    fesetround(FE_DOWNWARD);
}

inline double sigm(double x, bool up) {
    round(!up);
    double denom = 1 + exp(-x);
    round(up);
    return 1 / denom;
}

inline double tanh_my(double x, bool up) {
    round(up);
    return tanh(x);
}

// y belongs to [0, 1]
double sigm_tang_k(double y) {
    assert(0 <= y && y <= 1 && "y should be within range.");
    downward();
    return y * (1 - y);
}

// y belongs to [-1, 1]
double tanh_tang_k(double y) {
    assert(-1 <= y && y <= 1 && "y should be within range.");
    upward();
    double temp = y * y;
    downward();
    return 1 - temp;
}

// It x is greater - it should give upper bound. If x is smaller - lower bound.
void S_curve_tang_bound(double* k, double* b, double x, bool is_sigm) {
    check_round();
    ASRTF(x != 0, "x should equal zero.");
    bool up = x > 0;

    double y;
    if (is_sigm) {
        y = sigm(x, up);
        *k = sigm_tang_k(y);
    } else {
        y = tanh_my(x, up);
        *k = tanh_tang_k(y);
    }
    *k = max(0.0, *k - 0.00001);

    round(!up);
    double temp = (*k) * x;
    round(up);

    *b = y - temp;
    if (up) {
        *b += 0.00001;
    } else {
        *b -= 0.00001;
    }

    round(up);
    double correctness = *k * x + *b;

    if (up && correctness < y) {
        cout << setprecision(DIGITS) << "tang doing box for " << x << endl;
        *b = y + 0.00001;
        *k = 0;
    } else if (!up && correctness > y) {
        cout << setprecision(DIGITS) << "tang doing box for " << x << endl;
        *b = y - 0.00001;
        *k = 0;
    }

    fesetround(FE_TONEAREST);
}

void S_curve_chord_bound(double* k, double* b, double x_lb, double x_ub, bool is_sigm) {
    check_round();
    ASRTF(x_lb < x_ub, "x_lb < x_ub");
    ASRTF(x_ub <= 0 || 0 <= x_lb, "x_lb x_ub shouldn't split zero");
    bool up = x_lb < 0;
    double y_lb, y_ub;
    if (is_sigm) {
        y_lb = sigm(x_lb, up);
        y_ub = sigm(x_ub, up);
    } else {
        y_lb = tanh_my(x_lb, up);
        y_ub = tanh_my(x_ub, up);
    }
    if (x_ub - x_lb <= 0.00001) {
        *k = 0;
        if (up) {
            *b = y_ub;
            *b += 0.00001;
        } else {
            *b = y_lb;
            *b -= 0.00001;
        }
        return;
    }
    round(!up);
    double den = x_ub - x_lb;
    round(up);
    *k = (y_ub - y_lb) / den;
    if (up) {
        *k += 0.00001;
    } else {
        *k -= 0.00001;
    }

    round(!up);
    double temp = (*k) * x_lb;
    round(up);
    *b = y_lb - temp;

    if (up) {
        *b += 0.00001;
    } else {
        *b -= 0.00001;
    }

    round(up);
    double c1 = *k * x_lb + *b;
    double c2 = *k * x_ub + *b;
    if (up && (c1 < y_lb || c2 < y_ub)) {
        cout << setprecision(DIGITS) << "chord doing box for " << x_lb << " " << x_ub << endl;
        *b = y_ub + 0.00001;
        *k = 0;
    } else if (!up && (c1 > y_lb || c2 > y_ub)) {
        cout << setprecision(DIGITS) << "chord doing box for " << x_lb << " " << x_ub << endl;
        *b = y_lb - 0.00001;
        *k = 0;
    }

    fesetround(FE_TONEAREST);
}

void compute_curve_bounds(double x_bound, bool is_sigm,
                          double& k_lb, double& b_lb, double& k_ub, double& b_ub) {
    ASRTF(x_bound != 0, "x_bound cannot be zero.");

    check_round();
    fesetround(FE_UPWARD);
    ASRTF(fegetround() == FE_UPWARD, "Making sure that setting floating rounding mode works correctly.");

    if (x_bound > 0) {
        S_curve_chord_bound(&k_lb, &b_lb, 0, x_bound, is_sigm);
        S_curve_tang_bound(&k_ub, &b_ub, x_bound, is_sigm);
    } else {
        S_curve_chord_bound(&k_ub, &b_ub, x_bound, 0, is_sigm);
        S_curve_tang_bound(&k_lb, &b_lb, x_bound, is_sigm);
    }

    fesetround(FE_TONEAREST);
//
//    b_lb -= 1.0;
//    b_ub += 1.0;
}

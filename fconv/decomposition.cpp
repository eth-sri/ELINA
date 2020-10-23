#include <map>
#include <cassert>
#include <math.h>
#include "pdd.h"
#include "utils.h"
#include "fp_mat.h"
#include "S_curve.h"

using namespace std;

void lift_to_relu_y_branch(const int xi, PDD& pdd_dual, const Polarity polarity) {
    const int dim = pdd_dual.dim;
    pdd_dual.dim++;

    // Intentionally swapped because PDD here is given in the dual view.
    vector<double*>& V = pdd_dual.H;
    vector<double*>& H = pdd_dual.V;

    if (V.empty()) {
        // There are no vertices - the further convex hull with an empty vertex set
        // would not change the vertex set, thus I'm leaving it empty.
        assert(H.empty() && "Consistency checking that H is empty.");
        return;
    }
    vector<set_t>& incidence = pdd_dual.incidence;

    for (size_t i = 0; i < V.size(); i++) {
        V[i] = fp_arr_resize(dim + 1, dim, V[i]);
        if (polarity == PLUS) {
            V[i][dim] = V[i][xi + 1];
        }
    }
    for (size_t i = 0; i < H.size(); i++) {
        H[i] = fp_arr_resize(dim + 1, dim, H[i]);
    }
    H.resize(H.size() + 2);
    H[H.size() - 2] = fp_arr_create(dim + 1);
    H[H.size() - 1] = fp_arr_create(dim + 1);

    H[H.size() - 2][dim] = 1;
    H[H.size() - 1][dim] = -1;
    if (polarity == PLUS) {
        H[H.size() - 2][xi + 1] = -1;
        H[H.size() - 1][xi + 1] = 1;
    }

    // Incidence size equals previous size of H which is now increased by 2.
    incidence.resize(H.size());
    incidence[incidence.size() - 2] = set_create(V.size());
    incidence[incidence.size() - 1] = set_create(V.size());
    set_enable_all(incidence[incidence.size() - 2]);
    set_enable_all(incidence[incidence.size() - 1]);
}

void lift_to_tasi_y_branch(const int xi,
                              PDD& pdd_dual,
                              const double x_bound,
                              Activation activation) {
    ASRTF(x_bound != 0, "x_bound has to be either negative or positive.");
    const int dim = pdd_dual.dim;
    pdd_dual.dim++;

    // Intentionally swapped because PDD here is given in the dual view.
    vector<double*>& V = pdd_dual.H;
    vector<double*>& H = pdd_dual.V;

    if (V.empty()) {
        // There are no vertices - the further convex hull with an empty vertex set
        // would not change the vertex set, thus I'm leaving it empty.
        assert(H.empty() && "Consistency checking that H is empty.");
        return;
    }
    vector<set_t>& incidence = pdd_dual.incidence;

    double k_lb, b_lb, k_ub, b_ub;
    compute_curve_bounds(x_bound, activation==Sigm, k_lb, b_lb, k_ub, b_ub);

    vector<double*> V_new;
    V_new.reserve(V.size() * 2);
    vector<int> map_lb(V.size());
    vector<int> map_ub(V.size());

    for (size_t i = 0; i < V.size(); i++) {
        double* v_lb = fp_arr_resize(dim + 1, dim, V[i]);
        double x_cur = v_lb[xi + 1];
        double lb = k_lb * x_cur + b_lb;
        double ub = k_ub * x_cur + b_ub;
        ASRTF(lb <= ub, "Unsoundness detected.");
        v_lb[dim] = lb;
        map_lb[i] = V_new.size();
        V_new.push_back(v_lb);
        if (lb == ub) {
            map_ub[i] = map_lb[i];
        } else {
            double* v_ub = fp_arr_copy(dim + 1, v_lb);
            v_ub[dim] = ub;
            map_ub[i] = V_new.size();
            V_new.push_back(v_ub);
        }
    }

    for (size_t i = 0; i < H.size(); i++) {
        H[i] = fp_arr_resize(dim + 1, dim, H[i]);
    }
    H.resize(H.size() + 2);
    H[H.size() - 2] = fp_arr_create(dim + 1);
    H[H.size() - 1] = fp_arr_create(dim + 1);

    double* h_lb = H[H.size() - 2];
    double* h_ub = H[H.size() - 1];

    // In the minus branch:
    // y >= k * x + b equivalent -b - k * x + y >= 0
    h_lb[0] = -b_lb;
    h_lb[xi + 1] = -k_lb;
    h_lb[dim] = 1;

    // y <= k * x + b equivalent b + k * x - y >= 0
    h_ub[0] = b_ub;
    h_ub[xi + 1] = k_ub;
    h_ub[dim] = -1;

    vector<set_t> incidence_new = set_arr_create(H.size(), V_new.size());

    for (size_t h = 0; h < H.size() - 2; h++) {
        set_t inc = incidence[h];
        set_t inc_new = incidence_new[h];
        for (size_t v = 0; v < V.size(); v++) {
            if (set_test_bit(inc, v)) {
                // Note that they can map to the same vertex, but it's okay.
                set_enable_bit(inc_new, map_lb[v]);
                set_enable_bit(inc_new, map_ub[v]);
            }
        }
    }

    set_t inc_lb = incidence_new[H.size() - 2];
    set_t inc_ub = incidence_new[H.size() - 1];

    for (int v : map_lb) {
        set_enable_bit(inc_lb, v);
    }
    for (int v : map_ub) {
        set_enable_bit(inc_ub, v);
    }

    pdd_dual.H = V_new;
    // Since incidence is reference it is important
    // that I free it _before_ I update incidence in pdd.
    set_arr_free(incidence);
    pdd_dual.incidence = incidence_new;

}

PDD decomposition_recursive(Quadrant& quadrant, const map<Quadrant, PDD>& quadrant2pdd,
                            const int K, Activation activation,
                            const vector<double>& x_lb, const vector<double>& x_ub) {
    const int xi = (int) quadrant.size();
    if (xi == K) {
        PDD pdd = quadrant2pdd.at(quadrant);
        PDD_debug_consistency_check(pdd);
        return pdd;
    } else {
        quadrant.push_back(MINUS);
        PDD pdd_minus = decomposition_recursive(quadrant, quadrant2pdd, K, activation, x_lb, x_ub);
        quadrant.back() = PLUS;
        PDD pdd_plus = decomposition_recursive(quadrant, quadrant2pdd, K, activation, x_lb, x_ub);
        quadrant.pop_back();

        if (activation == Relu) {
            lift_to_relu_y_branch(xi, pdd_minus, MINUS);
            lift_to_relu_y_branch(xi, pdd_plus, PLUS);
        } else {
            lift_to_tasi_y_branch(xi, pdd_minus, x_lb[xi], activation);
            lift_to_tasi_y_branch(xi, pdd_plus, x_ub[xi], activation);
        }

        PDD_debug_consistency_check(pdd_minus);
        PDD_debug_consistency_check(pdd_plus);

        PDD res = PDD_intersect_two_PDDs(pdd_minus, pdd_plus);

        return res;
    }
}

// TODO[gleb] Add support for multiple passes.
vector<double*> decomposition(const int K, const map<Quadrant, PDD>& quadrant2pdd,
                              Activation activation,
                              const vector<double>& x_lb, const vector<double>& x_ub) {
    ASRTF(2 <= K && K <= 5, "Only 2 <= K <= 5 are currently supported.");
    ASRTF((int) quadrant2pdd.size() == POW2[K], "Sanity check - the number of quadrants should be 2^K.");

    Quadrant quadrant {};
    quadrant.reserve(K);
    PDD res = decomposition_recursive(quadrant, quadrant2pdd, K, activation, x_lb, x_ub);
    vector<double*>& H = res.V;
    vector<double*>& V = res.H;

    PDD_adjust_H_for_soundness_finite_polytope(2 * K + 1, H, V);

    // H is given in order (1, x1, ..., xk, yk, ..., y1) thus last k cols have to be reversed.
    // The desired order is (1, x1, ..., xk, y1, ..., yk).
    for (size_t hi = 0; hi < H.size(); hi++) {
        int first = K + 1;
        int last = 2 * K;
        double* h = H[hi];
        while (first < last) {
            swap(h[first], h[last]);
            first++;
            last--;
        }
    }

    fp_mat_free(V);
    set_arr_free(res.incidence);

    return H;
}

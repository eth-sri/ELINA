#include <map>
#include "pdd.h"
#include "utils.h"

using namespace std;

void project_to_relu_y_branch(const int xi, PDD& pdd_dual, const Polarity polarity) {
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
    vector<bset>& incidence = pdd_dual.incidence;

    for (size_t i = 0; i < V.size(); i++) {
        V[i] = (double*) realloc(V[i], sizeof(double) * (dim + 1));
        if (polarity == MINUS) {
            V[i][dim] = 0;
        } else {
            V[i][dim] = V[i][xi + 1];
        }
    }
    for (size_t i = 0; i < H.size(); i++) {
        H[i] = (double*) realloc(H[i], sizeof(double) * (dim + 1));
        H[i][dim] = 0;
    }
    H.resize(H.size() + 2);
    H[H.size() - 2] = (double*) calloc(dim + 1, sizeof(double));
    H[H.size() - 1] = (double*) calloc(dim + 1, sizeof(double));

    H[H.size() - 2][dim] = 1;
    H[H.size() - 1][dim] = -1;
    if (polarity == PLUS) {
        H[H.size() - 2][xi + 1] = -1;
        H[H.size() - 1][xi + 1] = 1;
    }

    incidence.reserve(V.size() + 2);
    incidence.emplace_back(V.size());
    incidence.back().set();
    incidence.emplace_back(V.size());
    incidence.back().set();
}

PDD decomposition_recursive(Quadrant& quadrant, const map<Quadrant, PDD>& quadrant2pdd, const int K) {
    if ((int) quadrant.size() == K) {
        return quadrant2pdd.at(quadrant);
    } else {
        quadrant.push_back(MINUS);
        PDD pdd_minus = decomposition_recursive(quadrant, quadrant2pdd, K);
        quadrant.back() = PLUS;
        PDD pdd_plus = decomposition_recursive(quadrant, quadrant2pdd, K);
        quadrant.pop_back();

        project_to_relu_y_branch(quadrant.size(), pdd_minus, MINUS);
        project_to_relu_y_branch(quadrant.size(), pdd_plus, PLUS);
        PDD res = PDD_intersect_two_PDDs(pdd_minus, pdd_plus);

        return res;
    }
}

// TODO[gleb] Add support for multiple passes.
vector<double*> decomposition(const int K, const map<Quadrant, PDD>& quadrant2pdd) {
    ASRTF(2 <= K && K <= 5, "Only 2 <= K <= 5 are currently supported.");
    ASRTF((int) quadrant2pdd.size() == POW2[K], "Sanity check - the number of quadrants should be 2^K.");

    Quadrant quadrant {};
    quadrant.reserve(K);
    PDD res = decomposition_recursive(quadrant, quadrant2pdd, K);
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

    free_mat(V);

    return H;
}

#include <map>
#include "pdd.h"
#include "utils.h"

using namespace std;

void project_to_relu_y_branch(const int xi, PDD& pdd_dual, Polarity polarity) {
    const int dim = pdd_dual.dim;
    pdd_dual.dim += 1;

    // Intentionally swapped because PDD here is given in the dual view.
    MatrixXd& V = pdd_dual.H;
    MatrixXd& H = pdd_dual.V;
    vector<bset>& incidence = pdd_dual.incidence;

    H.conservativeResize(H.rows() + 2, dim + 1);
    H.bottomRows(2) = MatrixXd::Zero(2, dim + 1);
    H.rightCols(1) = MatrixXd::Zero(H.rows(), 1);
    H(H.rows() - 2, dim) = 1;
    H(H.rows() - 1, dim) = -1;
    if (polarity == PLUS) {
        H(H.rows() - 2, xi + 1) = -1;
        H(H.rows() - 1, xi + 1) = 1;
    }

    V.conservativeResize(NoChange, dim + 1);
    if (polarity == MINUS) {
        V.rightCols(1) = MatrixXd::Zero(V.rows(), 1);
    } else {
        V.rightCols(1) = V.col(xi + 1);
    }

    incidence.emplace_back(V.rows());
    incidence.back().set();
    incidence.emplace_back(V.rows());
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

        return PDD_intersect_two_PDDs(pdd_minus, pdd_plus);
    }
}

// TODO[gleb] Add support for multiple passes.
MatrixXd decomposition(const map<Quadrant, PDD>& quadrant2pdd, const int K) {
    ASRTF(2 <= K && K <= 5, "Only 2 <= K <= 5 are currently supported.");
    ASRTF((int) quadrant2pdd.size() == POW2[K], "Sanity check - the number of quadrants should be 2^K.");

    Quadrant quadrant {};
    quadrant.reserve(K);
    PDD res = decomposition_recursive(quadrant, quadrant2pdd, K);
    MatrixXd& H = res.V;
    MatrixXd& V = res.H;

    PDD_adjust_H_for_soundness_finite_polytope(H, V);

    // H is given in order (1, x1, ..., xk, yk, ..., y1) thus last k rows have to be reversed.
    // For some reason if I don't do it through temporary variable it doesn't work.
    MatrixXd H_rightCols_reversed = H.rightCols(K).rowwise().reverse();
    H.rightCols(K) = H_rightCols_reversed;

    return H;
}

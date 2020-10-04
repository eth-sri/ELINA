#include <map>
#include <cassert>
#include <cstring>
#include "fkrelu.h"
#include "octahedron.h"
#include "split_in_quadrants.h"
#include "decomposition.h"
#include "pdd.h"
#include "mpq.h"
#include "utils.h"
#include "fp_mat.h"

using namespace std;

// Computation of relaxation for 1-relu easily done with analytical formula.
vector<double*> relu_1(double lb, double ub) {
    ASRTF(lb <= ub, "Unsoundness - lower bound should be <= then upper bound.");
    ASRTF(lb < 0 && 0 < ub, "Expecting non-trivial input where lb < 0 < ub.");

    double lmd = -lb * ub / (ub - lb);
    double mu = ub / (ub - lb);
    assert(lmd > 0 && "Expected lmd > 0.");
    assert(mu > 0 && "Expected mu > 0.");

    vector<double*> res = fp_mat_create(3, 3);

    // y >= 0
    res[0][0] = 0;
    res[0][1] = 0;
    res[0][2] = 1;

    // y >= x
    res[1][0] = 0;
    res[1][1] = -1;
    res[1][2] = 1;

    // y <= mu * x + lmd;
    res[2][0] = lmd;
    res[2][1] = mu;
    res[2][2] = -1;

    return res;
}

void verify_fkrelu_input(const int K, const vector<double*>& A) {
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
    ASRTF((int) A.size() == POW3[K] - 1, "Unexpected number of rows in the input.");
    const vector<vector<int>>& coefs = K2OCTAHEDRON_COEFS[K];
    for (int i = 0; i < (int) A.size(); i++) {
        for (int j = 0; j < K; j++) {
            ASRTF(A[i][j + 1] == coefs[i][j], "Input is not of correct format.");
        }
    }
}

vector<double*> fkrelu(const int K, const vector<double*>& A) {
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
    verify_fkrelu_input(K, A);
    if (K == 1) {
        return relu_1(-A[0][0], A[1][0]);
    }
    OctahedronV oct = compute_octahedron_V(K, A);

    // Split in quadrants takes care of memory management of input vertices.
    map<Quadrant, QuadrantInfo> quadrant2info = split_in_quadrants(oct.V,
                                                                   oct.incidence,
                                                                   oct.orthant_adjacencies,
                                                                   K);

    map<Quadrant, PDD> quadrant2pdd;
    for (auto& pair : quadrant2info) {
        const Quadrant& quadrant = pair.first;
        vector<mpq_t*>& V_mpq = pair.second.V;

        if (V_mpq.empty()) {
            // Input in the quadrant is the empty set.
            quadrant2pdd[quadrant] = {K + 1, {}, {}, {}};
            continue;
        }

        vector<double*> V = mpq_mat_to_fp(K + 1, V_mpq);
        mpq_mat_free(K + 1, V_mpq);

        const vector<set_t>& incidence_V_to_H = pair.second.V_to_H_incidence;
        assert(
                incidence_V_to_H.size() == V.size() &&
                "Incidence_V_to_H.size() should equal V.size()");
        vector<set_t> incidence_H_to_V_with_redundancy = set_arr_transpose(incidence_V_to_H);
        set_arr_free(incidence_V_to_H);
        assert(
                incidence_H_to_V_with_redundancy.size() == A.size() + K &&
                "Incidence_H_to_V_with_redundancy.size() should equal A.size() + K");
        vector<int> maximal_H = compute_maximal_indexes(incidence_H_to_V_with_redundancy);
        set_t is_maximal = set_create(incidence_H_to_V_with_redundancy.size());
        for (auto i : maximal_H) {
            set_enable_bit(is_maximal, i);
        }

        vector<double*> H(maximal_H.size());
        vector<set_t> incidence_H_to_V(maximal_H.size());

        int count = 0;
        for (size_t i = 0; i < incidence_H_to_V_with_redundancy.size(); i++) {
            if (!set_test_bit(is_maximal, i)) {
                set_free(incidence_H_to_V_with_redundancy[i]);
                continue;
            }
            double* h = (double*) calloc(K + 1, sizeof(double));
            H[count] = h;
            incidence_H_to_V[count] = incidence_H_to_V_with_redundancy[i];
            count++;
            if (i < A.size()) {
                memcpy(h, A[i], sizeof(double) * (K + 1));
            } else {
                int xi = i - (int) A.size();
                assert(0 <= xi && xi < K && "Sanity checking the range of xi.");
                if (quadrant[xi] == MINUS) {
                    h[xi + 1] = -1;
                }  else {
                    h[xi + 1] = 1;
                }
            }
        }
        assert(count == (int) maximal_H.size() && "count should equal maximal_H.size()");
        set_free(is_maximal);
        quadrant2pdd[quadrant] = {K + 1, V, H, incidence_H_to_V};
    }
    return decomposition(K, quadrant2pdd);
}

vector<double*> krelu_with_cdd(const int K, const vector<double*>& A) {
    // No need to verify since CDD can work with input in any format.
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
    map<Quadrant, QuadrantInfo> quadrant2info = compute_quadrants_with_cdd(K, A);

    size_t num_vertices = 0;
    for (auto& entry : quadrant2info) {
        num_vertices += entry.second.V.size();
    }

    dd_MatrixPtr vertices = dd_CreateMatrix(num_vertices, 2 * K + 1);
    vertices->representation = dd_Generator;
    size_t counter = 0;

    for (auto& entry : quadrant2info) {
        const auto& quadrant = entry.first;
        auto& V_quadrant = entry.second.V;

        for (const auto& v : V_quadrant) {
            mpq_t* v_projection = vertices->matrix[counter];
            counter++;
            for (int i = 0; i < K + 1; i++) {
                mpq_set(v_projection[i], v[i]);
            }
            for (int i = 0; i < K; i++) {
                if (quadrant[i] == PLUS) {
                    // Only need to set for the case of PLUS,
                    // because in case of MINUS there should be 0.
                    // And it is already there automatically.
                    mpq_set(v_projection[1 + i + K], v[1 + i]);
                }
            }
        }
        mpq_mat_free(K + 1, V_quadrant);
        set_arr_free(entry.second.V_to_H_incidence);
    }
    assert(counter == num_vertices && "Consistency checking that counter equals the number of vertices.");

    vector<double*> H = cdd_compute_inequalities_from_vertices(vertices);
    dd_FreeMatrix(vertices);

    return H;
}

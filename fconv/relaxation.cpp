#include <map>
#include <cassert>
#include <cstring>
#include <cmath>
#include "relaxation.h"
#include "octahedron.h"
#include "quadrants.h"
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

// Max pool in case of 1-variable is simply adding y = x.
vector<double*> pool_1(const vector<double*>& A) {
    vector<double*> res = fp_mat_create(A.size() + 2, 3);

    for (size_t i = 0; i < A.size(); i++) {
        res[i][0] = A[i][0];
        res[i][1] = A[i][1];
    }

    res[A.size()][1] = 1;
    res[A.size()][2] = -1;
    res[A.size() + 1][1] = -1;
    res[A.size() + 1][2] = 1;

    return res;
}

void verify_that_octahedron_and_all_xi_split_zero(const int K, const vector<double*>& A) {
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
    ASRTF((int) A.size() == POW3[K] - 1, "Unexpected number of rows in the input.");
    const vector<vector<int>>& coefs = K2OCTAHEDRON_COEFS[K];
    for (int i = 0; i < (int) A.size(); i++) {
        for (int j = 0; j < K; j++) {
            ASRTF(A[i][j + 1] == coefs[i][j], "Input is not of correct format.");
        }
    }
    for (int xi = 0; xi < K; xi++) {
        double lb = -A[LOWER_BOUND_INDEX[K][xi]][0];
        double ub = A[UPPER_BOUND_INDEX[K][xi]][0];
        ASRTF(lb < 0, "Lower bound should be negative.");
        ASRTF(0 < ub, "Upper bound should be positive.");
    }
}

vector<double*> cdd_compute_inequalities_from_vertices(dd_MatrixPtr vertices) {
    vertices->representation = dd_Generator;
    dd_PolyhedraPtr poly = cdd_Matrix_to_Poly(vertices);

    dd_MatrixPtr inequalities = dd_CopyInequalities(poly);
    // Note that linearities are quite rare.
    set_type linearities = inequalities->linset;

    const int num_ineq = inequalities->rowsize;
    const int dim = inequalities->colsize;
    vector<double*> H = fp_mat_create(num_ineq + set_card(linearities), dim);
    int counter = 0;
    for (int ineq = 0; ineq < num_ineq; ineq++) {
        for (int j = 0; j < dim; j++) {
            H[counter][j] = mpq_get_d(inequalities->matrix[ineq][j]);
        }
        counter++;
        if (set_member(ineq + 1, linearities)) {
            for (int j = 0; j < dim; j++) {
                H[counter][j] = -H[counter - 1][j];
            }
            counter++;
        }
    }
    assert(
            counter == num_ineq + set_card(linearities) &&
            "Counter should equal num_ineq + number of linearities");

    for (auto h : H) {
        double abs_coef = 0;
        for (int i = 0; i < dim; i++) {
            abs_coef = max(fabs(h[i]), abs_coef);
        }
        ASRTF(abs_coef != 0, "Inequality cannot consist fully of zeros.");
        for (int i = 0; i < dim; i++) {
            h[i] /= abs_coef;
        }
    }

    dd_FreePolyhedra(poly);
    dd_FreeMatrix(inequalities);
    return H;
}

vector<double*> fast_relaxation_through_decomposition(const int K,
                                                      const vector<double*>& A,
                                                      Activation activation) {
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
    ASRTF(activation == Relu || activation == Tanh || activation == Sigm,
          "Activation should be Relu, Tanh or Sigm.");
//    cout << "the input is" << endl;
//    fp_mat_print(K + 1, A);
    verify_that_octahedron_and_all_xi_split_zero(K, A);
    if (K == 1 && activation == Relu) {
        return relu_1(-A[0][0], A[1][0]);
    }
    if (K == 1) {
        return ktasi_with_cdd(K, A, activation);
    }
    OctahedronV oct = get_octahedron_V(K, A);

    // Split in quadrants takes care of memory management of input vertices.
    map<Quadrant, VInc_mpq> quadrant2vinc = split_in_quadrants(oct.V,
                                                               oct.incidence,
                                                               oct.orthant_adjacencies,
                                                               K);

    map<Quadrant, PDD> quadrant2pdd;
    for (auto& entry : quadrant2vinc) {
        const Quadrant& quadrant = entry.first;
        vector<mpq_t*>& V_mpq = entry.second.V;

        if (V_mpq.empty()) {
            // Input in the quadrant is the empty set.
            quadrant2pdd[quadrant] = {K + 1, {}, {}, {}};
            continue;
        }

        vector<double*> V = mpq_mat_to_fp(K + 1, V_mpq);
        mpq_mat_free(K + 1, V_mpq);

        const vector<set_t>& incidence_V_to_H = entry.second.V_to_H_incidence;
        assert(
                incidence_V_to_H.size() == V.size() &&
                "Incidence_V_to_H.size() should equal V.size()");
        vector<set_t> incidence_H_to_V = set_arr_transpose(incidence_V_to_H);
        set_arr_free(incidence_V_to_H);
        assert(
                incidence_H_to_V.size() == A.size() + K &&
                "Incidence_H_to_V.size() should equal A.size() + K");
        vector<int> irredund_idx = compute_maximal_indexes(incidence_H_to_V);
        set_t is_irredund = set_create(incidence_H_to_V.size());
        for (auto i : irredund_idx) {
            set_enable_bit(is_irredund, i);
        }

        vector<double*> H_irredund = fp_mat_create(irredund_idx.size(), K + 1);
        vector<set_t> incidence_irredund(irredund_idx.size());

        size_t count = 0;
        for (size_t i = 0; i < incidence_H_to_V.size(); i++) {
            if (!set_test_bit(is_irredund, i)) {
                set_free(incidence_H_to_V[i]);
                continue;
            }
            double* h = H_irredund[count];
            incidence_irredund[count] = incidence_H_to_V[i];
            count++;
            if (i < A.size()) {
                fp_arr_set(K + 1, h, A[i]);
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
        assert(count == irredund_idx.size() && "count should equal maximal_H.size()");
        set_free(is_irredund);
        quadrant2pdd[quadrant] = {K + 1, V, H_irredund, incidence_irredund};
    }

    // Lower and upper bounds are needed for decomposition of tanh and sigmoid functions.
    vector<double> x_lb(K);
    vector<double> x_ub(K);
    for (int xi = 0; xi < K; xi++) {
        x_lb[xi] = -A[LOWER_BOUND_INDEX[K][xi]][0];
        x_ub[xi] = A[UPPER_BOUND_INDEX[K][xi]][0];
    }

    vector<double*> res = decomposition(K, quadrant2pdd, activation, x_lb, x_ub);

    return res;
}

vector<double*> krelu_with_cdd(const int K, const vector<double*>& A) {
    // No need to verify since CDD can work with input in any format.
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
    map<Quadrant, VInc_mpq> quadrant2vinc = get_quadrants_cdd(K, A);

    size_t num_vertices = 0;
    for (auto& entry : quadrant2vinc) {
        num_vertices += entry.second.V.size();
    }

    dd_MatrixPtr vertices = dd_CreateMatrix(num_vertices, 2 * K + 1);
    size_t counter = 0;

    for (auto& entry : quadrant2vinc) {
        const auto& quadrant = entry.first;
        auto& V_quadrant = entry.second.V;

        for (const auto& v : V_quadrant) {
            mpq_t* v_projection = vertices->matrix[counter];
            counter++;
            mpq_arr_set(K + 1, v_projection, v);
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
    assert(counter == num_vertices && "Counter should equal the number of vertices.");

    vector<double*> H = cdd_compute_inequalities_from_vertices(vertices);
    dd_FreeMatrix(vertices);

    return H;
}

vector<double*> fkpool(const int K, const vector<double*>& A) {
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
    if (K == 1) {
        return pool_1(A);
    }

    vector<DD_mpq> quadrants = get_pool_quadrants(K, A);
    assert((int) quadrants.size() == K && "Number of quadrants should be K.");

    vector<PDD> quadrants_pdds;
    for (int xi = 0; xi < K; xi++) {
        vector<mpq_t*>& V_mpq = quadrants[xi].V;
        vector<double*>& H = quadrants[xi].H;

        if (V_mpq.empty()) {
            quadrants_pdds.push_back({K + 2, {}, {}, {}});
            continue;
        }

        vector<double*> V = mpq_mat_to_fp(K + 1, V_mpq);
        mpq_mat_free(K + 1, V_mpq);
        for (size_t i = 0; i < V.size(); i++) {
            V[i] = fp_arr_resize(K + 2, K + 1, V[i]);
            V[i][K + 1] = V[i][xi + 1];
        }

        const vector<set_t>& incidence_V_to_H = quadrants[xi].V_to_H_incidence;
        assert(
                incidence_V_to_H.size() == V.size() &&
                "Incidence_V_to_H.size() should equal V.size()");
        vector<set_t> incidence_H_to_V = set_arr_transpose(incidence_V_to_H);
        set_arr_free(incidence_V_to_H);
        assert(
                incidence_H_to_V.size() == H.size() &&
                "Incidence_H_to_V.size() should equal H.size()");

        vector<int> irredund_idx = compute_maximal_indexes(incidence_H_to_V);
        set_t is_irredund = set_create(incidence_H_to_V.size());
        for (auto i : irredund_idx) {
            set_enable_bit(is_irredund, i);
        }

        // irredund_idx.size() + 2 because there will also be new equality y = xi.
        vector<double*> H_irredund(irredund_idx.size() + 2);
        vector<set_t> incidence_irredund(irredund_idx.size() + 2);

        size_t count = 0;
        for (size_t i = 0; i < incidence_H_to_V.size(); i++) {
            if (!set_test_bit(is_irredund, i)) {
                free(H[i]);
                set_free(incidence_H_to_V[i]);
                continue;
            }
            H_irredund[count] = fp_arr_resize(K + 2, K + 1, H[i]);
            incidence_irredund[count] = incidence_H_to_V[i];
            count++;
        }
        assert(count == irredund_idx.size() && "count should equal irredund_idx.size()");
        set_free(is_irredund);

        // y = xi
        double* h1 = fp_arr_create(K + 2);
        double* h2 = fp_arr_create(K + 2);
        h1[xi + 1] = 1;
        h1[K + 1] = -1;
        h2[xi + 1] = -1;
        h2[K + 1] = 1;
        set_t inc1 = set_create(V.size());
        set_t inc2 = set_create(V.size());
        set_enable_all(inc1);
        set_enable_all(inc2);

        H_irredund[count] = h1;
        H_irredund[count + 1] = h2;
        incidence_irredund[count] = inc1;
        incidence_irredund[count + 1] = inc2;

        quadrants_pdds.push_back({K + 2, V, H_irredund, incidence_irredund});
    }

    PDD res = PDD_intersect_all(quadrants_pdds);
    vector<double*>& H = res.V;
    vector<double*>& V = res.H;
    PDD_adjust_H_for_soundness_finite_polytope(K + 2, H, V);

    fp_mat_free(V);
    set_arr_free(res.incidence);

    return H;
}

vector<double*> kpool_with_cdd(const int K, const vector<double*>& A) {
    // No need to verify since CDD can work with input in any format.
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");

    vector<DD_mpq> quadrants = get_pool_quadrants_cdd(K, A);
    assert((int) quadrants.size() == K && "Number of quadrants should be K.");

    size_t num_vertices = 0;
    for (auto& dd : quadrants) {
        num_vertices += dd.V.size();
    }

    dd_MatrixPtr vertices = dd_CreateMatrix(num_vertices, K + 2);
    size_t counter = 0;

    for (int xi = 0; xi < K; xi++) {
        auto& V_quadrant = quadrants[xi].V;
        for (mpq_t* v : V_quadrant) {
            mpq_t *v_projection = vertices->matrix[counter];
            counter++;
            mpq_arr_set(K + 1, v_projection, v);
            mpq_set(v_projection[K + 1], v[xi + 1]);
        }
        mpq_mat_free(K + 1, V_quadrant);
        set_arr_free(quadrants[xi].V_to_H_incidence);
    }
    assert(counter == num_vertices && "Counter should equal the number of vertices.");

    vector<double*> H = cdd_compute_inequalities_from_vertices(vertices);
    dd_FreeMatrix(vertices);

    return H;
}

vector<double*> ktasi_with_cdd(int K, const vector<double*>& A, Activation activation) {
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
    ASRTF(activation == Tanh || activation == Sigm, "Only Tanh and Sigm are supported.");
    verify_that_octahedron_and_all_xi_split_zero(K, A);

    map<Quadrant, vector<mpq_t*>> quadrant2V = get_tasi_quadrants_cdd(K, A, activation);

    size_t num_vertices = 0;
    for (auto& entry : quadrant2V) {
        num_vertices += entry.second.size();
    }

    dd_MatrixPtr vertices = dd_CreateMatrix(num_vertices, 2 * K + 1);
    size_t counter = 0;
    for (auto& entry : quadrant2V) {
        for (mpq_t* v : entry.second) {
            mpq_arr_set(2 * K + 1, vertices->matrix[counter], v);
            counter++;
        }
        mpq_mat_free(2 * K + 1, entry.second);
    }
    assert(counter == num_vertices && "Counter should equal the number of vertices.");

    vector<double*> H = cdd_compute_inequalities_from_vertices(vertices);
    dd_FreeMatrix(vertices);

    return H;
}

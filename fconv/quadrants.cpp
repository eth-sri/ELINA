#include <algorithm>
#include <cassert>
#include <math.h>
#include <limits>
#include <iomanip>
#include "quadrants.h"
#include "octahedron.h"
#include "fp_mat.h"
#include "mpq.h"
#include "S_curve.h"
#include "S_curve2.h"

constexpr int DIGITS = numeric_limits<double>::digits10;

using namespace std;

map<Quadrant, VInc_mpq> get_quadrants_cdd(const int K, const vector<double*>& A) {
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");

    const vector<Quadrant>& quadrants = K2QUADRANTS[K];

    const int NUM_H = (int) A.size();

    dd_MatrixPtr cdd_A = dd_CreateMatrix(NUM_H + K, K + 1);
    cdd_A->representation = dd_Inequality;

    for (int i = 0; i < NUM_H; i++) {
        mpq_arr_set_d(K + 1, cdd_A->matrix[i], A[i]);
    }

    map<Quadrant, VInc_mpq> quadrant2vinc;
    for (const auto& quadrant : quadrants) {
        for (int xi = 0; xi < K; xi++) {
            mpq_t* row = cdd_A->matrix[NUM_H + xi];
            mpq_arr_set_zero(K + 1, row);
            if (quadrant[xi] == MINUS) {
                mpq_set_si(row[xi + 1], -1, 1);
            } else {
                mpq_set_si(row[xi + 1], 1, 1);
            }
        }
        dd_PolyhedraPtr poly = cdd_Matrix_to_Poly(cdd_A);

        dd_MatrixPtr cdd_V = dd_CopyGenerators(poly);
        dd_SetFamilyPtr cdd_incidence = dd_CopyIncidence(poly);
        const size_t num_v = cdd_V->rowsize;
        ASRTF(
                cdd_incidence->famsize == (int) num_v && \
                cdd_incidence->setsize == NUM_H + K + 1,
                "Sanity checking cdd_incidence size.");

        vector<mpq_t*> V(num_v);
        vector<set_t> incidence(num_v);

        for (size_t v = 0; v < num_v; v++) {
            ASRTF(!mpq_cmp_si(cdd_V->matrix[v][0], 1, 1), "Polytope is not bounded.");
            V[v] = mpq_arr_copy(K + 1, cdd_V->matrix[v]);
            incidence[v] = set_from_cdd(cdd_incidence->set[v]);
        }
        dd_FreePolyhedra(poly);
        dd_FreeMatrix(cdd_V);
        dd_FreeSetFamily(cdd_incidence);

        quadrant2vinc[quadrant] = {K + 1, V, incidence};
    }

    dd_FreeMatrix(cdd_A);
    return quadrant2vinc;
}

vector<DD_mpq> get_pool_quadrants_cdd(const int K, const vector<double*>& A) {
    ASRTF(1 <= K && K <= 5, "K should be within allowed range.");

    vector<DD_mpq> quadrants;
    for (int xi = 0; xi < K; xi++) {
        vector<double*> A_quadrant = fp_mat_copy(K + 1, A);
        A_quadrant.reserve(A.size() + K - 1);
        for (int xj = 0; xj < K; xj++) {
            if (xi == xj) {
                continue;
            }
            double* row = fp_arr_create(K + 1);
            // xi - xj >= 0.
            row[xi + 1] = 1;
            row[xj + 1] = -1;
            A_quadrant.push_back(row);
        }
        assert(A_quadrant.size() == A.size() + K - 1 && "A_quadrant has unexpected size.");

        dd_MatrixPtr cdd_A = dd_CreateMatrix(A_quadrant.size(), K + 1);
        cdd_A->representation = dd_Inequality;
        for (size_t i = 0; i < A_quadrant.size(); i++) {
            mpq_arr_set_d(K + 1, cdd_A->matrix[i], A_quadrant[i]);
        }

        dd_PolyhedraPtr poly = cdd_Matrix_to_Poly(cdd_A);
        dd_MatrixPtr cdd_V = dd_CopyGenerators(poly);
        dd_SetFamilyPtr cdd_incidence = dd_CopyIncidence(poly);
        const size_t num_v = cdd_V->rowsize;
        ASRTF(
                cdd_incidence->famsize == (int) num_v && \
                cdd_incidence->setsize == (int) A.size() + K,
                "cdd_incidence has unexpected size.");

        vector<mpq_t*> V(num_v);
        vector<set_t> incidence(num_v);

        for (size_t v = 0; v < num_v; v++) {
            ASRTF(!mpq_cmp_si(cdd_V->matrix[v][0], 1, 1), "Polytope is not bounded.");
            V[v] = mpq_arr_copy(K + 1, cdd_V->matrix[v]);
            incidence[v] = set_from_cdd(cdd_incidence->set[v]);
        }

        dd_FreePolyhedra(poly);
        dd_FreeMatrix(cdd_A);
        dd_FreeMatrix(cdd_V);
        dd_FreeSetFamily(cdd_incidence);

        quadrants.push_back({K + 1, V, A_quadrant, incidence});
    }

    return quadrants;
}

vector<DD_mpq> get_pool_quadrants(const int K, const vector<double*>& A) {
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");

    vector<DD_mpq> quadrants;
    for (int xi = 0; xi < K; xi++) {
        vector<double*> A_quadrant = fp_mat_copy(K + 1, A);
        for (int xj = 0; xj < K; xj++) {
            if (xi == xj) {
                continue;
            }
            vector<int> constraint_coef(K);
            constraint_coef[xi] = 1;
            constraint_coef[xj] = -1;
            int idx = coef2index(constraint_coef);
            A_quadrant[idx][0] = min(A_quadrant[idx][0], 0.0);
        }
        OctahedronV octahedron = get_octahedron_V(K, A_quadrant);
        quadrants.push_back({K + 1, octahedron.V, A_quadrant, octahedron.incidence});
    }
    return quadrants;
}

struct SegmentCons {
    mpq_t k_lb, b_lb, k_ub, b_ub;

    void init() {
        mpq_inits(k_lb, b_lb, k_ub, b_ub, NULL);
    }

    void free() {
        mpq_clears(k_lb, b_lb, k_ub, b_ub, NULL);
    }

    void lb(const mpq_t& x, mpq_t& y) {
        mpq_mul(y, x, k_lb);
        mpq_add(y, y, b_lb);
    }

    void ub(const mpq_t& x, mpq_t& y) {
        mpq_mul(y, x, k_ub);
        mpq_add(y, y, b_ub);
    }
};

SegmentCons get_tasi_constraints(double x_bound, Activation activation) {
    ASRTF(x_bound != 0, "x_bound should be zero.");

    double k_lb, b_lb, k_ub, b_ub;
    compute_S_curve_bounds(min(x_bound, 0.0),
                           max(x_bound, 0.0),
                           activation==Sigm,
                           &k_lb, &b_lb, &k_ub, &b_ub);
//    compute_curve_bounds(x_bound, activation==Sigm, k_lb, b_lb, k_ub, b_ub);

    SegmentCons sc;
    sc.init();
    mpq_set_d(sc.k_lb, k_lb);
    mpq_set_d(sc.b_lb, b_lb);
    mpq_set_d(sc.k_ub, k_ub);
    mpq_set_d(sc.b_ub, b_ub);

    return sc;
}

SegmentCons** create_all_segment_cons(int K, const vector<double*>& A, Activation activation) {
    SegmentCons** segment_cons = (SegmentCons**) calloc(K, sizeof(SegmentCons*));
    for (int xi = 0; xi < K; xi++) {
        segment_cons[xi] = (SegmentCons*) calloc(2, sizeof(SegmentCons));
        segment_cons[xi][0] = get_tasi_constraints(-A[LOWER_BOUND_INDEX[K][xi]][0], activation);
        segment_cons[xi][1] = get_tasi_constraints(A[UPPER_BOUND_INDEX[K][xi]][0], activation);
    }
    return segment_cons;
}

void free_all_segment_cons(int K, SegmentCons** segment_cons) {
    for (int xi = 0; xi < K; xi++) {
        segment_cons[xi][0].free();
        segment_cons[xi][1].free();
        free(segment_cons[xi]);
    }
    free(segment_cons);
}

map<Quadrant, vector<mpq_t*>> get_tasi_quadrants_cdd(
        const int K,
        const vector<double*>& A,
        Activation activation) {
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");

    const vector<Quadrant>& quadrants = K2QUADRANTS[K];

    const int NUM_H = (int) A.size();
    dd_MatrixPtr cdd_A = dd_CreateMatrix(NUM_H + 3 * K, 2 * K + 1);
    cdd_A->representation = dd_Inequality;

    for (int i = 0; i < NUM_H; i++) {
        mpq_arr_set_d(K + 1, cdd_A->matrix[i], A[i]);
    }

    SegmentCons** segment_cons = create_all_segment_cons(K, A, activation);

    map<Quadrant, vector<mpq_t*>> quadrant2V;
    for (const auto& quadrant : quadrants) {
        for (int xi = 0; xi < K; xi++) {
            mpq_t* row = cdd_A->matrix[NUM_H + xi];
            mpq_arr_set_zero(2 * K + 1, row);
            if (quadrant[xi] == MINUS) {
                mpq_set_si(row[xi + 1], -1, 1);
            } else {
                mpq_set_si(row[xi + 1], 1, 1);
            }
        }
        for (int xi = 0; xi < K; xi++) {
            SegmentCons& sc = segment_cons[xi][quadrant[xi]];
            int base_i = NUM_H + K + 2 * xi;
            mpq_t* row_lb = cdd_A->matrix[base_i];
            mpq_t* row_ub = cdd_A->matrix[base_i + 1];
            mpq_arr_set_zero(2 * K + 1, row_lb);
            mpq_arr_set_zero(2 * K + 1, row_ub);

            // Lower bound y >= kx + b equivalent -b - kx + y >= 0
            mpq_neg(row_lb[0], sc.b_lb);
            mpq_neg(row_lb[xi + 1], sc.k_lb);
            mpq_set_si(row_lb[xi + 1 + K], 1, 1);

            // Upper bound y <= kx + b equivalent b + kx - y >= 0
            mpq_set(row_ub[0], sc.b_ub);
            mpq_set(row_ub[xi + 1], sc.k_ub);
            mpq_set_si(row_ub[xi + 1 + K], -1, 1);
        }

        dd_PolyhedraPtr poly = cdd_Matrix_to_Poly(cdd_A);
        dd_MatrixPtr cdd_V = dd_CopyGenerators(poly);
        const size_t num_v = cdd_V->rowsize;

        vector<mpq_t*> V(num_v);

        for (size_t v = 0; v < num_v; v++) {
            ASRTF(!mpq_cmp_si(cdd_V->matrix[v][0], 1, 1), "Polytope is not bounded.");
            V[v] = mpq_arr_copy(2 * K + 1, cdd_V->matrix[v]);
        }

        dd_FreePolyhedra(poly);
        dd_FreeMatrix(cdd_V);

        quadrant2V[quadrant] = V;
    }

    free_all_segment_cons(K, segment_cons);

    dd_FreeMatrix(cdd_A);

    return quadrant2V;
}

map<Quadrant, vector<mpq_t*>> get_tasi_quadrants_cdd_lift(
        const int K,
        const vector<double*>& A,
        Activation activation) {
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");

    const vector<Quadrant>& quadrants = K2QUADRANTS[K];

    const int NUM_H = (int) A.size();

    dd_MatrixPtr cdd_A = dd_CreateMatrix(NUM_H + K, K + 1);
    cdd_A->representation = dd_Inequality;

    for (int i = 0; i < NUM_H; i++) {
        mpq_arr_set_d(K + 1, cdd_A->matrix[i], A[i]);
    }

    SegmentCons** segment_cons = create_all_segment_cons(K, A, activation);

    mpq_t y_lb, y_ub;
    mpq_inits(y_lb, y_ub, NULL);

    map<Quadrant, vector<mpq_t*>> quadrant2V;
    for (const auto& quadrant : quadrants) {
        for (int xi = 0; xi < K; xi++) {
            mpq_t* row = cdd_A->matrix[NUM_H + xi];
            mpq_arr_set_zero(K + 1, row);
            if (quadrant[xi] == MINUS) {
                mpq_set_si(row[xi + 1], -1, 1);
            } else {
                mpq_set_si(row[xi + 1], 1, 1);
            }
        }
        dd_PolyhedraPtr poly = cdd_Matrix_to_Poly(cdd_A);

        dd_MatrixPtr cdd_V = dd_CopyGenerators(poly);
        const size_t num_v = cdd_V->rowsize;

        vector<mpq_t*> V;
        for (size_t vi = 0; vi < num_v; vi++) {
            const size_t base = V.size();
            mpq_t* v_base = cdd_V->matrix[vi];
            V.push_back(mpq_arr_create(2 * K + 1));
            mpq_arr_set(K + 1, V[base], v_base);
            for (int xi = 0; xi < K; xi++) {
                SegmentCons sc = segment_cons[xi][quadrant[xi]];
                sc.lb(v_base[xi + 1], y_lb);
                sc.ub(v_base[xi + 1], y_ub);
//                cout << setprecision(DIGITS) << "y_lb " << mpq_get_d(y_lb) << " y_ub "
//                << mpq_get_d(y_ub) << endl;
                int cmp = mpq_cmp(y_lb, y_ub);
                ASRTF(cmp <= 0, "Unsoundness detected.");
                if (cmp == 0) {
                    // In case x has extreme value both lower and upper bound of y
                    // have the same value - thus y can be set directly.
                    for (size_t i = base; i < V.size(); i++) {
                        mpq_set(V[i][xi + 1 + K], y_lb);
                    }
                    continue;
                }
                const size_t num = V.size() - base;
                for (size_t i = base; i < base + num; i++) {
                    mpq_t* v_lb = V[i];
                    mpq_t* v_ub = mpq_arr_copy(2 * K + 1, v_lb);
                    V.push_back(v_ub);
                    mpq_set(v_lb[xi + 1 + K], y_lb);
                    mpq_set(v_ub[xi + 1 + K], y_ub);
                }
                assert(V.size() - base == 2 * num &&
                       "The number of new vertices should've doubled.");
            }
        }

        dd_FreePolyhedra(poly);
        dd_FreeMatrix(cdd_V);

        quadrant2V[quadrant] = V;
    }

    mpq_clears(y_lb, y_ub, NULL);
    free_all_segment_cons(K, segment_cons);

    dd_FreeMatrix(cdd_A);

    return quadrant2V;
}

#include <cassert>
#include <string.h>
#include "octahedron.h"
#include "mpq.h"
#include "fkpool.h"
#include "fp_mat.h"
#include "pdd.h"

PDD intersect_all_PDDs_recursively(vector<PDD>& pdds) {
    ASRTF(!pdds.empty(), "Expected non-empty number of PDDs.");
    if (pdds.size() == 1) {
        return pdds[0];
    }
    vector<PDD> next_pdds;
    size_t i = 0;
    while (i < pdds.size()) {
        if (i < pdds.size() - 1) {
            next_pdds.push_back(PDD_intersect_two_PDDs(pdds[i], pdds[i + 1]));
            i += 2;
        } else {
            next_pdds.push_back(pdds[i]);
            i++;
        }
    }
    return intersect_all_PDDs_recursively(next_pdds);
}

// Computation of relaxation for 1-relu easily done with analytical formula.
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

vector<double*> fkpool(const int K, const vector<double*>& A) {
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
    if (K == 1) {
        return pool_1(A);
    }

    // TODO[gleb] Implement specialized fast compute_max_pool_quadrants.
    vector<QuadrantInfo> quadrant_infos = compute_max_pool_quadrants_with_cdd(K, A);
    assert((int) quadrant_infos.size() == K && "Number of quadrants should be K.");

    vector<PDD> quadrant_pdds;
    for (int xi = 0; xi < K; xi++) {
        vector<mpq_t*>& V_mpq = quadrant_infos[xi].V;

        if (V_mpq.empty()) {
            quadrant_pdds.push_back({K + 2, {}, {}, {}});
            continue;
        }

        vector<double*> V = mpq_mat_to_fp(K + 1, V_mpq);
        mpq_mat_free(K + 1, V_mpq);
        for (size_t i = 0; i < V.size(); i++) {
            V[i] = (double*) realloc(V[i], (K + 2) * sizeof(double));
            V[i][K + 1] = V[i][xi + 1];
        }

        const vector<set_t>& incidence_V_to_H = quadrant_infos[xi].V_to_H_incidence;
        assert(
                incidence_V_to_H.size() == V.size() &&
                "Incidence_V_to_H.size() should equal V.size()");
        vector<set_t> incidence_H_to_V_with_redundancy = set_arr_transpose(incidence_V_to_H);
        set_arr_free(incidence_V_to_H);
        assert(
                incidence_H_to_V_with_redundancy.size() == A.size() + K - 1 &&
                "Incidence_H_to_V_with_redundancy.size() should equal A.size() + K - 1");

        vector<int> maximal_H = compute_maximal_indexes(incidence_H_to_V_with_redundancy);
        set_t is_maximal = set_create(incidence_H_to_V_with_redundancy.size());
        for (auto i : maximal_H) {
            set_enable_bit(is_maximal, i);
        }

        // maximal_H.size() + 2 because there will also be new equality y = xi.
        vector<double*> H = fp_mat_create(maximal_H.size() + 2, K + 2);
        vector<set_t> incidence_H_to_V(maximal_H.size() + 2);

        int count = 0;
        for (size_t i = 0; i < incidence_H_to_V_with_redundancy.size(); i++) {
            if (!set_test_bit(is_maximal, i)) {
                set_free(incidence_H_to_V_with_redundancy[i]);
                continue;
            }
            double* h = H[count];
            incidence_H_to_V[count] = incidence_H_to_V_with_redundancy[i];
            count++;
            if (i < A.size()) {
                memcpy(h, A[i], sizeof(double) * (K + 1));
            } else {
                int xj = i - (int) A.size();
                if (xj >= xi) {
                    xj++;
                }
                assert(0 <= xj && xj < K && "Sanity checking the range of xj.");
                // xi >= xj
                h[xi + 1] = 1;
                h[xj + 1] = -1;
            }
        }
        assert(count == (int) maximal_H.size() && "count should equal maximal_H.size()");
        set_free(is_maximal);

        // y = xi
        H[count][xi + 1] = 1;
        H[count][K + 1] = -1;
        H[count + 1][xi + 1] = -1;
        H[count + 1][K + 1] = 1;

        incidence_H_to_V[count] = set_create(V.size());
        set_enable_all(incidence_H_to_V[count]);
        incidence_H_to_V[count + 1] = set_create(V.size());
        set_enable_all(incidence_H_to_V[count + 1]);

        quadrant_pdds.push_back({K + 2, V, H, incidence_H_to_V});
    }

    PDD res = intersect_all_PDDs_recursively(quadrant_pdds);
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

    vector<QuadrantInfo> quadrant_infos = compute_max_pool_quadrants_with_cdd(K, A);
    assert((int) quadrant_infos.size() == K && "Number of quadrants should be K.");

    size_t num_vertices = 0;
    for (auto& info : quadrant_infos) {
        num_vertices += info.V.size();
    }

    dd_MatrixPtr vertices = dd_CreateMatrix(num_vertices, K + 2);
    vertices->representation = dd_Generator;
    size_t counter = 0;

    for (int xi = 0; xi < K; xi++) {
        auto& V_quadrant = quadrant_infos[xi].V;
        for (mpq_t* v : V_quadrant) {
            mpq_t *v_projection = vertices->matrix[counter];
            counter++;
            for (int i = 0; i < K + 1; i++) {
                mpq_set(v_projection[i], v[i]);
            }
            mpq_set(v_projection[K + 1], v[xi + 1]);
        }
        mpq_mat_free(K + 1, V_quadrant);
        set_arr_free(quadrant_infos[xi].V_to_H_incidence);
    }
    assert(counter == num_vertices && "Consistency checking that counter equals the number of vertices.");

    vector<double*> H = cdd_compute_inequalities_from_vertices(vertices);
    dd_FreeMatrix(vertices);

    return H;
}

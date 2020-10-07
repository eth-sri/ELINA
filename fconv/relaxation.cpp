#include <map>
#include <cassert>
#include <cstring>
#include "relaxation.h"
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
          abs_coef = max(abs(h[i]), abs_coef);
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

vector<double *> fkrelu(const int K, const vector<double *> &A) {
  ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
  verify_that_octahedron_and_all_xi_split_zero(K, A);
  if (K == 1) {
    return relu_1(-A[0][0], A[1][0]);
  }
  OctahedronV oct = compute_octahedron_V(K, A);

  // Split in quadrants takes care of memory management of input vertices.
  map<Quadrant, QuadrantInfo> quadrant2info =
      split_in_quadrants(oct.V, oct.incidence, oct.orthant_adjacencies, K);

  map<Quadrant, PDD> quadrant2pdd;
  for (auto &pair : quadrant2info) {
    const Quadrant &quadrant = pair.first;
    vector<mpq_t *> &V_mpq = pair.second.V;

    if (V_mpq.empty()) {
      // Input in the quadrant is the empty set.
      quadrant2pdd[quadrant] = {K + 1, {}, {}, {}};
      continue;
    }

    vector<double *> V = mpq_mat_to_fp(K + 1, V_mpq);
    mpq_mat_free(K + 1, V_mpq);

    const vector<set_t> &incidence_V_to_H = pair.second.V_to_H_incidence;
    assert(incidence_V_to_H.size() == V.size() &&
           "Incidence_V_to_H.size() should equal V.size()");
    vector<set_t> incidence_H_to_V_with_redundancy =
        set_arr_transpose(incidence_V_to_H);
    set_arr_free(incidence_V_to_H);
    assert(incidence_H_to_V_with_redundancy.size() == A.size() + K &&
           "Incidence_H_to_V_with_redundancy.size() should equal A.size() + K");
    vector<int> maximal_H =
        compute_maximal_indexes(incidence_H_to_V_with_redundancy);
    set_t is_maximal = set_create(incidence_H_to_V_with_redundancy.size());
    for (auto i : maximal_H) {
      set_enable_bit(is_maximal, i);
    }

    vector<double *> H(maximal_H.size());
    vector<set_t> incidence_H_to_V(maximal_H.size());

    int count = 0;
    for (size_t i = 0; i < incidence_H_to_V_with_redundancy.size(); i++) {
      if (!set_test_bit(is_maximal, i)) {
        set_free(incidence_H_to_V_with_redundancy[i]);
        continue;
      }
      double *h = (double *)calloc(K + 1, sizeof(double));
      H[count] = h;
      incidence_H_to_V[count] = incidence_H_to_V_with_redundancy[i];
      count++;
      if (i < A.size()) {
        memcpy(h, A[i], sizeof(double) * (K + 1));
      } else {
        int xi = i - (int)A.size();
        assert(0 <= xi && xi < K && "Sanity checking the range of xi.");
        if (quadrant[xi] == MINUS) {
          h[xi + 1] = -1;
        } else {
          h[xi + 1] = 1;
        }
      }
    }
    assert(count == (int)maximal_H.size() &&
           "count should equal maximal_H.size()");
    set_free(is_maximal);
    quadrant2pdd[quadrant] = {K + 1, V, H, incidence_H_to_V};
  }
  return decomposition(K, quadrant2pdd);
}

vector<double*> krelu_with_cdd(const int K, const vector<double*>& A) {
    // No need to verify since CDD can work with input in any format.
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
    map<Quadrant, QuadrantInfo> quadrant2info =
        compute_quadrants_with_cdd(K, A);

    size_t num_vertices = 0;
    for (auto &entry : quadrant2info) {
      num_vertices += entry.second.V.size();
    }

    dd_MatrixPtr vertices = dd_CreateMatrix(num_vertices, 2 * K + 1);
    size_t counter = 0;

    for (auto &entry : quadrant2info) {
      const auto &quadrant = entry.first;
      auto &V_quadrant = entry.second.V;

      for (const auto &v : V_quadrant) {
        mpq_t *v_projection = vertices->matrix[counter];
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

    // TODO[gleb] Implement specialized fast compute_max_pool_quadrants.
    vector<QuadrantInfo> quadrant_infos =
        compute_max_pool_quadrants_with_cdd(K, A);
    assert((int)quadrant_infos.size() == K &&
           "Number of quadrants should be K.");

    vector<PDD> quadrant_pdds;
    for (int xi = 0; xi < K; xi++) {
      vector<mpq_t *> &V_mpq = quadrant_infos[xi].V;

      if (V_mpq.empty()) {
        quadrant_pdds.push_back({K + 2, {}, {}, {}});
        continue;
        }

        vector<double*> V = mpq_mat_to_fp(K + 1, V_mpq);
        mpq_mat_free(K + 1, V_mpq);
        for (size_t i = 0; i < V.size(); i++) {
          V[i] = (double *)realloc(V[i], (K + 2) * sizeof(double));
          V[i][K + 1] = V[i][xi + 1];
        }

        const vector<set_t> &incidence_V_to_H =
            quadrant_infos[xi].V_to_H_incidence;
        assert(
                incidence_V_to_H.size() == V.size() &&
                "Incidence_V_to_H.size() should equal V.size()");
        vector<set_t> incidence_H_to_V_with_redundancy =
            set_arr_transpose(incidence_V_to_H);
        set_arr_free(incidence_V_to_H);
        assert(incidence_H_to_V_with_redundancy.size() == A.size() + K - 1 &&
               "Incidence_H_to_V_with_redundancy.size() should equal A.size() "
               "+ K - 1");

        vector<int> maximal_H =
            compute_maximal_indexes(incidence_H_to_V_with_redundancy);
        set_t is_maximal = set_create(incidence_H_to_V_with_redundancy.size());
        for (auto i : maximal_H) {
          set_enable_bit(is_maximal, i);
        }

        // maximal_H.size() + 2 because there will also be new equality y = xi.
        vector<double *> H = fp_mat_create(maximal_H.size() + 2, K + 2);
        vector<set_t> incidence_H_to_V(maximal_H.size() + 2);

        int count = 0;
        for (size_t i = 0; i < incidence_H_to_V_with_redundancy.size(); i++) {
          if (!set_test_bit(is_maximal, i)) {
            set_free(incidence_H_to_V_with_redundancy[i]);
            continue;
          }
          double *h = H[count];
          incidence_H_to_V[count] = incidence_H_to_V_with_redundancy[i];
          count++;
          if (i < A.size()) {
            memcpy(h, A[i], sizeof(double) * (K + 1));
          } else {
            int xj = i - (int)A.size();
            if (xj >= xi) {
              xj++;
            }
            assert(0 <= xj && xj < K && "Sanity checking the range of xj.");
            // xi >= xj
            h[xi + 1] = 1;
            h[xj + 1] = -1;
          }
        }
        assert(count == (int)maximal_H.size() &&
               "count should equal maximal_H.size()");
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

    PDD res = PDD_intersect_all(quadrant_pdds);
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

    vector<QuadrantInfo> quadrant_infos =
        compute_max_pool_quadrants_with_cdd(K, A);
    assert((int)quadrant_infos.size() == K &&
           "Number of quadrants should be K.");

    size_t num_vertices = 0;
    for (auto &info : quadrant_infos) {
      num_vertices += info.V.size();
    }

    dd_MatrixPtr vertices = dd_CreateMatrix(num_vertices, K + 2);
    size_t counter = 0;

    for (int xi = 0; xi < K; xi++) {
      auto &V_quadrant = quadrant_infos[xi].V;
      for (mpq_t *v : V_quadrant) {
        mpq_t *v_projection = vertices->matrix[counter];
        counter++;
        mpq_arr_set(K + 1, v_projection, v);
        mpq_set(v_projection[K + 1], v[xi + 1]);
        }
        mpq_mat_free(K + 1, V_quadrant);
        set_arr_free(quadrant_infos[xi].V_to_H_incidence);
    }
    assert(counter == num_vertices && "Counter should equal the number of vertices.");

    vector<double*> H = cdd_compute_inequalities_from_vertices(vertices);
    dd_FreeMatrix(vertices);

    return H;
}

vector<double *> ktanh_with_cdd(int K, const vector<double *> &A) {
  ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
  verify_that_octahedron_and_all_xi_split_zero(K, A);

  map<Quadrant, vector<mpq_t *>> quadrant2vertices =
      compute_tanh_quadrants_with_cdd(K, A);

  size_t num_vertices = 0;
  for (auto &entry : quadrant2vertices) {
    num_vertices += entry.second.size();
  }

  dd_MatrixPtr vertices = dd_CreateMatrix(num_vertices, 2 * K + 1);
  size_t counter = 0;
  for (auto &entry : quadrant2vertices) {
    for (mpq_t *v : entry.second) {
      mpq_arr_set(2 * K + 1, vertices->matrix[counter], v);
      counter++;
    }
    mpq_mat_free(K + 1, entry.second);
  }
  assert(counter == num_vertices &&
         "Counter should equal the number of vertices.");

  vector<double *> H = cdd_compute_inequalities_from_vertices(vertices);
  dd_FreeMatrix(vertices);

  return H;
}

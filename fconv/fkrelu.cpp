#include "fkrelu.h"
#include "decomposition.h"
#include "mpq.h"
#include "octahedron.h"
#include "pdd.h"
#include "split_in_quadrants.h"
#include "utils.h"
#include <map>

vector<bset> transpose_incidence(const vector<bset> &input) {
  ASRTF(!input.empty(), "The function doesn't support empty input.");

  const size_t num_rows = input.size();
  const size_t num_cols = input[0].size();

  vector<bset> output(num_cols, bset(num_rows));
  for (size_t row = 0; row < num_rows; row++) {
    const bset &in_row = input[row];
    assert(in_row.size() == num_cols &&
           "All incidences should be of the same size.");
    for (size_t col = 0; col < num_cols; col++) {
      if (!in_row.test(col)) {
        continue;
      }
      output[col].set(row);
    }
  }

  return output;
}

// Computation of relaxation for 1-relu easily done with analytical formula.
vector<double *> relu_1(double lb, double ub) {
  ASRTF(lb <= ub, "Unsoundness - lower bound should be <= then upper bound.");
  ASRTF(lb < 0 && 0 < ub, "Expecting non-trivial input where lb < 0 < ub.");

  double lmd = -lb * ub / (ub - lb);
  double mu = ub / (ub - lb);
  assert(lmd > 0 && "Expected lmd > 0.");
  assert(mu > 0 && "Expected mu > 0.");

  vector<double *> res = create_mat(3, 3);

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

void verify_fkrelu_input(const int K, const vector<double *> &A) {
  ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
  ASRTF((int)A.size() == POW3[K] - 1,
        "Unexpected number of rows in the input.");
  const vector<vector<int>> &coefs = K2OCTAHEDRON_COEFS[K];
  for (int i = 0; i < (int)A.size(); i++) {
    for (int j = 0; j < K; j++) {
      ASRTF(A[i][j + 1] == coefs[i][j], "Input is not of correct format.");
    }
  }
}

vector<double *> fkrelu(const int K, const vector<double *> &A) {
  ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
  //    cout << "In fkrelu K is "  << K << endl;
  verify_fkrelu_input(K, A);
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

    vector<double *> V = mpq_to_double(K + 1, V_mpq);
    mpq_free_array_vector(K + 1, V_mpq);

    const vector<bset> &incidence_V_to_H = pair.second.V_to_H_incidence;
    assert(incidence_V_to_H.size() == V.size() &&
           "Incidence_V_to_H size should equal V.size()");
    vector<bset> incidence_H_to_V_with_redundancy =
        transpose_incidence(incidence_V_to_H);
    assert(incidence_H_to_V_with_redundancy.size() == A.size() + K &&
           "Incidence_H_to_V_with_redundancy size should equal A.size() + K");
    vector<int> maximal_H =
        compute_maximal_indexes(incidence_H_to_V_with_redundancy);

    vector<double *> H(maximal_H.size());
    vector<bset> incidence_H_to_V(maximal_H.size());

    for (size_t i = 0; i < maximal_H.size(); i++) {
      double *h = (double *)calloc(K + 1, sizeof(double));
      H[i] = h;
      int maximal = maximal_H[i];
      incidence_H_to_V[i] = incidence_H_to_V_with_redundancy[maximal];
      if (maximal < (int)A.size()) {
        for (int i = 0; i < K + 1; i++) {
          h[i] = A[maximal][i];
        }
      } else {
        int xi = maximal - (int)A.size();
        assert(0 <= xi && xi < K && "Sanity checking the range of xi.");
        if (quadrant[xi] == MINUS) {
          h[xi + 1] = -1;
        } else {
          h[xi + 1] = 1;
        }
      }
    }

    quadrant2pdd[quadrant] = {K + 1, V, H, incidence_H_to_V};
  }
  return decomposition(K, quadrant2pdd);
}

vector<double *> krelu_with_cdd(const int K, const vector<double *> &A) {
  // No need to verify since CDD can work with input in any format.
  ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
  map<Quadrant, QuadrantInfo> quadrant2info = compute_quadrants_with_cdd(K, A);

  size_t num_vertices = 0;
  for (auto &entry : quadrant2info) {
    num_vertices += entry.second.V.size();
  }

  dd_MatrixPtr vertices = dd_CreateMatrix(num_vertices, 2 * K + 1);
  vertices->representation = dd_Generator;
  size_t counter = 0;

  // Now I will compute the final V. This will allow me to verify that produced
  // constraints do not violate any of the original vertices and thus I produce
  // a sound overapproximation.
  vector<mpq_t *> V;
  for (auto &entry : quadrant2info) {
    const auto &quadrant = entry.first;
    auto &V_quadrant = entry.second.V;

    for (const auto &v : V_quadrant) {
      mpq_t *v_projection = vertices->matrix[counter];
      counter++;
      for (int i = 0; i < K + 1; i++) {
        mpq_set(v_projection[i], v[i]);
      }
      for (int i = 0; i < K; i++) {
        if (quadrant[i] == PLUS) {
          // Only need to set for the case of PLUS, because in case of MINUS
          // there should be 0. And it is already there automatically.
          mpq_set(v_projection[1 + i + K], v[1 + i]);
        }
      }
    }

    mpq_free_array_vector(K + 1, V_quadrant);
  }
  assert(counter == num_vertices &&
         "Consistency checking that counter equals the number of vertices.");

  dd_ErrorType err = dd_NoError;
  dd_PolyhedraPtr poly = dd_DDMatrix2Poly(vertices, &err);
  ASRTF(err == dd_NoError,
        "Converting matrix to polytope failed with error " + to_string(err));

  dd_MatrixPtr inequalities = dd_CopyInequalities(poly);

  vector<double *> H = cdd2double(inequalities);
  for (auto h : H) {
    double abs_coef = abs(h[0]);
    for (int i = 1; i < 2 * K + 1; i++) {
      abs_coef = min(abs(h[i]), abs_coef);
    }
    for (int i = 0; i < 2 * K + 1; i++) {
      h[i] /= abs_coef;
    }
  }

  dd_FreePolyhedra(poly);
  dd_FreeMatrix(vertices);
  dd_FreeMatrix(inequalities);

  return H;
}

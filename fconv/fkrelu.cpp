#include "fkrelu.h"
#include "decomposition.h"
#include "mpq.h"
#include "octahedron.h"
#include "pdd.h"
#include "split_in_quadrants.h"
#include "utils.h"
#include <Eigen/Dense>
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
MatrixXd relu_1(double lb, double ub) {
  ASRTF(lb <= ub, "Unsoundness - lower bound should be <= then upper bound.");
  ASRTF(lb < 0 && 0 < ub, "Expecting non-trivial input where lb < 0 < ub.");

  double lmd = -lb * ub / (ub - lb);
  double mu = ub / (ub - lb);
  assert(lmd > 0 && "Expected lmd > 0.");
  assert(mu > 0 && "Expected mu > 0.");

  MatrixXd res(3, 3);

  res << 0, 0, 1,  // y >= 0
      0, -1, 1,    // y >= x
      lmd, mu, -1; // y <= mu * x + lmd;

  return res;
}

void verify_fkrelu_input(const MatrixXd &A) {
  const int K = A.cols() - 1;
  ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
  ASRTF(A.rows() == POW3[K] - 1, "Unexpected number of rows in the input.");
  const vector<vector<int>> &coefs = K2OCTAHEDRON_COEFS[K];
  for (int i = 0; i < A.rows(); i++) {
    const RowXd &row = A.row(i);
    const vector<int> &coef_row = coefs[i];
    for (int j = 0; j < K; j++) {
      ASRTF(row(0, j + 1) == coef_row[j], "Input is not of correct format.");
    }
  }
}

MatrixXd fkrelu(const MatrixXd &A) {
  const int K = A.cols() - 1;
  ASRTF(1 <= K && K <= 4, "K should be within allowed range.");
  verify_fkrelu_input(A);
  if (K == 1) {
    return relu_1(-A(0, 0), A(1, 0));
  }
  OctahedronV oct = compute_octahedron_V(A);

  // Split in quadrants takes care of memory management of input vertices.
  map<Quadrant, QuadrantInfo> quadrant2info =
      split_in_quadrants(oct.V, oct.incidence, oct.orthant_adjacencies, K);

  map<Quadrant, PDD> quadrant2pdd;
  for (auto &pair : quadrant2info) {
    const Quadrant &quadrant = pair.first;
    vector<mpq_t *> &V_mpq = pair.second.V;

    if (V_mpq.empty()) {
      // Input in the quadrant is the empty set.
      MatrixXd empty(0, K + 1);
      quadrant2pdd[quadrant] = {K + 1, empty, empty, {}};
      continue;
    }

    MatrixXd V = mpq_convert2eigen(K + 1, V_mpq);
    mpq_free_array_vector(K + 1, V_mpq);

    const vector<bset> &incidence_V_to_H = pair.second.V_to_H_incidence;
    assert((int)incidence_V_to_H.size() == V.rows() &&
           "Incidence_V_to_H size should equal V.rows()");
    vector<bset> incidence_H_to_V_with_redundancy =
        transpose_incidence(incidence_V_to_H);
    assert((int)incidence_H_to_V_with_redundancy.size() == A.rows() + K &&
           "Incidence_H_to_V_with_redundancy size should equal A.rows() + K");
    vector<int> maximal_H =
        compute_maximal_indexes(incidence_H_to_V_with_redundancy);

    MatrixXd H(maximal_H.size(), K + 1);
    vector<bset> incidence_H_to_V(maximal_H.size());

    for (int i = 0; i < (int)maximal_H.size(); i++) {
      int maximal = maximal_H[i];
      incidence_H_to_V[i] = incidence_H_to_V_with_redundancy[maximal];
      if (maximal < A.rows()) {
        H.row(i) = A.row(maximal);
      } else {
        H.row(i) = MatrixXd::Zero(1, K + 1);
        int xi = maximal - A.rows();
        assert(0 <= xi && xi < K && "Sanity checking the range of xi.");
        if (quadrant[xi] == MINUS) {
          H(i, xi + 1) = -1;
        } else {
          H(i, xi + 1) = 1;
        }
      }
    }

    quadrant2pdd[quadrant] = {K + 1, V, H, incidence_H_to_V};
  }
  return decomposition(quadrant2pdd, K);
}
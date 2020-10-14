#include <algorithm>
#include <cassert>
#include <math.h>
#include "quadrants.h"
#include "octahedron.h"
#include "fp_mat.h"
#include "mpq.h"

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

void precompute_tanh(double x_fp, mpq_t *x, mpq_t *y, mpq_t *k_lb, mpq_t *b_lb,
                     mpq_t *k_ub, mpq_t *b_ub) {
  // This computations are important to do in exact precision,
  // because in that case we can practically demonstrate dimension trick.
  assert(x_fp != 0 && "x should be non-zero.");
  mpq_t one;
  mpq_init(one);
  mpq_set_si(one, 1, 1);

  mpq_t *dy;
  mpq_t *b;
  mpq_t *k;

  if (x_fp < 0) {
    dy = k_lb;
    b = b_lb;
    k = k_ub;
  } else {
    dy = k_ub;
    b = b_ub;
    k = k_lb;
  }

  mpq_set_d(*x, x_fp);
  mpq_set_d(*y, tanh(x_fp));

  // One of the lines is x * dy + b
  // dy is the derivative of tanh: dy = 1 - y^2
  mpq_mul(*dy, *y, *y);
  mpq_sub(*dy, one, *dy);
  mpq_mul(*b, *x, *dy);
  mpq_sub(*b, *y, *b);

  // Another line has coefficient y / x
  mpq_div(*k, *y, *x);

  mpq_clear(one);
}

void precompute_sigm(double x_fp, mpq_t *x, mpq_t *y, mpq_t *k_lb, mpq_t *b_lb,
                     mpq_t *k_ub, mpq_t *b_ub) {
  // This computations are important to do in exact precision,
  // because in that case we can practically demonstrate dimension trick.
  assert(x_fp != 0 && "x should be non-zero.");
  mpq_t one;
  mpq_init(one);
  mpq_set_si(one, 1, 1);

  mpq_t *dy;
  mpq_t *dy_b;
  mpq_t *k;
  mpq_t *k_b;

  if (x_fp < 0) {
    dy = k_lb;
    dy_b = b_lb;
    k = k_ub;
    k_b = b_ub;
  } else {
    dy = k_ub;
    dy_b = b_ub;
    k = k_lb;
    k_b = b_lb;
  }

  mpq_set_d(*x, x_fp);
  // y = 1 / (1 + e^(-x))
  mpq_set_d(*y, exp(-x_fp));
  mpq_add(*y, one, *y);
  mpq_div(*y, one, *y);

  // One of the lines is x * dy + b
  // dy is the derivative of sigm: dy = y * (1 - y)
  mpq_sub(*dy, one, *y);
  mpq_mul(*dy, *y, *dy);
  mpq_mul(*dy_b, *x, *dy);
  mpq_sub(*dy_b, *y, *dy_b);

  // Another line is k * x + 1/2
  mpq_set_si(*k_b, 1, 2);
  mpq_sub(*k, *y, *k_b);
  mpq_div(*k, *k, *x);

  mpq_clear(one);
}

void precompute_arr_tasi(const int K, const vector<double *> &A,
                         Activation activation, vector<mpq_t *> &x,
                         vector<mpq_t *> &y, vector<mpq_t *> &k_lb,
                         vector<mpq_t *> &b_lb, vector<mpq_t *> &k_ub,
                         vector<mpq_t *> &b_ub) {
  x = mpq_mat_create(2, K);
  y = mpq_mat_create(2, K);
  k_lb = mpq_mat_create(2, K);
  b_lb = mpq_mat_create(2, K);
  k_ub = mpq_mat_create(2, K);
  b_ub = mpq_mat_create(2, K);

  for (int p = 0; p <= 1; p++) {
    for (int xi = 0; xi < K; xi++) {
      double x_fp = (p == 0) ? -A[LOWER_BOUND_INDEX[K][xi]][0]
                             : A[UPPER_BOUND_INDEX[K][xi]][0];
      if (activation == Tanh) {
        precompute_tanh(x_fp, &(x[p][xi]), &(y[p][xi]), &(k_lb[p][xi]),
                        &(b_lb[p][xi]), &(k_ub[p][xi]), &(b_ub[p][xi]));
      } else {
        precompute_sigm(x_fp, &(x[p][xi]), &(y[p][xi]), &(k_lb[p][xi]),
                        &(b_lb[p][xi]), &(k_ub[p][xi]), &(b_ub[p][xi]));
      }
    }
  }
}

map<Quadrant, vector<mpq_t*>> get_tasi_quadrants_cdd(
        const int K,
        const vector<double*>& A,
        Activation activation) {
    ASRTF(1 <= K && K <= 4, "K should be within allowed range.");

    const vector<Quadrant>& quadrants = K2QUADRANTS[K];

    const int NUM_H = (int) A.size();

    // Total NUM_H + 3 * K inequalities:
    // - NUM_H original inequalities
    // - K inequalities for a quadrant
    // - 2 * K inequalities for y (2 per quadrant)
    dd_MatrixPtr cdd_A = dd_CreateMatrix(NUM_H + 3 * K, 2 * K + 1);
    cdd_A->representation = dd_Inequality;

    for (int i = 0; i < NUM_H; i++) {
        mpq_arr_set_d(K + 1, cdd_A->matrix[i], A[i]);
    }

    vector<mpq_t *> x, y, k_lb, b_lb, k_ub, b_ub;
    precompute_arr_tasi(K, A, activation, x, y, k_lb, b_lb, k_ub, b_ub);

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
          mpq_t *row_lb = cdd_A->matrix[NUM_H + K + 2 * xi];
          mpq_t *row_ub = cdd_A->matrix[NUM_H + K + 2 * xi + 1];
          mpq_arr_set_zero(2 * K + 1, row_lb);
          mpq_arr_set_zero(2 * K + 1, row_ub);

          Polarity p = quadrant[xi];

          // Lower bound y >= kx + b equivalent -b - kx + y >= 0
          mpq_neg(row_lb[0], b_lb[p][xi]);
          mpq_neg(row_lb[xi + 1], k_lb[p][xi]);
          mpq_set_si(row_lb[xi + 1 + K], 1, 1);

          // Upper bound y <= kx + b equivalent b + kx - y >= 0
          mpq_set(row_ub[0], b_ub[p][xi]);
          mpq_set(row_ub[xi + 1], k_ub[p][xi]);
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

    mpq_mat_free(K, x);
    mpq_mat_free(K, y);
    mpq_mat_free(K, k_lb);
    mpq_mat_free(K, b_lb);
    mpq_mat_free(K, k_ub);
    mpq_mat_free(K, b_ub);

    dd_FreeMatrix(cdd_A);

    return quadrant2V;
}

map<Quadrant, vector<mpq_t *>>
get_tasi_quadrants_cdd_dim(const int K, const vector<double *> &A,
                           Activation activation) {
  ASRTF(1 <= K && K <= 4, "K should be within allowed range.");

  const vector<Quadrant> &quadrants = K2QUADRANTS[K];

  const int NUM_H = (int)A.size();

  dd_MatrixPtr cdd_A = dd_CreateMatrix(NUM_H + K, K + 1);
  cdd_A->representation = dd_Inequality;

  for (int i = 0; i < NUM_H; i++) {
    mpq_arr_set_d(K + 1, cdd_A->matrix[i], A[i]);
  }

  vector<mpq_t *> x, y, k_lb, b_lb, k_ub, b_ub;
  precompute_arr_tasi(K, A, activation, x, y, k_lb, b_lb, k_ub, b_ub);

  map<Quadrant, vector<mpq_t *>> quadrant2V;
  for (const auto &quadrant : quadrants) {
    for (int xi = 0; xi < K; xi++) {
      mpq_t *row = cdd_A->matrix[NUM_H + xi];
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

    vector<mpq_t *> V;
    for (size_t vi = 0; vi < num_v; vi++) {
      const size_t base = V.size();
      mpq_t *v_base = cdd_V->matrix[vi];
      V.push_back(mpq_arr_create(2 * K + 1));
      mpq_arr_set(K + 1, V[base], v_base);
      for (int xi = 0; xi < K; xi++) {
        Polarity p = quadrant[xi];
        if (!mpq_cmp(v_base[xi + 1], x[p][xi])) {
          // In case x has extreme value both lower and upper bound of y
          // have the same value - thus y can be set directly.
          for (size_t i = base; i < V.size(); i++) {
            mpq_set(V[i][xi + 1 + K], y[p][xi]);
          }
          continue;
        }
        const size_t num = V.size() - base;
        for (size_t i = base; i < base + num; i++) {
          mpq_t *v_lb = V[i];
          mpq_t *v_ub = mpq_arr_copy(2 * K + 1, v_lb);
          V.push_back(v_ub);

          mpq_mul(v_lb[xi + 1 + K], v_lb[xi + 1], k_lb[p][xi]);
          mpq_add(v_lb[xi + 1 + K], v_lb[xi + 1 + K], b_lb[p][xi]);

          mpq_mul(v_ub[xi + 1 + K], v_ub[xi + 1], k_ub[p][xi]);
          mpq_add(v_ub[xi + 1 + K], v_ub[xi + 1 + K], b_ub[p][xi]);
        }
        assert(V.size() - base == 2 * num &&
               "The number of new vertices should've doubled.");
      }
    }

    dd_FreePolyhedra(poly);
    dd_FreeMatrix(cdd_V);

    quadrant2V[quadrant] = V;
  }

  mpq_mat_free(K, x);
  mpq_mat_free(K, y);
  mpq_mat_free(K, k_lb);
  mpq_mat_free(K, b_lb);
  mpq_mat_free(K, k_ub);
  mpq_mat_free(K, b_ub);

  dd_FreeMatrix(cdd_A);

  return quadrant2V;
}

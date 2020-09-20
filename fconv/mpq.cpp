#include "mpq.h"
#include "asrt.h"
#include <Eigen/Dense>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

mpq_t *mpq_create_array(const int n) {
  auto arr = (mpq_t *)calloc(n, sizeof(mpq_t));
  ASRTF(arr != nullptr, "Calloc has failed.");
  for (int i = 0; i < n; i++) {
    mpq_init(arr[i]);
  }
  return arr;
}

mpq_t *mpq_create_and_copy_array(int n, const mpq_t *src) {
  mpq_t *arr = mpq_create_array(n);
  for (int i = 0; i < n; i++) {
    mpq_set(arr[i], src[i]);
  }
  return arr;
}

void mpq_free_array(int n, mpq_t *arr) {
  for (int i = 0; i < n; i++) {
    mpq_clear(arr[i]);
  }
  free(arr);
}

bool mpq_arrays_are_equal(const int n, mpq_t *first, mpq_t *second) {
  for (int i = 0; i < n; i++) {
    if (mpq_cmp(first[i], second[i]) != 0) {
      return false;
    }
  }
  return true;
}

Eigen::MatrixXd mpq_convert2eigen(const int n,
                                  const std::vector<mpq_t *> &mpq_A) {
  Eigen::MatrixXd eig_A(mpq_A.size(), n);
  for (size_t i = 0; i < mpq_A.size(); i++) {
    mpq_t *row = mpq_A[i];
    for (int j = 0; j < n; j++) {
      eig_A(i, j) = mpq_get_d(row[j]);
    }
  }
  return eig_A;
}

vector<mpq_t *> mpq_from_eigen(const Eigen::MatrixXd &A) {
  vector<mpq_t *> out(A.rows());

  const int num_cols = A.cols();
  for (int i = 0; i < (int)A.rows(); i++) {
    const auto &row_fp = A.row(i);
    mpq_t *row_mpq = mpq_create_array(num_cols);
    out[i] = row_mpq;
    for (int j = 0; j < num_cols; j++) {
      mpq_set_d(row_mpq[j], row_fp(j));
    }
  }

  return out;
}

vector<mpq_t *> mpq_copy_array_vector(int n,
                                      const vector<mpq_t *> &src_vector) {
  vector<mpq_t *> copy(src_vector.size());

  for (size_t i = 0; i < src_vector.size(); i++) {
    copy[i] = mpq_create_and_copy_array(n, src_vector[i]);
  }

  return copy;
}

void mpq_free_array_vector(int n, vector<mpq_t *> &vec) {
  for (const auto &arr : vec) {
    mpq_free_array(n, arr);
  }
}

void mpq_print_array(int n, const mpq_t *v) {
  for (int i = 0; i < n; i++) {
    cout << mpq_get_d(v[i]) << " ";
  }
  cout << endl;
}

vector<mpq_t *> mpq_matrix_mul(const int n, const vector<mpq_t *> &A,
                               const vector<mpq_t *> &B) {
  // This function does A * B' multiplication.
  ASRTF(!A.empty() && !B.empty(), "A and B should be non-empty");

  mpq_t prod, prod2, cur, cur2;
  mpq_inits(prod, prod2, cur, cur2, NULL);

  vector<mpq_t *> res(A.size());
  for (size_t i = 0; i < A.size(); i++) {
    auto &res_row = res[i];
    auto &A_row = A[i];

    res_row = mpq_create_array(B.size());

    for (size_t j = 0; j < B.size(); j++) {
      auto &B_row = B[j];
      mpq_set_si(cur, 0, 1);

      for (int d = 0; d < n; d++) {
        mpq_mul(prod, A_row[d], B_row[d]);
        mpq_add(cur, cur, prod);
      }

      mpq_set(res_row[j], cur);
    }
  }

  mpq_clears(prod, prod2, cur, cur2, NULL);

  return res;
}

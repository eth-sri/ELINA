#include "mpq.h"
#include "asrt.h"
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

mpq_t *mpq_create_and_copy_array(const int n, const mpq_t *src) {
  mpq_t *arr = mpq_create_array(n);
  for (int i = 0; i < n; i++) {
    mpq_set(arr[i], src[i]);
  }
  return arr;
}

void mpq_free_array(const int n, mpq_t *arr) {
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

vector<double *> mpq_to_double(const int n, const vector<mpq_t *> &mpq_A) {
  vector<double *> A(mpq_A.size());
  for (size_t i = 0; i < mpq_A.size(); i++) {
    double *row = (double *)calloc(n, sizeof(double));
    A[i] = row;
    const mpq_t *mpq_row = mpq_A[i];
    for (int j = 0; j < n; j++) {
      row[j] = mpq_get_d(mpq_row[j]);
    }
  }
  return A;
}

vector<mpq_t *> mpq_from_double(const int n, const vector<double *> &A) {
  vector<mpq_t *> mpq_A(A.size());
  for (size_t i = 0; i < A.size(); i++) {
    mpq_t *row_mpq = mpq_create_array(n);
    mpq_A[i] = row_mpq;
    const double *row = A[i];
    for (int j = 0; j < n; j++) {
      mpq_set_d(row_mpq[j], row[j]);
    }
  }
  return mpq_A;
}

vector<mpq_t *> mpq_copy_array_vector(const int n,
                                      const vector<mpq_t *> &src_vector) {
  vector<mpq_t *> copy(src_vector.size());

  for (size_t i = 0; i < src_vector.size(); i++) {
    copy[i] = mpq_create_and_copy_array(n, src_vector[i]);
  }

  return copy;
}

void mpq_free_array_vector(const int n, vector<mpq_t *> &vec) {
  for (const auto &arr : vec) {
    mpq_free_array(n, arr);
  }
}

void mpq_print_array(const int n, const mpq_t *v) {
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

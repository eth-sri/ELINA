#include <cassert>
#include <cstdlib>
#include <iostream>
#include "mpq.h"
#include "asrt.h"
#include <vector>

using namespace std;

mpq_t* mpq_arr_create(const int n) {
    assert(n >= 0 && "n should be non-negative.");
    auto arr = (mpq_t*) calloc(n, sizeof(mpq_t));
    for (int i = 0; i < n; i++) {
        mpq_init(arr[i]);
    }
    return arr;
}

mpq_t *mpq_arr_copy(const int n, const mpq_t *src) {
  assert(n > 0 && "n should be positive.");
  mpq_t *arr = mpq_arr_create(n);
  for (int i = 0; i < n; i++) {
    mpq_set(arr[i], src[i]);
  }
  return arr;
}

mpq_t* mpq_arr_resize(const int new_n, const int old_n, mpq_t* arr) {
    assert(old_n >= 0 && new_n >= 0 && "new_n and old_n should be non-negative.");
    if (new_n == old_n) {
        return arr;
    }
    // If new_n is smaller than old_n, this will free all removed mpq_t.
    for (int i = new_n; i < old_n; i++) {
        mpq_clear(arr[i]);
    }
    mpq_t* new_arr = (mpq_t*) realloc(arr, new_n * sizeof(mpq_t));
    // If new_n is greater than old_n, this will allocate all added mpq_t.
    for (int i = old_n; i < new_n; i++) {
        mpq_init(new_arr[i]);
    }
    return new_arr;
}

void mpq_arr_free(const int n, mpq_t* arr) {
    assert(n >= 0 && "n should be non-negative.");
    for (int i = 0; i < n; i++) {
        mpq_clear(arr[i]);
    }
    free(arr);
}

void mpq_arr_set_zero(const int n, mpq_t* arr) {
    assert(n >= 0 && "n should be non-negative.");
    for (int i = 0; i < n; i++) {
        mpq_set_si(arr[i], 0, 1);
    }
}

bool mpq_arr_equal(const int n, mpq_t* first, mpq_t* second) {
    assert(n >= 0 && "n should be non-negative.");
    for (int i = 0; i < n; i++) {
        if (mpq_cmp(first[i], second[i]) != 0) {
            return false;
        }
    }
    return true;
}

void mpq_arr_print(const int n, const mpq_t *arr) {
  for (int i = 0; i < n; i++) {
    cout << mpq_get_d(arr[i]) << " ";
  }
  cout << endl;
}

vector<mpq_t*> mpq_mat_copy(const int n, const vector<mpq_t*>& src) {
    assert(n >= 0 && "n should be non-negative.");
    vector<mpq_t*> copy(src.size());

    for (size_t i = 0; i < src.size(); i++) {
        copy[i] = mpq_arr_copy(n, src[i]);
    }

    return copy;
}

void mpq_mat_free(const int n, const vector<mpq_t*>& mat) {
    for (auto row : mat) {
        mpq_arr_free(n, row);
    }
}

vector<mpq_t*> mpq_mat_mul_with_transpose(const int n,
                                         const vector<mpq_t*>& A,
                                         const vector<mpq_t*>& B) {

    mpq_t prod, accum;
    mpq_init(prod);
    mpq_init(accum);

    vector<mpq_t*> Res(A.size());
    for (size_t i = 0; i < A.size(); i++) {
        mpq_t* a = A[i];
        mpq_t* res = mpq_arr_create(B.size());
        Res[i] = res;
        for (size_t j = 0; j < B.size(); j++) {
            mpq_t* b = B[j];
            mpq_set_si(accum, 0, 1);
            for (int d = 0; d < n; d++) {
                mpq_mul(prod, a[d], b[d]);
                mpq_add(accum, accum, prod);
            }
            mpq_set(res[j], accum);
        }
    }

    mpq_clear(prod);
    mpq_clear(accum);

    return Res;
}

vector<double*> mpq_mat_to_fp(const int n, const vector<mpq_t*>& mpq_A) {
    vector<double*> A(mpq_A.size());
    for (size_t i = 0; i < mpq_A.size(); i++) {
        double* row = (double*) calloc(n, sizeof(double));
        A[i] = row;
        const mpq_t* mpq_row = mpq_A[i];
        for (int j = 0; j < n; j++) {
            row[j] = mpq_get_d(mpq_row[j]);
        }
    }
    return A;
}

vector<mpq_t*> mpq_mat_from_fp(const int n, const vector<double*>& A) {
    vector<mpq_t*> mpq_A(A.size());
    for (size_t i = 0; i < A.size(); i++) {
        mpq_t* row_mpq = mpq_arr_create(n);
        mpq_A[i] = row_mpq;
        const double* row = A[i];
        for (int j = 0; j < n; j++) {
            mpq_set_d(row_mpq[j], row[j]);
        }
    }
    return mpq_A;
}

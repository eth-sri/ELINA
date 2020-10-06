#pragma once

#include <gmp.h>
#include <vector>

using namespace std;

mpq_t* mpq_arr_create(int n);

mpq_t* mpq_arr_copy(int n, mpq_t* src);

void mpq_arr_set(int n, mpq_t* dst, mpq_t* src);

void mpq_arr_set_d(int n, mpq_t* dst, double* src);

mpq_t* mpq_arr_resize(int new_n, int old_n, mpq_t* arr);

void mpq_arr_free(int n, mpq_t* arr);

void mpq_arr_set_zero(int n, mpq_t* arr);

bool mpq_arr_equal(int n, mpq_t* first, mpq_t* second);

void mpq_arr_print(int n, mpq_t* arr);

vector<mpq_t*> mpq_mat_copy(int n, const vector<mpq_t*>& src);

void mpq_mat_free(int n, const vector<mpq_t*>& mat);

vector<mpq_t*> mpq_mat_mul_with_transpose(int n, const vector<mpq_t*>& A, const vector<mpq_t*>& B);

vector<double*> mpq_mat_to_fp(int n, const vector<mpq_t*>& mpq_A);

vector<mpq_t*> mpq_mat_from_fp(int n, const vector<double*>& A);

#pragma once

#include <vector>
#include "setoper.h"
#include "cdd.h"

using namespace std;

double* fp_arr_create(int n);

double* fp_arr_copy(int n, double* src);

void fp_arr_set(int n, double* dst, double* src);

double* fp_arr_resize(int new_n, int old_n, double* arr);

vector<double*> fp_mat_create(int rows, int cols);

vector<double*> fp_mat_copy(int cols, const vector<double*>& src);

void fp_mat_free(const vector<double*>& mat);

vector<double*> fp_mat_mul_with_transpose(int n,
                                       const vector<double*>& mat1,
                                       const vector<double*>& mat2);

vector<double*> fp_mat_read(int cols, const string& path);

void fp_mat_print(int cols, const vector<double*>& mat);

dd_MatrixPtr fp_mat_to_cdd(int cols, const vector<double*>& A);

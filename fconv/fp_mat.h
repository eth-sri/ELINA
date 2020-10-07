#pragma once

#include <vector>
#include "setoper.h"
#include "cdd.h"

using namespace std;

double* fp_arr_copy(const int n, double* src);

vector<double*> fp_mat_create(int rows, int cols);

void fp_mat_free(const vector<double*>& mat);

vector<double*> fp_mat_mul_with_transpose(int n,
                                       const vector<double*>& mat1,
                                       const vector<double*>& mat2);

vector<double*> fp_mat_read(int cols, const string& path);

void fp_mat_print(int n, const vector<double*>& mat);

dd_MatrixPtr fp_mat_to_cdd(int n, const vector<double*>& A);

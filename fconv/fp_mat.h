#pragma once

#include <vector>
#include "setoper.h"
#include "cdd.h"

using namespace std;

vector<double*> fp_mat_create(const int rows, const int cols);

void fp_mat_free(const vector<double*>& mat);

vector<double*> fp_mat_mul_with_transpose(int n,
                                       const vector<double*>& mat1,
                                       const vector<double*>& mat2);

vector<double*> fp_mat_read(int cols, const string& path);

dd_MatrixPtr fp_mat_to_cdd(const int n, const vector<double*>& A);

vector<double*> fp_mat_from_cdd(dd_MatrixPtr cdd_A);

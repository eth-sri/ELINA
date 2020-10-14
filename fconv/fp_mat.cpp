#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <string.h>
#include "fp_mat.h"
#include "asrt.h"
#include "setoper.h"
#include "cdd.h"

using namespace std;

double* fp_arr_create(const int n) {
    return (double*) calloc(n, sizeof(double));
}

double* fp_arr_copy(const int n, double* src) {
    auto arr = (double*) calloc(n, sizeof(double));
    memcpy(arr, src, n * sizeof(double));
    return arr;
}

void fp_arr_set(const int n, double* dst, double* src) {
    memcpy(dst, src, n * sizeof(double));
}

double* fp_arr_resize(const int new_n, const int old_n, double* arr) {
    if (new_n == old_n) {
        return arr;
    }
    arr = (double*) realloc(arr, new_n * sizeof(double));
    for (int i = old_n; i < new_n; i++) {
        arr[i] = 0;
    }
    return arr;
}

vector<double*> fp_mat_create(const int rows, const int cols) {
    vector<double*> mat(rows);
    for (int i = 0; i < rows; i++) {
        mat[i] = (double*) calloc(cols, sizeof(double));
    }
    return mat;
}

vector<double*> fp_mat_copy(const int cols, const vector<double*>& src) {
    vector<double*> mat(src.size());
    for (size_t i = 0; i < src.size(); i++) {
        mat[i] = fp_arr_copy(cols, src[i]);
    }
    return mat;
}

void fp_mat_free(const vector<double*>& mat) {
    for (auto v : mat) {
        free(v);
    }
}

// This can be optimized if there will be a need for it.
vector<double*> fp_mat_mul_with_transpose(const int n,
                                         const vector<double*>& A,
                                         const vector<double*>& B) {
    vector<double*> Res = fp_mat_create(A.size(), B.size());
    for (size_t i = 0; i < A.size(); i++) {
        double* a = A[i];
        double* res = Res[i];
        for (size_t j = 0; j < B.size(); j++) {
            double* b = B[j];
            double accum = 0;
            for (int d = 0; d < n; d++) {
                accum += a[d] * b[d];
            }
            res[j] = accum;
        }
    }
    return Res;
}

vector<double*> fp_mat_read(const int cols, const string& path) {
    ifstream file(path);
    ASRTF(file.is_open(), "File could not be opened.");

    string row_string;
    vector<double*> mat;
    while (getline(file, row_string)) {
        istringstream row_stream(row_string);
        double* row = (double*) calloc(cols, sizeof(double));
        int i = 0;
        double number;
        while (row_stream >> number) {
            ASRTF(i < cols, "The number of doubles in the row exceed cols.");
            row[i] = number;
            i++;
        }
        mat.push_back(row);
    }
    return mat;
}

void fp_mat_print(const int cols, const vector<double*>& mat) {
    for (size_t i = 0; i < mat.size(); i++) {
        for (int j = 0; j < cols; j++) {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
}

dd_MatrixPtr fp_mat_to_cdd(const int cols, const vector<double*>& A) {
    dd_MatrixPtr cdd_A = dd_CreateMatrix(A.size(), cols);
    ASRTF(cdd_A != nullptr, "Failed to create cdd_A.");
    for (int i = 0; i < (int) A.size(); i++) {
        for (int j = 0; j < cols; j++) {
            dd_set_d(cdd_A->matrix[i][j], A[i][j]);
        }
    }
    return cdd_A;
}

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "fp_mat.h"
#include "asrt.h"
#include "setoper.h"
#include "cdd.h"

using namespace std;

vector<double*> fp_mat_create(const int rows, const int cols) {
    vector<double*> mat(rows);
    for (int i = 0; i < rows; i++) {
        mat[i] = (double*) calloc(cols, sizeof(double));
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

dd_MatrixPtr fp_mat_to_cdd(const int n, const vector<double*>& A) {
    dd_MatrixPtr cdd_A = dd_CreateMatrix(A.size(), n);
    ASRTF(cdd_A != nullptr, "Failed to create cdd_A.");
    for (int i = 0; i < (int) A.size(); i++) {
        for (int j = 0; j < n; j++) {
            dd_set_d(cdd_A->matrix[i][j], A[i][j]);
        }
    }
    return cdd_A;
}

vector<double*> fp_mat_from_cdd(dd_MatrixPtr cdd_A) {
    const int rows = cdd_A->rowsize;
    const int cols = cdd_A->colsize;

    vector<double*> A = fp_mat_create(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            A[i][j] = mpq_get_d(cdd_A->matrix[i][j]);
        }
    }
    return A;
}

#include <string.h>
#include "utils.h"
#include "fconv.h"
#include "relaxation.h"
#include "sparse_cover.h"
#include "fp_mat.h"

vector<double*> mat_external_to_internal_format(const MatDouble &cmat) {
    vector<double*> A = fp_mat_create(cmat.rows, cmat.cols);
    for (int i = 0; i < cmat.rows; i++) {
        memcpy(A[i], &(cmat.data[i * cmat.cols]), cmat.cols * sizeof(double));
    }
    return A;
}

MatDouble mat_internal_to_external_format(const int n, const vector<double*>& A) {
    auto mat = (double *) calloc(A.size() * n, sizeof(double));
    for (int i = 0; i < (int) A.size(); i++) {
        memcpy(&(mat[i * n]), A[i], n * sizeof(double));
    }
    return {(int) A.size(), n, mat};
}

MatDouble new_MatDouble(int rows, int cols, const double *data) {
    ASRTF(rows > 0, "Only rows > 0 is allowed");
    ASRTF(cols > 0, "Only cols > 0 is allowed");
    auto data_copy = (double *) calloc(rows * cols, sizeof(double));
    for (int i = 0; i < rows * cols; i++) {
        data_copy[i] = data[i];
    }
    return {rows, cols, data_copy};
}

void free_MatDouble(MatDouble cmat) {
    free((double *) cmat.data);
}

void free_MatInt(MatInt cmat) {
    free((int *) cmat.data);
}

enum RelaxationType {
    ReluFast,
    ReluCDD,
    PoolFast,
    PoolCDD
};

MatDouble compute_relaxation(MatDouble input_hrep,
                             RelaxationType type) {
    dd_set_global_constants();
    const int K = input_hrep.cols - 1;
    vector<double*> A = mat_external_to_internal_format(input_hrep);

    vector<double*> H;
    switch (type) {
        case ReluFast:
            H = fkrelu(K, A);
            break;
        case ReluCDD:
            H = krelu_with_cdd(K, A);
            break;
        case PoolFast:
            H = fkpool(K, A);
            break;
        case PoolCDD:
            H = kpool_with_cdd(K, A);
            break;
        default:
            throw runtime_error("Unknown relaxation type.");
    }

    MatDouble out = (type == PoolFast || type == PoolCDD) ?
            mat_internal_to_external_format(K + 2, H) :
            mat_internal_to_external_format(2 * K + 1, H);

    fp_mat_free(A);
    fp_mat_free(H);
    dd_free_global_constants();
    return out;
}

MatDouble fkrelu(MatDouble input_hrep) {
    return compute_relaxation(input_hrep, ReluFast);
}

MatDouble krelu_with_cdd(MatDouble input_hrep) {
    return compute_relaxation(input_hrep, ReluCDD);
}

MatDouble fkpool(MatDouble input_hrep) {
    return compute_relaxation(input_hrep, PoolFast);
}

MatDouble kpool_with_cdd(MatDouble input_hrep) {
    return compute_relaxation(input_hrep, PoolCDD);
}

MatInt generate_sparse_cover(const int N, const int K) {
    vector<vector<int>> cover = sparse_cover(N, K);
    // I'm not sure how to combine std::vector and ctypes thus converting to plain array format.
    auto mat = (int *) calloc(cover.size() * K, sizeof(int));

    int cur = 0;
    for (const auto& comb : cover) {
        for (int i = 0; i < K; i++) {
            mat[cur] = comb[i];
            cur++;
        }
    }

    return {(int) cover.size(), K, mat};
}

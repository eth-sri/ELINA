#include <string.h>
#include "utils.h"
#include "fconv.h"
#include "relaxation.h"
#include "sparse_cover.h"
#include "fp_mat.h"

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

enum Version { Fast, CDD };

MatDouble compute_relaxation(MatDouble input_hrep,
                             Activation activation,
                             Version version) {
    dd_set_global_constants();
    const int K = input_hrep.cols - 1;
    vector<double*> A = mat_external_to_internal_format(input_hrep);

    vector<double*> H;
    if (activation == Relu && version == Fast) {
      H = fast_relaxation_through_decomposition(K, A, Relu);
    } else if (activation == Relu && version == CDD) {
      H = krelu_with_cdd(K, A);
    } else if (activation == Pool && version == Fast) {
      H = fkpool(K, A);
    } else if (activation == Pool && version == CDD) {
      H = kpool_with_cdd(K, A);
    } else if (activation == Tanh && version == Fast) {
      H = fast_relaxation_through_decomposition(K, A, Tanh);
    } else if (activation == Tanh && version == CDD) {
      H = ktanh_with_cdd(K, A);
    } else {
      throw runtime_error("Unknown activation function and version.");
    }

    MatDouble out = (activation == Pool) ?
            mat_internal_to_external_format(K + 2, H) :
            mat_internal_to_external_format(2 * K + 1, H);

    fp_mat_free(A);
    fp_mat_free(H);
    dd_free_global_constants();
    return out;
}

MatDouble fkrelu(MatDouble input_hrep) {
    return compute_relaxation(input_hrep, Relu, Fast);
}

MatDouble krelu_with_cdd(MatDouble input_hrep) {
    return compute_relaxation(input_hrep, Relu, CDD);
}

MatDouble fkpool(MatDouble input_hrep) {
    return compute_relaxation(input_hrep, Pool, Fast);
}

MatDouble kpool_with_cdd(MatDouble input_hrep) {
    return compute_relaxation(input_hrep, Pool, CDD);
}

MatDouble fktanh(MatDouble input_hrep) {
    return compute_relaxation(input_hrep, Tanh, Fast);
}

MatDouble ktanh_with_cdd(MatDouble input_hrep) {
    return compute_relaxation(input_hrep, Tanh, CDD);
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

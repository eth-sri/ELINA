#include <iostream>
#include <algorithm>
#include <cassert>
#include <string.h>
#include "utils.h"
#include "setoper.h"
#include "cdd.h"
#include "dynamic_bitset.h"
#include "fp_mat.h"
#include "mpq.h"

using namespace std;

Timer::Timer() {
    start = chrono::high_resolution_clock::now();
}

int Timer::micros() {
    auto now = chrono::high_resolution_clock::now();
    return chrono::duration_cast<chrono::microseconds>(now - start).count();
}

// It also does not remove members that have full incidence.
vector<int> compute_maximal_indexes(const vector<set_t>& incidence) {
    if (incidence.empty()) {
        return {};
    }
    size_t num_members = incidence.size();
    int num_containers = set_size(incidence[0]);
    if (num_containers == 0) {
        return {};
    }

    vector<int> indexes_by_cardinality(num_members);
    vector<int> cardinalities(num_members);
    for (size_t i = 0; i < num_members; i++) {
        indexes_by_cardinality[i] = i;
        cardinalities[i] = set_count(incidence[i]);
    }

    sort(indexes_by_cardinality.begin(), indexes_by_cardinality.end(),
         [&cardinalities](int i, int j){return cardinalities[i] > cardinalities[j];});

    vector<int> maximal;
    int count_full_incidence = 0;
    for (int i : indexes_by_cardinality) {
        if (cardinalities[i] == num_containers) {
            count_full_incidence++;
            maximal.push_back(i);
            continue;
        }
        bool is_subset = false;
        const set_t inc_i = incidence[i];
        for (int j = count_full_incidence; j < (int) maximal.size(); j++) {
            if (set_is_subset_of(inc_i, incidence[maximal[j]])) {
                is_subset = true;
                break;
            }
        }
        if (!is_subset) {
            maximal.push_back(i);
        }
    }

    sort(maximal.begin(), maximal.end());
    return maximal;
}

// A simple wrapper that has error control.
dd_PolyhedraPtr cdd_Matrix_to_Poly(dd_MatrixPtr A) {
    dd_ErrorType err = dd_NoError;
    dd_PolyhedraPtr poly = dd_DDMatrix2Poly(A, &err);
    ASRTF(err == dd_NoError, "Converting matrix to polytope failed with error " + to_string(err));
    return poly;
}

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

int coef2index(const vector<int>& coef) {
    const int K = (int) coef.size();
    int index = 0;
    int offset = 1;
    for (int i = K - 1; i >= 0; i--) {
        if (coef[i] == 0) {
            index += offset;
        } else if (coef[i] == -1) {
            index += 2 * offset;
        }
        offset *= 3;
    }
    assert(index != offset / 2 && "Index cannot be in the middle because it means coefs are all zero.");
    if (index > offset / 2) {
        index--;
    }
    return index;
}

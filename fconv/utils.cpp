#include <iostream>
#include <algorithm>
#include <cassert>
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

vector<double *> cdd_compute_inequalities_from_vertices(dd_MatrixPtr vertices) {
  dd_ErrorType err = dd_NoError;
  dd_PolyhedraPtr poly = dd_DDMatrix2Poly(vertices, &err);
  ASRTF(err == dd_NoError,
        "Converting matrix to polytope failed with error " + to_string(err));

  dd_MatrixPtr inequalities = dd_CopyInequalities(poly);
  // Note that linearities are quite rare.
  set_type linearities = inequalities->linset;

  const int num_ineq = inequalities->rowsize;
  const int dim = inequalities->colsize;
  vector<double *> H = fp_mat_create(num_ineq + set_card(linearities), dim);
  int counter = 0;
  for (int ineq = 0; ineq < num_ineq; ineq++) {
    for (int j = 0; j < dim; j++) {
      H[counter][j] = mpq_get_d(inequalities->matrix[ineq][j]);
    }
    counter++;
    if (set_member(ineq + 1, linearities)) {
      for (int j = 0; j < dim; j++) {
        H[counter][j] = -H[counter - 1][j];
      }
      counter++;
    }
  }
  assert(counter == num_ineq + set_card(linearities) &&
         "Counter should equal num_ineq + number of linearities");

  for (auto h : H) {
    double abs_coef = 0;
    for (int i = 0; i < dim; i++) {
      abs_coef = max(abs(h[i]), abs_coef);
    }
    ASRTF(abs_coef != 0, "Inequality cannot consist fully of zeros.");
    for (int i = 0; i < dim; i++) {
      h[i] /= abs_coef;
    }
  }

  dd_FreePolyhedra(poly);
  dd_FreeMatrix(inequalities);
  return H;
}

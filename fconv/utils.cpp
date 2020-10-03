#include <iostream>
#include <algorithm>
#include "utils.h"
#include "setoper.h"
#include "cdd.h"
#include "dynamic_bitset.h"

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

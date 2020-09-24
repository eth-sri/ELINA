#include <vector>
#include <algorithm>
#include "sparse_cover.h"

#include "asrt.h"

using namespace std;

bool intersects_too_much(const int K, const vector<int>& first, const vector<int>& second) {
    int intersection_size = 0;

    int i = 0;
    int j = 0;

    while (i < K && j < K) {
        if (first[i] == second[j]) {
            i++;
            j++;
            intersection_size++;
        } else if (first[i] < second[j]) {
            i++;
        } else {
            j++;
        }
    }

    return intersection_size >= K - 1;
}

vector<vector<int>> sparse_cover(const int N, const int K) {
    ASRTF(3 <= K && K <= 5, "K is not within allowed range.");
    ASRTF(K <= N, "N is not within allowed range.");

    vector<vector<int>> all_selected_combs;
    vector<bool> v(N);
    fill(v.begin(), v.begin() + K, true);

    // This heuristic has to run fast, the main bottleneck is iterating through
    // the array of combinations to see if there is already a combination that has big
    // intersection with combination that we are trying to add.
    // The combinations are added in the alphabetical order, thus we can use indexing
    // to first traverse the part of the array where we already no that first elements match -
    // that area has a higher probability of rejecting new combination.
    vector<int> indexing_on_first(N, -1);

    // Improved heuristic that allows indexing based on second indexes.
    vector<vector<int>> indexing_on_second_all(N, vector<int>(N, -1));

    do {
        vector<int> new_comb(K);

        int cur_i = 0;
        for (int i = 0; i < N; ++i) {
            if (v[i]) {
                new_comb[cur_i] = i;
                cur_i++;
            }
        }

        bool to_add = true;

        const int first_elem = new_comb[0];
        if (indexing_on_first[first_elem] != -1) {
            int selected_comb_i = indexing_on_first[first_elem];
            while (selected_comb_i < (int) all_selected_combs.size()) {
                const auto &selected_comb = all_selected_combs[selected_comb_i];
                // We have left the promising region that allows to do fast rejection -
                // now we have to do the full traversal.
                if (selected_comb[0] != first_elem) {
                    break;
                }

                if (intersects_too_much(K, new_comb, selected_comb)) {
                    to_add = false;
                    break;
                }

                selected_comb_i++;
            }
        }

        // The combination was rejected based on first-element indexing heuristic.
        if (!to_add) {
            continue;
        }

        const int second_elem = new_comb[1];
        vector<int> &indexing_on_second = indexing_on_second_all[second_elem];
        // Doing a copy here intentionally.
        for (int selected_comb_i : indexing_on_second) {
            if (selected_comb_i == -1) {
                continue;
            }
            while (selected_comb_i < (int) all_selected_combs.size()) {
                const auto &selected_comb = all_selected_combs[selected_comb_i];
                // We have left the promising region that allows to do fast rejection -
                // now we have to do the full traversal.
                if (selected_comb[1] != second_elem) {
                    break;
                }

                if (intersects_too_much(K, new_comb, selected_comb)) {
                    to_add = false;
                    break;
                }

                selected_comb_i++;
            }
        }

        // The combination was rejected based on second-element indexing heuristic.
        if (!to_add) {
            continue;
        }

        // Doing full traversal.
        for (const auto &selected_comb : all_selected_combs) {
            // Both new_comb and selected_comb are sorted.
            // Now I compute their intersection size.
            if (intersects_too_much(K, new_comb, selected_comb)) {
                to_add = false;
                break;
            }
        }

        if (to_add) {
            if (indexing_on_first[first_elem] == -1) {
                indexing_on_first[first_elem] = (int) all_selected_combs.size();
            }
            // Update indexing on second element.
            // Note that having first_elem in the indexing is not a mistake.
            if (indexing_on_second[first_elem] == -1) {
                indexing_on_second[first_elem] = (int) all_selected_combs.size();
            }
            all_selected_combs.push_back(new_comb);
        }
    } while (prev_permutation(v.begin(), v.end()));

    return all_selected_combs;
}

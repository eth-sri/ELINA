#include <set>
#include "split_in_quadrants.h"
#include "utils.h"
#include "mpq.h"

vector<Adj> compute_adjacency_from_incidence(const vector<bset>& incidence) {
    size_t n = incidence.size();

    vector<Adj> potential_adjacencies;
    vector<bset> adjacencies_incidence;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            // A potential optimization here is analogue of Chernikova criteria.
            // Pair can create an adjacency, only if it's incident more to a certain number of constraints.
            potential_adjacencies.emplace_back(i, j);
            adjacencies_incidence.push_back(incidence[i] & incidence[j]);
        }
    }
    vector<int> maximal = compute_maximal_indexes(adjacencies_incidence);

    vector<Adj> maximal_adjacencies(maximal.size());
    for (size_t i = 0; i < maximal.size(); i++) {
        maximal_adjacencies[i] = potential_adjacencies[maximal[i]];
    }

    return maximal_adjacencies;
}

void split_in_quadrants_recursive(
        Quadrant& quadrant,
        const vector<Adj>& adjacencies,
        // Incidence here is V to H.
        vector<bset>& incidence,
        vector<mpq_t*>& vrep,
        map<Quadrant, QuadrantInfo>& quadrant2info,
        const int K,
        const int NUM_H) {
    const int xi = (int) quadrant.size();
    assert(incidence.size() == vrep.size() && "The size of incidence should equal the size of vrep.");

    if (xi == K) {
        quadrant2info[quadrant] = {K + 1, vrep, incidence};
        return;
    }

    const size_t num_v = vrep.size();

    vector<int> V_to_sign(num_v, 0); // -1 means xi < 0; 0 means xi == 0; 1 means xi > 0
    vector<int> V_to_minus(num_v, -1);
    vector<int> V_to_plus(num_v, -1);
    vector<int> V_zero;
    int count_minus = 0, count_plus = 0;

    for (size_t vi = 0; vi < num_v; vi++) {
        const mpq_t* v = vrep[vi];
        int sign = mpq_sgn(v[xi + 1]);
        if (sign <= 0) {
            V_to_minus[vi] = count_minus;
            count_minus++;
            if (sign < 0) {
                V_to_sign[vi] = -1;
            }
        }
        if (sign >= 0) {
            V_to_plus[vi] = count_plus;
            count_plus++;
            if (sign > 0) {
                V_to_sign[vi] = 1;
            }
        }
        if (sign == 0) {
            V_zero.push_back(vi);
            V_to_sign[vi] = 0;
            // It's okay to modify incidence because from it I will be creating two copies
            // that will be passed down to children recursively.
            assert(!incidence[vi].test(NUM_H + xi) &&
                  "Consistency check, the bit with this index could not be set.");
            incidence[vi].set(NUM_H + xi);
        }
    }

    // adj_split is a list of adjacencies that are split by this orthant - i.e. one vertex of the
    // adjacency is on one side, and another vertex is on the another side.
    vector<Adj> adj_split;
    vector<Adj> adj_minus;
    vector<Adj> adj_plus;
    // The procedure might lead to creating the duplicate adjacencies between existing vertices that are
    // incident to the new orthant. In order to avoid this problem, I record the adjacencies between such
    // vertices to ensure there is no duplication.
    set<Adj> adj_zero;
    for (const auto& adj : adjacencies) {
        assert((0 <= adj.first && adj.first < adj.second && adj.second < (int) num_v) &&
                "Checking the consistency of adj.");

        int sign1 = V_to_sign[adj.first];
        int sign2 = V_to_sign[adj.second];

        if (sign1 != 0 && sign2 != 0 && sign1 != sign2) {
            // This adjacency is guaranteed to create new vertex.
            adj_split.push_back(adj);
        }
        if (sign1 <= 0 && sign2 <= 0) {
            int new_first = V_to_minus[adj.first];
            int new_second = V_to_minus[adj.second];
            assert(new_first != -1 && new_second != -1 && "Consistency check for generating minus adj.");
            adj_minus.emplace_back(new_first, new_second);
        }
        // In that case adjacency goes to plus.
        if (sign1 >= 0 && sign2 >= 0) {
            int new_first = V_to_plus[adj.first];
            int new_second = V_to_plus[adj.second];
            assert(new_first != -1 && new_second != -1 && "Consistency check for generating plus adj.");
            adj_plus.emplace_back(new_first, new_second);
        }
        if (sign1 == 0 && sign2 == 0) {
            adj_zero.insert({adj.first, adj.second});
        }
    }

    const size_t num_new_v = adj_split.size();
    vector<mpq_t*> new_vrep(num_new_v);
    // Because I would need to compute adjacency between newly created vertices and vertices
    // that are incident to the new orthant, I create a list of incidences of these vertices.
    vector<bset> incidence_of_V_new_and_zero(num_new_v + V_zero.size(), bset(NUM_H + K));

    mpq_t delta, a_coef, b_coef, temp1, temp2;
    mpq_inits(delta, a_coef, b_coef, temp1, temp2, NULL);
    for (int adj_i = 0; adj_i < (int) adj_split.size(); adj_i++) {
        const Adj& adj = adj_split[adj_i];
        incidence_of_V_new_and_zero[adj_i] = incidence[adj.first] & incidence[adj.second];
        incidence_of_V_new_and_zero[adj_i].set(NUM_H + xi);
        const mpq_t* a = vrep[adj.first];
        const mpq_t* b = vrep[adj.second];
        new_vrep[adj_i] = mpq_create_array(K + 1);
        mpq_t* v = new_vrep[adj_i];

        int a_sign = mpq_sgn(a[xi + 1]);
        int b_sign = mpq_sgn(b[xi + 1]);
        ASRTF(((a_sign > 0 && b_sign < 0) || (a_sign < 0 && b_sign > 0)),
              "Checking consistency of adjacencies chosen for processing.");

        mpq_sub(delta, b[xi + 1], a[xi + 1]);
        assert(mpq_sgn(delta) != 0 && "Consistency check - delta can't be 0.");
        mpq_div(a_coef, b[xi + 1], delta);
        mpq_div(b_coef, a[xi + 1], delta);

        for (int j = 0; j < K + 1; j++) {
            mpq_mul(temp1, a_coef, a[j]);
            mpq_mul(temp2, b_coef, b[j]);
            mpq_sub(v[j], temp1, temp2);
        }

        // Basic correctness checks.
        assert(mpq_sgn(v[xi + 1]) == 0 && "The cur element of the new vertex have to be 0.");
        assert(mpq_cmp_si(v[0], 1, 1) == 0 && "First element of the new vertex have to be 1.");

        if (a_sign < 0) {
            // adj.first -> in negative branch, adj.second -> in positive.
            // count_minus + adj_i -> will be a shifted index of a generated point.
            assert(V_to_minus[adj.first] != -1 && "Consistency check when creating new adj for minus.");
            assert(V_to_plus[adj.second] != -1 && "Consistency check when creating new adj for plus.");
            adj_minus.emplace_back(V_to_minus[adj.first], count_minus + adj_i);
            adj_plus.emplace_back(V_to_plus[adj.second], count_plus + adj_i);
        } else {
            // adj.first -> in positive branch, adj.second -> in negative.
            assert(V_to_minus[adj.second] != -1 && "Consistency check when creating new adj for minus.");
            assert(V_to_plus[adj.first] != -1 && "Consistency check when creating new adj for plus.");
            adj_minus.emplace_back(V_to_minus[adj.second], count_minus + adj_i);
            adj_plus.emplace_back(V_to_plus[adj.first], count_plus + adj_i);
        }
    }
    mpq_clears(delta, a_coef, b_coef, temp1, temp2, NULL);

    for (size_t i = 0; i < V_zero.size(); i++) {
        incidence_of_V_new_and_zero[num_new_v + i] = incidence[V_zero[i]];
    }

    vector<Adj> new_adjacencies = compute_adjacency_from_incidence(incidence_of_V_new_and_zero);

    for (const auto& adj : new_adjacencies) {
        int minus_first = -1, minus_second = -1, plus_first = -1, plus_second = -1;
        mpq_t* v_first = NULL;
        mpq_t* v_second = NULL;
        if (adj.first >= (int) num_new_v && adj.second >= (int) num_new_v) {
            // If both vertices of the adjacency are old vertices that are incident to the new
            // orthant I need to check if I have already enumerated their adjacency.
            // If that's the case, adjacency has to be skipped to avoid duplication.
            int vi_1 = V_zero[adj.first - (int) adj_split.size()];
            int vi_2 = V_zero[adj.second - (int) adj_split.size()];
            if (adj_zero.count({vi_1, vi_2}) > 0) {
                // This adjacency has already been recorded.
                // Thus to avoid enumerating it two times, we skip it.
                continue;
            }
        }

        if (adj.first < (int) num_new_v) {
            minus_first = count_minus + adj.first;
            plus_first = count_plus + adj.first;
            v_first = new_vrep[adj.first];
        } else {
            int vi = V_zero[adj.first - (int) num_new_v];
            minus_first = V_to_minus[vi];
            plus_first = V_to_plus[vi];
            v_first = vrep[vi];
        }
        if (adj.second < (int) num_new_v) {
            minus_second = count_minus + adj.second;
            plus_second = count_plus + adj.second;
            v_second = new_vrep[adj.second];
        } else {
            int vi = V_zero[adj.second - (int) num_new_v];
            minus_second = V_to_minus[vi];
            plus_second = V_to_plus[vi];
            v_second = vrep[vi];
        }

        bool useful = false;
        for (int i = xi + 1; i < K; i++) {
            int sign_h1 = mpq_sgn(v_first[i + 1]);
            int sign_h2 = mpq_sgn(v_second[i + 1]);
            if (sign_h1 != 0 && sign_h2 != 0 && sign_h1 != sign_h2) {
                // Adjacency is useful only if it is to be split by one of the future orthants.
                useful = true;
                break;
            }
        }
        if (!useful) {
            continue;
        }

        adj_minus.emplace_back(minus_first, minus_second);
        adj_plus.emplace_back(plus_first, plus_second);
    }

    vector<mpq_t*> vrep_minus(count_minus + num_new_v);
    vector<mpq_t*> vrep_plus(count_plus + num_new_v);
    vector<bset> incidence_minus(count_minus + num_new_v);
    vector<bset> incidence_plus(count_plus + num_new_v);
    for (size_t vi = 0; vi < num_v; vi++) {
        int sign = V_to_sign[vi];
        if (sign <= 0) {
            int minus_i = V_to_minus[vi];
            // The minus branch assumes ownership - note that no copy is happening.
            vrep_minus[minus_i] = vrep[vi];
            incidence_minus[minus_i] = incidence[vi];
        }
        if (sign >= 0) {
            int plus_i = V_to_plus[vi];
            if (sign > 0) {
                // If vertex goes exclusively to the plus branch - then plus branch assumes ownership.
                // And no copy happens.
                vrep_plus[plus_i] = vrep[vi];
            } else {
                vrep_plus[plus_i] = mpq_create_and_copy_array(K + 1, vrep[vi]);
            }
            incidence_plus[plus_i] = incidence[vi];
        }
    }

    for (size_t vi = 0; vi < num_new_v; vi++) {
        // Minus branch takes ownership.
        vrep_minus[count_minus + vi] = new_vrep[vi];
        incidence_minus[count_minus + vi] = incidence_of_V_new_and_zero[vi];
        // Plus branch copies.
        vrep_plus[count_plus + vi] = mpq_create_and_copy_array(K + 1, new_vrep[vi]);
        incidence_plus[count_plus + vi] = incidence_of_V_new_and_zero[vi];
    }

    quadrant.push_back(MINUS);
    split_in_quadrants_recursive(quadrant, adj_minus, incidence_minus, vrep_minus, quadrant2info, K, NUM_H);
    quadrant.back() = PLUS;
    split_in_quadrants_recursive(quadrant, adj_plus, incidence_plus, vrep_plus, quadrant2info, K, NUM_H);
    quadrant.pop_back();
}

map<Quadrant, QuadrantInfo> split_in_quadrants(vector<mpq_t*>& V,
                                               vector<bset>& incidence,
                                               const vector<Adj>& orthant_adjacencies,
                                               const int K) {
    ASRTF(2 <= K && K <= 5, "K should be within allowed range.");
    ASRTF(!V.empty(), "Assuming that the input is not empty.");
    ASRTF(V.size() == incidence.size(), "Number of vertexes should equal number of incidences.");

    const int NUM_H = (int) incidence[0].size();
    for (auto& inc : incidence) {
        inc.resize(NUM_H + K);
    }

    Quadrant quadrant {};
    quadrant.reserve(K);
    map<Quadrant, QuadrantInfo> quadrant2info;
    split_in_quadrants_recursive(quadrant, orthant_adjacencies, incidence, V, quadrant2info, K, NUM_H);

    return quadrant2info;
}

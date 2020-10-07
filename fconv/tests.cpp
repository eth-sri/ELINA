#include <algorithm>
#include <string>
#include <iostream>
#include <set>
#include <map>
#include <unistd.h>
#include "octahedron.h"
#include "utils.h"
#include "mpq.h"
#include "split_in_quadrants.h"
#include "relaxation.h"
#include "sparse_cover.h"
#include "fp_mat.h"

// Temporary disabled the last test for k=4 before I add proper support
// for inputs that split zero.
constexpr int K2NUM_TESTS[6] = {0, 2, 4, 5, 3, 4};

enum Activation {
    Relu,
    Pool,
    Tanh
};

// Creates bijective mapping between the rows of two matrices.
// Will throw an exception if such mapping does not exist.
vector<int> get_bijective_mapping_for_matrix_rows(
        const int dim,
        const vector<mpq_t*>& first,
        const vector<mpq_t*>& second) {
    ASRTF(first.size() == second.size(), "Num rows should match.");
    const int num_v = (int) first.size();
    vector<int> mapping_first_to_second(num_v, -1);
    vector<int> mapping_second_to_first(num_v, -1);
    for (int fi = 0; fi < (int) num_v; fi++) {
        mpq_t* vertex_f = first[fi];
        int mapped_v = -1;
        for (int si = 0; si < num_v; si++) {
            if (mpq_arr_equal(dim, vertex_f, second[si])) {
                ASRTF(mapped_v == -1, "Vertex can only be mapped to one vertex.");
                ASRTF(mapping_second_to_first[si] == -1, "Has to be one-to-one correspondence.");
                mapped_v = si;
            }
        }
        ASRTF(mapped_v != -1, "Vertex should be mapped to at least one vertex.");
        mapping_first_to_second[fi] = mapped_v;
        mapping_second_to_first[mapped_v] = fi;
    }
    return mapping_first_to_second;
}

void run_octahedron_test(const int K, const string& path) {
    cout << "running octahedron test: " << path << endl;

    vector<double*> A = fp_mat_read(K + 1, path);

    Timer t_fast;
    OctahedronV fast = compute_octahedron_V(K, A);
    int micros_fast = t_fast.micros();
    Timer t_slow;
    OctahedronV slow = compute_V_with_cdd(K, A);
    int micros_slow = t_slow.micros();

    cout << "\tAcceleration is " << (double) micros_slow / micros_fast << endl;
    cout << "\tFrom " << micros_slow / 1000 << " ms to " << micros_fast / 1000 << " ms" << endl;

    vector<int> mapping = get_bijective_mapping_for_matrix_rows(K + 1, fast.V, slow.V);

    set<Adj> adjacencies_set(fast.orthant_adjacencies.begin(), fast.orthant_adjacencies.end());

    // Comparing adjacencies
    ASRTF(fast.orthant_adjacencies.size() == slow.orthant_adjacencies.size(), "The size of adjacencies should match.");
    vector<Adj> adjacencies_mapped;
    for (const auto& adj : fast.orthant_adjacencies) {
        int first = min(mapping[adj.first], mapping[adj.second]);
        int second = max(mapping[adj.first], mapping[adj.second]);
        adjacencies_mapped.emplace_back(first, second);
    }
    // Both slow_V and fast_V return adjacencies in the sorted order. So comparing them is easy.
    sort(adjacencies_mapped.begin(), adjacencies_mapped.end());

    for (size_t i = 0; i < fast.orthant_adjacencies.size(); i++) {
        Adj adj_fast = adjacencies_mapped[i];
        Adj adj_slow = slow.orthant_adjacencies[i];

        ASRTF(adj_fast.first == adj_slow.first && adj_fast.second == adj_slow.second, "Adjacency should equal.");
    }

    // Comparing incidence
    ASRTF(fast.incidence.size() == slow.incidence.size(),
          "The size of incidence should match.");

    for (size_t i = 0; i < fast.incidence.size(); i++) {
        // I intentionally do a copy here so I can modify it.
        const set_t inc_fast = fast.incidence[i];
        const set_t inc_slow = slow.incidence[mapping[i]];
        ASRTF(set_size(inc_fast) == set_size(inc_slow), "Sizes of incidence should match.");
        ASRTF(set_equal(inc_fast, inc_slow), "Incidence should match.");
    }

    set_arr_free(fast.incidence);
    set_arr_free(slow.incidence);
    mpq_mat_free(K + 1, fast.V);
    mpq_mat_free(K + 1, slow.V);
    fp_mat_free(A);

    cout << "\tpassed" << endl;
}

void run_split_in_quadrants_test(const int K, const string& path) {
    cout << "running split in quadrants test: " << path << endl;

    vector<double*> A = fp_mat_read(K + 1, path);
    const int num_h = (int) A.size();

    Timer t_fast;
    OctahedronV oct_fast = compute_octahedron_V(K, A);
    map<Quadrant, QuadrantInfo> fast = split_in_quadrants(
            oct_fast.V,
            oct_fast.incidence,
            oct_fast.orthant_adjacencies,
            K);
    int micros_fast = t_fast.micros();

    Timer t_slow;
    map<Quadrant, QuadrantInfo> slow = compute_quadrants_with_cdd(K, A);
    int micros_slow = t_slow.micros();

    cout << "\tAcceleration is " << (double) micros_slow / micros_fast << endl;
    cout << "\tFrom " << micros_slow / 1000 << " ms to " << micros_fast / 1000 << " ms" << endl;

    for (auto& p : slow) {
        QuadrantInfo mine = fast.find(p.first)->second;
        QuadrantInfo gt = p.second;

        ASRTF(mine.V.size() == gt.V.size(), "Sizes should match.");

        // If the function completes - bijective mapping exists - set of vertices is equal.
        vector<int> mapping = get_bijective_mapping_for_matrix_rows(K + 1, mine.V, gt.V);

        const size_t num_v = mine.V.size();
        const vector<set_t>& gt_incidence = gt.V_to_H_incidence;
        const vector<set_t>& mine_incidence = mine.V_to_H_incidence;
        ASRTF(gt_incidence.size() == num_v, "The size of incidence should equal num_v.");
        ASRTF(mine_incidence.size() == num_v, "The size of incidence should equal num_v.");

        for (size_t i = 0; i < num_v; i++) {
            const set_t inc_fast = mine_incidence[i];
            const set_t inc_slow = gt_incidence[mapping[i]];
            ASRTF(set_size(inc_fast) == num_h + K, "Size of inc_fast should be num_h + K.");
            ASRTF(set_size(inc_slow) == num_h + K, "Size of inc_slow should be num_h + K.");
            ASRTF(set_equal(inc_fast, inc_slow), "Incidence should be equal.");
        }
    }

    for (auto& entry : fast) {
        mpq_mat_free(K + 1, entry.second.V);
        set_arr_free(entry.second.V_to_H_incidence);
    }
    for (auto& entry : slow) {
        mpq_mat_free(K + 1, entry.second.V);
        set_arr_free(entry.second.V_to_H_incidence);
    }
    fp_mat_free(A);

    cout << "\tpassed" << endl;
}

void run_fkrelu_test(const int K, const string& path) {
    cout << "running fkrelu test: " << path << endl;

    vector<double*> A = fp_mat_read(K + 1, path);

    Timer t;
    vector<double*> H_double = fkrelu(K, A);
    int micros = t.micros();

    cout << "\tK = " << K << " took " << micros / 1000 << \
        " ms and discovered " << H_double.size() << " constraints" << endl;

    vector<mpq_t*> H = mpq_mat_from_fp(2 * K + 1, H_double);

    map<Quadrant, QuadrantInfo> quadrant2info = compute_quadrants_with_cdd(K, A);

    // Now I will compute the final V. This will allow me to verify that produced constraints do not
    // violate any of the original vertices and thus I produce a sound overapproximation.
    vector<mpq_t*> V;
    for (auto& entry : quadrant2info) {
        const auto& quadrant = entry.first;
        auto& V_quadrant = entry.second.V;

        for (auto v : V_quadrant) {
            v = mpq_arr_resize(2 * K + 1, K + 1, v);
            for (int i = 0; i < K; i++) {
                if (quadrant[i] == PLUS) {
                    // Only need to set for the case of PLUS,
                    // because in case of MINUS there should be 0.
                    // And it is already there automatically.
                    mpq_set(v[1 + i + K], v[1 + i]);
                }
            }
            V.push_back(v);
        }
        set_arr_free(entry.second.V_to_H_incidence);
    }

    vector<mpq_t*> V_x_H = mpq_mat_mul_with_transpose(2 * K + 1, V, H);
    for (const auto& row : V_x_H) {
        for (size_t i = 0; i < H.size(); i++) {
            ASRTF(mpq_sgn(row[i]) >= 0, "All discovered constraints should be sound with respect to V");
        }
    }

    mpq_mat_free(2 * K + 1, H);
    mpq_mat_free(2 * K + 1, V);
    mpq_mat_free(H.size(), V_x_H);
    fp_mat_free(A);
    fp_mat_free(H_double);

    cout << "\tpassed" << endl;
}

void run_fkpool_test(const int K, const string& path) {
    cout << "running fkpool test: " << path << endl;

    vector<double*> A = fp_mat_read(K + 1, path);

    Timer t;
    vector<double*> H_double = fkpool(K, A);
    int micros = t.micros();

    cout << "\tK = " << K << " took " << micros / 1000 << \
        " ms and discovered " << H_double.size() << " constraints" << endl;

    vector<mpq_t*> H = mpq_mat_from_fp(K + 2, H_double);

    vector<QuadrantInfo> quadrant_infos = compute_max_pool_quadrants_with_cdd(K, A);

    // Now I will compute the final V. This will allow me to verify that produced constraints do not
    // violate any of the original vertices and thus I produce a sound overapproximation.
    vector<mpq_t*> V;
    for (int xi = 0; xi < K; xi++) {
        auto& V_quadrant = quadrant_infos[xi].V;
        for (auto v : V_quadrant) {
            v = mpq_arr_resize(K + 2, K + 1, v);
            mpq_set(v[K + 1], v[xi + 1]);
            V.push_back(v);
        }
        set_arr_free(quadrant_infos[xi].V_to_H_incidence);
    }

    vector<mpq_t*> V_x_H = mpq_mat_mul_with_transpose(K + 2, V, H);
    for (const auto& row : V_x_H) {
        for (size_t i = 0; i < H.size(); i++) {
            ASRTF(mpq_sgn(row[i]) >= 0, "All discovered constraints should be sound with respect to V");
        }
    }

    mpq_mat_free(K + 2, H);
    mpq_mat_free(K + 2, V);
    mpq_mat_free(H.size(), V_x_H);
    fp_mat_free(A);
    fp_mat_free(H_double);

    cout << "\tpassed" << endl;
}

// The test doesn't check the correctness, it showcases the number of constraints and runtime.
void run_relaxation_with_cdd_test(const int K, const string& path, Activation activation) {
    cout << "running relaxation with cdd test: " << path << endl;

    vector<double*> A = fp_mat_read(K + 1, path);

    Timer t;
    vector<double*> H;
    switch (activation) {
        case Relu:
            cout << "\tRelu" << endl;
            H = krelu_with_cdd(K, A);
            break;
        case Pool:
            cout << "\tPool" << endl;
            H = kpool_with_cdd(K, A);
            break;
        case Tanh:
            cout << "\tTanh" << endl;
            H = ktanh_with_cdd(K, A);
            break;
        default:
            throw runtime_error("Unknown activation.");
    }
    int micros = t.micros();

    cout << "\tK = " << K << " took " << micros / 1000 << \
        " ms and discovered " << H.size() << " constraints" << endl;

    fp_mat_free(A);
    fp_mat_free(H);

    cout << "\tpassed" << endl;
}

void run_1relu_test() {
    cout << "running 1-relu test:" << endl;
    vector<double*> inp = fp_mat_create(2, 2);

    // x >= -2
    inp[0][0] = 2;
    inp[0][1] = 1;

    // x <= 2
    inp[1][0] = 2;
    inp[1][1] = -1;

    vector<double*> expected = fp_mat_create(3, 3);

    // y >= 0
    expected[0][0] = 0;
    expected[0][1] = 0;
    expected[0][2] = 1;

    // y >= x
    expected[1][0] = 0;
    expected[1][1] = -1;
    expected[1][2] = 1;

    // y <= 1 + 1/2 * x
    expected[2][0] = 1;
    expected[2][1] = 0.5;
    expected[2][2] = -1;

    vector<double*> out = fkrelu(1, inp);

    ASRTF(out.size() == 3, "Expected that output has 3 rows.");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ASRTF(out[i][j] == expected[i][j], "Actual and expected elements match.");
        }
    }
    fp_mat_free(inp);
    fp_mat_free(expected);
    fp_mat_free(out);

    cout << "\tpassed" << endl;
}

void run_sparse_cover_test(const int N, const int K) {
    cout << "running sparse cover test: N " << N << " K " << K << endl;

    Timer t;
    vector<vector<int>> cover = sparse_cover(N, K);
    int micros = t.micros();

    cout << "\ttook " << micros / 1000 \
        << " ms and generated " << cover.size() \
        << " combinations" << endl;

    for (const auto& comb : cover) {
        ASRTF((int) comb.size() == K, "Combination size should equal K.");
        for (auto elem : comb) {
            ASRTF(0 <= elem && elem < N, "Every element of combination should be within range.");
        }
        for (int i = 0; i < K - 1; i++) {
            ASRTF(comb[i] < comb[i + 1], "Element of combination should be in ascending order.");
        }
    }

    cout << "\tpassed" << endl;
}

void run_tanh_quadrants_with_cdd_dimension_trick_test(const int K, const string& path) {
    cout << "running tanh quadrants with cdd dimension trick test: " << path << endl;

    vector<double*> A = fp_mat_read(K + 1, path);

    Timer t_fast;
    map<Quadrant, vector<mpq_t*>> quadrant2vertices_fast =
            compute_tanh_quadrants_with_cdd_dimension_trick(K, A);
    int micros_fast = t_fast.micros();

    Timer t_slow;
    map<Quadrant, vector<mpq_t*>> quadrant2vertices_slow =
            compute_tanh_quadrants_with_cdd(K, A);
    int micros_slow = t_slow.micros();

    cout << "\tAcceleration is " << (double) micros_slow / micros_fast << endl;
    cout << "\tFrom " << micros_slow / 1000 << " ms to " << micros_fast / 1000 << " ms" << endl;

    const vector<Quadrant>& quadrants = K2QUADRANTS[K];

    for (const auto& q : quadrants) {
        const vector<mpq_t*>& V_fast = quadrant2vertices_fast[q];
        const vector<mpq_t*>& V_slow = quadrant2vertices_slow[q];
        get_bijective_mapping_for_matrix_rows(2 * K + 1, V_fast, V_slow);
        mpq_mat_free(2 * K + 1, V_fast);
        mpq_mat_free(2 * K + 1, V_slow);
    }

    cout << "\tpassed" << endl;
}


void run_all_octahedron_tests() {
    cout << "Running all fast V octahedron tests" << endl;
    for (int k = 2; k <= 4; k++) {
        for (int i = 1; i <= K2NUM_TESTS[k]; i++) {
            run_octahedron_test(
                    k,
                    "octahedron_hrep/k" + to_string(k) + "/" + to_string(i) + ".txt");
        }
    }
}

void run_all_split_in_quadrants_tests() {
    cout << "Running all split in quadrant tests" << endl;
    for (int k = 2; k <= 4; k++) {
        for (int i = 1; i <= K2NUM_TESTS[k]; i++) {
            run_split_in_quadrants_test(
                    k,
                    "octahedron_hrep/k" + to_string(k) + "/" + to_string(i) + ".txt");
        }
    }
}

void run_all_fkrelu_tests() {
    cout << "Running all fkrelu tests" << endl;
    for (int k = 1; k <= 4; k++) {
        for (int i = 1; i <= K2NUM_TESTS[k]; i++) {
            run_fkrelu_test(
                    k,
                    "octahedron_hrep/k" + to_string(k) + "/" + to_string(i) + ".txt");
        }
    }
}

void run_all_fkpool_tests() {
    cout << "Running all fkpool tests" << endl;
    for (int k = 1; k <= 4; k++) {
        for (int i = 1; i <= K2NUM_TESTS[k]; i++) {
            run_fkpool_test(
                    k,
                    "octahedron_hrep/k" + to_string(k) + "/" + to_string(i) + ".txt");
        }
    }
}

void run_all_krelu_with_cdd_tests() {
    cout << "Running all krelu with cdd tests" << endl;
    for (int k = 1; k <= 3; k++) {
        for (int i = 1; i <= K2NUM_TESTS[k]; i++) {
            run_relaxation_with_cdd_test(
                    k,
                    "octahedron_hrep/k" + to_string(k) + "/" + to_string(i) + ".txt",
                    Relu);
        }
    }
}

void run_all_kpool_with_cdd_tests() {
    cout << "Running all kpool with cdd tests" << endl;
    // kpool with cdd doesn't scale to k=4 (takes 1-2 minutes).
    for (int k = 1; k <= 3; k++) {
        for (int i = 1; i <= K2NUM_TESTS[k]; i++) {
            run_relaxation_with_cdd_test(
                    k,
                    "octahedron_hrep/k" + to_string(k) + "/" + to_string(i) + ".txt",
                    Pool);
        }
    }
}

void run_all_ktanh_with_cdd_tests() {
    cout << "Running all ktanh with cdd tests" << endl;
    // kpool with cdd doesn't scale to k=3 (takes 1-2 minutes).
    for (int k = 1; k <= 2; k++) {
        for (int i = 1; i <= K2NUM_TESTS[k]; i++) {
            run_relaxation_with_cdd_test(
                    k,
                    "octahedron_hrep/k" + to_string(k) + "/" + to_string(i) + ".txt",
                    Tanh);
        }
    }
}

void run_all_sparse_cover_tests() {
    cout << "Running all sparse cover tests" << endl;
    run_sparse_cover_test(50, 3);
    run_sparse_cover_test(100, 3);
    run_sparse_cover_test(25, 4);
    run_sparse_cover_test(20, 5);
}

void run_all_tanh_quadrants_with_cdd_dimension_trick_tests() {
    cout << "Running all tanh quadrant with cdd dimension trick tests" << endl;
    // k = 4 is somewhat slow (several seconds), so doing tests until k = 3.
    for (int k = 1; k <= 3; k++) {
        for (int i = 1; i <= K2NUM_TESTS[k]; i++) {
            run_tanh_quadrants_with_cdd_dimension_trick_test(
                    k,
                    "octahedron_hrep/k" + to_string(k) + "/" + to_string(i) + ".txt");
        }
    }
}


void handler(int sig) {
    void *array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

int main() {
    signal(SIGSEGV, handler);
    dd_set_global_constants();

    run_all_octahedron_tests();
    run_all_split_in_quadrants_tests();
    run_all_fkrelu_tests();
    run_all_krelu_with_cdd_tests();
    run_1relu_test();
    run_all_fkpool_tests();
    run_all_kpool_with_cdd_tests();
    run_all_ktanh_with_cdd_tests();
    run_all_tanh_quadrants_with_cdd_dimension_trick_tests();
    run_all_sparse_cover_tests();

    dd_free_global_constants();
    return 0;
}

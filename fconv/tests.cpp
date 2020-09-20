#include "fkrelu.h"
#include "mpq.h"
#include "octahedron.h"
#include "split_in_quadrants.h"
#include "utils.h"
#include <iostream>
#include <map>
#include <set>
#include <string>

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
          if (mpq_arrays_are_equal(dim, vertex_f, second[si])) {
            ASRTF(mapped_v == -1, "Vertex can only be mapped to one vertex.");
            ASRTF(mapping_second_to_first[si] == -1,
                  "Has to be one-to-one correspondence.");
            mapped_v = si;
          }
        }
        ASRTF(mapped_v != -1, "Vertex should be mapped to at least one vertex.");
        mapping_first_to_second[fi] = mapped_v;
        mapping_second_to_first[mapped_v] = fi;
    }
    return mapping_first_to_second;
}

void run_octahedron_test(const string &path) {
  cout << "running octahedron test: " << path << endl;
  MatrixXd A = read_matrix(path);
  const int K = (int)A.cols() - 1;
  Timer t_fast;
  OctahedronV fast = compute_octahedron_V(A);
  int micros_fast = t_fast.micros();
  Timer t_slow;
  OctahedronV slow = compute_V_with_cdd(A);
  int micros_slow = t_slow.micros();

  cout << "\tAcceleration is " << (double)micros_slow / micros_fast << endl;
  cout << "\tFrom " << micros_slow / 1000 << " ms to " << micros_fast / 1000
       << " ms" << endl;

  vector<int> mapping =
      get_bijective_mapping_for_matrix_rows(K + 1, fast.V, slow.V);

  set<Adj> adjacencies_set(fast.orthant_adjacencies.begin(),
                           fast.orthant_adjacencies.end());

  // Comparing adjacencies
  ASRTF(fast.orthant_adjacencies.size() == slow.orthant_adjacencies.size(),
        "The size of adjacencies should match.");
  vector<Adj> adjacencies_mapped;
  for (const auto &adj : fast.orthant_adjacencies) {
    int first = min(mapping[adj.first], mapping[adj.second]);
    int second = max(mapping[adj.first], mapping[adj.second]);
    adjacencies_mapped.emplace_back(first, second);
  }
  // Both slow_V and fast_V return adjacencies in the sorted order. So comparing
  // them is easy.
  sort(adjacencies_mapped.begin(), adjacencies_mapped.end());

  for (size_t i = 0; i < fast.orthant_adjacencies.size(); i++) {
    Adj adj_fast = adjacencies_mapped[i];
    Adj adj_slow = slow.orthant_adjacencies[i];

    ASRTF(adj_fast.first == adj_slow.first &&
              adj_fast.second == adj_slow.second,
          "Adjacency should equal.");
  }

  // Comparing incidence
  ASRTF(fast.incidence.size() == slow.incidence.size(),
        "The size of incidence should match.");

  for (size_t i = 0; i < fast.incidence.size(); i++) {
    // I intentionally do a copy here so I can modify it.
    const bset &inc_fast = fast.incidence[i];
    const bset &inc_slow = slow.incidence[mapping[i]];
    ASRTF(inc_fast.size() == inc_slow.size(),
          "Sizes of incidence should match.");
    ASRTF(inc_fast == inc_slow, "Incidence should match.");
  }

  mpq_free_array_vector(K + 1, fast.V);
  mpq_free_array_vector(K + 1, slow.V);

  cout << "\tpassed" << endl;
}

void run_split_in_quadrants_test(const string &path) {
  cout << "running split in quadrants test: " << path << endl;

  MatrixXd A = read_matrix(path);
  const int K = (int)A.cols() - 1;
  const int num_h = (int)A.rows();

  Timer t_fast;
  OctahedronV oct_fast = compute_octahedron_V(A);
  map<Quadrant, QuadrantInfo> fast = split_in_quadrants(
      oct_fast.V, oct_fast.incidence, oct_fast.orthant_adjacencies, K);
  int micros_fast = t_fast.micros();

  Timer t_slow;
  map<Quadrant, QuadrantInfo> slow = compute_quadrants_with_cdd(A);
  int micros_slow = t_slow.micros();

  cout << "\tAcceleration is " << (double)micros_slow / micros_fast << endl;
  cout << "\tFrom " << micros_slow / 1000 << " ms to " << micros_fast / 1000
       << " ms" << endl;

  for (auto &p : slow) {
    QuadrantInfo mine = fast.find(p.first)->second;
    QuadrantInfo gt = p.second;

    ASRTF(mine.V.size() == gt.V.size(), "Sizes should match.");

    // If the function completes - bijective mapping exists - set of vertices is
    // equal.
    vector<int> mapping =
        get_bijective_mapping_for_matrix_rows(K + 1, mine.V, gt.V);

    const size_t num_v = mine.V.size();
    const vector<bset> &gt_incidence = gt.V_to_H_incidence;
    const vector<bset> &mine_incidence = mine.V_to_H_incidence;
    ASRTF(gt_incidence.size() == num_v,
          "The size of incidence should equal num_v.");
    ASRTF(mine_incidence.size() == num_v,
          "The size of incidence should equal num_v.");

    for (size_t i = 0; i < num_v; i++) {
      // I intentionally do a copy here so I can modify it.
      const bset &inc_fast = mine_incidence[i];
      const bset &inc_slow = gt_incidence[mapping[i]];
      ASRTF((int)inc_fast.size() == num_h + K,
            "Size of inc_fast should be num_h + K.");
      ASRTF((int)inc_slow.size() == num_h + K,
            "Size of inc_slow should be num_h + K.");
      ASRTF(inc_fast == inc_slow, "Incidence should be equal.");
    }
  }

  for (auto &entry : fast) {
    mpq_free_array_vector(K + 1, entry.second.V);
  }
  for (auto &entry : slow) {
    mpq_free_array_vector(K + 1, entry.second.V);
  }

  cout << "\tpassed" << endl;
}

void run_fkrelu_test(const string &path) {
  cout << "running fkrelu test: " << path << endl;

  MatrixXd A = read_matrix(path);

  Timer t;
  MatrixXd H_double = fkrelu(A);
  int micros = t.micros();

  const int K = (int)A.cols() - 1;
  cout << "\tK = " << K << " took " << micros / 1000 << " ms" << endl;

  vector<mpq_t *> H = mpq_from_eigen(H_double);

  map<Quadrant, QuadrantInfo> quadrant2info = compute_quadrants_with_cdd(A);

  // Now I will compute the final V. This will allow me to verify that produced
  // constraints do not violate any of the original vertices and thus I produce
  // a sound overapproximation.
  vector<mpq_t *> V;
  for (auto &entry : quadrant2info) {
    const auto &quadrant = entry.first;
    auto &V_quadrant = entry.second.V;

    for (const auto &v : V_quadrant) {
      mpq_t *v_projection = mpq_create_array(2 * K + 1);
      for (int i = 0; i < K + 1; i++) {
        mpq_set(v_projection[i], v[i]);
      }
      for (int i = 0; i < K; i++) {
        if (quadrant[i] == PLUS) {
          // Only need to set for the case of PLUS, because in case of MINUS
          // there should be 0. And it is already there automatically.
          mpq_set(v_projection[1 + i + K], v[1 + i]);
        }
      }
      V.push_back(v_projection);
    }

    mpq_free_array_vector(K + 1, V_quadrant);
  }

  vector<mpq_t *> V_x_H = mpq_matrix_mul(2 * K + 1, V, H);
  for (const auto &row : V_x_H) {
    for (size_t i = 0; i < H.size(); i++) {
      ASRTF(mpq_sgn(row[i]) >= 0,
            "All discovered constraints should be sound with respect to V");
    }
  }

  mpq_free_array_vector(2 * K + 1, H);
  mpq_free_array_vector(2 * K + 1, V);
  mpq_free_array_vector(H.size(), V_x_H);

  cout << "\tpassed" << endl;
}

void run_all_octahedron_tests() {
    cout << "Running all fast V octahedron tests" << endl;
    for (int k = 2; k <= 4; k++) {
      int num_tests = k == 3 ? 5 : 4;
      for (int i = 1; i <= num_tests; i++) {
        run_octahedron_test("octahedron_hrep/k" + to_string(k) + "/" +
                            to_string(i) + ".txt");
      }
    }
}

void run_all_split_in_quadrants_tests() {
    cout << "Running all split in quadrant tests" << endl;
    for (int k = 2; k <= 4; k++) {
      int num_tests = k == 3 ? 5 : 4;
      for (int i = 1; i <= num_tests; i++) {
        run_split_in_quadrants_test("octahedron_hrep/k" + to_string(k) + "/" +
                                    to_string(i) + ".txt");
      }
    }
}

void run_all_fkrelu_tests() {
    cout << "Running all fkrelu tests" << endl;
    for (int k = 2; k <= 4; k++) {
      int num_tests = k == 3 ? 5 : 4;
      for (int i = 1; i <= num_tests; i++) {
        run_fkrelu_test("octahedron_hrep/k" + to_string(k) + "/" +
                        to_string(i) + ".txt");
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

    dd_free_global_constants();
    return 0;
}

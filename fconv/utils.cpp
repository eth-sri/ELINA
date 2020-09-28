#include "utils.h"
#include "cdd.h"
#include "setoper.h"
#include <iostream>
#include <numeric>

using namespace std;

vector<double *> create_mat(const int rows, const int cols) {
  vector<double *> mat(rows);
  for (int i = 0; i < rows; i++) {
    mat[i] = (double *)calloc(cols, sizeof(double));
  }
  return mat;
}

void free_mat(const vector<double *> &mat) {
  for (auto v : mat) {
    free(v);
  }
}

Timer::Timer() {
    start = chrono::high_resolution_clock::now();
}

int Timer::micros() {
    auto now = chrono::high_resolution_clock::now();
    return chrono::duration_cast<chrono::microseconds>(now - start).count();
}

// It also does not remove members that have full incidence.
vector<int> compute_maximal_indexes(const vector<bset> &incidence) {
  if (incidence.empty()) {
    return {};
  }
  size_t num_members = incidence.size();
  size_t num_containers = incidence[0].size();
  if (num_containers == 0) {
    return {};
  }
  vector<int> indexes_by_cardinality(num_members);
  iota(indexes_by_cardinality.begin(), indexes_by_cardinality.end(), 0);

  // TODO[gleb] If the count() function is expensive - there is a potential
  // optimization.
  sort(indexes_by_cardinality.begin(), indexes_by_cardinality.end(),
       [&incidence](int i, int j) {
         return incidence[i].count() > incidence[j].count();
       });

  vector<int> maximal;
  int count_full_incidence = 0;
  for (int i : indexes_by_cardinality) {
    const bset &inc_i = incidence[i];
    if (inc_i.count() == num_containers) {
      count_full_incidence++;
      maximal.push_back(i);
      continue;
    }
    bool is_subset = false;
    for (int j = count_full_incidence; j < (int)maximal.size(); j++) {
      if (inc_i.is_subset_of(incidence[maximal[j]])) {
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

vector<double *> read_matrix(int &num_cols, const string &path) {
  ifstream file(path);
  asrt(file.is_open(), "File could not be opened.");

  string row_string;
  vector<vector<double>> rows;
  while (getline(file, row_string)) {
    istringstream row_stream(row_string);
    vector<double> row;
    double number;
    while (row_stream >> number) {
      row.push_back(number);
    }
    rows.push_back(row);
  }

  int num_rows = (int)rows.size();
  ASRTF(num_rows > 0, "Expected to read at least one row.");
  num_cols = (int)rows[0].size();

  for (auto &row : rows) {
    ASRTF((int)row.size() == num_cols,
          "The number of doubles in every row should match.");
  }

  vector<double *> mat = create_mat(num_rows, num_cols);
  for (int i = 0; i < num_rows; i++) {
    for (int j = 0; j < num_cols; j++) {
      mat[i][j] = rows[i][j];
    }
  }

  return mat;
}

dd_MatrixPtr double2cdd(const int n, const vector<double *> &A) {
  dd_MatrixPtr cdd_A = dd_CreateMatrix(A.size(), n);
  ASRTF(cdd_A != nullptr, "Failed to create cdd_A.");
  for (int i = 0; i < (int)A.size(); i++) {
    for (int j = 0; j < n; j++) {
      dd_set_d(cdd_A->matrix[i][j], A[i][j]);
    }
  }
  return cdd_A;
}

vector<double *> cdd2double(dd_MatrixPtr cdd_A) {
  const int rows = cdd_A->rowsize;
  const int cols = cdd_A->colsize;

  vector<double *> A = create_mat(rows, cols);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      A[i][j] = mpq_get_d(cdd_A->matrix[i][j]);
    }
  }
  return A;
}

void print_vertices(const int dim, const vector<mpq_t *> &vertices) {
  cout << "printing " << vertices.size() << " vertices" << endl;
  for (size_t vi = 0; vi < vertices.size(); vi++) {
    const auto &v = vertices[vi];
    cout << "i = " << vi << ": ";
    for (int i = 0; i < dim; i++) {
      cout << mpq_get_d(v[i]) << " ";
    }
    cout << endl;
  }
}

#include <iostream>
#include "utils.h"
#include "setoper.h"
#include "cdd.h"
#include <Eigen/Dense>
#include <numeric>

using namespace std;

Timer::Timer() {
    start = chrono::high_resolution_clock::now();
}

int Timer::micros() {
    auto now = chrono::high_resolution_clock::now();
    return chrono::duration_cast<chrono::microseconds>(now - start).count();
}

// It also does not remove members that have full incidence.
vector<int> compute_maximal_indexes(const vector<bset>& incidence) {
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

    sort(indexes_by_cardinality.begin(), indexes_by_cardinality.end(),
         [&incidence](int i, int j){return incidence[i].count() > incidence[j].count();});

    vector<int> maximal;
    int count_full_incidence = 0;
    for (int i : indexes_by_cardinality) {
        const bset& inc_i = incidence[i];
        if (inc_i.count() == num_containers) {
            count_full_incidence++;
            maximal.push_back(i);
            continue;
        }
        bool is_subset = false;
        for (int j = count_full_incidence; j < (int) maximal.size(); j++) {
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

MatrixXd read_matrix(const string& path) {
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

    size_t num_rows = rows.size();
    ASRTF(num_rows > 0, "Expected to read at least one row.");
    size_t num_cols = rows[0].size();

    for (auto& row : rows) {
        ASRTF(row.size() == num_cols, "The number of doubles in every row should match.");
    }

    MatrixXd mat(num_rows, num_cols);
    for (size_t i = 0; i < num_rows; i++) {
        for (size_t j = 0; j < num_cols; j++) {
            mat(i, j) = rows[i][j];
        }
    }

    return mat;
}

dd_MatrixPtr eigen2cdd(const MatrixXd& A) {
    dd_MatrixPtr cdd_A = dd_CreateMatrix(A.rows(), A.cols());
    ASRTF(cdd_A != nullptr, "Failed to create cdd_A.");
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            dd_set_d(cdd_A->matrix[i][j], A(i, j));
        }
    }
    return cdd_A;
}

MatrixXd cdd2eigen(dd_MatrixPtr cdd_A) {
    const int rows = cdd_A->rowsize;
    const int cols = cdd_A->colsize;

    MatrixXd A(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            A(i, j) = mpq_get_d(cdd_A->matrix[i][j]);
        }
    }
    return A;
}

void print_vertices(const int dim, const vector<mpq_t*>& vertices) {
    cout << "printing " << vertices.size() << " vertices" << endl;
    for (size_t vi = 0; vi < vertices.size(); vi++) {
        const auto& v = vertices[vi];
        cout << "i = " << vi << ": ";
        for (int i = 0; i < dim; i++) {
            cout << mpq_get_d(v[i]) << " ";
        }
        cout << endl;
    }
}

#pragma once

#include <unordered_set>
#include "setoper.h"
#include "cdd.h"
#include <Eigen/Dense>
#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include "asrt.h"

// I consider to be a zero everything that is abs(x) <= EPS.
constexpr double EPS = 1.0E-7;

using namespace std;
using namespace Eigen;

using bset = boost::dynamic_bitset<>;

class Timer {
    chrono::time_point<chrono::high_resolution_clock> start;

public:
    Timer();

    int micros();
};

vector<int> compute_maximal_indexes(const vector<bset>& incidence);

MatrixXd read_matrix(const string& path);

dd_MatrixPtr eigen2cdd(const MatrixXd& A);

void print_vertices(const int dim, const vector<mpq_t*>& vertices);

using Adj = pair<int, int>;

using RowXd = Matrix<double, 1, Eigen::Dynamic>;

using ColXd = Matrix<double, Eigen::Dynamic, 1>;

enum Polarity {
    MINUS,
    PLUS
};

using Quadrant = vector<Polarity>;

struct QuadrantInfo {
    int dim;
    vector<mpq_t*> V;
    vector<bset> V_to_H_incidence;
};

constexpr int POW2[11] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

constexpr int POW3[11] = {1, 3, 9, 27, 81, 243, 729, 2187, 6561, 19683, 59049};

constexpr int MAX_V2 = 8;
constexpr int MAX_V3 = 48;
constexpr int MAX_V4 = 384;
constexpr int MAX_V5 = 3840;

constexpr int NUM_H2 = POW3[2] - 1;
constexpr int NUM_H3 = POW3[3] - 1;
constexpr int NUM_H4 = POW3[4] - 1;
constexpr int NUM_H5 = POW3[5] - 1;

constexpr int K2NUM_H[6] = {-1, -1, NUM_H2, NUM_H3, NUM_H4, NUM_H5};
constexpr int K2MAX_V[6] = {-1, -1, MAX_V2, MAX_V3, MAX_V4, MAX_V5};

#pragma once

#include "asrt.h"
#include "cdd.h"
#include "setoper.h"
#include <array>
#include <boost/dynamic_bitset.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>

// I consider to be a zero everything that is abs(x) <= EPS.
constexpr double EPS = 1.0E-7;

using namespace std;

using bset = boost::dynamic_bitset<>;

vector<double *> create_mat(int rows, int cols);

void free_mat(const vector<double *> &mat);

class Timer {
    chrono::time_point<chrono::high_resolution_clock> start;

public:
    Timer();

    int micros();
};

vector<int> compute_maximal_indexes(const vector<bset> &incidence);

vector<double *> read_matrix(int &cols, const string &path);

dd_MatrixPtr double2cdd(int n, const vector<double *> &A);

vector<double *> cdd2double(dd_MatrixPtr cdd_A);

void print_vertices(const int dim, const vector<mpq_t *> &vertices);

using Adj = pair<int, int>;

enum Polarity { MINUS, PLUS };

using Quadrant = vector<Polarity>;

struct QuadrantInfo {
  int dim;
  vector<mpq_t *> V;
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

const vector<vector<int>> OCTAHEDRON_COEFS_1 = {{1}, {-1}};

const vector<vector<int>> OCTAHEDRON_COEFS_2 = {
        {1, 1},
        {1, 0},
        {1, -1},
        {0, 1},
        {0, -1},
        {-1, 1},
        {-1, 0},
        {-1, -1},
};

const vector<vector<int>> OCTAHEDRON_COEFS_3 = {
        {1, 1, 1},
        {1, 1, 0},
        {1, 1, -1},
        {1, 0, 1},
        {1, 0, 0},
        {1, 0, -1},
        {1, -1, 1},
        {1, -1, 0},
        {1, -1, -1},
        {0, 1, 1},
        {0, 1, 0},
        {0, 1, -1},
        {0, 0, 1},
        {0, 0, -1},
        {0, -1, 1},
        {0, -1, 0},
        {0, -1, -1},
        {-1, 1, 1},
        {-1, 1, 0},
        {-1, 1, -1},
        {-1, 0, 1},
        {-1, 0, 0},
        {-1, 0, -1},
        {-1, -1, 1},
        {-1, -1, 0},
        {-1, -1, -1},
};

const vector<vector<int>> OCTAHEDRON_COEFS_4 = {
        {1, 1, 1, 1},
        {1, 1, 1, 0},
        {1, 1, 1, -1},
        {1, 1, 0, 1},
        {1, 1, 0, 0},
        {1, 1, 0, -1},
        {1, 1, -1, 1},
        {1, 1, -1, 0},
        {1, 1, -1, -1},
        {1, 0, 1, 1},
        {1, 0, 1, 0},
        {1, 0, 1, -1},
        {1, 0, 0, 1},
        {1, 0, 0, 0},
        {1, 0, 0, -1},
        {1, 0, -1, 1},
        {1, 0, -1, 0},
        {1, 0, -1, -1},
        {1, -1, 1, 1},
        {1, -1, 1, 0},
        {1, -1, 1, -1},
        {1, -1, 0, 1},
        {1, -1, 0, 0},
        {1, -1, 0, -1},
        {1, -1, -1, 1},
        {1, -1, -1, 0},
        {1, -1, -1, -1},
        {0, 1, 1, 1},
        {0, 1, 1, 0},
        {0, 1, 1, -1},
        {0, 1, 0, 1},
        {0, 1, 0, 0},
        {0, 1, 0, -1},
        {0, 1, -1, 1},
        {0, 1, -1, 0},
        {0, 1, -1, -1},
        {0, 0, 1, 1},
        {0, 0, 1, 0},
        {0, 0, 1, -1},
        {0, 0, 0, 1},
        {0, 0, 0, -1},
        {0, 0, -1, 1},
        {0, 0, -1, 0},
        {0, 0, -1, -1},
        {0, -1, 1, 1},
        {0, -1, 1, 0},
        {0, -1, 1, -1},
        {0, -1, 0, 1},
        {0, -1, 0, 0},
        {0, -1, 0, -1},
        {0, -1, -1, 1},
        {0, -1, -1, 0},
        {0, -1, -1, -1},
        {-1, 1, 1, 1},
        {-1, 1, 1, 0},
        {-1, 1, 1, -1},
        {-1, 1, 0, 1},
        {-1, 1, 0, 0},
        {-1, 1, 0, -1},
        {-1, 1, -1, 1},
        {-1, 1, -1, 0},
        {-1, 1, -1, -1},
        {-1, 0, 1, 1},
        {-1, 0, 1, 0},
        {-1, 0, 1, -1},
        {-1, 0, 0, 1},
        {-1, 0, 0, 0},
        {-1, 0, 0, -1},
        {-1, 0, -1, 1},
        {-1, 0, -1, 0},
        {-1, 0, -1, -1},
        {-1, -1, 1, 1},
        {-1, -1, 1, 0},
        {-1, -1, 1, -1},
        {-1, -1, 0, 1},
        {-1, -1, 0, 0},
        {-1, -1, 0, -1},
        {-1, -1, -1, 1},
        {-1, -1, -1, 0},
        {-1, -1, -1, -1},
};


const vector<vector<vector<int>>> K2OCTAHEDRON_COEFS = {
        {},
        OCTAHEDRON_COEFS_1,
        OCTAHEDRON_COEFS_2,
        OCTAHEDRON_COEFS_3,
        OCTAHEDRON_COEFS_4};

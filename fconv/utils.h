#pragma once

#include <unordered_set>
#include "setoper.h"
#include "cdd.h"
#include <array>
#include <chrono>
#include <iostream>
#include "asrt.h"
#include "dynamic_bitset.h"

// I consider to be a zero everything that is abs(x) <= EPS.
constexpr double EPS = 1.0E-7;

using namespace std;

class Timer {
    chrono::time_point<chrono::high_resolution_clock> start;

public:
    Timer();

    int micros();
};

using Adj = pair<int, int>;

enum Polarity { MINUS, PLUS };

using Quadrant = vector<Polarity>;

struct QuadrantInfo {
  int dim;
  vector<mpq_t *> V;
  vector<set_t> V_to_H_incidence;
};

vector<int> compute_maximal_indexes(const vector<set_t>& incidence);

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

// The coefficient should be taken with the minus sign to get an upper bound.
constexpr double LOWER_BOUND_INDEX[5][4] = {
        {-1, -1, -1, -1},
        {0, -1, -1, -1},
        {1, 3, -1, -1},
        {4, 10, 12, -1},
        {13, 31, 37, 39}
};

constexpr double UPPER_BOUND_INDEX[5][4] = {
        {-1, -1, -1, -1},
        {1, -1, -1, -1},
        {6, 4, -1, -1},
        {21, 15, 13, -1},
        {66, 48, 42, 40}
};

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

const vector<Quadrant> QUADRANTS_1 = {
        {MINUS},
        {PLUS}
};

const vector<Quadrant> QUADRANTS_2 = {
        {MINUS, MINUS},
        {MINUS, PLUS},
        {PLUS, MINUS},
        {PLUS, PLUS}
};

const vector<Quadrant> QUADRANTS_3 = {
        {MINUS, MINUS, MINUS},
        {MINUS, MINUS, PLUS},
        {MINUS, PLUS, MINUS},
        {MINUS, PLUS, PLUS},
        {PLUS, MINUS, MINUS},
        {PLUS, MINUS, PLUS},
        {PLUS, PLUS, MINUS},
        {PLUS, PLUS, PLUS}
};

const vector<Quadrant> QUADRANTS_4 = {
        {MINUS, MINUS, MINUS, MINUS},
        {MINUS, MINUS, MINUS, PLUS},
        {MINUS, MINUS, PLUS, MINUS},
        {MINUS, MINUS, PLUS, PLUS},
        {MINUS, PLUS, MINUS, MINUS},
        {MINUS, PLUS, MINUS, PLUS},
        {MINUS, PLUS, PLUS, MINUS},
        {MINUS, PLUS, PLUS, PLUS},
        {PLUS, MINUS, MINUS, MINUS},
        {PLUS, MINUS, MINUS, PLUS},
        {PLUS, MINUS, PLUS, MINUS},
        {PLUS, MINUS, PLUS, PLUS},
        {PLUS, PLUS, MINUS, MINUS},
        {PLUS, PLUS, MINUS, PLUS},
        {PLUS, PLUS, PLUS, MINUS},
        {PLUS, PLUS, PLUS, PLUS}
};

const vector<vector<Quadrant>> K2QUADRANTS = {
        {},
        QUADRANTS_1,
        QUADRANTS_2,
        QUADRANTS_3,
        QUADRANTS_4};

#pragma once
#include <execinfo.h>
#include <signal.h>
#include <map>
#include "gmp.h"
#include "utils.h"
#include "dynamic_bitset.h"

using namespace std;

struct OctahedronV {
    int dim;

    vector<mpq_t*> V;

    // This is the limited set of adjacencies - I only store the adjacencies that are on different
    // sides of at least one orthant. These are the adjacencies that will produce a new vertex
    // when I will be splitting octahedron by the orthant.
    vector<Adj> orthant_adjacencies;

    // Incidence of vertices to input octahedron constraints.
    vector<set_t> incidence;
};

OctahedronV compute_octahedron_V(int K, const vector<double *> &A);

OctahedronV compute_V_with_cdd(int K, const vector<double *> &A);

map<Quadrant, QuadrantInfo>
compute_quadrants_with_cdd(int K, const vector<double *> &A);

vector<QuadrantInfo>
compute_max_pool_quadrants_with_cdd(const int K, const vector<double *> &A);

// Tanh/Sigmoid
map<Quadrant, vector<mpq_t *>>
compute_tasi_quadrants_with_cdd(const int K, const vector<double *> &A,
                                Activation activation);

// Tanh/Sigmoid
map<Quadrant, vector<mpq_t *>>
compute_tasi_quadrants_with_cdd_dim(const int K, const vector<double *> &A,
                                    Activation activation);

#pragma once

#include <map>
#include "utils.h"

map<Quadrant, VInc_mpq> get_quadrants_cdd(int K, const vector<double*>& A);

vector<DD_mpq> get_pool_quadrants_cdd(int K, const vector<double*>& A);

vector<DD_mpq> get_pool_quadrants(int K, const vector<double*>& A);

// Tanh/Sigmoid
map<Quadrant, vector<mpq_t*>> get_tasi_quadrants_cdd(int K,
                                                     const vector<double*>& A,
                                                     Activation activation);

// Tanh/Sigmoid
map<Quadrant, vector<mpq_t*>> get_tasi_quadrants_cdd_dim(int K,
                                                         const vector<double*>& A,
                                                         Activation activation);

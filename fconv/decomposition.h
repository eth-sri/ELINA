#pragma once

#include "pdd.h"

/* The function take ownership of all memory in the input. */
vector<double*> decomposition(const int K, const map<Quadrant, PDD>& quadrant2pdd, Activation activation,
                              const vector<double>& x_lb, const vector<double>& x_ub,
                              const vector<double>& orthants);

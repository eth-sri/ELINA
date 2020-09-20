#pragma once
#include <map>
#include "utils.h"

map<Quadrant, QuadrantInfo>
split_in_quadrants(vector<mpq_t *> &V, vector<bset> &incidence,
                   const vector<Adj> &orthant_adjacencies, const int K);

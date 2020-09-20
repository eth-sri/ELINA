#pragma once

#include "pdd.h"

MatrixXd decomposition(const map<Quadrant, PDD> &quadrant2pdd, const int K);

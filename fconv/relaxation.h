#pragma once

#include <vector>
#include "utils.h"

using namespace std;

vector<double*> fast_relaxation_through_decomposition(int K,
                                                      const vector<double*>& A,
                                                      Activation activation);

vector<double*> krelu_with_cdd(int K, const vector<double*>& A);

vector<double*> fkpool(int K, const vector<double*>& A);

vector<double*> kpool_with_cdd(int K, const vector<double*>& A);

vector<double*> ktanh_with_cdd(int K, const vector<double*>& A);

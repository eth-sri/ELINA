#pragma once

#include <vector>

using namespace std;

vector<double *> fkpool(int K, const vector<double *> &A);

vector<double *> kpool_with_cdd(int K, const vector<double *> &A);

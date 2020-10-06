#pragma once

#include <vector>

using namespace std;

vector<double*> fkrelu(int K, const vector<double*>& A);

vector<double*> krelu_with_cdd(int K, const vector<double*>& A);

vector<double*> fkpool(int K, const vector<double*>& A);

vector<double*> kpool_with_cdd(int K, const vector<double*>& A);

vector<double*> ktanh_with_cdd(int K, const vector<double*>& A);

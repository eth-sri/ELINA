#pragma once

std::vector<double*> fkrelu(int K, const std::vector<double*>& A);

std::vector<double*> krelu_with_cdd(int K, const std::vector<double*>& A);

#pragma once

#include <gmp.h>
#include <vector>

using namespace std;

mpq_t* mpq_create_array(int n);

mpq_t* mpq_create_and_copy_array(int n, const mpq_t* src);

void mpq_free_array(int n, mpq_t* arr);

bool mpq_arrays_are_equal(int n, mpq_t* first, mpq_t* second);

vector<double*> mpq_to_double(int n, const vector<mpq_t*>& mpq_A);

vector<mpq_t*> mpq_from_double(int n, const vector<double*>& A);

vector<mpq_t*> mpq_copy_array_vector(int n, const vector<mpq_t*>& src_vector);

void mpq_free_array_vector(int n, vector<mpq_t*>& vec);

void mpq_print_array(int n, const mpq_t* v);

vector<mpq_t*> mpq_matrix_mul(int n, const vector<mpq_t*>& A, const vector<mpq_t*>& B);

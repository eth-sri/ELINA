#pragma once
#include "S_curve.h"

typedef struct {
    int rows;
    int cols;
    double *data;
} MatDouble;

typedef struct {
    int rows;
    int cols;
    int *data;
} MatInt;

#ifdef __cplusplus
extern "C" {
#endif

MatDouble new_MatDouble(
        int rows,
        int cols,
        const double *data);

void free_MatDouble(MatDouble mat);

void free_MatInt(MatInt mat);

MatDouble fkrelu(MatDouble input_hrep);

MatDouble krelu_with_cdd(MatDouble input_hrep);

MatDouble fkpool(MatDouble input_hrep);

MatDouble kpool_with_cdd(MatDouble input_hrep);

MatDouble fktanh(MatDouble input_hrep);

MatDouble ktanh_with_cdd(MatDouble input_hrep);

MatDouble fksigm(MatDouble input_hrep);

MatDouble ksigm_with_cdd(MatDouble input_hrep);

MatDouble ftanh_orthant(MatDouble input_hrep);

MatDouble fsigm_orthant(MatDouble input_hrep);

MatInt generate_sparse_cover(int N, int K);

void S_curve_chord_bound(double* k, double* b, double x_lb, double x_ub, bool is_sigm);

void S_curve_tang_bound(double* k, double* b, double x, bool is_sigm);

#ifdef __cplusplus
}
#endif

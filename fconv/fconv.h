#pragma once

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

MatInt generate_sparse_cover(int N, int K);

#ifdef __cplusplus
}
#endif

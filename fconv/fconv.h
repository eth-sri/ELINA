#pragma once

typedef struct {
  const int rows;
  const int cols;
  const double *data;
} MatDouble;

typedef struct {
  const int rows;
  const int cols;
  const int *data;
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

MatInt generate_sparse_cover(int N, int K);

#ifdef __cplusplus
}
#endif

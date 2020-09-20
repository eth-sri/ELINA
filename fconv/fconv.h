#pragma once

typedef struct {
    const int rows;
    const int cols;
    const double *data;
} MatDouble;

#ifdef __cplusplus
extern "C" {
#endif

MatDouble new_MatDouble(
        int rows,
        int cols,
        const double *data);

void free_MatDouble(MatDouble mat);

MatDouble fkrelu(MatDouble input_hrep);

#ifdef __cplusplus
}
#endif

#include "fconv.h"
#include "fkrelu.h"
#include "utils.h"
#include <Eigen/Dense>

using namespace Eigen;

MatrixXd cmat2eigen(const MatDouble &cmat) {
  MatrixXd A(cmat.rows, cmat.cols);
  for (int i = 0; i < cmat.rows; i++) {
    for (int j = 0; j < cmat.cols; j++) {
      A(i, j) = cmat.data[i * cmat.cols + j];
    }
  }
  return A;
}

MatDouble eigen2cmat(const MatrixXd &A) {
  int rows = A.rows();
  int cols = A.cols();
  auto mat = (double *)calloc(rows * cols, sizeof(double));

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      mat[i * cols + j] = A(i, j);
    }
  }
  return {rows, cols, mat};
}

MatDouble new_MatDouble(int rows, int cols, const double *data) {
    ASRTF(rows > 0, "Only rows > 0 is allowed");
    ASRTF(cols > 0, "Only cols > 0 is allowed");
    auto data_copy = (double *) calloc(rows * cols, sizeof(double));
    for (int i = 0; i < rows * cols; i++) {
        data_copy[i] = data[i];
    }
    return {rows, cols, data_copy};
}

void free_MatDouble(MatDouble cmat) {
    free((double *) cmat.data);
}

MatDouble fkrelu(MatDouble input_hrep) {
  dd_set_global_constants();
  MatrixXd A = cmat2eigen(input_hrep);
  MatrixXd H = fkrelu(A);
  dd_free_global_constants();
  return eigen2cmat(H);
}

MatDouble krelu_with_cdd(MatDouble input_hrep) {
  dd_set_global_constants();
  MatrixXd A = cmat2eigen(input_hrep);
  MatrixXd H = krelu_with_cdd(A);
  dd_free_global_constants();
  return eigen2cmat(H);
}

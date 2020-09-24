#pragma once

#include <Eigen/Dense>

Eigen::MatrixXd fkrelu(const Eigen::MatrixXd& A);

Eigen::MatrixXd krelu_with_cdd(const Eigen::MatrixXd& A);

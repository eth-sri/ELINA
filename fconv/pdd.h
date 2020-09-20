#pragma once

#include "utils.h"
#include <Eigen/Dense>
#include <boost/dynamic_bitset.hpp>

using namespace std;

struct PDD {
    int dim;
    MatrixXd H;
    MatrixXd V;
    vector<bset> incidence; // V to H incidence.
};

PDD PDD_from_H_and_V(const MatrixXd &H, const MatrixXd &V);

vector<bset> PDD_compute_incidence(const MatrixXd &H, const MatrixXd &V);

PDD PDD_batch_intersect(const PDD &pdd, const MatrixXd &H_new);

PDD PDD_intersect_two_PDDs(const PDD &pdd1, const PDD &pdd2);

void PDD_adjust_H_for_soundness_finite_polytope(MatrixXd &H, const MatrixXd &V);

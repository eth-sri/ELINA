#include "pdd.h"

struct PDD_VInc {
    MatrixXd V;
    vector<bset> incidence;
};

MatrixXd vert_concat2(const MatrixXd& top, const MatrixXd& bottom) {
    ASRTF(top.cols() == bottom.cols(), "To vertically concatenate matrices number of columns should match.");
    MatrixXd res(top.rows() + bottom.rows(), top.cols());
    res.topRows(top.rows()) = top;
    res.bottomRows(bottom.rows()) = bottom;
    return res;
}

// Add normalization of H.
// Add checking that V is finite.
// Add checking that H doesn't contain null vector.
// TODO[gleb] I have not properly worked with this function yet, so it is not debugged well yet.
PDD PDD_from_H_and_V(const MatrixXd& H, const MatrixXd& V) {
    ASRTF(H.cols() == V.cols(), "Dimensionality of H and V should match.");
    ASRTF(H.cols() > 0, "Dimensionality should be at least 1.");
    ASRTF(H.rows() > 0, "There should be at least one constraint.");
    ASRTF(V.rows() > 0, "There should be at least one vertex.");

    ASRTF(H.cwiseAbs().rowwise().maxCoeff().minCoeff() >= EPS, "There should be no null constraints.");
    vector<bset> incidence = PDD_compute_incidence(H, V);

    ASRTF(incidence.size() == (size_t) V.rows(), "PDD_from_H_and_V consistency 1");
    ASRTF(incidence[0].size() == (size_t) H.rows(), "PDD_from_H_and_V consistency 2");

    return {(int) H.cols(), H, V, incidence};
}

vector<bset> PDD_compute_incidence(const MatrixXd& H, const MatrixXd& V) {
    ASRTF(H.cols() == V.cols(), "Dimensionality of H and V should match.");
    MatrixXd H_x_V = H * V.transpose();
    ASRTF(H_x_V.minCoeff() >= -EPS, "H and V do not form a PDD.");

    const int num_H = H.rows();
    const int num_V = V.rows();
    vector<bset> incidence(num_V, bset(num_H));

    for (int vi = 0; vi < num_V; vi++) {
        bset& inc = incidence[vi];
        const ColXd& H_x_v = H_x_V.col(vi);
        for (int hi = 0; hi < num_H; hi++) {
            if (H_x_v(hi, 0) <= EPS) {
                inc.set(hi);
            }
        }
    }

    return incidence;
}

PDD_VInc PDD_make_irredundant(const MatrixXd& V, const vector<bset>& incidence) {
    vector<int> maximal = compute_maximal_indexes(incidence);
    MatrixXd V_res(maximal.size(), V.cols());
    vector<bset> incidence_res(maximal.size());
    for (size_t i = 0; i < maximal.size(); i++) {
        V_res.row(i) = V.row(maximal[i]);
        incidence_res[i] = incidence[maximal[i]];
    }
    return {V_res, incidence_res};
}

PDD_VInc PDD_batch_intersect_helper(const PDD& pdd, const MatrixXd& H_new) {
    const int dim = pdd.dim;

    ASRTF(H_new.cols() == dim, "Dimensionality of new hyperplanes should equal PDD dimensionality.");
    const int num_H_base = (int) pdd.H.rows();
    const int num_H_new = (int) H_new.rows();
    ASRTF(num_H_new > 0, "There should be non-empty number of new hyperplanes.");

    assert(
            pdd.V.rows() == (int) pdd.incidence.size() &&
            "Consistency checking number of vertices and size of incidence.");

    const MatrixXd& V_base = pdd.V;
    const Eigen::MatrixXd H_x_V = H_new * V_base.transpose();

    // Points that violate at least one of the constraints.
    vector<int> outs;
    vector<vector<int>> outs_violation_vec;
    vector<bset> outs_violation;
    vector<bset> outs_incidence;
    vector<bset> outs_incidence_new;
    // Points that don't violate any constraints.
    vector<int> ins;
    vector<bset> ins_incidence;
    vector<bset> ins_incidence_new;

    for (int vi = 0; vi < (int) V_base.rows(); vi++) {
        vector<int> vio_vec;
        bset vio(num_H_new);
        bset inc_new(num_H_new);
        bset inc = pdd.incidence[vi]; // Intentionally doing a copy.
        assert((int) inc.size() == num_H_base && "Sanity checking the size of incidence.");
        inc.resize(num_H_base + num_H_new);
        const Eigen::Matrix<double, Eigen::Dynamic, 1>& H_x_v = H_x_V.col(vi);

        for (int hi = 0; hi < num_H_new; hi++) {
            double val = H_x_v(hi, 0);
            if (val < -EPS) {
                vio_vec.push_back(hi);
                vio.set(hi);
            } else if (val <= EPS) {
                inc_new.set(hi);
                inc.set(num_H_base + hi);
            }
        }
        assert(vio_vec.size() == vio.count() && "Sanity check violation sizes should match.");
        if (!vio_vec.empty()) {
            outs.push_back(vi);
            outs_violation_vec.push_back(vio_vec);
            outs_violation.push_back(vio);
            outs_incidence.push_back(inc);
            outs_incidence_new.push_back(inc_new);
        } else {
            ins.push_back(vi);
            ins_incidence.push_back(inc);
            ins_incidence_new.push_back(inc_new);
        }
    }

    size_t max_count = outs.size() * ins.size() + ins.size();
    MatrixXd V_res(max_count, dim);
    vector<bset> V_res_incidence(max_count);

    size_t count = 0;
    for (int i_out = 0; i_out < (int) outs.size(); i_out++) {
        const int out = outs[i_out];
        const ColXd& H_x_out = H_x_V.col(out);
        const RowXd& V_out = V_base.row(out);

        const vector<int>& vio_vec = outs_violation_vec[i_out];
        const bset& vio = outs_violation[i_out];
        const bset& out_inc = outs_incidence[i_out];

        for (int i_in = 0; i_in < (int) ins.size(); i_in++) {
            const bset& in_inc_new = ins_incidence_new[i_in];

            if ((vio & in_inc_new).any()) {
                continue;
            }

            const int in = ins[i_in];
            const bset& in_inc = ins_incidence[i_in];

            // This code can be accelerated using the analog of the Chernikova criteria.
            // That is based on the cardinality of intersection of in_inc and out_inc.

            const Matrix<double, Eigen::Dynamic, 1>& H_x_in = H_x_V.col(in);

            int argmin_h = -1;
            double min_ratio_h = -1;
            for (int hi : vio_vec) {
                double cur_ratio_h = -H_x_in(hi, 0) / H_x_out(hi, 0);
                assert(cur_ratio_h > 0 && "Sanity check that cur_ratio is positive.");
                assert(H_x_in(hi, 0) > 0 && "Sanity check H_x_in(hi, 0) > 0");
                assert(H_x_out(hi, 0) < 0 && "Sanity check H_x_out(hi, 0) < 0");
                if (cur_ratio_h < min_ratio_h || argmin_h == -1) {
                    argmin_h = hi;
                    min_ratio_h = cur_ratio_h;
                }
            }
            const RowXd& V_in = V_base.row(in);

            RowXd v = H_x_in(argmin_h, 0) * V_out - H_x_out(argmin_h, 0) * V_in;
            double abs_coef = v.array().abs().maxCoeff();

            if (abs_coef <= EPS) {
                // I.e. got the null vector, that is possible. Just ignoring it.
                continue;
            }

            // Scaling for numerical stability.
            V_res.row(count) = v / abs_coef;
            V_res_incidence[count] = in_inc & out_inc;
            // Note that theoretically ray-shooting might hit multiple hyperplanes at a time.
            // So the incidence potentially can be slightly bigger for this vertex, and it is very
            // easy to implement it in exact precision. However, with fp its more tricky and
            // since we are anyway approximating, I'm keeping it simple.
            V_res_incidence[count].set(num_H_base + argmin_h);
            count++;
        }
    }

    // Also can be useful to do the negative ray-shooting.
    // But for now I have disabled it since it is not very useful for our problem.
//    for (int i_o1 = 0; i_o1 < (int) outs.size() - 1; i_o1++) {
//        const int o1 = outs[i_o1];
//        const ColXd& H_x_o1 = H_x_V.col(o1);
//        const RowXd& V_o1 = V_base.row(o1);
//
//        const vector<int>& vio_vec1 = outs_violation_vec[i_o1];
//        const bset& vio1 = outs_violation[i_o1];
//        const bset& inc1 = outs_incidence[i_o1];
//        const bset& inc_new1 = outs_incidence_new[i_o1];
//
//        for (int i_o2 = i_o1 + 1; i_o2 < (int) outs.size(); i_o2++) {
//            const bset& vio2 = outs_violation[i_o2];
//            const bset& inc_new2 = outs_incidence_new[i_o2];
//
//            if ((vio1 & vio2).any() || (vio1 & inc_new2).any() || (vio2 & inc_new1).any()) {
//                continue;
//            }
//
//            const vector<int>& vio_vec2 = outs_violation_vec[i_o2];
//            const int o2 = outs[i_o2];
//            const bset& inc2 = outs_incidence[i_o2];
//
//            // This code can be accelerated using the analog of the Chernikova criteria.
//            // That is based on the cardinality of intersection of in_inc and out_inc.
//
//            const ColXd& H_x_o2 = H_x_V.col(o2);
//            const RowXd& V_o2 = V_base.row(o2);
//
//            int argmin_h1, argmin_h2;
//            double min_ratio_h1, min_ratio_h2;
//            shoot_ray_to_find_closest_hyperplane(argmin_h1, min_ratio_h1, vio_vec1, H_x_o2, H_x_o1);
//            shoot_ray_to_find_closest_hyperplane(argmin_h2, min_ratio_h2, vio_vec2, H_x_o1, H_x_o2);
//
//            if (min_ratio_h1 < (1 / min_ratio_h2)) {
//                // If that's the case - there is no sound region between two out vertices.
//                continue;
//            }
//
//            RowXd v1 = H_x_o2(argmin_h1, 0) * V_o1 - H_x_o1(argmin_h1, 0) * V_o2;
//            bset inc = inc1 & inc2;
//            inc.set(num_H_base + argmin_h1);
//            maybe_add_new_vertex(v1, inc);
//
//            RowXd v2 = H_x_o1(argmin_h2, 0) * V_o2 - H_x_o2(argmin_h2, 0) * V_o1;
//            inc.set(num_H_base + argmin_h1, false);
//            inc.set(num_H_base + argmin_h2);
//            maybe_add_new_vertex(v2, inc);
//        }
//    }

    for (size_t i = 0; i < ins.size(); i++) {
        V_res.row(count) = V_base.row(ins[i]);
        V_res_incidence[count] = ins_incidence[i];
        count++;
    }

    assert(count <= max_count && "Sanity checking the count size.");

    V_res.conservativeResize(count, NoChange);
    V_res_incidence.resize(count);

    return {V_res, V_res_incidence};
}


PDD PDD_batch_intersect(const PDD& pdd, const MatrixXd& H_new) {
    PDD_VInc VInc = PDD_batch_intersect_helper(pdd, H_new);
    VInc = PDD_make_irredundant(VInc.V, VInc.incidence);
    MatrixXd H = vert_concat2(pdd.H, H_new);
    return {pdd.dim, H, VInc.V, VInc.incidence};
}

PDD PDD_intersect_two_PDDs(const PDD& pdd1, const PDD& pdd2) {
    ASRTF(pdd1.dim == pdd2.dim, "Dimensions of two PDDs should match.");

    const MatrixXd& H1 = pdd1.H;
    const MatrixXd& H2 = pdd2.H;

    const int num_H1 = (int) H1.rows();
    const int num_H2 = (int) H2.rows();

    if (num_H2 == 0) {
        return pdd1;
    } else if (num_H1 == 0) {
        return pdd2;
    }

    PDD_VInc VInc1 = PDD_batch_intersect_helper(pdd1, H2);
    PDD_VInc VInc2 = PDD_batch_intersect_helper(pdd2, H1);

    MatrixXd V = vert_concat2(VInc1.V, VInc2.V);

    vector<bset> V_incidence = VInc1.incidence; // Intentionally doing a copy.
    V_incidence.reserve(V.rows());

    const vector<bset>& V2_incidence = VInc2.incidence;
    for (auto& inc_H2_H1 : V2_incidence) {
        bset inc(num_H1 + num_H2);
        // Probably there is a much more efficient way to do this operation,
        // but I doubt that it will be a bottleneck. Should double check.
        for (int i = 0; i < num_H1; i++) {
            inc.set(i, inc_H2_H1.test(num_H2 + i));
        }
        for (int i = 0; i < num_H2; i++) {
            inc.set(num_H1 + i, inc_H2_H1.test(i));
        }
        V_incidence.push_back(inc);
    }

    PDD_VInc VInc = PDD_make_irredundant(V, V_incidence);

    return {pdd1.dim, vert_concat2(H1, H2), VInc.V, VInc.incidence};
}

void PDD_adjust_H_for_soundness_finite_polytope(MatrixXd& H, const MatrixXd& V) {
    ASRTF(H.cols() == V.cols(), "Dimensions of H and V should match.");

    for (int vi = 0; vi < V.rows(); vi++) {
        ASRTF(V(vi, 0) == 1, "The adjustment strategy only works with vertices, i.e. no rays allowed.");
    }

    MatrixXd H_x_V = H * V.transpose();

    MatrixXd H_abs = H.cwiseAbs();
    MatrixXd V_abs = V.cwiseAbs();

    MatrixXd H_x_V_err = H_abs * V_abs.transpose();
    constexpr double REL_ERR = 8.8817842e-16; // 1 / 2^50.
    H_x_V_err *= REL_ERR;

    // Sound lower bound of the scalar product.
    MatrixXd H_x_V_lb = H_x_V - H_x_V_err;

    RowXd Adjustment = H_x_V_lb.rowwise().minCoeff();
    ASRTF(Adjustment.size() == H.rows(), "Sanity check of the size of the adjustment vector.");

    for (int hi = 0; hi < (int) H.rows(); hi++) {
        H(hi, 0) -= Adjustment(hi);
    }
}

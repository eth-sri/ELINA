#include "pdd.h"

// This can be optimized if there will be a need for it.
vector<double *> mat_mul_with_transpose(const int n,
                                        const vector<double *> &mat1,
                                        const vector<double *> &mat2) {
  vector<double *> res = create_mat(mat1.size(), mat2.size());
  for (size_t i1 = 0; i1 < mat1.size(); i1++) {
    const double *r1 = mat1[i1];
    double *res_row = res[i1];
    for (size_t i2 = 0; i2 < mat2.size(); i2++) {
      const double *r2 = mat2[i2];
      double accum = 0;
      for (int j = 0; j < n; j++) {
        accum += r1[j] * r2[j];
      }
      res_row[i2] = accum;
    }
  }
  return res;
}

struct PDD_VInc {
    vector<double*> V;
    vector<bset> incidence;
};

/*
 * The function takes ownership of the input's memory.
 */
PDD_VInc PDD_make_irredundant(const vector<double *> &V,
                              const vector<bset> &incidence) {
  assert(V.size() == incidence.size() &&
         "Consistency checking that size of V equals size of incidence.");
  vector<int> maximal = compute_maximal_indexes(incidence);
  vector<double *> V_res(maximal.size());
  vector<bset> incidence_res(maximal.size());
  bset is_maximal(V.size());
  for (int i : maximal) {
    is_maximal.set(i);
  }
  size_t count = 0;
  for (size_t i = 0; i < V.size(); i++) {
    if (is_maximal.test(i)) {
      V_res[count] = V[i];
      incidence_res[count] = incidence[i];
      count++;
    } else {
      free(V[i]);
    }
  }
  assert(count == maximal.size() &&
         "Consistency checking that count == maximal.size().");

  return {V_res, incidence_res};
}

/*
 * Note that function doesn't take ownership of H_new memory and transfers ownership of
 * pdd.V and pdd.incidence to PDD_VInc.
 * It leaves pdd in inconsistent state - and the function is only to be used internally.
 */
PDD_VInc PDD_batch_intersect_helper(const PDD &pdd,
                                    const vector<double *> &H_new) {
  const int dim = pdd.dim;

  ASRTF(!H_new.empty(), "There should be non-empty number of new hyperplanes.");

  assert(pdd.V.size() == pdd.incidence.size() &&
         "Consistency checking number of vertices and size of incidence.");

  vector<double *> V_x_H = mat_mul_with_transpose(dim, pdd.V, H_new);

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

  for (size_t vi = 0; vi < pdd.V.size(); vi++) {
    vector<int> vio_vec;
    bset vio(H_new.size());
    bset inc_new(H_new.size());
    bset inc = pdd.incidence[vi]; // Intentionally doing a copy.
    assert(inc.size() == pdd.H.size() &&
           "Sanity checking the size of incidence.");
    inc.resize(pdd.H.size() + H_new.size());
    const double *v_x_H = V_x_H[vi];

    for (size_t hi = 0; hi < H_new.size(); hi++) {
      double val = v_x_H[hi];
      if (val < -EPS) {
        vio_vec.push_back(hi);
        vio.set(hi);
      } else if (val <= EPS) {
        inc_new.set(hi);
        inc.set(pdd.H.size() + hi);
      }
    }
    assert(vio_vec.size() == vio.count() &&
           "Sanity check violation sizes should match.");
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
  vector<double *> V_res(max_count);
  vector<bset> V_res_incidence(max_count);

  size_t count = 0;
  for (size_t i_out = 0; i_out < outs.size(); i_out++) {
    const int out = outs[i_out];
    const double *out_x_H = V_x_H[out];
    const double *V_out = pdd.V[out];

    const vector<int> &vio_vec = outs_violation_vec[i_out];
    const bset &vio = outs_violation[i_out];
    const bset &out_inc = outs_incidence[i_out];

    for (size_t i_in = 0; i_in < ins.size(); i_in++) {
      const bset &in_inc_new = ins_incidence_new[i_in];

      if ((vio & in_inc_new).any()) {
        continue;
      }

      const int in = ins[i_in];
      const bset &in_inc = ins_incidence[i_in];

      // This code can be accelerated using the analog of the Chernikova
      // criteria. That is based on the cardinality of intersection of in_inc
      // and out_inc.

      const double *in_x_H = V_x_H[in];

      int argmin_h = -1;
      double min_ratio_h = -1;
      for (int hi : vio_vec) {
        double cur_ratio_h = -in_x_H[hi] / out_x_H[hi];
        assert(cur_ratio_h > 0 && "Sanity check that cur_ratio is positive.");
        assert(in_x_H[hi] > 0 && "Sanity check in_x_H[hi] > 0");
        assert(out_x_H[hi] < 0 && "Sanity check out_x_H[hi] > 0");
        if (cur_ratio_h < min_ratio_h || argmin_h == -1) {
          argmin_h = hi;
          min_ratio_h = cur_ratio_h;
        }
      }
      const double *V_in = pdd.V[in];

      double *v = (double *)calloc(dim, sizeof(double));
      double abs_coef = 0;
      double out_coef = in_x_H[argmin_h];
      double in_coef = -out_x_H[argmin_h];
      for (int i = 0; i < dim; i++) {
        double val = out_coef * V_out[i] + in_coef * V_in[i];
        abs_coef = max(abs(val), abs_coef);
        v[i] = val;
      }
      assert(abs_coef >= 0 && "abs_coef can't be negative.");
      if (abs_coef <= EPS) {
        // I.e. got the null vector, that is possible. Just ignoring it.
        continue;
      }

      // Scaling for numerical stability.
      for (int i = 0; i < dim; i++) {
        v[i] /= abs_coef;
      }

      V_res[count] = v;
      V_res_incidence[count] = in_inc & out_inc;

      // Note that theoretically ray-shooting might hit multiple hyperplanes at
      // a time. So the incidence potentially can be slightly bigger for this
      // vertex, and it is very easy to implement it in exact precision.
      // However, with fp its more tricky and since we are anyway approximating,
      // I'm keeping it simple.
      V_res_incidence[count].set(pdd.H.size() + argmin_h);
      count++;
    }
  }

  // It is also possible to do negative-negative ray-shooting.
  // But for now I have removed it since it didn't prove to be very useful for
  // our problem.

  for (size_t i = 0; i < ins.size(); i++) {
    V_res[count] = pdd.V[ins[i]];
    V_res_incidence[count] = ins_incidence[i];
    count++;
  }

  assert(count <= max_count && "Sanity checking the count size.");
  V_res.resize(count);
  V_res_incidence.resize(count);

  for (int out : outs) {
    free(pdd.V[out]);
  }

  free_mat(V_x_H);

  return {V_res, V_res_incidence};
}

void PDD_batch_intersect(PDD& pdd, vector<double*>& H_new) {
    PDD_VInc VInc = PDD_batch_intersect_helper(pdd, H_new);
    VInc = PDD_make_irredundant(VInc.V, VInc.incidence);
    pdd.H.reserve(pdd.H.size() + H_new.size());
    for (auto h : H_new) {
        pdd.H.push_back(h);
    }
    pdd.V = VInc.V;
    pdd.incidence = VInc.incidence;
}

PDD PDD_intersect_two_PDDs(PDD& pdd1, PDD& pdd2) {
    ASRTF(pdd1.dim == pdd2.dim, "Dimensions of two PDDs should match.");

    const vector<double*>& H1 = pdd1.H;
    const vector<double*>& H2 = pdd2.H;

    if (H1.empty()) {
        assert(pdd1.V.empty() && "Consistency checking that pdd1.V is empty.");
        return pdd2;
    } else if (H2.empty()) {
        assert(pdd2.V.empty() && "Consistency checking that pdd2.V is empty.");
        return pdd1;
    }

    PDD_VInc VInc1 = PDD_batch_intersect_helper(pdd1, H2);
    PDD_VInc VInc2 = PDD_batch_intersect_helper(pdd2, H1);

    vector<double*>& V = VInc1.V;
    V.reserve(V.size() + VInc2.V.size());
    for (auto v : VInc2.V) {
        V.push_back(v);
    }

    vector<bset> &V_incidence = VInc1.incidence;
    V_incidence.reserve(V.size());

    const size_t num_H1 = H1.size();
    const size_t num_H2 = H2.size();
    for (auto& inc_H2_H1 : VInc2.incidence) {
      bset inc(num_H1 + num_H2);
      // Probably there is a much more efficient way to do this operation,
      // but I doubt that it will be a bottleneck. Should double check.
      for (size_t i = 0; i < num_H1; i++) {
        inc.set(i, inc_H2_H1.test(num_H2 + i));
        }
        for (size_t i = 0; i < num_H2; i++) {
          inc.set(num_H1 + i, inc_H2_H1.test(i));
        }
        V_incidence.push_back(inc);
    }

    PDD_VInc VInc = PDD_make_irredundant(V, V_incidence);

    vector<double*>& H = pdd1.H;
    H.reserve(H1.size() + H2.size());
    for (auto h : H2) {
        H.push_back(h);
    }

    return {pdd1.dim, H, VInc.V, VInc.incidence};
}

void PDD_adjust_H_for_soundness_finite_polytope(const int dim,
                                                vector<double*>& H,
                                                const vector<double*>& V) {
    ASRTF(dim >= 2, "Expected that dimension >= 2.");
    for (auto v : V) {
        ASRTF(v[0] == 1, "The adjustment strategy only works with vertices, i.e. no rays allowed.");
    }
    for (size_t hi = 0; hi < H.size(); hi++) {
        double* h = H[hi];
        double min_adjustment = 0;
        for (size_t vi = 0; vi < V.size(); vi++) {
            double* v = V[vi];
            double accum = 0;
            double accum_abs = 0;
            for (int i = 0; i < dim; i++) {
                accum += h[i] * v[i];
                accum_abs += abs(h[i] * v[i]);
            }
            constexpr double REL_ERR = 8.8817842e-16; // 1 / 2^50.
            double adjustment = accum - REL_ERR * accum_abs;
            min_adjustment = min(adjustment, min_adjustment);
        }
        h[0] -= min_adjustment;
    }
}

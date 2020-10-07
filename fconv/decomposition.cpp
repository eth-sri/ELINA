#include <map>
#include <cassert>
#include <math.h>
#include "pdd.h"
#include "utils.h"
#include "fp_mat.h"

using namespace std;

void project_to_relu_y_branch(const int xi, PDD &pdd_dual,
                              const Polarity polarity) {
  const int dim = pdd_dual.dim;
  pdd_dual.dim++;

  // Intentionally swapped because PDD here is given in the dual view.
  vector<double *> &V = pdd_dual.H;
  vector<double *> &H = pdd_dual.V;

  if (V.empty()) {
    // There are no vertices - the further convex hull with an empty vertex set
    // would not change the vertex set, thus I'm leaving it empty.
    assert(H.empty() && "Consistency checking that H is empty.");
    return;
  }
  vector<set_t> &incidence = pdd_dual.incidence;

  for (size_t i = 0; i < V.size(); i++) {
    V[i] = (double *)realloc(V[i], sizeof(double) * (dim + 1));
    if (polarity == MINUS) {
      V[i][dim] = 0;
    } else {
      V[i][dim] = V[i][xi + 1];
    }
  }
  for (size_t i = 0; i < H.size(); i++) {
    H[i] = (double *)realloc(H[i], sizeof(double) * (dim + 1));
    H[i][dim] = 0;
  }
  H.resize(H.size() + 2);
  H[H.size() - 2] = (double *)calloc(dim + 1, sizeof(double));
  H[H.size() - 1] = (double *)calloc(dim + 1, sizeof(double));

  H[H.size() - 2][dim] = 1;
  H[H.size() - 1][dim] = -1;
  if (polarity == PLUS) {
    H[H.size() - 2][xi + 1] = -1;
    H[H.size() - 1][xi + 1] = 1;
  }

  // Incidence size equals previous size of H which is now increased by 2.
  incidence.resize(H.size());
  incidence[incidence.size() - 2] = set_create(V.size());
  incidence[incidence.size() - 1] = set_create(V.size());
  set_enable_all(incidence[incidence.size() - 2]);
  set_enable_all(incidence[incidence.size() - 1]);
}

void project_to_tanh_y_branch(const int xi, PDD &pdd_dual,
                              const double x_bound) {
  ASRTF(x_bound != 0, "x_bound has to be either negative or positive.");
  const int dim = pdd_dual.dim;
  pdd_dual.dim++;

  // Intentionally swapped because PDD here is given in the dual view.
  vector<double *> &V = pdd_dual.H;
  vector<double *> &H = pdd_dual.V;

  if (V.empty()) {
    // There are no vertices - the further convex hull with an empty vertex set
    // would not change the vertex set, thus I'm leaving it empty.
    assert(H.empty() && "Consistency checking that H is empty.");
    return;
  }
  vector<set_t> &incidence = pdd_dual.incidence;

  double y_bound = tanh(x_bound);
  double dy = 1 - y_bound * y_bound;
  double b = y_bound - x_bound * dy;
  double k = abs(y_bound / x_bound);

  vector<double *> V_new;
  V_new.reserve(V.size() * 2);
  vector<int> map1(V.size());
  vector<int> map2(V.size());

  for (size_t i = 0; i < V.size(); i++) {
    double *v = (double *)realloc(V[i], sizeof(double) * (dim + 1));
    double x_cur = v[xi + 1];
    if (v[xi + 1] == x_bound) {
      // Both lower and upper bound would map to the same y.
      v[dim] = y_bound;
      map1[i] = V_new.size();
      map2[i] = V_new.size();
      V_new.push_back(v);
      continue;
    }
    double *v2 = fp_arr_copy(dim + 1, v);
    v[dim] = dy * x_cur + b;
    v2[dim] = k * x_cur;

    map1[i] = V_new.size();
    V_new.push_back(v);
    map2[i] = V_new.size();
    V_new.push_back(v2);
  }

  for (size_t i = 0; i < H.size(); i++) {
    H[i] = (double *)realloc(H[i], sizeof(double) * (dim + 1));
    H[i][dim] = 0;
  }
  H.resize(H.size() + 2);
  H[H.size() - 2] = (double *)calloc(dim + 1, sizeof(double));
  H[H.size() - 1] = (double *)calloc(dim + 1, sizeof(double));

  double *h1 = H[H.size() - 2];
  double *h2 = H[H.size() - 1];

  // Note that order h1 and h2 matters because according
  // to it the incidence information will be adjusted.
  if (x_bound < 0) {
    // In the minus branch:
    // y >= dy * x + b equivalent -b - dy * x + y >= 0
    h1[0] = -b;
    h1[xi + 1] = -dy;
    h1[dim] = 1;
    // y <= k * x equivalent k * x - y >= 0
    h2[xi + 1] = k;
    h2[dim] = -1;
  } else {
    // In the plus branch:
    // y <= dy * x + b equivalent b + dy * x - y >= 0
    h1[0] = b;
    h1[xi + 1] = dy;
    h1[dim] = -1;
    // y >= k * x equivalent -k * x + y >= 0
    h2[xi + 1] = -k;
    h2[dim] = 1;
  }

  vector<set_t> incidence_new = set_arr_create(H.size(), V_new.size());

  for (size_t h = 0; h < H.size() - 2; h++) {
    set_t inc = incidence[h];
    set_t inc_new = incidence_new[h];
    for (size_t v = 0; v < V.size(); v++) {
      if (set_test_bit(inc, v)) {
        set_enable_bit(inc_new, map1[v]);
        set_enable_bit(inc_new, map2[v]);
      }
    }
  }

  set_t inc1 = incidence_new[H.size() - 2];
  set_t inc2 = incidence_new[H.size() - 1];

  for (int v : map1) {
    set_enable_bit(inc1, v);
  }
  for (int v : map2) {
    set_enable_bit(inc2, v);
  }

  pdd_dual.H = V_new;
  // Since incidence is reference it is important
  // that I free it _before_ I update incidence in pdd.
  set_arr_free(incidence);
  pdd_dual.incidence = incidence_new;
}

PDD decomposition_recursive(Quadrant &quadrant,
                            const map<Quadrant, PDD> &quadrant2pdd, const int K,
                            Activation activation, const vector<double> &x_lb,
                            const vector<double> &x_ub) {
  const int xi = (int)quadrant.size();
  if (xi == K) {
    PDD pdd = quadrant2pdd.at(quadrant);
    PDD_debug_consistency_check(pdd);
    return pdd;
  } else {
    quadrant.push_back(MINUS);
    PDD pdd_minus = decomposition_recursive(quadrant, quadrant2pdd, K,
                                            activation, x_lb, x_ub);
    quadrant.back() = PLUS;
    PDD pdd_plus = decomposition_recursive(quadrant, quadrant2pdd, K,
                                           activation, x_lb, x_ub);
    quadrant.pop_back();

    if (activation == Relu) {
      project_to_relu_y_branch(xi, pdd_minus, MINUS);
      project_to_relu_y_branch(xi, pdd_plus, PLUS);
    } else {
      project_to_tanh_y_branch(xi, pdd_minus, x_lb[xi]);
      project_to_tanh_y_branch(xi, pdd_plus, x_ub[xi]);
    }

    PDD_debug_consistency_check(pdd_minus);
    PDD_debug_consistency_check(pdd_plus);

    PDD res = PDD_intersect_two_PDDs(pdd_minus, pdd_plus);

    return res;
  }
}

// TODO[gleb] Add support for multiple passes.
vector<double *> decomposition(const int K,
                               const map<Quadrant, PDD> &quadrant2pdd,
                               Activation activation,
                               const vector<double> &x_lb,
                               const vector<double> &x_ub) {
  ASRTF(2 <= K && K <= 5, "Only 2 <= K <= 5 are currently supported.");
  ASRTF((int)quadrant2pdd.size() == POW2[K],
        "Sanity check - the number of quadrants should be 2^K.");
  ASRTF(activation == Relu || activation == Tanh,
        "Activation should be Relu or Tanh.");

  Quadrant quadrant{};
  quadrant.reserve(K);
  PDD res = decomposition_recursive(quadrant, quadrant2pdd, K, activation, x_lb,
                                    x_ub);
  vector<double *> &H = res.V;
  vector<double *> &V = res.H;

  PDD_adjust_H_for_soundness_finite_polytope(2 * K + 1, H, V);

  // H is given in order (1, x1, ..., xk, yk, ..., y1) thus last k cols have to
  // be reversed. The desired order is (1, x1, ..., xk, y1, ..., yk).
  for (size_t hi = 0; hi < H.size(); hi++) {
    int first = K + 1;
    int last = 2 * K;
    double *h = H[hi];
    while (first < last) {
      swap(h[first], h[last]);
      first++;
      last--;
    }
  }

  fp_mat_free(V);
  set_arr_free(res.incidence);

  return H;
}

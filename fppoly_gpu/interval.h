#pragma once

__device__ void interval_mul(float_type &a_inf, float_type &a_sup,
                             const float_type b_inf, const float_type b_sup,
                             const float_type c_inf, const float_type c_sup) {
  float_type inf_inf = b_inf * c_inf;
  float_type inf_sup = b_inf * c_sup;
  float_type sup_inf = b_sup * c_inf;
  float_type sup_sup = b_sup * c_sup;

  if (inf_inf < inf_sup) {
    if (sup_inf < sup_sup) {
      if (inf_inf < sup_inf) {
        a_inf = inf_inf;
      } else {
        a_inf = sup_inf;
      }

      if (inf_sup < sup_sup) {
        a_sup = sup_sup;
      } else {
        a_sup = inf_sup;
      }
    } else {
      if (inf_inf < sup_sup) {
        a_inf = inf_inf;
      } else {
        a_inf = sup_sup;
      }

      if (inf_sup < sup_inf) {
        a_sup = sup_inf;
      } else {
        a_sup = inf_sup;
      }
    }
  } else {
    if (sup_inf < sup_sup) {
      if (inf_sup < sup_inf) {
        a_inf = inf_sup;
      } else {
        a_inf = sup_inf;
      }

      if (inf_inf < sup_sup) {
        a_sup = sup_sup;
      } else {
        a_sup = inf_inf;
      }
    } else {
      if (inf_sup < sup_sup) {
        a_inf = inf_sup;
      } else {
        a_inf = sup_sup;
      }

      if (inf_inf < sup_inf) {
        a_sup = sup_inf;
      } else {
        a_sup = inf_inf;
      }
    }
  }
}

__device__ void interval_mul_symmetric_c(float_type &a_inf, float_type &a_sup,
                                         const float_type b_inf,
                                         const float_type b_sup,
                                         const float_type c) {
  if (-b_inf < b_sup) {
    a_inf = b_sup * -c;
    a_sup = b_sup * c;
  } else {
    a_inf = b_inf * c;
    a_sup = b_inf * -c;
  }
}

__device__ void interval_mul_const_c(float_type &a_inf, float_type &a_sup,
                                     const float_type b_inf,
                                     const float_type b_sup,
                                     const float_type c) {
  if (c >= 0) {
    a_inf = b_inf * c;
    a_sup = b_sup * c;
  } else {
    a_inf = b_sup * c;
    a_sup = b_inf * c;
  }
}

__device__ void
interval_mul_expr_coeff(float_type &res_inf, float_type &res_sup,
                        const float_type inf, const float_type sup,
                        const float_type inf_expr, const float_type sup_expr) {
  interval_mul(res_inf, res_sup, inf, sup, inf_expr, sup_expr);

  const float_type maxA = max(abs(inf_expr), abs(sup_expr));
  float_type tmp1, tmp2;

  interval_mul_symmetric_c(tmp1, tmp2, inf, sup, maxA * ulp);

  res_inf += tmp1;
  res_sup += tmp2;
}

__device__ void interval_mul_cst_coeff(float_type &res_inf, float_type &res_sup,
                                       const float_type inf,
                                       const float_type sup,
                                       const float_type inf_expr,
                                       const float_type sup_expr) {
  interval_mul_expr_coeff(res_inf, res_sup, inf, sup, inf_expr, sup_expr);

  res_inf -= min_denormal;
  res_sup += min_denormal;
}

__device__ void interval_mul_expr_coeff_const_expr(float_type &res_inf,
                                                   float_type &res_sup,
                                                   const float_type inf,
                                                   const float_type sup,
                                                   const float_type expr) {
  interval_mul_const_c(res_inf, res_sup, inf, sup, expr);

  float_type tmp1, tmp2;

  interval_mul_symmetric_c(tmp1, tmp2, inf, sup, abs(expr) * ulp);

  res_inf += tmp1;
  res_sup += tmp2;
}

__device__ void interval_mul_cst_coeff_const_expr(float_type &res_inf,
                                                  float_type &res_sup,
                                                  const float_type inf,
                                                  const float_type sup,
                                                  const float_type expr) {
  interval_mul_expr_coeff_const_expr(res_inf, res_sup, inf, sup, expr);

  res_inf -= min_denormal;
  res_sup += min_denormal;
}

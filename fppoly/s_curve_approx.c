#include "s_curve_approx.h"

void compute_chord_slope(double *slope_inf, double *slope_sup, double f_sup_l, double f_sup_u, 
			 double f_inf_l, double f_inf_u, double inf_l, double inf_u, double sup_l, double sup_u){
	double num_l =  f_sup_l + f_inf_u;
	double num_u =  f_sup_u + f_inf_l;

	double den_l = sup_l + inf_u; 
	double den_u = sup_u + inf_l;

        elina_double_interval_div(slope_inf, slope_sup, num_l, num_u, den_l, den_u);
}

void compute_derivative(double *slope_inf, double *slope_sup, double s_curve_l, double s_curve_u, double sq_l, double sq_u, bool is_sigmoid){
    double sq_den_sup_l, sq_den_sup_u;
    elina_double_interval_mul(&sq_den_sup_l, &sq_den_sup_u, sq_l, sq_u, sq_l, sq_u);
    if(is_sigmoid){
            elina_double_interval_div(slope_inf, slope_sup, s_curve_l, s_curve_u, sq_den_sup_l, sq_den_sup_u);
    }
    else{
            *slope_inf = -1 + sq_den_sup_u;
            *slope_sup = 1 + sq_den_sup_l;
    }

}

void compute_slope_and_intercept_s_curve_lexpr(
    fppoly_internal_t *pr, double *slope_inf, double *slope_sup,
    double *intercept_inf, double *intercept_sup, double inf_coeff,
    double sup_coeff, double lb, double ub, bool is_sigmoid, bool *boxify) {

  fesetround(FE_DOWNWARD);
  double e_sup_l = is_sigmoid ? -exp(ub) : -tanh(ub);
  double e_inf_l = is_sigmoid ? -exp(-lb) : -tanh(-lb);
  fesetround(FE_UPWARD);
  double e_sup_u = is_sigmoid ? exp(ub) : tanh(ub);
  double e_inf_u = is_sigmoid ? exp(-lb) : tanh(-lb);
  double f_sup_l, f_sup_u;
  double f_inf_l, f_inf_u;
  double den_sup_l, den_sup_u;
  double den_inf_l, den_inf_u;
  if (is_sigmoid) {
    den_sup_l = -1 + e_sup_l;
    den_sup_u = 1 + e_sup_u;
    den_inf_l = -1 + e_inf_l;
    den_inf_u = 1 + e_inf_u;
    elina_double_interval_div(&f_sup_l, &f_sup_u, e_sup_l, e_sup_u, den_sup_l,
                              den_sup_u);
    elina_double_interval_div(&f_inf_l, &f_inf_u, e_inf_l, e_inf_u, den_inf_l,
                              den_inf_u);
  } else {
    f_inf_l = e_inf_l;
    f_inf_u = e_inf_u;
    f_sup_l = e_sup_l;
    f_sup_u = e_sup_u;
    den_inf_l = e_inf_l;
    den_inf_u = e_inf_u;
    den_sup_l = e_sup_l;
    den_sup_u = e_sup_u;
  }

  if ((-lb == ub) || (-f_inf_l == f_sup_u)) {
    *slope_inf = 0.0;
    *slope_sup = 0.0;
    //
    double tmp1, tmp2;
    elina_double_interval_mul(&tmp1, &tmp2, inf_coeff, sup_coeff, f_inf_l,
                              f_sup_u);
    *intercept_inf = tmp1;
    *intercept_sup = tmp2;
    *boxify = true;
  } else if (sup_coeff < 0 || inf_coeff < 0) {
    double add_inf, add_sup;
    double mul_inf, mul_sup;
    double x_l, x_u;
    double f_x_l, f_x_u;
    if (sup_coeff < 0) {
      if (ub < 0) {
        compute_chord_slope(slope_inf, slope_sup, f_sup_l, f_sup_u, f_inf_l,
                            f_inf_u, lb, -lb, -ub, ub);
        x_l = ub;
        x_u = -ub;
        f_x_l = f_sup_l;
        f_x_u = f_sup_u;

        elina_double_interval_mul_cst_coeff(pr, &add_inf, &add_sup, x_l, x_u,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, intercept_inf, intercept_sup,
                                            f_x_l, f_x_u, add_inf, add_sup);
        double tmp1, tmp2, tmp3, tmp4;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, lb, -lb,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, &tmp3, &tmp4, *intercept_inf,
                                            *intercept_sup, tmp1, tmp2);
        if (tmp4 < f_inf_u) {

          *boxify = true;
        }

      } else if (lb <= 0) {

        compute_derivative(slope_inf, slope_sup, e_sup_l, e_sup_u, den_sup_l,
                           den_sup_u, is_sigmoid);
        x_l = ub;
        x_u = -ub;
        f_x_l = f_sup_l;
        f_x_u = f_sup_u;
        elina_double_interval_mul_cst_coeff(pr, &add_inf, &add_sup, x_l, x_u,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, intercept_inf, intercept_sup,
                                            f_x_l, f_x_u, add_inf, add_sup);
        double tmp1, tmp2, tmp3, tmp4;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, lb, -lb,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, &tmp3, &tmp4, *intercept_inf,
                                            *intercept_sup, tmp1, tmp2);
        if (-tmp3 < f_inf_u) {

          *boxify = true;
        }
      } else {

        if (lb <= ub) {
          // double slope_inf1, slope_sup1;
          // double slope_inf2, slope_sup2;
          compute_derivative(slope_inf, slope_sup, e_sup_l, e_sup_u, den_sup_l,
                             den_sup_u, is_sigmoid);

        } else {
          compute_derivative(slope_inf, slope_sup, e_inf_l, e_inf_u, den_inf_l,
                             den_inf_u, is_sigmoid);
        }
        x_l = ub;
        x_u = -ub;
        f_x_l = f_sup_l;
        f_x_u = f_sup_u;
        elina_double_interval_mul_cst_coeff(pr, &add_inf, &add_sup, x_l, x_u,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, intercept_inf, intercept_sup,
                                            f_x_l, f_x_u, add_inf, add_sup);
        double tmp1, tmp2, tmp3, tmp4;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, lb, -lb,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, &tmp3, &tmp4, *intercept_inf,
                                            *intercept_sup, tmp1, tmp2);
        if (-tmp3 < f_inf_u) {
          *boxify = true;
        }
      }
    } else {
      if (ub < 0) {
        compute_derivative(slope_inf, slope_sup, e_inf_l, e_inf_u, den_inf_l,
                           den_inf_u, is_sigmoid);
        x_l = -lb;
        x_u = lb;
        f_x_l = f_inf_l;
        f_x_u = f_inf_u;
        elina_double_interval_mul_cst_coeff(pr, &add_inf, &add_sup, x_l, x_u,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, intercept_inf, intercept_sup,
                                            f_x_l, f_x_u, add_inf, add_sup);
        double tmp1, tmp2, tmp3, tmp4;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, -ub, ub,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, &tmp3, &tmp4, *intercept_inf,
                                            *intercept_sup, tmp1, tmp2);
        if (tmp4 > -f_sup_l) {
          *boxify = true;
        }
      } else if (lb <= 0) {
        compute_chord_slope(slope_inf, slope_sup, f_sup_l, f_sup_u, f_inf_l,
                            f_inf_u, lb, -lb, -ub, ub);

        x_l = -lb;
        x_u = lb;
        f_x_l = f_inf_l;
        f_x_u = f_inf_u;
        elina_double_interval_mul(&add_inf, &add_sup, x_l, x_u, *slope_inf,
                                  *slope_sup);
        elina_double_interval_add_cst_coeff(pr, intercept_inf, intercept_sup,
                                            f_x_l, f_x_u, add_inf, add_sup);
        double tmp1, tmp2, tmp3, tmp4;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, -ub, ub,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, &tmp3, &tmp4, *intercept_inf,
                                            *intercept_sup, tmp1, tmp2);

        if (-tmp3 > f_sup_u) {
          *boxify = true;
        }
      } else {
        if (lb <= ub) {
          compute_derivative(slope_inf, slope_sup, e_sup_l, e_sup_u, den_sup_l,
                             den_sup_u, is_sigmoid);
        } else {
          compute_derivative(slope_inf, slope_sup, e_inf_l, e_inf_u, den_inf_l,
                             den_inf_u, is_sigmoid);
        }
        x_l = -lb;
        x_u = lb;
        f_x_l = f_inf_l;
        f_x_u = f_inf_u;
        elina_double_interval_mul_cst_coeff(pr, &add_inf, &add_sup, x_l, x_u,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, intercept_inf, intercept_sup,
                                            f_x_l, f_x_u, add_inf, add_sup);
        double tmp1, tmp2, tmp3, tmp4;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, -ub, ub,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, &tmp3, &tmp4, *intercept_inf,
                                            *intercept_sup, tmp1, tmp2);
        if (tmp4 > -f_sup_l) {

          *boxify = true;
        }
      }
    }
    if (*boxify) {
      *slope_inf = 0.0;
      *slope_sup = 0.0;
      double tmp1, tmp2;
      elina_double_interval_mul(&tmp1, &tmp2, inf_coeff, sup_coeff, f_inf_l,
                                f_sup_u);
      *intercept_inf = tmp1;
      *intercept_sup = tmp2;
    }
  } else {

    *slope_inf = 0.0;
    *slope_sup = 0.0;
    double tmp1, tmp2;
    elina_double_interval_mul(&tmp1, &tmp2, inf_coeff, sup_coeff, f_inf_l,
                              f_sup_u);
    *intercept_inf = tmp1;
    *intercept_sup = tmp2;
    *boxify = true;
  }

  return;
}

void compute_slope_and_intercept_s_curve_uexpr(
    fppoly_internal_t *pr, double *slope_inf, double *slope_sup,
    double *intercept_inf, double *intercept_sup, double inf_coeff,
    double sup_coeff, double lb, double ub, bool is_sigmoid, bool *boxify) {
  fesetround(FE_DOWNWARD);
  double e_sup_l = is_sigmoid ? -exp(ub) : -tanh(ub);
  double e_inf_l = is_sigmoid ? -exp(-lb) : -tanh(-lb);

  fesetround(FE_UPWARD);
  double e_sup_u = is_sigmoid ? exp(ub) : tanh(ub);
  double e_inf_u = is_sigmoid ? exp(-lb) : tanh(-lb);

  double f_sup_l, f_sup_u;
  double f_inf_l, f_inf_u;
  double den_sup_l, den_sup_u;
  double den_inf_l, den_inf_u;
  double connecting_slope_l, connecting_slope_u;

  if (is_sigmoid) {
    den_sup_l = -1 + e_sup_l;
    den_sup_u = 1 + e_sup_u;
    den_inf_l = -1 + e_inf_l;
    den_inf_u = 1 + e_inf_u;
    elina_double_interval_div(&f_sup_l, &f_sup_u, e_sup_l, e_sup_u, den_sup_l,
                              den_sup_u);
    elina_double_interval_div(&f_inf_l, &f_inf_u, e_inf_l, e_inf_u, den_inf_l,
                              den_inf_u);
  } else {
    f_inf_l = e_inf_l;
    f_inf_u = e_inf_u;
    f_sup_l = e_sup_l;
    f_sup_u = e_sup_u;
    den_inf_l = e_inf_l;
    den_inf_u = e_inf_u;
    den_sup_l = e_sup_l;
    den_sup_u = e_sup_u;
  }

  if ((-lb == ub) || (-f_inf_l == f_sup_u)) {
    *slope_inf = 0.0;
    *slope_sup = 0.0;
    *boxify = true;
    double tmp1, tmp2;
    elina_double_interval_mul(&tmp1, &tmp2, inf_coeff, sup_coeff, f_inf_l,
                              f_sup_u);
    *intercept_inf = tmp1;
    *intercept_sup = tmp2;

  }

  else if (sup_coeff < 0 || inf_coeff < 0) {
    double add_inf, add_sup;
    double mul_inf, mul_sup;
    double x_l, x_u;
    double f_x_l, f_x_u;

    if (sup_coeff < 0) {
      if (ub < 0) {

        compute_derivative(slope_inf, slope_sup, e_inf_l, e_inf_u, den_inf_l,
                           den_inf_u, is_sigmoid);
        x_l = -lb;
        x_u = lb;
        f_x_l = f_inf_l;
        f_x_u = f_inf_u;
        elina_double_interval_mul_cst_coeff(pr, &add_inf, &add_sup, x_l, x_u,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, intercept_inf, intercept_sup,
                                            f_x_l, f_x_u, add_inf, add_sup);
        double tmp1, tmp2, tmp3, tmp4;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, -ub, ub,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, &tmp3, &tmp4, *intercept_inf,
                                            *intercept_sup, tmp1, tmp2);
        if (tmp4 > -f_sup_l) {
          *boxify = true;
        }
      } else if (lb <= 0) {

        compute_chord_slope(slope_inf, slope_sup, f_sup_l, f_sup_u, f_inf_l,
                            f_inf_u, lb, -lb, -ub, ub);

        x_l = -lb;
        x_u = lb;
        f_x_l = f_inf_l;
        f_x_u = f_inf_u;
        elina_double_interval_mul_cst_coeff(pr, &add_inf, &add_sup, x_l, x_u,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, intercept_inf, intercept_sup,
                                            f_x_l, f_x_u, add_inf, add_sup);
        double tmp1, tmp2, tmp3, tmp4;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, -ub, ub,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, &tmp3, &tmp4, *intercept_inf,
                                            *intercept_sup, tmp1, tmp2);
        if (-tmp3 > f_sup_u) {
          *boxify = true;
        }
        //}
      } else {
        // double slope_inf1, slope_sup1;
        // double slope_inf2, slope_sup2;
        if (lb <= ub) {
          compute_derivative(slope_inf, slope_sup, e_sup_l, e_sup_u, den_sup_l,
                             den_sup_u, is_sigmoid);

        } else {
          compute_derivative(slope_inf, slope_sup, e_inf_l, e_inf_u, den_inf_l,
                             den_inf_u, is_sigmoid);
        }

        x_l = -lb;
        x_u = lb;
        f_x_l = f_inf_l;
        f_x_u = f_inf_u;
        elina_double_interval_mul_cst_coeff(pr, &add_inf, &add_sup, x_l, x_u,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, intercept_inf, intercept_sup,
                                            f_x_l, f_x_u, add_inf, add_sup);
        double tmp1, tmp2, tmp3, tmp4;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, -ub, ub,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, &tmp3, &tmp4, *intercept_inf,
                                            *intercept_sup, tmp1, tmp2);
        if (tmp4 > -f_sup_l) {
          *boxify = true;
        }
      }

    } else {
      if (ub < 0) {

        compute_chord_slope(slope_inf, slope_sup, f_sup_l, f_sup_u, f_inf_l,
                            f_inf_u, lb, -lb, -ub, ub);

        x_l = ub;
        x_u = -ub;
        f_x_l = f_sup_l;
        f_x_u = f_sup_u;
        elina_double_interval_mul_cst_coeff(pr, &add_inf, &add_sup, x_l, x_u,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, intercept_inf, intercept_sup,
                                            f_x_l, f_x_u, add_inf, add_sup);
        double tmp1, tmp2, tmp3, tmp4;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, lb, -lb,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, &tmp3, &tmp4, *intercept_inf,
                                            *intercept_sup, tmp1, tmp2);

        if (tmp4 < f_inf_u) {
          *boxify = true;
        }
      } else if (lb <= 0) {

        compute_derivative(slope_inf, slope_sup, e_sup_l, e_sup_u, den_sup_l,
                           den_sup_u, is_sigmoid);

        x_l = ub;
        x_u = -ub;
        f_x_l = f_sup_l;
        f_x_u = f_sup_u;
        elina_double_interval_mul_cst_coeff(pr, &add_inf, &add_sup, x_l, x_u,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, intercept_inf, intercept_sup,
                                            f_x_l, f_x_u, add_inf, add_sup);
        double tmp1, tmp2, tmp3, tmp4;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, lb, -lb,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, &tmp3, &tmp4, *intercept_inf,
                                            *intercept_sup, tmp1, tmp2);

        if (-tmp3 < f_inf_u) {
          *boxify = true;
        }

      } else {

        if (lb <= ub) {
          compute_derivative(slope_inf, slope_sup, e_sup_l, e_sup_u, den_sup_l,
                             den_sup_u, is_sigmoid);
        } else {
          compute_derivative(slope_inf, slope_sup, e_inf_l, e_inf_u, den_inf_l,
                             den_inf_u, is_sigmoid);
        }

        x_l = ub;
        x_u = -ub;
        f_x_l = f_sup_l;
        f_x_u = f_sup_u;
        elina_double_interval_mul_cst_coeff(pr, &add_inf, &add_sup, x_l, x_u,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, intercept_inf, intercept_sup,
                                            f_x_l, f_x_u, add_inf, add_sup);
        double tmp1, tmp2, tmp3, tmp4;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, lb, -lb,
                                            *slope_inf, *slope_sup);
        elina_double_interval_add_cst_coeff(pr, &tmp3, &tmp4, *intercept_inf,
                                            *intercept_sup, tmp1, tmp2);
        if (-tmp3 < f_inf_u) {
          *boxify = true;
        }
      }
    }

    if (*boxify) {
      *slope_inf = 0.0;
      *slope_sup = 0.0;
      double tmp1, tmp2;
      elina_double_interval_mul(&tmp1, &tmp2, inf_coeff, sup_coeff, f_inf_l,
                                f_sup_u);
      *intercept_inf = tmp1;
      *intercept_sup = tmp2;
    }

  }

  else {

    *slope_inf = 0.0;
    *slope_sup = 0.0;
    double tmp1, tmp2;
    elina_double_interval_mul(&tmp1, &tmp2, inf_coeff, sup_coeff, f_inf_l,
                              f_sup_u);
    *intercept_inf = tmp1;
    *intercept_sup = tmp2;
    *boxify = true;
  }
}

expr_t *lexpr_replace_s_curve_bounds(fppoly_internal_t *pr, expr_t *expr,
                                     neuron_t **neurons, bool is_sigmoid) {
  size_t num_neurons = expr->size;
  size_t i, k;
  expr_t *res = alloc_expr();
  res->inf_coeff = (double *)malloc(num_neurons * sizeof(double));
  res->sup_coeff = (double *)malloc(num_neurons * sizeof(double));
  res->inf_cst = expr->inf_cst;
  res->sup_cst = expr->sup_cst;
  res->type = expr->type;
  res->size = num_neurons;

  for (i = 0; i < num_neurons; i++) {
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }
    neuron_t *neuron_k = neurons[k];
    double lb = neurons[k]->lb;
    double ub = neurons[k]->ub;
    // if(expr->sup_coeff[i]<0 || expr->inf_coeff[i] < 0){
    double slope_inf, slope_sup;
    double intercept_inf, intercept_sup;
    double mul_inf, mul_sup;
    bool boxify = false;
    compute_slope_and_intercept_s_curve_lexpr(
        pr, &slope_inf, &slope_sup, &intercept_inf, &intercept_sup,
        expr->inf_coeff[i], expr->sup_coeff[i], lb, ub, is_sigmoid, &boxify);
    if (boxify) {
      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;
      res->inf_cst = res->inf_cst + intercept_inf;
      res->sup_cst = res->sup_cst + intercept_sup;
    } else {
      elina_double_interval_mul_expr_coeff(
          pr, &res->inf_coeff[i], &res->sup_coeff[i], slope_inf, slope_sup,
          expr->inf_coeff[i], expr->sup_coeff[i]);
      elina_double_interval_mul_cst_coeff(pr, &mul_inf, &mul_sup, intercept_inf,
                                          intercept_sup, expr->inf_coeff[i],
                                          expr->sup_coeff[i]);
      elina_double_interval_add_cst_coeff(pr, &res->inf_cst, &res->sup_cst,
                                          mul_inf, mul_sup, res->inf_cst,
                                          res->sup_cst);
    }
  }
  if (expr->type == SPARSE) {
    res->dim = (size_t *)malloc(num_neurons * sizeof(size_t));
    for (i = 0; i < num_neurons; i++) {
      res->dim[i] = expr->dim[i];
    }
  }
  return res;
}

expr_t *uexpr_replace_s_curve_bounds(fppoly_internal_t *pr, expr_t *expr,
                                     neuron_t **neurons, bool is_sigmoid) {
  size_t num_neurons = expr->size;
  size_t i, k;
  expr_t *res = alloc_expr();
  res->inf_coeff = (double *)malloc(num_neurons * sizeof(double));
  res->sup_coeff = (double *)malloc(num_neurons * sizeof(double));
  res->inf_cst = expr->inf_cst;
  res->sup_cst = expr->sup_cst;
  res->type = expr->type;
  res->size = num_neurons;

  for (i = 0; i < num_neurons; i++) {
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }
    neuron_t *neuron_k = neurons[k];
    double lb = neurons[k]->lb;
    double ub = neurons[k]->ub;
    double slope_inf, slope_sup;
    double intercept_inf, intercept_sup;
    double mul_inf, mul_sup;
    bool boxify = false;
    compute_slope_and_intercept_s_curve_uexpr(
        pr, &slope_inf, &slope_sup, &intercept_inf, &intercept_sup,
        expr->inf_coeff[i], expr->sup_coeff[i], lb, ub, is_sigmoid, &boxify);
    if (boxify) {
      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;
      res->inf_cst = res->inf_cst + intercept_inf;
      res->sup_cst = res->sup_cst + intercept_sup;
    } else {
      elina_double_interval_mul_expr_coeff(
          pr, &res->inf_coeff[i], &res->sup_coeff[i], slope_inf, slope_sup,
          expr->inf_coeff[i], expr->sup_coeff[i]);
      elina_double_interval_mul_cst_coeff(pr, &mul_inf, &mul_sup, intercept_inf,
                                          intercept_sup, expr->inf_coeff[i],
                                          expr->sup_coeff[i]);
      elina_double_interval_add_cst_coeff(pr, &res->inf_cst, &res->sup_cst,
                                          mul_inf, mul_sup, res->inf_cst,
                                          res->sup_cst);
    }
  }
  if (expr->type == SPARSE) {
    res->dim = (size_t *)malloc(num_neurons * sizeof(size_t));
    for (i = 0; i < num_neurons; i++) {
      res->dim[i] = expr->dim[i];
    }
  }
  return res;
}

expr_t *uexpr_replace_sigmoid_bounds(fppoly_internal_t *pr, expr_t *expr,
                                     neuron_t **neurons) {
  return uexpr_replace_s_curve_bounds(pr, expr, neurons, true);
}

expr_t *uexpr_replace_tanh_bounds(fppoly_internal_t *pr, expr_t *expr,
                                  neuron_t **neurons) {
  return uexpr_replace_s_curve_bounds(pr, expr, neurons, false);
}

expr_t *lexpr_replace_sigmoid_bounds(fppoly_internal_t *pr, expr_t *expr,
                                     neuron_t **neurons) {
  return lexpr_replace_s_curve_bounds(pr, expr, neurons, true);
}

expr_t *lexpr_replace_tanh_bounds(fppoly_internal_t *pr, expr_t *expr,
                                  neuron_t **neurons) {
  return lexpr_replace_s_curve_bounds(pr, expr, neurons, false);
}

double apply_s_curve_lexpr(fppoly_internal_t *pr, expr_t **lexpr_p,
                           neuron_t *neuron, bool is_sigmoid) {
  expr_t *lexpr = *lexpr_p;
  size_t i;
  size_t size = lexpr->size;
  double lb = neuron->lb;
  double ub = neuron->ub;
  bool boxify = false;
  double slope_inf, slope_sup;

  double intercept_inf, intercept_sup;
  compute_slope_and_intercept_s_curve_lexpr(pr, &slope_inf, &slope_sup,
                                            &intercept_inf, &intercept_sup, -1,
                                            1, lb, ub, is_sigmoid, &boxify);
  fesetround(FE_DOWNWARD);
  double e_inf_l = is_sigmoid ? -exp(-lb) : -tanh(-lb);
  fesetround(FE_UPWARD);
  double e_inf_u = is_sigmoid ? exp(-lb) : tanh(-lb);
  double f_inf_l, f_inf_u;
  double den_inf_l, den_inf_u;
  if (is_sigmoid) {
    den_inf_l = -1 + e_inf_l;
    den_inf_u = 1 + e_inf_u;
    elina_double_interval_div(&f_inf_l, &f_inf_u, e_inf_l, e_inf_u, den_inf_l,
                              den_inf_u);
  } else {
    f_inf_l = e_inf_l;
    f_inf_u = e_inf_u;
  }
  if (boxify) {
    for (i = 0; i < size; i++) {
      lexpr->inf_coeff[i] = 0.0;
      lexpr->sup_coeff[i] = 0.0;
    }
    lexpr->inf_cst = lexpr->inf_cst + intercept_inf;
    lexpr->sup_cst = lexpr->sup_cst + intercept_sup;
  } else {
    double mul_inf, mul_sup;
    for (i = 0; i < size; i++) {
      elina_double_interval_mul_expr_coeff(
          pr, &lexpr->inf_coeff[i], &lexpr->sup_coeff[i], slope_inf, slope_sup,
          lexpr->inf_coeff[i], lexpr->sup_coeff[i]);
    }
    elina_double_interval_mul_cst_coeff(pr, &lexpr->inf_cst, &lexpr->sup_cst,
                                        slope_inf, slope_sup, lexpr->inf_cst,
                                        lexpr->sup_cst);
    elina_double_interval_add_cst_coeff(pr, &lexpr->inf_cst, &lexpr->sup_cst,
                                        intercept_inf, intercept_sup,
                                        lexpr->inf_cst, lexpr->sup_cst);
  }
  return f_inf_l;
}

double apply_s_curve_uexpr(fppoly_internal_t *pr, expr_t **uexpr_p,
                           neuron_t *neuron, bool is_sigmoid) {
  expr_t *uexpr = *uexpr_p;
  size_t i;
  size_t size = uexpr->size;
  double lb = neuron->lb;
  double ub = neuron->ub;
  bool boxify = false;
  double slope_inf, slope_sup;
  double intercept_inf, intercept_sup;
  compute_slope_and_intercept_s_curve_uexpr(pr, &slope_inf, &slope_sup,
                                            &intercept_inf, &intercept_sup, -1,
                                            1, lb, ub, is_sigmoid, &boxify);

  fesetround(FE_DOWNWARD);
  double e_sup_l = is_sigmoid ? -exp(ub) : -tanh(ub);
  fesetround(FE_UPWARD);
  double e_sup_u = is_sigmoid ? exp(ub) : tanh(ub);
  double f_sup_l, f_sup_u;
  double den_sup_l, den_sup_u;
  if (is_sigmoid) {
    den_sup_l = -1 + e_sup_l;
    den_sup_u = 1 + e_sup_u;
    elina_double_interval_div(&f_sup_l, &f_sup_u, e_sup_l, e_sup_u, den_sup_l,
                              den_sup_u);
  } else {
    f_sup_l = e_sup_l;
    f_sup_u = e_sup_u;
  }

  if (boxify) {
    for (i = 0; i < size; i++) {
      uexpr->inf_coeff[i] = 0.0;
      uexpr->sup_coeff[i] = 0.0;
    }
    uexpr->inf_cst = uexpr->inf_cst + intercept_inf;
    uexpr->sup_cst = uexpr->sup_cst + intercept_sup;
  } else {
    double mul_inf, mul_sup;
    for (i = 0; i < size; i++) {
      elina_double_interval_mul_expr_coeff(
          pr, &uexpr->inf_coeff[i], &uexpr->sup_coeff[i], slope_inf, slope_sup,
          uexpr->inf_coeff[i], uexpr->sup_coeff[i]);
    }
    elina_double_interval_mul_cst_coeff(pr, &uexpr->inf_cst, &uexpr->sup_cst,
                                        slope_inf, slope_sup, uexpr->inf_cst,
                                        uexpr->sup_cst);
    elina_double_interval_add_cst_coeff(pr, &uexpr->inf_cst, &uexpr->sup_cst,
                                        intercept_inf, intercept_sup,
                                        uexpr->inf_cst, uexpr->sup_cst);
  }
  return f_sup_u;
}

double apply_sigmoid_lexpr(fppoly_internal_t *pr, expr_t **lexpr_p, neuron_t * neuron){
  return apply_s_curve_lexpr(pr, lexpr_p, neuron, true);
}

double apply_tanh_lexpr(fppoly_internal_t *pr, expr_t **lexpr_p, neuron_t * neuron){
  return apply_s_curve_lexpr(pr, lexpr_p, neuron, false);
}

double apply_sigmoid_uexpr(fppoly_internal_t *pr, expr_t **uexpr_p, neuron_t * neuron){
  return apply_s_curve_uexpr(pr, uexpr_p, neuron, true);
}

double apply_tanh_uexpr(fppoly_internal_t *pr, expr_t **uexpr_p, neuron_t * neuron){
  return apply_s_curve_uexpr(pr, uexpr_p, neuron, false);
}

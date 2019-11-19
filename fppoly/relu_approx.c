#include "relu_approx.h"

expr_t *lexpr_replace_relu_bounds(fppoly_internal_t *pr, expr_t *expr,
                                  neuron_t **neurons, bool use_area_heuristic) {
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
    double width = ub + lb;
    double lambda_inf = -ub / width;
    double lambda_sup = ub / width;
    if ((expr->sup_coeff[i] == 0) && (expr->inf_coeff[i] == 0)) {
      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;
      continue;
    } else if (neuron_k->ub <= 0) {
      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;
      continue;
    } else if (neuron_k->lb < 0) {
      res->inf_coeff[i] = expr->inf_coeff[i];
      res->sup_coeff[i] = expr->sup_coeff[i];
    } else if (expr->sup_coeff[i] < 0) {

      double mu_inf = lambda_inf * neurons[k]->lb;
      double mu_sup = lambda_sup * neurons[k]->lb;
      // res->coeff[i] = lambda*expr->coeff[i];
      // res->cst = res->cst + expr->coeff[i]*mu;
      elina_double_interval_mul_expr_coeff(
          pr, &res->inf_coeff[i], &res->sup_coeff[i], lambda_inf, lambda_sup,
          expr->inf_coeff[i], expr->sup_coeff[i]);
      double tmp1, tmp2;
      elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, mu_inf, mu_sup,
                                          expr->inf_coeff[i],
                                          expr->sup_coeff[i]);
      res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
      res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
    } else if (expr->inf_coeff[i] < 0) {

      double area1 = lb * ub;
      double area2 = 0.5 * ub * width;
      double area3 = 0.5 * lb * width;
      // if((area1 < area2) && (area1 < area3)){
      // if(1){
      // res->coeff[i] = lambda*expr->coeff[i];
      //	elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],lambda_inf,lambda_sup,expr->inf_coeff[i],expr->sup_coeff[i]);

      //}
      if (use_area_heuristic) {
        if ((area2 < area1) && (area2 < area3)) {
          res->inf_coeff[i] = 0.0;
          res->sup_coeff[i] = 0.0;
        } else {
          res->inf_coeff[i] = expr->inf_coeff[i];
          res->sup_coeff[i] = expr->sup_coeff[i];
        }
      } else {
        res->inf_coeff[i] = 0.0;
        res->sup_coeff[i] = 0.0;
      }
    } else {

      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;
      double tmp1, tmp2;
      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i], 0, ub);
      res->inf_cst = res->inf_cst + tmp1;
      res->sup_cst = res->sup_cst + tmp2;
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

expr_t *uexpr_replace_relu_bounds(fppoly_internal_t *pr, expr_t *expr,
                                  neuron_t **neurons, bool use_area_heuristic) {
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
    double width = ub + lb;
    double lambda_inf = -ub / width;
    double lambda_sup = ub / width;
    if ((expr->sup_coeff[i] == 0) && (expr->inf_coeff[i] == 0)) {
      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;
      continue;
    } else if (neuron_k->ub <= 0) {
      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;
      continue;
    } else if (neuron_k->lb < 0) {
      res->inf_coeff[i] = expr->inf_coeff[i];
      res->sup_coeff[i] = expr->sup_coeff[i];
    } else if (expr->inf_coeff[i] < 0) {

      double mu_inf = lambda_inf * neurons[k]->lb;
      double mu_sup = lambda_sup * neurons[k]->lb;
      // res->coeff[i] = lambda*expr->coeff[i];
      // res->cst = res->cst + expr->coeff[i]*mu;
      elina_double_interval_mul_expr_coeff(
          pr, &res->inf_coeff[i], &res->sup_coeff[i], lambda_inf, lambda_sup,
          expr->inf_coeff[i], expr->sup_coeff[i]);
      double tmp1, tmp2;
      elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, mu_inf, mu_sup,
                                          expr->inf_coeff[i],
                                          expr->sup_coeff[i]);
      res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
      res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
    } else if (expr->sup_coeff[i] < 0) {

      double area1 = lb * ub;
      double area2 = 0.5 * ub * width;
      double area3 = 0.5 * lb * width;
      // if((area1 < area2) && (area1 < area3)){
      // if(1){
      // res->coeff[i] = lambda*expr->coeff[i];
      //	elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],lambda_inf,lambda_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
      //}
      if (use_area_heuristic) {
        if ((area2 < area1) && (area2 < area3)) {
          res->inf_coeff[i] = 0.0;
          res->sup_coeff[i] = 0.0;
        } else {
          res->inf_coeff[i] = expr->inf_coeff[i];
          res->sup_coeff[i] = expr->sup_coeff[i];
        }
      } else {
        res->inf_coeff[i] = 0.0;
        res->sup_coeff[i] = 0.0;
      }
      //
    } else {

      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;
      double tmp1, tmp2;
      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i], 0, ub);
      res->inf_cst = res->inf_cst + tmp1;
      res->sup_cst = res->sup_cst + tmp2;
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

double apply_relu_lexpr(fppoly_internal_t *pr, expr_t **lexpr_p,
                        neuron_t *neuron) {
  expr_t *lexpr = *lexpr_p;
  size_t i;
  size_t size = lexpr->size;
  double lb = neuron->lb;
  double ub = neuron->ub;
  double width = lb + ub;
  if (ub < 0) {
    free_expr(*lexpr_p);
    *lexpr_p = NULL;
    return 0;
  }
  if (lb < 0) {
    return lb;
  }
  double area1 = lb * ub;
  double area2 = 0.5 * ub * width;
  if (area1 < area2) {
    double lambda_inf = -ub / width;
    double lambda_sup = ub / width;
    for (i = 0; i < size; i++) {
      // lexpr->coeff[i] = lexpr->coeff[i]*lambda;
      elina_double_interval_mul_expr_coeff(
          pr, &lexpr->inf_coeff[i], &lexpr->sup_coeff[i], lambda_inf,
          lambda_sup, lexpr->inf_coeff[i], lexpr->sup_coeff[i]);
    }
    // lexpr->cst = lexpr->cst*lambda;
    elina_double_interval_mul_cst_coeff(pr, &lexpr->inf_cst, &lexpr->sup_cst,
                                        lambda_inf, lambda_sup, lexpr->inf_cst,
                                        lexpr->sup_cst);
    // double res, res1;

    return -(lambda_inf * lb);
  } else {
    free_expr(*lexpr_p);
    *lexpr_p = NULL;
    return 0;
  }
}

double apply_relu_uexpr(fppoly_internal_t *pr, expr_t **uexpr_p,
                        neuron_t *neuron) {
  expr_t *uexpr = *uexpr_p;
  size_t i;
  size_t size = uexpr->size;
  double lb = neuron->lb;
  double ub = neuron->ub;
  double width = lb + ub;
  if (ub < 0) {
    free_expr(*uexpr_p);
    *uexpr_p = NULL;
    return 0;
  }
  if (lb < 0) {
    return ub;
  }
  double lambda_inf = -ub / width;
  double lambda_sup = ub / width;
  for (i = 0; i < size; i++) {
    // uexpr->coeff[i] = uexpr->coeff[i]*lambda;
    elina_double_interval_mul_expr_coeff(
        pr, &uexpr->inf_coeff[i], &uexpr->sup_coeff[i], lambda_inf, lambda_sup,
        uexpr->inf_coeff[i], uexpr->sup_coeff[i]);
  }
  elina_double_interval_mul_cst_coeff(pr, &uexpr->inf_cst, &uexpr->sup_cst,
                                      lambda_inf, lambda_sup, uexpr->inf_cst,
                                      uexpr->sup_cst);
  double mu_inf = lambda_inf * lb;
  double mu_sup = lambda_sup * lb;
  // uexpr->cst = uexpr->cst*lambda;
  // uexpr->cst = uexpr->cst + lambda*lb;
  uexpr->inf_cst += mu_inf;
  uexpr->sup_cst += mu_sup;
  return ub;
}

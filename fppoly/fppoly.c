/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
 *  This software is distributed under GNU Lesser General Public License Version 3.0.
 *  For more information, see the ELINA project website at:
 *  http://elina.ethz.ch
 *
 *  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
 *  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
 *  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
 *  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 *  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY     
 *  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
 *  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
 *  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
 *  CONTRACT, TORT OR OTHERWISE).
 *
 */

#include "fppoly.h"

fppoly_t* fppoly_of_abstract0(elina_abstract0_t* a)
{
  return (fppoly_t*)a->value;
}

elina_abstract0_t* abstract0_of_fppoly(elina_manager_t* man, fppoly_t* fp)
{
  elina_abstract0_t* r = malloc(sizeof(elina_abstract0_t));
  assert(r);
  r->value = fp;
  r->man = elina_manager_copy(man);
  return r;
}


static inline void fppoly_internal_free(fppoly_internal_t* pr)
{
    if (pr) {
	pr->funid = ELINA_FUNID_UNKNOWN;
	free(pr);
	pr = NULL;
    }
}


static inline fppoly_internal_t* fppoly_internal_alloc(void)
{
    fppoly_internal_t* pr = (fppoly_internal_t*)malloc(sizeof(fppoly_internal_t));
    pr->funid = ELINA_FUNID_UNKNOWN;
    pr->man = NULL;
    pr->funopt = NULL; 
    pr->min_denormal = ldexpl(1.0,-1074);
    pr->ulp = ldexpl(1.0,-52);
    return pr;
}


/* back pointer to our internal structure from the manager */
fppoly_internal_t* fppoly_init_from_manager(elina_manager_t* man, elina_funid_t funid)
{
	
    fppoly_internal_t* pr = (fppoly_internal_t*)man->internal;
    pr->funid = funid;
	
    if (!(pr->man)) pr->man = man;
	
    return pr;
}


elina_manager_t * fppoly_manager_alloc(void){
	void** funptr;
	fesetround(FE_UPWARD);
	fppoly_internal_t *pr = fppoly_internal_alloc();
	elina_manager_t *man = elina_manager_alloc("fppoly",/* Library name */
			"1.0", /* version */
			pr, /* internal structure */
			(void (*)(void*))fppoly_internal_free /* free function for internal */
			);
	funptr = man->funptr;
	funptr[ELINA_FUNID_FREE] = &fppoly_free;
	/* 3.Printing */
	funptr[ELINA_FUNID_FPRINT] = &fppoly_fprint;
	return man;
}

void expr_fprint(FILE *stream, expr_t *expr) {
  if ((expr->inf_coeff == NULL) || (expr->sup_coeff == NULL)) {
    fprintf(stdout, "+ [%g, %g]\n", -expr->inf_cst, expr->sup_cst);
    return;
  }
  size_t size = expr->size;
  size_t i;
  for (i = 0; i < size; i++) {
    if (i == 0) {
      if (expr->type == DENSE) {
        fprintf(stream, "[%g, %g]x0 ", -expr->inf_coeff[0], expr->sup_coeff[0]);
      } else {
        fprintf(stream, "[%g, %g]x%zu ", -expr->inf_coeff[0],
                expr->sup_coeff[0], expr->dim[0]);
      }
    }

    else {
      if (expr->type == DENSE) {
        fprintf(stream, "+ [%g, %g]x%zu ", -expr->inf_coeff[i],
                expr->sup_coeff[i], i);
      } else {
        fprintf(stream, "+ [%g, %g]x%zu ", -expr->inf_coeff[i],
                expr->sup_coeff[i], expr->dim[i]);
      }
    }
  }

  fprintf(stdout, "+ [%g, %g]\n", -expr->inf_cst, expr->sup_cst);
}

void elina_double_interval_mul_expr_coeff(fppoly_internal_t *pr,
                                          double *res_inf, double *res_sup,
                                          double inf, double sup,
                                          double inf_expr, double sup_expr) {
  elina_double_interval_mul(res_inf, res_sup, inf, sup, inf_expr, sup_expr);
  double maxA = fmax(fabs(inf_expr), fabs(sup_expr));
  double tmp1, tmp2;
  elina_double_interval_mul(&tmp1, &tmp2, inf, sup, maxA * pr->ulp,
                            maxA * pr->ulp);
  *res_inf += tmp1;
  *res_sup += tmp2;
}

void elina_double_interval_mul_cst_coeff(fppoly_internal_t *pr, double *res_inf,
                                         double *res_sup, double inf,
                                         double sup, double inf_expr,
                                         double sup_expr) {
  elina_double_interval_mul_expr_coeff(pr, res_inf, res_sup, inf, sup, inf_expr,
                                       sup_expr);
  *res_inf += pr->min_denormal;
  *res_sup += pr->min_denormal;
}

void expr_print(expr_t *expr) { expr_fprint(stdout, expr); }

expr_t *alloc_expr(void) {
  expr_t *expr = (expr_t *)malloc(sizeof(expr_t));
  expr->inf_coeff = NULL;
  expr->sup_coeff = NULL;
  expr->dim = NULL;
  return expr;
}

expr_t *create_dense_expr(double *coeff, double cst, size_t size) {
  expr_t *expr = (expr_t *)malloc(sizeof(expr_t));
  expr->inf_coeff = (double *)malloc(size * sizeof(double));
  expr->sup_coeff = (double *)malloc(size * sizeof(double));
  expr->dim = NULL;
  size_t i;
  expr->size = size;
  expr->inf_cst = -cst;
  expr->sup_cst = cst;
  expr->type = DENSE;
  for (i = 0; i < size; i++) {
    expr->inf_coeff[i] = -coeff[i];
    expr->sup_coeff[i] = coeff[i];
  }
  return expr;
}

expr_t *create_cst_expr(double l, double u) {
  expr_t *expr = (expr_t *)malloc(sizeof(expr_t));
  expr->inf_coeff = NULL;
  expr->sup_coeff = NULL;
  expr->dim = NULL;
  expr->type = SPARSE;
  expr->size = 0;
  expr->inf_cst = l;
  expr->sup_cst = u;
  return expr;
}

expr_t *create_sparse_expr(double *coeff, double cst, size_t *dim,
                           size_t size) {
  expr_t *expr = (expr_t *)malloc(sizeof(expr_t));
  if (size > 0) {
    expr->inf_coeff = (double *)malloc(size * sizeof(double));
    expr->sup_coeff = (double *)malloc(size * sizeof(double));
    expr->dim = (size_t *)malloc(size * sizeof(size_t));
  } else {
    expr->inf_coeff = NULL;
    expr->sup_coeff = NULL;
    expr->dim = NULL;
  }
  size_t i;
  expr->size = size;
  expr->inf_cst = -cst;
  expr->sup_cst = cst;
  expr->type = SPARSE;
  for (i = 0; i < size; i++) {
    expr->inf_coeff[i] = -coeff[i];
    expr->sup_coeff[i] = coeff[i];
    expr->dim[i] = dim[i];
  }
  return expr;
}

void free_expr(expr_t *expr) {
  if (expr->inf_coeff) {
    free(expr->inf_coeff);
    expr->inf_coeff = NULL;
  }
  if (expr->sup_coeff) {
    free(expr->sup_coeff);
    expr->sup_coeff = NULL;
  }
  if (expr->type == SPARSE && expr->dim) {
    free(expr->dim);
  }
  expr->dim = NULL;
  free(expr);
  expr = NULL;
}

expr_t *copy_cst_expr(expr_t *src) {
  expr_t *dst = (expr_t *)malloc(sizeof(expr_t));
  dst->inf_coeff = NULL;
  dst->sup_coeff = NULL;
  dst->inf_cst = src->inf_cst;
  dst->sup_cst = src->sup_cst;
  dst->type = src->type;
  dst->dim = NULL;
  dst->size = src->size;
  return dst;
}

expr_t *copy_expr(expr_t *src) {
  expr_t *dst = (expr_t *)malloc(sizeof(expr_t));
  dst->inf_coeff = (double *)malloc(src->size * sizeof(double));
  dst->sup_coeff = (double *)malloc(src->size * sizeof(double));

  size_t i;
  dst->inf_cst = src->inf_cst;
  dst->sup_cst = src->sup_cst;
  dst->type = src->type;
  for (i = 0; i < src->size; i++) {
    dst->inf_coeff[i] = src->inf_coeff[i];
    dst->sup_coeff[i] = src->sup_coeff[i];
  }
  if (src->type == SPARSE) {
    dst->dim = (size_t *)malloc(src->size * sizeof(size_t));
    for (i = 0; i < src->size; i++) {
      dst->dim[i] = src->dim[i];
    }
  }
  dst->size = src->size;
  return dst;
}

expr_t *concretize_dense_sub_expr(fppoly_internal_t *pr, expr_t *expr,
                                  double *inf, double *sup, size_t start,
                                  size_t size) {
  expr_t *res = (expr_t *)malloc(sizeof(expr_t));
  res->inf_coeff = (double *)malloc(start * sizeof(double));
  res->sup_coeff = (double *)malloc(start * sizeof(double));
  size_t i;
  res->inf_cst = expr->inf_cst;
  res->sup_cst = expr->sup_cst;
  res->type = expr->type;
  for (i = 0; i < start; i++) {
    res->inf_coeff[i] = expr->inf_coeff[i];
    res->sup_coeff[i] = expr->sup_coeff[i];
  }
  for (i = start; i < size; i++) {
    double tmp1, tmp2;
    elina_double_interval_mul_expr_coeff(pr, &tmp1, &tmp2, inf[i - start],
                                         sup[i - start], expr->inf_coeff[i],
                                         expr->sup_coeff[i]);
    res->inf_cst += tmp1;
    res->sup_cst += tmp2;
  }
  res->size = start;
  return res;
}

void merge_sparse_expr(expr_t *expr, size_t l, size_t m, size_t r) {
  int i, j, k;
  int n1 = m - l + 1;
  int n2 = r - m;

  /* create temp arrays */
  size_t *L = (size_t *)malloc(n1 * sizeof(size_t));
  size_t *R = (size_t *)malloc(n2 * sizeof(size_t));
  double *L2 = (double *)malloc(n1 * sizeof(double));
  double *R2 = (double *)malloc(n2 * sizeof(double));
  double *L3 = (double *)malloc(n1 * sizeof(double));
  double *R3 = (double *)malloc(n2 * sizeof(double));

  /* Copy data to temp arrays L[] and R[] */
  for (i = 0; i < n1; i++) {
    L[i] = expr->dim[l + i];
    L2[i] = expr->inf_coeff[l + i];
    L3[i] = expr->sup_coeff[l + i];
  }
  for (j = 0; j < n2; j++) {
    R[j] = expr->dim[m + 1 + j];
    R2[j] = expr->inf_coeff[m + 1 + j];
    R3[j] = expr->sup_coeff[m + 1 + j];
  }

  /* Merge the temp arrays back into arr[l..r]*/
  i = 0; // Initial index of first subarray
  j = 0; // Initial index of second subarray
  k = l; // Initial index of merged subarray
  while (i < n1 && j < n2) {
    if (L[i] <= R[j]) {
      expr->dim[k] = L[i];
      expr->inf_coeff[k] = L2[i];
      expr->sup_coeff[k] = L3[i];
      i++;
    } else {
      expr->dim[k] = R[j];
      expr->inf_coeff[k] = R2[j];
      expr->sup_coeff[k] = R3[j];
      j++;
    }
    k++;
  }

  /* Copy the remaining elements of L[], if there
     are any */
  while (i < n1) {
    expr->dim[k] = L[i];
    expr->inf_coeff[k] = L2[i];
    expr->sup_coeff[k] = L3[i];
    i++;
    k++;
  }

  /* Copy the remaining elements of R[], if there
     are any */
  while (j < n2) {
    expr->dim[k] = R[j];
    expr->inf_coeff[k] = R2[j];
    expr->sup_coeff[k] = R3[j];
    j++;
    k++;
  }
  free(L);
  free(R);
  free(L2);
  free(R2);
  free(L3);
  free(R3);
}

/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void merge_sort_sparse_expr(expr_t *expr, size_t l, size_t r) {
  if (l < r) {
    // Same as (l+r)/2, but avoids overflow for
    // large l and h
    size_t m = l + (r - l) / 2;

    // Sort first and second halves
    merge_sort_sparse_expr(expr, l, m);
    merge_sort_sparse_expr(expr, m + 1, r);

    merge_sparse_expr(expr, l, m, r);
  }
}

void sort_sparse_expr(expr_t *expr) {
  merge_sort_sparse_expr(expr, 0, expr->size - 1);
}

neuron_t *neuron_alloc(void){
	neuron_t *res =  (neuron_t *)malloc(sizeof(neuron_t));
        res->expr = NULL;
        res->lb = -INFINITY;
	res->ub = INFINITY;
	res->lexpr = NULL;
	res->uexpr = NULL;
	return res;
}

layer_t *create_layer(size_t size, layertype_t type,
                      activation_type_t activation) {
  layer_t *layer = (layer_t *)malloc(sizeof(layer_t));
  layer->dims = size;
  layer->type = type;
  layer->activation = activation;
  layer->neurons = (neuron_t **)malloc(size * sizeof(neuron_t *));
  size_t i;
  for (i = 0; i < size; i++) {
    layer->neurons[i] = neuron_alloc();
  }
  layer->h_t_inf = NULL;
  layer->h_t_sup = NULL;
  layer->c_t_inf = NULL;
  layer->c_t_sup = NULL;
  return layer;
}

void fppoly_from_network_input_box(fppoly_t *res, size_t intdim, size_t realdim, double *inf_array, double *sup_array){
	
	res->layers = NULL;
	res->numlayers = 0;
	res->lstm_index = 0;
	size_t num_pixels = intdim + realdim;
	res->input_inf = (double *)malloc(num_pixels*sizeof(double));
	res->input_sup = (double *)malloc(num_pixels*sizeof(double));
	res->input_lexpr = NULL;
	res->input_uexpr = NULL;
	size_t i;
	for(i=0; i < num_pixels; i++){
		res->input_inf[i] = -inf_array[i];
		res->input_sup[i] = sup_array[i];
	}
	res->num_pixels = num_pixels;
        res->out = NULL;
}


elina_abstract0_t * fppoly_from_network_input(elina_manager_t *man, size_t intdim, size_t realdim, double *inf_array, double *sup_array){
	fppoly_t * res = (fppoly_t *)malloc(sizeof(fppoly_t));
	fppoly_from_network_input_box(res, intdim, realdim, inf_array, sup_array);
	return abstract0_of_fppoly(man,res);
}

void fppoly_set_network_input_box(elina_manager_t *man, elina_abstract0_t* element, size_t intdim, size_t realdim, double *inf_array, double * sup_array){
    fppoly_t * res = fppoly_of_abstract0(element);
    size_t num_pixels = intdim + realdim;
    res->numlayers = 0;
    size_t i;
    for(i=0; i < num_pixels; i++){
        res->input_inf[i] = -inf_array[i];
        res->input_sup[i] = sup_array[i];
    }
}

elina_abstract0_t *fppoly_from_network_input_poly(
    elina_manager_t *man, size_t intdim, size_t realdim, double *inf_array,
    double *sup_array, double *lexpr_weights, double *lexpr_cst,
    size_t *lexpr_dim, double *uexpr_weights, double *uexpr_cst,
    size_t *uexpr_dim, size_t expr_size) {
  fppoly_t *res = (fppoly_t *)malloc(sizeof(fppoly_t));

  fppoly_from_network_input_box(res, intdim, realdim, inf_array, sup_array);
  size_t num_pixels = intdim + realdim;
  res->input_lexpr = (expr_t **)malloc(num_pixels * sizeof(expr_t *));
  res->input_uexpr = (expr_t **)malloc(num_pixels * sizeof(expr_t *));

  size_t i;
  double *tmp_weights = (double *)malloc(expr_size * sizeof(double));
  size_t *tmp_dim = (size_t *)malloc(expr_size * sizeof(size_t));

  for (i = 0; i < num_pixels; i++) {

    size_t j;
    for (j = 0; j < expr_size; j++) {
      tmp_weights[j] = lexpr_weights[i * expr_size + j];
      tmp_dim[j] = lexpr_dim[i * expr_size + j];
    }
    res->input_lexpr[i] =
        create_sparse_expr(tmp_weights, lexpr_cst[i], tmp_dim, expr_size);
    sort_sparse_expr(res->input_lexpr[i]);
    // printf("w: %p %g %g %g cst: %g dim: %p %zu %zu
    // %zu\n",lexpr_weights[i],lexpr_weights[i][0],lexpr_weights[i][1],
    // lexpr_weights[i][2],lexpr_cst[i],lexpr_dim[i],lexpr_dim[i][0],lexpr_dim[i][1],
    // lexpr_dim[i][2]); expr_print(res->input_lexpr[i]); fflush(stdout);
    for (j = 0; j < expr_size; j++) {
      tmp_weights[j] = uexpr_weights[i * expr_size + j];
      tmp_dim[j] = uexpr_dim[i * expr_size + j];
    }
    res->input_uexpr[i] =
        create_sparse_expr(tmp_weights, uexpr_cst[i], tmp_dim, expr_size);
    sort_sparse_expr(res->input_uexpr[i]);
    //	expr_print(res->input_uexpr[i]);
    //	fflush(stdout);
  }
  free(tmp_weights);
  free(tmp_dim);
  return abstract0_of_fppoly(man, res);
}

void fppoly_alloc_first_layer(fppoly_t *fp, size_t size, layertype_t type,
                              activation_type_t activation) {
  layer_t *layer = create_layer(size, type, activation);
  fp->layers = (layer_t **)malloc(2000 * sizeof(layer_t *));
  fp->layers[0] = layer;
  fp->numlayers = 1;
  return;
}

void fppoly_add_new_layer(fppoly_t *fp, size_t size, layertype_t type,
                          activation_type_t activation) {
  size_t numlayers = fp->numlayers;
  fp->layers[numlayers] = create_layer(size, type, activation);
  fp->numlayers++;
  return;
}

void elina_double_interval_add_expr_coeff(fppoly_internal_t *pr,
                                          double *res_inf, double *res_sup,
                                          double inf, double sup,
                                          double inf_expr, double sup_expr) {
  *res_inf = inf + inf_expr;
  *res_sup = sup + sup_expr;
  double maxA = fmax(fabs(inf_expr), fabs(sup_expr));
  double tmp1, tmp2;
  elina_double_interval_mul(&tmp1, &tmp2, inf, sup, maxA * pr->ulp,
                            maxA * pr->ulp);
  *res_inf += tmp1;
  *res_sup += tmp2;
}

void elina_double_interval_add_cst_coeff(fppoly_internal_t *pr, double *res_inf,
                                         double *res_sup, double inf,
                                         double sup, double inf_expr,
                                         double sup_expr) {
  elina_double_interval_add_expr_coeff(pr, res_inf, res_sup, inf, sup, inf_expr,
                                       sup_expr);
  *res_inf += pr->min_denormal;
  *res_sup += pr->min_denormal;
}

expr_t *multiply_expr(fppoly_internal_t *pr, expr_t *expr, double mul_inf,
                      double mul_sup) {
  expr_t *res = alloc_expr();
  if (expr->size > 0) {
    res->inf_coeff = malloc(expr->size * sizeof(double));
    res->sup_coeff = malloc(expr->size * sizeof(double));
  } else {
    res->inf_coeff = NULL;
    res->sup_coeff = NULL;
  }
  res->type = expr->type;
  size_t i;
  for (i = 0; i < expr->size; i++) {
    // res->coeff[i] = mul_coeff*expr->coeff[i];
    elina_double_interval_mul_expr_coeff(
        pr, &res->inf_coeff[i], &res->sup_coeff[i], mul_inf, mul_sup,
        expr->inf_coeff[i], expr->sup_coeff[i]);
  }
  if (expr->type == SPARSE) {
    if (expr->size > 0) {
      res->dim = (size_t *)malloc(expr->size * sizeof(size_t));
      for (i = 0; i < expr->size; i++) {
        res->dim[i] = expr->dim[i];
      }
    } else {
      res->dim = NULL;
    }
  }
  res->size = expr->size;

  elina_double_interval_mul_cst_coeff(pr, &res->inf_cst, &res->sup_cst, mul_inf,
                                      mul_sup, expr->inf_cst, expr->sup_cst);

  // res->cst = mul_coeff*expr->cst;
  return res;
}

expr_t *multiply_cst_expr(fppoly_internal_t *pr, expr_t *expr, double mul_inf,
                          double mul_sup) {
  expr_t *res = alloc_expr();
  res->inf_coeff = NULL;
  res->sup_coeff = NULL;
  res->dim = NULL;
  res->type = expr->type;
  res->size = expr->size;
  elina_double_interval_mul_cst_coeff(pr, &res->inf_cst, &res->sup_cst, mul_inf,
                                      mul_sup, expr->inf_cst, expr->sup_cst);
  // res->cst = mul_coeff*expr->cst;
  return res;
}

void add_cst_expr(fppoly_internal_t *pr, expr_t *exprA, expr_t *exprB) {
  double maxA = fmax(fabs(exprA->inf_cst), fabs(exprA->sup_cst));
  double maxB = fmax(fabs(exprB->inf_cst), fabs(exprB->sup_cst));
  exprA->inf_cst = exprA->inf_cst + exprB->inf_cst + (maxA + maxB) * pr->ulp +
                   pr->min_denormal;
  exprA->sup_cst = exprA->sup_cst + exprB->sup_cst + (maxA + maxB) * pr->ulp +
                   pr->min_denormal;
  return;
}

// A = A + B
void add_expr(fppoly_internal_t *pr, expr_t *exprA, expr_t *exprB) {
  //
  size_t sizeB = exprB->size;
  if (sizeB == 0) {
    double maxA = fmax(fabs(exprA->inf_cst), fabs(exprA->sup_cst));
    double maxB = fmax(fabs(exprB->inf_cst), fabs(exprB->sup_cst));
    exprA->inf_cst = exprA->inf_cst + exprB->inf_cst + (maxA + maxB) * pr->ulp +
                     pr->min_denormal;
    exprA->sup_cst = exprA->sup_cst + exprB->sup_cst + (maxA + maxB) * pr->ulp +
                     pr->min_denormal;
    return;
  }
  size_t i;
  if (exprA->size == 0) {

    exprA->size = exprB->size;
    double maxA = fmax(fabs(exprA->inf_cst), fabs(exprA->sup_cst));
    double maxB = fmax(fabs(exprB->inf_cst), fabs(exprB->sup_cst));
    exprA->inf_cst +=
        exprB->inf_cst + (maxA + maxB) * pr->ulp + pr->min_denormal;
    exprA->sup_cst +=
        exprB->sup_cst + (maxA + maxB) * pr->ulp + pr->min_denormal;
    exprA->inf_coeff = (double *)malloc(sizeB * sizeof(double));
    exprA->sup_coeff = (double *)malloc(sizeB * sizeof(double));
    for (i = 0; i < sizeB; i++) {
      exprA->inf_coeff[i] = exprB->inf_coeff[i];
      exprA->sup_coeff[i] = exprB->sup_coeff[i];
    }
    exprA->type = exprB->type;
    if (exprA->type == SPARSE) {
      exprA->dim = (size_t *)malloc(sizeB * sizeof(size_t));
      for (i = 0; i < sizeB; i++) {
        exprA->dim[i] = exprB->dim[i];
      }
    }

    return;
  } else {
    size_t sizeA = exprA->size;
    assert(sizeA == sizeB);
    double maxA = fmax(fabs(exprA->inf_cst), fabs(exprA->sup_cst));
    double maxB = fmax(fabs(exprB->inf_cst), fabs(exprB->sup_cst));
    exprA->inf_cst +=
        exprB->inf_cst + (maxA + maxB) * pr->ulp + pr->min_denormal;
    exprA->sup_cst +=
        exprB->sup_cst + (maxA + maxB) * pr->ulp + pr->min_denormal;
    if (exprA->type == DENSE) {
      if (exprB->type == DENSE) {
        for (i = 0; i < sizeB; i++) {
          maxA = fmax(fabs(exprA->inf_coeff[i]), fabs(exprA->sup_coeff[i]));
          maxB = fmax(fabs(exprB->inf_coeff[i]), fabs(exprB->sup_coeff[i]));
          exprA->inf_coeff[i] = exprA->inf_coeff[i] + exprB->inf_coeff[i] +
                                (maxA + maxB) * pr->ulp;
          exprA->sup_coeff[i] = exprA->sup_coeff[i] + exprB->sup_coeff[i] +
                                (maxA + maxB) * pr->ulp;
        }
      } else {

        size_t k = 0;
        for (i = 0; i < sizeA; i++) {
          if (k < sizeB && exprB->dim[k] == i) {
            maxA = fmax(fabs(exprA->inf_coeff[i]), fabs(exprA->sup_coeff[i]));
            maxB = fmax(fabs(exprB->inf_coeff[k]), fabs(exprB->sup_coeff[k]));
            exprA->inf_coeff[i] = exprA->inf_coeff[i] + exprB->inf_coeff[k] +
                                  (maxA + maxB) * pr->ulp;
            exprA->sup_coeff[i] = exprA->sup_coeff[i] + exprB->sup_coeff[k] +
                                  (maxA + maxB) * pr->ulp;
            k++;
          }
        }
      }
    } else {
      size_t sizeB = exprB->size;
      size_t k;
      double *new_inf_coeff;
      double *new_sup_coeff;
      if (exprB->type == DENSE) {

        i = 0;
        new_inf_coeff = (double *)malloc(sizeB * sizeof(double));
        new_sup_coeff = (double *)malloc(sizeB * sizeof(double));
        for (k = 0; k < sizeB; k++) {
          if (i < sizeA && exprA->dim[i] == k) {
            maxA = fmax(fabs(exprA->inf_coeff[i]), fabs(exprA->sup_coeff[i]));
            maxB = fmax(fabs(exprB->inf_coeff[k]), fabs(exprB->sup_coeff[k]));
            new_inf_coeff[k] = exprA->inf_coeff[i] + exprB->inf_coeff[k] +
                               (maxA + maxB) * pr->ulp;
            new_sup_coeff[k] = exprA->sup_coeff[i] + exprB->sup_coeff[k] +
                               (maxA + maxB) * pr->ulp;
            i++;
          } else {
            new_inf_coeff[k] = exprB->inf_coeff[k];
            new_sup_coeff[k] = exprB->sup_coeff[k];
          }
        }
        exprA->type = DENSE;
        exprA->size = sizeB;
        free(exprA->dim);
        exprA->dim = NULL;
      } else {
        i = 0;
        k = 0;
        size_t l = 0;
        // printf("before sort\n");
        // expr_print(exprA);
        // fflush(stdout);
        // if(exprA->size>0){
        //	sort_sparse_expr(exprA);
        //}
        // printf("after sort\n");
        // expr_print(exprA);
        // fflush(stdout);

        new_inf_coeff = (double *)malloc((sizeA + sizeB) * sizeof(double));
        new_sup_coeff = (double *)malloc((sizeA + sizeB) * sizeof(double));
        size_t *new_dim = (size_t *)malloc((sizeA + sizeB) * sizeof(size_t));
        while (i < sizeA && k < sizeB) {
          if (exprA->dim[i] < exprB->dim[k]) {
            new_inf_coeff[l] = exprA->inf_coeff[i];
            new_sup_coeff[l] = exprA->sup_coeff[i];
            new_dim[l] = exprA->dim[i];
            i++;

          } else if (exprB->dim[k] < exprA->dim[i]) {
            new_inf_coeff[l] = exprB->inf_coeff[k];
            new_sup_coeff[l] = exprB->sup_coeff[k];
            new_dim[l] = exprB->dim[k];
            k++;
          } else {
            maxA = fmax(fabs(exprA->inf_coeff[i]), fabs(exprA->sup_coeff[i]));
            maxB = fmax(fabs(exprB->inf_coeff[k]), fabs(exprB->sup_coeff[k]));
            new_inf_coeff[l] = exprA->inf_coeff[i] + exprB->inf_coeff[k] +
                               (maxA + maxB) * pr->ulp;
            new_sup_coeff[l] = exprA->sup_coeff[i] + exprB->sup_coeff[k] +
                               (maxA + maxB) * pr->ulp;
            new_dim[l] = exprA->dim[i];
            i++;
            k++;
          }
          l++;
        }
        while (i < sizeA) {
          new_inf_coeff[l] = exprA->inf_coeff[i];
          new_sup_coeff[l] = exprA->sup_coeff[i];
          new_dim[l] = exprA->dim[i];
          i++;
          l++;
        }
        while (k < sizeB) {
          new_inf_coeff[l] = exprB->inf_coeff[k];
          new_sup_coeff[l] = exprB->sup_coeff[k];
          new_dim[l] = exprB->dim[k];
          k++;
          l++;
        }

        new_inf_coeff = (double *)realloc(new_inf_coeff, l * sizeof(double));
        new_sup_coeff = (double *)realloc(new_sup_coeff, l * sizeof(double));
        free(exprA->dim);
        exprA->dim = NULL;
        new_dim = (size_t *)realloc(new_dim, l * sizeof(size_t));
        exprA->dim = new_dim;
        exprA->size = l;
      }
      if (exprA->inf_coeff) {
        free(exprA->inf_coeff);
        exprA->inf_coeff = NULL;
      }
      if (exprA->sup_coeff) {
        free(exprA->sup_coeff);
        exprA->sup_coeff = NULL;
      }
      exprA->inf_coeff = new_inf_coeff;
      exprA->sup_coeff = new_sup_coeff;
    }
  }
}

expr_t *replace_input_poly_cons_in_lexpr(fppoly_internal_t *pr, expr_t *expr,
                                         fppoly_t *fp) {
  size_t dims = expr->size;
  size_t i, k;
  double tmp1, tmp2;
  expr_t *res;
  if (expr->type == DENSE) {
    k = 0;
  } else {
    k = expr->dim[0];
  }
  expr_t *mul_expr = NULL;

  if (expr->sup_coeff[0] < 0) {
    mul_expr = fp->input_uexpr[k];
  } else if (expr->inf_coeff[0] < 0) {
    mul_expr = fp->input_lexpr[k];
  }

  if (mul_expr != NULL) {
    if (mul_expr->size == 0) {
      res = multiply_cst_expr(pr, mul_expr, expr->inf_coeff[0],
                              expr->sup_coeff[0]);
    } else {
      res = multiply_expr(pr, mul_expr, expr->inf_coeff[0], expr->sup_coeff[0]);
    }
  }

  else {
    elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[0],
                              expr->sup_coeff[0], fp->input_inf[k],
                              fp->input_sup[k]);
    res = create_cst_expr(tmp1, -tmp1);
  }
  for (i = 1; i < dims; i++) {
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }

    expr_t *mul_expr = NULL;
    expr_t *sum_expr = NULL;
    if (expr->sup_coeff[i] < 0) {
      mul_expr = fp->input_uexpr[k];
    } else if (expr->inf_coeff[i] < 0) {
      mul_expr = fp->input_lexpr[k];
    }

    if (mul_expr != NULL) {
      if (mul_expr->size == 0) {
        sum_expr = multiply_cst_expr(pr, mul_expr, expr->inf_coeff[i],
                                     expr->sup_coeff[i]);
        add_cst_expr(pr, res, sum_expr);
      } else if (expr->inf_coeff[i] != 0 && expr->sup_coeff[i] != 0) {
        sum_expr =
            multiply_expr(pr, mul_expr, expr->inf_coeff[i], expr->sup_coeff[i]);
        add_expr(pr, res, sum_expr);
      }
      // free_expr(mul_expr);
      if (sum_expr != NULL) {
        free_expr(sum_expr);
      }
    } else {
      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i], fp->input_inf[k],
                                fp->input_sup[k]);
      res->inf_cst = res->inf_cst + tmp1;
      res->sup_cst = res->sup_cst - tmp1;
    }
  }

  res->inf_cst = res->inf_cst + expr->inf_cst;
  res->sup_cst = res->sup_cst + expr->sup_cst;
  return res;
}

expr_t *replace_input_poly_cons_in_uexpr(fppoly_internal_t *pr, expr_t *expr,
                                         fppoly_t *fp) {
  size_t dims = expr->size;
  size_t i, k;
  double tmp1, tmp2;
  expr_t *res;
  if (expr->type == DENSE) {
    k = 0;
  } else {
    k = expr->dim[0];
  }
  expr_t *mul_expr = NULL;

  if (expr->sup_coeff[0] < 0) {
    mul_expr = fp->input_lexpr[k];
  } else if (expr->inf_coeff[0] < 0) {
    mul_expr = fp->input_uexpr[k];
  }

  if (mul_expr != NULL) {
    if (mul_expr->size == 0) {
      res = multiply_cst_expr(pr, mul_expr, expr->inf_coeff[0],
                              expr->sup_coeff[0]);
    } else {
      res = multiply_expr(pr, mul_expr, expr->inf_coeff[0], expr->sup_coeff[0]);
    }
  } else {
    elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[0],
                              expr->sup_coeff[0], fp->input_inf[k],
                              fp->input_sup[k]);
    res = create_cst_expr(-tmp2, tmp2);
  }
  // printf("finish\n");
  // fflush(stdout);
  for (i = 1; i < dims; i++) {
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }
    expr_t *mul_expr = NULL;
    expr_t *sum_expr = NULL;
    if (expr->sup_coeff[i] < 0) {
      mul_expr = fp->input_lexpr[k];
    } else if (expr->inf_coeff[i] < 0) {
      mul_expr = fp->input_uexpr[k];
    }

    if (mul_expr != NULL) {
      if (mul_expr->size == 0) {
        sum_expr = multiply_cst_expr(pr, mul_expr, expr->inf_coeff[i],
                                     expr->sup_coeff[i]);
        add_cst_expr(pr, res, sum_expr);
      } else if (expr->inf_coeff[i] != 0 && expr->sup_coeff[i] != 0) {
        sum_expr =
            multiply_expr(pr, mul_expr, expr->inf_coeff[i], expr->sup_coeff[i]);
        add_expr(pr, res, sum_expr);
      }
      // free_expr(mul_expr);
      if (sum_expr != NULL) {
        free_expr(sum_expr);
      }
    } else {
      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i], fp->input_inf[k],
                                fp->input_sup[k]);
      res->inf_cst = res->inf_cst - tmp2;
      res->sup_cst = res->sup_cst + tmp2;
    }
  }
  res->inf_cst = res->inf_cst + expr->inf_cst;
  res->sup_cst = res->sup_cst + expr->sup_cst;
  return res;
}

double compute_lb_from_expr(fppoly_internal_t *pr, expr_t *expr, fppoly_t *fp,
                            int layerno) {
  size_t i, k;
  double tmp1, tmp2;
  // printf("start\n");
  // fflush(stdout);
  if ((fp->input_lexpr != NULL) && (fp->input_uexpr != NULL)) {
    expr = replace_input_poly_cons_in_lexpr(pr, expr, fp);
  }
  // expr_print(expr);
  // fflush(stdout);
  size_t dims = expr->size;
  double res_inf = expr->inf_cst;
  if (expr->inf_coeff == NULL || expr->sup_coeff == NULL) {
    return 0;
  }
  for (i = 0; i < dims; i++) {
    // if(expr->inf_coeff[i]<0){
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }
    if (layerno == -1) {
      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i], fp->input_inf[k],
                                fp->input_sup[k]);
    } else {

      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i],
                                fp->layers[layerno]->neurons[k]->lb,
                                fp->layers[layerno]->neurons[k]->ub);
    }
    // printf("tmp1: %g\n",tmp1);
    res_inf = res_inf + tmp1;
  }
  //	printf("inf: %g\n",-res_inf);
  //	fflush(stdout);
  if (fp->input_lexpr != NULL && fp->input_uexpr != NULL) {
    free_expr(expr);
  }
  // printf("finish\n");
  // fflush(stdout);
  return res_inf;
}

double compute_ub_from_expr(fppoly_internal_t *pr, expr_t *expr, fppoly_t *fp,
                            int layerno) {
  size_t i, k;
  double tmp1, tmp2;

  if ((fp->input_lexpr != NULL) && (fp->input_uexpr != NULL)) {
    expr = replace_input_poly_cons_in_uexpr(pr, expr, fp);
  }

  size_t dims = expr->size;
  double res_sup = expr->sup_cst;
  if (expr->inf_coeff == NULL || expr->sup_coeff == NULL) {
    return 0;
  }
  for (i = 0; i < dims; i++) {
    // if(expr->inf_coeff[i]<0){
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }
    if (layerno == -1) {
      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i], fp->input_inf[k],
                                fp->input_sup[k]);
    } else {
      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i],
                                fp->layers[layerno]->neurons[k]->lb,
                                fp->layers[layerno]->neurons[k]->ub);
    }
    res_sup = res_sup + tmp2;
  }
  // printf("sup: %g\n",res_sup);
  // fflush(stdout);
  if (fp->input_lexpr != NULL && fp->input_uexpr != NULL) {
    free_expr(expr);
  }
  return res_sup;
}

void ffn_handle_first_layer(elina_manager_t *man, elina_abstract0_t *abs,
                            double **weights, double *bias, size_t size,
                            size_t num_pixels, size_t *predecessors,
                            activation_type_t activation, bool alloc) {
  // printf("start \n");
  // fflush(stdout);
  fppoly_t *res = fppoly_of_abstract0(abs);
  // printf("coming here\n");
  // fflush(stdout);

  if (alloc) {
    fppoly_alloc_first_layer(res, size, FFN, activation);
  }
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  size_t i, j;

  // for(i=0; i < num_pixels; i++){
  //	elina_interval_print(itv[i]);
  //	printf("\n");
  //}
  // fflush(stdout);
  neuron_t **neurons = res->layers[0]->neurons;
  res->layers[0]->predecessors = predecessors;

  for (i = 0; i < size; i++) {
    neuron_t *neuron = neurons[i];
    double *weight_i = weights[i];
    double bias_i = bias[i];
    neuron->expr = create_dense_expr(weight_i, bias_i, num_pixels);
    neuron->lb = compute_lb_from_expr(pr, neuron->expr, res, -1);
    neuron->ub = compute_ub_from_expr(pr, neuron->expr, res, -1);
  }

  // printf("return here\n");
  // fppoly_fprint(stdout,man,res,NULL);
  // fflush(stdout);
  return;
}

void ffn_handle_first_relu_layer(elina_manager_t *man, elina_abstract0_t *abs,
                                 double **weights, double *bias, size_t size,
                                 size_t num_pixels, size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, RELU, true);
}

void ffn_handle_first_sigmoid_layer(elina_manager_t *man,
                                    elina_abstract0_t *abs, double **weights,
                                    double *bias, size_t size,
                                    size_t num_pixels, size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, SIGMOID, true);
}

void ffn_handle_first_tanh_layer(elina_manager_t *man, elina_abstract0_t *abs,
                                 double **weights, double *bias, size_t size,
                                 size_t num_pixels, size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, TANH, true);
}

void ffn_handle_first_parabola_layer(elina_manager_t *man,
                                     elina_abstract0_t *abs, double **weights,
                                     double *bias, size_t size,
                                     size_t num_pixels, size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, PARABOLA, true);
}

void ffn_handle_first_log_layer(elina_manager_t *man, elina_abstract0_t *abs,
                                double **weights, double *bias, size_t size,
                                size_t num_pixels, size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, LOG, true);
}

void ffn_handle_first_relu_layer_no_alloc(elina_manager_t *man,
                                          elina_abstract0_t *abs,
                                          double **weights, double *bias,
                                          size_t size, size_t num_pixels,
                                          size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, RELU, false);
}

void ffn_handle_first_sigmoid_layer_no_alloc(elina_manager_t *man,
                                             elina_abstract0_t *abs,
                                             double **weights, double *bias,
                                             size_t size, size_t num_pixels,
                                             size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, SIGMOID, false);
}

void ffn_handle_first_tanh_layer_no_alloc(elina_manager_t *man,
                                          elina_abstract0_t *abs,
                                          double **weights, double *bias,
                                          size_t size, size_t num_pixels,
                                          size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, TANH, false);
}

void ffn_handle_first_parabola_layer_no_alloc(elina_manager_t *man,
                                              elina_abstract0_t *abs,
                                              double **weights, double *bias,
                                              size_t size, size_t num_pixels,
                                              size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, PARABOLA, false);
}

void ffn_handle_first_log_layer_no_alloc(elina_manager_t *man,
                                         elina_abstract0_t *abs,
                                         double **weights, double *bias,
                                         size_t size, size_t num_pixels,
                                         size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, LOG, false);
}

expr_t *lexpr_replace_parabola_bounds(fppoly_internal_t *pr, expr_t *expr,
                                      neuron_t **neurons) {
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
    res->inf_coeff[i] = 0.0;
    res->sup_coeff[i] = 0.0;
    if ((expr->sup_coeff[i] == 0) && (expr->inf_coeff[i] == 0)) {
      continue;
    }
    double u_plus_l_sup = (ub - lb);
    double u_plus_l_inf = -(ub - lb);

    //= (u_plus_l)*(u_plus_l)
    if (expr->sup_coeff[i] < 0) {

      double lu_sup;
      double lu_inf;
      elina_double_interval_mul_cst_coeff(pr, &lu_inf, &lu_sup, lb, -lb, -ub,
                                          ub);
      // res->coeff[i] = lambda*expr->coeff[i];
      // res->cst = res->cst + expr->coeff[i]*mu;
      elina_double_interval_mul_expr_coeff(
          pr, &res->inf_coeff[i], &res->sup_coeff[i], u_plus_l_inf,
          u_plus_l_sup, expr->inf_coeff[i], expr->sup_coeff[i]);
      double tmp1, tmp2;
      elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, lu_inf, lu_sup,
                                          expr->inf_coeff[i],
                                          expr->sup_coeff[i]);
      res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
      res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
    } else if (expr->inf_coeff[i] < 0) {
      res->inf_coeff[i] = 0;
      res->sup_coeff[i] = 0;
      // double u_plus_l_sq_inf, u_plus_l_sq_sup;
      // u_plus_l_sq_inf = u_plus_l_inf/2;
      // u_plus_l_sq_sup = u_plus_l_sup/2;
      // elina_double_interval_mul_cst_coeff(pr,&u_plus_l_sq_inf,&u_plus_l_sq_sup,u_plus_l_inf/2,u_plus_l_sup/2,u_plus_l_inf/2,u_plus_l_sup/2);
      // elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],u_plus_l_inf,u_plus_l_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
      // double tmp1, tmp2;
      // elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,-u_plus_l_sq_inf,-u_plus_l_sq_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
      double tmp1, tmp2;
      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i], -lb * lb, lb * lb);
      res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
      res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
    } else {
      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;

      double tmp1, tmp2;
      if (lb * lb < ub * ub) {
        elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                  expr->sup_coeff[i], -lb * lb, ub * ub);
      } else {
        elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                  expr->sup_coeff[i], -ub * ub, lb * lb);
      }
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

expr_t *uexpr_replace_parabola_bounds(fppoly_internal_t *pr, expr_t *expr,
                                      neuron_t **neurons) {
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
    res->inf_coeff[i] = 0.0;
    res->sup_coeff[i] = 0.0;
    if ((expr->sup_coeff[i] == 0) && (expr->inf_coeff[i] == 0)) {
      continue;
    }
    double u_plus_l_sup = (ub - lb);
    double u_plus_l_inf = -(ub - lb);

    //= (u_plus_l)*(u_plus_l)
    if (expr->inf_coeff[i] < 0) {

      double lu_sup;
      double lu_inf;
      elina_double_interval_mul_cst_coeff(pr, &lu_inf, &lu_sup, lb, -lb, -ub,
                                          ub);
      // res->coeff[i] = lambda*expr->coeff[i];
      // res->cst = res->cst + expr->coeff[i]*mu;
      elina_double_interval_mul_expr_coeff(
          pr, &res->inf_coeff[i], &res->sup_coeff[i], u_plus_l_inf,
          u_plus_l_sup, expr->inf_coeff[i], expr->sup_coeff[i]);
      double tmp1, tmp2;
      elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, lu_inf, lu_sup,
                                          expr->inf_coeff[i],
                                          expr->sup_coeff[i]);
      res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
      res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
    } else if (expr->sup_coeff[i] < 0) {
      res->inf_coeff[i] = 0;
      res->sup_coeff[i] = 0;
      // double u_plus_l_sq_inf, u_plus_l_sq_sup;
      // u_plus_l_sq_inf = u_plus_l_inf/2;
      // u_plus_l_sq_sup = u_plus_l_sup/2;
      // elina_double_interval_mul_cst_coeff(pr,&u_plus_l_sq_inf,&u_plus_l_sq_sup,u_plus_l_inf/2,u_plus_l_sup/2,u_plus_l_inf/2,u_plus_l_sup/2);
      // elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],u_plus_l_inf,u_plus_l_sup,expr->inf_coeff[i],expr->sup_coeff[i]);

      double tmp1, tmp2;
      elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, -lb * lb, lb * lb,
                                          expr->inf_coeff[i],
                                          expr->sup_coeff[i]);
      // double tmp =lb*lb;
      res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
      res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
    } else {
      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;
      double tmp1, tmp2;
      if (lb * lb < ub * ub) {
        elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                  expr->sup_coeff[i], -lb * lb, ub * ub);
      } else {
        elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                  expr->sup_coeff[i], -ub * ub, lb * lb);
      }
      res->inf_cst = res->inf_cst + tmp2;
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

expr_t *lexpr_replace_log_bounds(fppoly_internal_t *pr, expr_t *expr,
                                 neuron_t **neurons) {
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
    res->inf_coeff[i] = 0.0;
    res->sup_coeff[i] = 0.0;
    if ((expr->sup_coeff[i] == 0) && (expr->inf_coeff[i] == 0)) {
      continue;
    }
    double u_plus_l_inf = -(lb + ub);
    double u_plus_l_sup = -u_plus_l_inf;

    //= (u_plus_l)*(u_plus_l)
    if (expr->sup_coeff[i] < 0) {
      double one_inf = 1;
      double one_sup = -1;
      double u_plus_l_sup = ub - lb;
      double u_plus_l_inf = -(ub - lb);

      double lambda_sup = -2 / u_plus_l_inf;
      double lambda_inf = -2 / u_plus_l_sup;
      elina_double_interval_mul_expr_coeff(
          pr, &res->inf_coeff[i], &res->sup_coeff[i], lambda_inf, lambda_sup,
          expr->inf_coeff[i], expr->sup_coeff[i]);

      fesetround(FE_DOWNWARD);
      double mu_inf = -log(-u_plus_l_inf / 2);
      fesetround(FE_UPWARD);
      double mu_sup = log(u_plus_l_sup / 2);
      elina_double_interval_add_cst_coeff(pr, &mu_inf, &mu_sup, one_inf,
                                          one_sup, mu_inf, mu_sup);
      double tmp1, tmp2;
      elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, mu_inf, mu_sup,
                                          expr->inf_coeff[i],
                                          expr->sup_coeff[i]);
      res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
      res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
    } else if (expr->inf_coeff[i] < 0) {
      double u_minus_l_sup = ub + lb;
      double u_minus_l_inf = -(ub + lb);

      if (u_minus_l_sup < 1e-9) {
        double tmp1, tmp2;
        double log_lb = -log(-lb), log_ub = log(-lb);
        if (log_lb > 0) {
          log_lb = -log_lb;
          log_ub = -log_ub;
        }
        elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                  expr->sup_coeff[i], -log(-lb), log(-lb));
        res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
        res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
      } else {
        double inv_u_by_l_sup = -1 / u_minus_l_inf;
        double inv_u_by_l_inf = -1 / u_minus_l_sup;

        double u_by_l_sup = -ub / lb;
        double u_by_l_inf = ub / lb;

        fesetround(FE_DOWNWARD);
        double log_u_by_l_inf = -log(-u_by_l_inf);
        double log_l_inf = -log(-lb);

        fesetround(FE_UPWARD);
        double log_u_by_l_sup = log(u_by_l_sup);
        double log_l_sup = log(-lb);

        double lambda_inf, lambda_sup;
        elina_double_interval_mul_cst_coeff(pr, &lambda_inf, &lambda_sup,
                                            log_u_by_l_inf, log_u_by_l_sup,
                                            inv_u_by_l_inf, inv_u_by_l_sup);
        elina_double_interval_mul_expr_coeff(
            pr, &res->inf_coeff[i], &res->sup_coeff[i], lambda_inf, lambda_sup,
            expr->inf_coeff[i], expr->sup_coeff[i]);

        double mu_inf, mu_sup;
        elina_double_interval_mul_cst_coeff(pr, &mu_inf, &mu_sup, -lb, lb,
                                            lambda_inf, lambda_sup);
        elina_double_interval_add_cst_coeff(pr, &mu_inf, &mu_sup, log_l_inf,
                                            log_l_sup, mu_inf, mu_sup);

        double tmp1, tmp2;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, mu_inf, mu_sup,
                                            expr->inf_coeff[i],
                                            expr->sup_coeff[i]);
        res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
        res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
      }
    } else {

      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;
      double log_lb = -log(-lb);
      double log_ub = log(ub);
      double tmp1, tmp2;
      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i], log_lb, log_ub);
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

expr_t *uexpr_replace_log_bounds(fppoly_internal_t *pr, expr_t *expr,
                                 neuron_t **neurons) {
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
    res->inf_coeff[i] = 0.0;
    res->sup_coeff[i] = 0.0;
    if ((expr->sup_coeff[i] == 0) && (expr->inf_coeff[i] == 0)) {
      continue;
    }
    double u_plus_l_inf = -(lb + ub);
    double u_plus_l_sup = -u_plus_l_inf;

    //= (u_plus_l)*(u_plus_l)
    if (expr->inf_coeff[i] < 0) {
      double one_inf = 1;
      double one_sup = -1;
      double u_plus_l_sup = ub - lb;
      double u_plus_l_inf = -(ub - lb);

      double lambda_sup = -2 / u_plus_l_inf;
      double lambda_inf = -2 / u_plus_l_sup;
      elina_double_interval_mul_expr_coeff(
          pr, &res->inf_coeff[i], &res->sup_coeff[i], lambda_inf, lambda_sup,
          expr->inf_coeff[i], expr->sup_coeff[i]);

      fesetround(FE_DOWNWARD);
      double mu_inf = -log(-u_plus_l_inf / 2);
      fesetround(FE_UPWARD);
      double mu_sup = log(u_plus_l_sup / 2);
      elina_double_interval_add_cst_coeff(pr, &mu_inf, &mu_sup, one_inf,
                                          one_sup, mu_inf, mu_sup);
      double tmp1, tmp2;
      elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, mu_inf, mu_sup,
                                          expr->inf_coeff[i],
                                          expr->sup_coeff[i]);
      res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
      res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
    } else if (expr->sup_coeff[i] < 0) {
      double u_minus_l_sup = ub + lb;
      double u_minus_l_inf = -(ub + lb);

      if (u_minus_l_sup < 1e-9) {
        double tmp1, tmp2;
        double log_lb = -log(-lb), log_ub = log(-lb);
        if (log_lb > 0) {
          log_lb = -log_lb;
          log_ub = -log_ub;
        }
        elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                  expr->sup_coeff[i], -log(-lb), log(-lb));
        res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
        res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
      } else {
        double inv_u_by_l_sup = -1 / u_minus_l_inf;
        double inv_u_by_l_inf = -1 / u_minus_l_sup;

        double u_by_l_sup = -ub / lb;
        double u_by_l_inf = ub / lb;

        fesetround(FE_DOWNWARD);
        double log_u_by_l_inf = -log(-u_by_l_inf);
        double log_l_inf = -log(-lb);

        fesetround(FE_UPWARD);
        double log_u_by_l_sup = log(u_by_l_sup);
        double log_l_sup = log(-lb);

        double lambda_inf, lambda_sup;
        elina_double_interval_mul_cst_coeff(pr, &lambda_inf, &lambda_sup,
                                            log_u_by_l_inf, log_u_by_l_sup,
                                            inv_u_by_l_inf, inv_u_by_l_sup);
        elina_double_interval_mul_expr_coeff(
            pr, &res->inf_coeff[i], &res->sup_coeff[i], lambda_inf, lambda_sup,
            expr->inf_coeff[i], expr->sup_coeff[i]);

        double mu_inf, mu_sup;
        elina_double_interval_mul_cst_coeff(pr, &mu_inf, &mu_sup, -lb, lb,
                                            lambda_inf, lambda_sup);
        elina_double_interval_add_cst_coeff(pr, &mu_inf, &mu_sup, log_l_inf,
                                            log_l_sup, mu_inf, mu_sup);

        double tmp1, tmp2;
        elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, mu_inf, mu_sup,
                                            expr->inf_coeff[i],
                                            expr->sup_coeff[i]);
        res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
        res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
      }
    } else {

      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;
      double log_lb = -log(-lb);
      double log_ub = log(ub);
      double tmp1, tmp2;
      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i], log_lb, log_ub);
      res->inf_cst = res->inf_cst + tmp1;
      res->sup_cst = res->sup_cst + tmp2;
      // res->inf_cst = -INFINITY;
      // res->sup_cst = INFINITY;
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

void compute_chord_slope(double *slope_inf, double *slope_sup, double f_sup_l,
                         double f_sup_u, double f_inf_l, double f_inf_u,
                         double inf_l, double inf_u, double sup_l,
                         double sup_u) {
  double num_l = f_sup_l + f_inf_u;
  double num_u = f_sup_u + f_inf_l;

  double den_l = sup_l + inf_u;
  double den_u = sup_u + inf_l;

  elina_double_interval_div(slope_inf, slope_sup, num_l, num_u, den_l, den_u);
}

void compute_derivative(double *slope_inf, double *slope_sup, double s_curve_l,
                        double s_curve_u, double sq_l, double sq_u,
                        bool is_sigmoid) {
  double sq_den_sup_l, sq_den_sup_u;
  elina_double_interval_mul(&sq_den_sup_l, &sq_den_sup_u, sq_l, sq_u, sq_l,
                            sq_u);
  if (is_sigmoid) {
    elina_double_interval_div(slope_inf, slope_sup, s_curve_l, s_curve_u,
                              sq_den_sup_l, sq_den_sup_u);
  } else {
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

expr_t *lexpr_replace_maxpool_or_lstm_bounds(fppoly_internal_t *pr,
                                             expr_t *expr, neuron_t **neurons) {
  // printf("begin\n");
  // fflush(stdout);
  size_t num_neurons = expr->size;
  size_t i, k;
  expr_t *res;
  if (expr->type == DENSE) {
    k = 0;
  } else {
    k = expr->dim[0];
  }
  neuron_t *neuron_k = neurons[k];
  if (expr->sup_coeff[0] < 0) {
    // expr_print(neuron_k->uexpr);
    if (neuron_k->uexpr == NULL) {
      res = (expr_t *)malloc(sizeof(expr_t));
      res->inf_coeff = res->sup_coeff = NULL;
      res->dim = NULL;
      res->size = 0;
      res->type = SPARSE;
      elina_double_interval_mul_cst_coeff(
          pr, &res->inf_cst, &res->sup_cst, neuron_k->lb, neuron_k->ub,
          expr->inf_coeff[0], expr->sup_coeff[0]);
    } else {
      res = multiply_expr(pr, neuron_k->uexpr, expr->inf_coeff[0],
                          expr->sup_coeff[0]);
    }
    // printf("multiply end %zu \n",k);
    // expr_print(res);
    // fflush(stdout);
  } else if (expr->inf_coeff[0] < 0) {
    // expr_print(neuron_k->lexpr);
    if (neuron_k->lexpr == NULL) {
      res = (expr_t *)malloc(sizeof(expr_t));
      res->inf_coeff = res->sup_coeff = NULL;
      res->dim = NULL;
      res->size = 0;
      res->type = SPARSE;
      elina_double_interval_mul_cst_coeff(
          pr, &res->inf_cst, &res->sup_cst, neuron_k->lb, neuron_k->ub,
          expr->inf_coeff[0], expr->sup_coeff[0]);
    } else {
      res = multiply_expr(pr, neuron_k->lexpr, expr->inf_coeff[0],
                          expr->sup_coeff[0]);
    }
    // printf("multiply end %zu \n",k);
    // expr_print(res);
    // fflush(stdout);
  } else {
    // printf("WTF1\n");
    // fflush(stdout);
    double tmp1, tmp2;
    elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, neuron_k->lb,
                                        neuron_k->ub, expr->inf_coeff[0],
                                        expr->sup_coeff[0]);
    double coeff[1];
    size_t dim[1];
    coeff[0] = 0;
    dim[0] = 0;
    res = create_sparse_expr(coeff, -tmp1, dim, 1);
  }
  // printf("middle\n");
  // fflush(stdout);
  expr_t *mul_expr;
  for (i = 1; i < num_neurons; i++) {
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }
    neuron_t *neuron_k = neurons[k];
    if (expr->sup_coeff[i] < 0) {
      // expr_print(neuron_k->uexpr);
      // printf("add start %zu %zu\n",k,i);

      // expr_print(res);

      if (neuron_k->uexpr == NULL) {
        mul_expr = (expr_t *)malloc(sizeof(expr_t));
        mul_expr->inf_coeff = mul_expr->sup_coeff = NULL;
        mul_expr->dim = NULL;
        mul_expr->size = 0;
        mul_expr->type = SPARSE;
        // printf("lb: %g %g\n");
        elina_double_interval_mul_cst_coeff(
            pr, &mul_expr->inf_cst, &mul_expr->sup_cst, neuron_k->lb,
            neuron_k->ub, expr->inf_coeff[i], expr->sup_coeff[i]);
        res->inf_cst += mul_expr->inf_cst;
        res->sup_cst += mul_expr->sup_cst;
      } else {
        mul_expr = multiply_expr(pr, neuron_k->uexpr, expr->inf_coeff[i],
                                 expr->sup_coeff[i]);

        add_expr(pr, res, mul_expr);
      }
      // expr_print(mul_expr);
      // fflush(stdout);
      // printf("add finish\n");
      // expr_print(res);
      // fflush(stdout);
      free_expr(mul_expr);
    } else if (expr->inf_coeff[i] < 0) {
      // expr_print(neuron_k->lexpr);
      // printf("add start %zu %zu\n",k,i);

      // expr_print(res);

      if (neuron_k->lexpr == NULL) {
        mul_expr = (expr_t *)malloc(sizeof(expr_t));
        mul_expr->inf_coeff = mul_expr->sup_coeff = NULL;
        mul_expr->dim = NULL;
        mul_expr->size = 0;
        mul_expr->type = SPARSE;
        elina_double_interval_mul_cst_coeff(
            pr, &mul_expr->inf_cst, &mul_expr->sup_cst, neuron_k->lb,
            neuron_k->ub, expr->inf_coeff[i], expr->sup_coeff[i]);
        res->inf_cst += mul_expr->inf_cst;
        res->sup_cst += mul_expr->sup_cst;
      } else {
        mul_expr = multiply_expr(pr, neuron_k->lexpr, expr->inf_coeff[i],
                                 expr->sup_coeff[i]);
        // printf("add start1 %zu %zu\n",k,i);
        // expr_print(res);
        // expr_print(mul_expr);
        // fflush(stdout);
        add_expr(pr, res, mul_expr);
      }
      // expr_print(mul_expr);
      //	fflush(stdout);
      // printf("add finish1\n");
      // expr_print(res);
      // fflush(stdout);
      free_expr(mul_expr);
    } else {
      // printf("WTF2\n");
      // fflush(stdout);
      double tmp1, tmp2;
      elina_double_interval_mul_expr_coeff(pr, &tmp1, &tmp2, neuron_k->lb,
                                           neuron_k->ub, expr->inf_coeff[i],
                                           expr->sup_coeff[i]);
      res->inf_cst = res->inf_cst + tmp1;
      res->sup_cst = res->sup_cst - tmp1;
    }
  }
  // printf("finish\n");
  // fflush(stdout);
  res->inf_cst = res->inf_cst + expr->inf_cst;
  res->sup_cst = res->sup_cst + expr->sup_cst;
  return res;
}

expr_t *uexpr_replace_maxpool_or_lstm_bounds(fppoly_internal_t *pr,
                                             expr_t *expr, neuron_t **neurons) {
  size_t num_neurons = expr->size;
  size_t i, k;
  expr_t *res;
  if (expr->type == DENSE) {
    k = 0;
  } else {
    k = expr->dim[0];
  }
  neuron_t *neuron_k = neurons[k];
  if (expr->sup_coeff[0] < 0) {
    res = multiply_expr(pr, neuron_k->lexpr, expr->inf_coeff[0],
                        expr->sup_coeff[0]);
  } else if (expr->inf_coeff[0] < 0) {
    res = multiply_expr(pr, neuron_k->uexpr, expr->inf_coeff[0],
                        expr->sup_coeff[0]);
  } else {

    double tmp1, tmp2;
    elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, neuron_k->lb,
                                        neuron_k->ub, expr->inf_coeff[0],
                                        expr->sup_coeff[0]);
    double coeff[1];
    size_t dim[1];
    coeff[0] = 0;
    dim[0] = 0;
    res = create_sparse_expr(coeff, tmp2, dim, 1);
  }

  for (i = 1; i < num_neurons; i++) {
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }
    neuron_t *neuron_k = neurons[k];
    if (expr->sup_coeff[i] < 0) {
      expr_t *mul_expr = multiply_expr(pr, neuron_k->lexpr, expr->inf_coeff[i],
                                       expr->sup_coeff[i]);
      add_expr(pr, res, mul_expr);
      free_expr(mul_expr);
    } else if (expr->inf_coeff[i] < 0) {
      expr_t *mul_expr = multiply_expr(pr, neuron_k->uexpr, expr->inf_coeff[i],
                                       expr->sup_coeff[i]);
      add_expr(pr, res, mul_expr);
      free_expr(mul_expr);
    } else {

      double tmp1, tmp2;
      elina_double_interval_mul_cst_coeff(pr, &tmp1, &tmp2, neuron_k->lb,
                                          neuron_k->ub, expr->inf_coeff[i],
                                          expr->sup_coeff[i]);
      res->inf_cst = res->inf_cst - tmp2;
      res->sup_cst = res->sup_cst + tmp2;
    }
  }
  res->inf_cst = res->inf_cst + expr->inf_cst;
  res->sup_cst = res->sup_cst + expr->sup_cst;
  return res;
}

expr_t *expr_from_previous_layer(fppoly_internal_t *pr, expr_t *expr,
                                 layer_t *prev_layer) {
  if (expr->size == 0) {
    return copy_cst_expr(expr);
  }
  if (expr->inf_coeff == NULL || expr->sup_coeff == NULL) {
    return alloc_expr();
  }
  // printf("coming here %zu\n",expr->size);
  //	fflush(stdout);
  neuron_t **prev_neurons = prev_layer->neurons;
  size_t out_num_neurons = prev_layer->dims;
  size_t in_num_neurons = expr->size;
  size_t i, k;
  expr_t *res;

  if (expr->type == DENSE) {
    k = 0;
  } else {
    k = expr->dim[0];
  }

  // printf("start2 %p %lu %p\n",expr->dim,k,prev_neurons[k]->expr);
  // fflush(stdout);
  // if(prev_layer->type==MAXPOOL){
  // printf("i: %zu k: %zu\n",i,k);
  // expr_print(prev_neurons[k]->expr);
  // fflush(stdout);
  //}
  if (prev_neurons[k]->expr->size == 0) {

    res = multiply_cst_expr(pr, prev_neurons[k]->expr, expr->inf_coeff[0],
                            expr->sup_coeff[0]);

  } else {

    res = multiply_expr(pr, prev_neurons[k]->expr, expr->inf_coeff[0],
                        expr->sup_coeff[0]);
  }
  // printf("debug\n");
  // fflush(stdout);
  for (i = 1; i < in_num_neurons; i++) {
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }
    expr_t *mul_expr;
    // printf("i: %zu k: %zu\n",i,k);
    // expr_print(res);
    // expr_print(prev_neurons[k]->expr);
    // fflush(stdout);
    // if(prev_layer->type==MAXPOOL){

    //}
    if (prev_neurons[k]->expr->size == 0) {
      mul_expr = multiply_cst_expr(pr, prev_neurons[k]->expr,
                                   expr->inf_coeff[i], expr->sup_coeff[i]);
      add_cst_expr(pr, res, mul_expr);
      free_expr(mul_expr);
    } else if (expr->inf_coeff[i] != 0 || expr->sup_coeff[i] != 0) {
      mul_expr = multiply_expr(pr, prev_neurons[k]->expr, expr->inf_coeff[i],
                               expr->sup_coeff[i]);
      // printf("start\n");
      // fflush(stdout);
      add_expr(pr, res, mul_expr);
      free_expr(mul_expr);
      // printf("finish\n");
      // fflush(stdout);
    }
    // if(prev_layer->type==MAXPOOL){
    // printf("i: %zu k: %zu\n",i,k);
    // expr_print(mul_expr);
    // printf("output\n");
    // expr_print(res);
    // fflush(stdout);
    //}
  }
  res->inf_cst = res->inf_cst + expr->inf_cst;
  res->sup_cst = res->sup_cst + expr->sup_cst;
  // printf("finish\n");
  // fflush(stdout);
  return res;
}

expr_t *lexpr_unroll_lstm_layer(fppoly_internal_t *pr, expr_t *expr,
                                neuron_t **neurons) {
  return NULL;
}

void update_state_using_predecessor_layer(fppoly_internal_t *pr, fppoly_t *fp,
                                          expr_t **lexpr_ptr,
                                          expr_t **uexpr_ptr, size_t k,
                                          bool use_area_heuristic) {
  expr_t *lexpr = *lexpr_ptr;
  expr_t *uexpr = *uexpr_ptr;
  expr_t *tmp_l = lexpr;
  expr_t *tmp_u = uexpr;
  neuron_t **aux_neurons = fp->layers[k]->neurons;
  if (fp->layers[k]->type == FFN || fp->layers[k]->type == CONV) {
    if (fp->layers[k]->activation == RELU) {

      lexpr =
          lexpr_replace_relu_bounds(pr, lexpr, aux_neurons, use_area_heuristic);

      uexpr =
          uexpr_replace_relu_bounds(pr, uexpr, aux_neurons, use_area_heuristic);
    } else if (fp->layers[k]->activation == SIGMOID) {

      lexpr = lexpr_replace_sigmoid_bounds(pr, lexpr, aux_neurons);

      uexpr = uexpr_replace_sigmoid_bounds(pr, uexpr, aux_neurons);
    } else if (fp->layers[k]->activation == TANH) {
      lexpr = lexpr_replace_tanh_bounds(pr, lexpr, aux_neurons);
      uexpr = uexpr_replace_tanh_bounds(pr, uexpr, aux_neurons);
    } else if (fp->layers[k]->activation == PARABOLA) {
      lexpr = lexpr_replace_parabola_bounds(pr, lexpr, aux_neurons);
      uexpr = uexpr_replace_parabola_bounds(pr, uexpr, aux_neurons);
    } else if (fp->layers[k]->activation == LOG) {
      lexpr = lexpr_replace_log_bounds(pr, lexpr, aux_neurons);
      uexpr = uexpr_replace_log_bounds(pr, uexpr, aux_neurons);
    }

    if (fp->layers[k]->activation != NONE) {
      free_expr(tmp_l);

      free_expr(tmp_u);
    }

    tmp_l = lexpr;
    tmp_u = uexpr;

    *lexpr_ptr = expr_from_previous_layer(pr, lexpr, fp->layers[k]);

    *uexpr_ptr = expr_from_previous_layer(pr, uexpr, fp->layers[k]);

    free_expr(tmp_l);
    free_expr(tmp_u);
  } else if (fp->layers[k]->type == MAXPOOL || fp->layers[k]->type == LSTM) {
    expr_t *tmp_l = lexpr;
    expr_t *tmp_u = uexpr;
    *lexpr_ptr = lexpr_replace_maxpool_or_lstm_bounds(pr, lexpr, aux_neurons);

    *uexpr_ptr = uexpr_replace_maxpool_or_lstm_bounds(pr, uexpr, aux_neurons);
    free_expr(tmp_l);

    free_expr(tmp_u);
  }
  //}
  else {
    expr_t *tmp_l = lexpr;
    expr_t *tmp_u = uexpr;
    *lexpr_ptr = lexpr_unroll_lstm_layer(pr, lexpr, aux_neurons);

    // uexpr = uexpr_unroll_lstm_layer(pr, uexpr, aux_neurons);

    free_expr(tmp_l);
    free_expr(tmp_u);
  }
}

void *update_state_using_previous_layers(void *args) {
  nn_thread_t *data = (nn_thread_t *)args;
  elina_manager_t *man = data->man;
  fppoly_t *fp = data->fp;
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  size_t layerno = data->layerno;
  size_t idx_start = data->start;
  size_t idx_end = data->end;
  bool use_area_heuristic = data->use_area_heuristic;
  size_t i;
  int k;

  neuron_t **out_neurons = fp->layers[layerno]->neurons;
  size_t num_out_neurons = fp->layers[layerno]->dims;
  // printf("idx: %zu %zu\n",idx_start,idx_end);
  // fflush(stdout);
  for (i = idx_start; i < idx_end; i++) {
    bool already_computed = false;
    expr_t *lexpr = copy_expr(out_neurons[i]->expr);
    expr_t *uexpr = copy_expr(out_neurons[i]->expr);
    if (fp->layers[layerno]->type == RESIDUAL) {
      expr_t *lexpr_copy = copy_expr(lexpr);
      lexpr_copy->inf_cst = 0;
      lexpr_copy->sup_cst = 0;
      expr_t *uexpr_copy = copy_expr(uexpr);
      uexpr_copy->inf_cst = 0;
      uexpr_copy->sup_cst = 0;
      size_t predecessor1 = fp->layers[layerno]->predecessors[0] - 1;
      size_t predecessor2 = fp->layers[layerno]->predecessors[1] - 1;
      char *predecessor_map = (char *)calloc(layerno, sizeof(char));
      // Assume no nested residual layers
      int iter = predecessor1;
      while (iter >= 0) {
        predecessor_map[iter] = 1;
        iter = fp->layers[iter]->predecessors[0] - 1;
      }
      iter = predecessor2;
      int common_predecessor = 0;
      while (iter >= 0) {
        if (predecessor_map[iter] == 1) {
          common_predecessor = iter;
          break;
        }
        iter = fp->layers[iter]->predecessors[0] - 1;
      }

      iter = predecessor1;
      while (iter != common_predecessor) {
        update_state_using_predecessor_layer(pr, fp, &lexpr, &uexpr, iter,
                                             use_area_heuristic);

        iter = fp->layers[iter]->predecessors[0] - 1;
      }
      iter = predecessor2;
      while (iter != common_predecessor) {

        update_state_using_predecessor_layer(pr, fp, &lexpr_copy, &uexpr_copy,
                                             iter, use_area_heuristic);

        iter = fp->layers[iter]->predecessors[0] - 1;
      }
      free(predecessor_map);
      add_expr(pr, lexpr, lexpr_copy);
      add_expr(pr, uexpr, uexpr_copy);
      free_expr(lexpr_copy);
      free_expr(uexpr_copy);
      // Assume at least one non-residual layer between two residual layers

      k = common_predecessor;

    } else {

      k = fp->layers[layerno]->predecessors[0] - 1;
    }
    // printf("k: %d layerno: %zu\n",k,layerno);
    while (k >= 0) {
      // for(k=fp->layers[layerno]->predecessors[0]-1; k >=0; k =
      // fp->layers[k]->predecessors[0]-1){
      neuron_t **aux_neurons = fp->layers[k]->neurons;

      if (fp->layers[k]->type == RESIDUAL) {
        if (fp->layers[k]->activation == RELU) {
          neuron_t **aux_neurons = fp->layers[k]->neurons;
          expr_t *tmp_l = lexpr;
          expr_t *tmp_u = uexpr;
          lexpr = lexpr_replace_relu_bounds(pr, lexpr, aux_neurons,
                                            use_area_heuristic);
          uexpr = uexpr_replace_relu_bounds(pr, uexpr, aux_neurons,
                                            use_area_heuristic);
          free_expr(tmp_l);
          free_expr(tmp_u);
        }
        expr_t *lexpr_copy = copy_expr(lexpr);
        lexpr_copy->inf_cst = 0;
        lexpr_copy->sup_cst = 0;
        expr_t *uexpr_copy = copy_expr(uexpr);
        uexpr_copy->inf_cst = 0;
        uexpr_copy->sup_cst = 0;
        size_t predecessor1 = fp->layers[k]->predecessors[0] - 1;
        size_t predecessor2 = fp->layers[k]->predecessors[1] - 1;

        char *predecessor_map = (char *)calloc(k, sizeof(char));
        // Assume no nested residual layers
        int iter = predecessor1;
        while (iter >= 0) {
          predecessor_map[iter] = 1;
          iter = fp->layers[iter]->predecessors[0] - 1;
        }
        iter = predecessor2;
        int common_predecessor = 0;
        while (iter >= 0) {
          if (predecessor_map[iter] == 1) {
            common_predecessor = iter;
            break;
          }
          iter = fp->layers[iter]->predecessors[0] - 1;
        }

        iter = predecessor1;
        while (iter != common_predecessor) {

          update_state_using_predecessor_layer(pr, fp, &lexpr, &uexpr, iter,
                                               use_area_heuristic);

          iter = fp->layers[iter]->predecessors[0] - 1;
        }
        iter = predecessor2;
        while (iter != common_predecessor) {

          update_state_using_predecessor_layer(pr, fp, &lexpr_copy, &uexpr_copy,
                                               iter, use_area_heuristic);

          iter = fp->layers[iter]->predecessors[0] - 1;
        }
        free(predecessor_map);
        add_expr(pr, lexpr, lexpr_copy);
        add_expr(pr, uexpr, uexpr_copy);
        free_expr(lexpr_copy);
        free_expr(uexpr_copy);
        // Assume at least one non-residual layer between two residual layers

        k = common_predecessor;
        if (fp->layers[k]->activation == RELU) {
          expr_t *tmp_l = lexpr;
          expr_t *tmp_u = uexpr;
          lexpr = lexpr_replace_relu_bounds(pr, lexpr, fp->layers[k]->neurons,
                                            use_area_heuristic);
          uexpr = uexpr_replace_relu_bounds(pr, uexpr, fp->layers[k]->neurons,
                                            use_area_heuristic);
          free_expr(tmp_l);
          free_expr(tmp_u);
        }
        out_neurons[i]->lb = compute_lb_from_expr(pr, lexpr, fp, k);
        out_neurons[i]->ub = compute_ub_from_expr(pr, uexpr, fp, k);
        already_computed = true;
        break;
        // continue;
      } else {

        update_state_using_predecessor_layer(pr, fp, &lexpr, &uexpr, k,
                                             use_area_heuristic);

        k = fp->layers[k]->predecessors[0] - 1;
      }

      // printf("k %d\n",k);
    }

    if (!already_computed) {
      out_neurons[i]->lb = compute_lb_from_expr(pr, lexpr, fp, -1);
      //- bias_i;
      out_neurons[i]->ub = compute_ub_from_expr(pr, uexpr, fp, -1); //+ bias_i;
    }
    if (fp->out != NULL) {

      fp->out->lexpr[i] = lexpr;
      fp->out->uexpr[i] = uexpr;
    } else {
      free_expr(lexpr);
      free_expr(uexpr);
    }
  }
  // printf("thread finish\n");
  // fflush(stdout);
  return NULL;
}

void update_state_using_previous_layers_parallel(elina_manager_t *man,
                                                 fppoly_t *fp, size_t layerno,
                                                 bool use_area_heuristic) {
  // size_t NUM_THREADS = get_nprocs();
  size_t NUM_THREADS = sysconf(_SC_NPROCESSORS_ONLN);
  nn_thread_t args[NUM_THREADS];
  pthread_t threads[NUM_THREADS];
  size_t num_out_neurons = fp->layers[layerno]->dims;
  size_t i;
  // printf("layerno %zu %zu\n",layerno,fp->layers[layerno]->predecessors[0]-1);
  // fflush(stdout);
  if (num_out_neurons < NUM_THREADS) {
    for (i = 0; i < num_out_neurons; i++) {
      args[i].start = i;
      args[i].end = i + 1;
      args[i].man = man;
      args[i].fp = fp;
      args[i].layerno = layerno;
      args[i].use_area_heuristic = use_area_heuristic;
      pthread_create(&threads[i], NULL, update_state_using_previous_layers,
                     (void *)&args[i]);
    }
    for (i = 0; i < num_out_neurons; i = i + 1) {
      pthread_join(threads[i], NULL);
    }
  } else {
    size_t idx_start = 0;
    size_t idx_n = num_out_neurons / NUM_THREADS;
    size_t idx_end = idx_start + idx_n;

    for (i = 0; i < NUM_THREADS; i++) {
      args[i].start = idx_start;
      args[i].end = idx_end;
      args[i].man = man;
      args[i].fp = fp;
      args[i].layerno = layerno;
      args[i].use_area_heuristic = use_area_heuristic;
      pthread_create(&threads[i], NULL, update_state_using_previous_layers,
                     (void *)&args[i]);
      idx_start = idx_end;
      idx_end = idx_start + idx_n;
      if (idx_end > num_out_neurons) {
        idx_end = num_out_neurons;
      }
      if ((i == NUM_THREADS - 2)) {
        idx_end = num_out_neurons;
      }
    }
    // idx_start = idx_end;
    // idx_end = num_out_neurons;
    // args[i].start = idx_start;
    // args[i].end = idx_end;
    // args[i].man = man;
    // args[i].fp = fp;
    // args[i].layerno = layerno;
    // pthread_create(&threads[i], NULL,update_state_using_previous_layers,
    // (void*)&args[i]);
    for (i = 0; i < NUM_THREADS; i = i + 1) {
      pthread_join(threads[i], NULL);
    }
  }
  // printf("end\n");
  // fflush(stdout);
}

void ffn_handle_intermediate_layer(elina_manager_t *man,
                                   elina_abstract0_t *element, double **weights,
                                   double *bias, size_t num_out_neurons,
                                   size_t num_in_neurons, size_t *predecessors,
                                   activation_type_t activation, bool alloc,
                                   bool use_area_heuristic) {
  // printf("ReLU start here %zu %zu\n",num_in_neurons,num_out_neurons);
  // fflush(stdout);
  fppoly_t *fp = fppoly_of_abstract0(element);
  size_t numlayers = fp->numlayers;
  if (alloc) {
    fppoly_add_new_layer(fp, num_out_neurons, FFN, activation);
  }
  neuron_t **out_neurons = fp->layers[numlayers]->neurons;
  fp->layers[numlayers]->predecessors = predecessors;

  size_t i;
  for (i = 0; i < num_out_neurons; i++) {
    double *weight_i = weights[i];
    double bias_i = bias[i];

    out_neurons[i]->expr = create_dense_expr(weight_i, bias_i, num_in_neurons);
  }
  update_state_using_previous_layers_parallel(man, fp, numlayers,
                                              use_area_heuristic);

  // printf("return here2\n");
  // fppoly_fprint(stdout,man,fp,NULL);
  // fflush(stdout);
  return;
}

void ffn_handle_intermediate_affine_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, NONE, true,
                                use_area_heuristic);
}

void ffn_handle_intermediate_relu_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, RELU, true,
                                use_area_heuristic);
}

void ffn_handle_intermediate_sigmoid_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, SIGMOID, true,
                                use_area_heuristic);
}

void ffn_handle_intermediate_tanh_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, TANH, true,
                                use_area_heuristic);
}

void ffn_handle_intermediate_parabola_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {

  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, PARABOLA, true,
                                use_area_heuristic);
}

void ffn_handle_intermediate_log_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, LOG, true,
                                use_area_heuristic);
}

void ffn_handle_intermediate_affine_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, NONE, false,
                                use_area_heuristic);
}

void ffn_handle_intermediate_relu_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, RELU, false,
                                use_area_heuristic);
}

void ffn_handle_intermediate_sigmoid_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, SIGMOID, false,
                                use_area_heuristic);
}

void ffn_handle_intermediate_tanh_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, TANH, false,
                                use_area_heuristic);
}

void ffn_handle_intermediate_parabola_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {

  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, PARABOLA, false,
                                use_area_heuristic);
}

void ffn_handle_intermediate_log_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, LOG, false,
                                use_area_heuristic);
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

void handle_final_relu_layer(fppoly_internal_t *pr, output_abstract_t *out,
                             neuron_t **neurons, size_t size, bool has_relu) {
  size_t i;
  if (has_relu) {
    for (i = 0; i < size; i++) {
      out->output_inf[i] = apply_relu_lexpr(pr, &out->lexpr[i], neurons[i]);
      out->output_sup[i] = apply_relu_uexpr(pr, &out->uexpr[i], neurons[i]);
    }
  } else {
    for (i = 0; i < size; i++) {
      out->output_inf[i] = neurons[i]->lb;
      out->output_sup[i] = neurons[i]->ub;
    }
  }
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

double apply_sigmoid_lexpr(fppoly_internal_t *pr, expr_t **lexpr_p,
                           neuron_t *neuron) {
  return apply_s_curve_lexpr(pr, lexpr_p, neuron, true);
}

double apply_tanh_lexpr(fppoly_internal_t *pr, expr_t **lexpr_p,
                        neuron_t *neuron) {
  return apply_s_curve_lexpr(pr, lexpr_p, neuron, false);
}

double apply_sigmoid_uexpr(fppoly_internal_t *pr, expr_t **uexpr_p,
                           neuron_t *neuron) {
  return apply_s_curve_uexpr(pr, uexpr_p, neuron, true);
}

double apply_tanh_uexpr(fppoly_internal_t *pr, expr_t **uexpr_p,
                        neuron_t *neuron) {
  return apply_s_curve_uexpr(pr, uexpr_p, neuron, false);
}

void handle_final_sigmoid_layer(fppoly_internal_t *pr, output_abstract_t *out,
                                neuron_t **neurons, size_t size,
                                bool has_sigmoid) {
  size_t i;
  if (has_sigmoid) {
    for (i = 0; i < size; i++) {
      out->output_inf[i] = apply_sigmoid_lexpr(pr, &out->lexpr[i], neurons[i]);
      out->output_sup[i] = apply_sigmoid_uexpr(pr, &out->uexpr[i], neurons[i]);
    }
  } else {
    for (i = 0; i < size; i++) {
      out->output_inf[i] = neurons[i]->lb;
      out->output_sup[i] = neurons[i]->ub;
    }
  }
}

void handle_final_tanh_layer(fppoly_internal_t *pr, output_abstract_t *out,
                             neuron_t **neurons, size_t size, bool has_tanh) {
  size_t i;
  if (has_tanh) {
    for (i = 0; i < size; i++) {
      out->output_inf[i] = apply_tanh_lexpr(pr, &out->lexpr[i], neurons[i]);
      out->output_sup[i] = apply_tanh_uexpr(pr, &out->uexpr[i], neurons[i]);
    }
  } else {
    for (i = 0; i < size; i++) {
      out->output_inf[i] = neurons[i]->lb;
      out->output_sup[i] = neurons[i]->ub;
    }
  }
}

void neuron_fprint(FILE *stream, neuron_t *neuron, char **name_of_dim) {
  // expr_fprint(stream,neuron->expr);
  fprintf(stream, "[%g, %g]\n", -neuron->lb, neuron->ub);
}

void layer_fprint(FILE *stream, layer_t *layer, char **name_of_dim) {
  size_t dims = layer->dims;
  size_t i;
  for (i = 0; i < dims; i++) {
    fprintf(stream, "neuron: %zu ", i);
    neuron_fprint(stream, layer->neurons[i], name_of_dim);
  }
}

void ffn_handle_last_layer(elina_manager_t *man, elina_abstract0_t *element,
                           double **weights, double *bias,
                           size_t num_out_neurons, size_t num_in_neurons,
                           size_t *predecessors, bool has_activation,
                           activation_type_t activation, bool alloc,
                           bool use_area_heuristic) {
  // printf("last\n");
  // fflush(stdout);
  fppoly_t *fp = fppoly_of_abstract0(element);
  size_t numlayers = fp->numlayers;
  if (alloc) {
    if (has_activation) {
      fppoly_add_new_layer(fp, num_out_neurons, FFN, activation);
    } else {
      fppoly_add_new_layer(fp, num_out_neurons, FFN, NONE);
    }
  }
  output_abstract_t *out =
      (output_abstract_t *)malloc(sizeof(output_abstract_t));
  out->output_inf = (double *)malloc(num_out_neurons * sizeof(double));
  out->output_sup = (double *)malloc(num_out_neurons * sizeof(double));
  out->lexpr = (expr_t **)malloc(num_out_neurons * sizeof(expr_t *));
  out->uexpr = (expr_t **)malloc(num_out_neurons * sizeof(expr_t *));
  fp->out = out;
  neuron_t **out_neurons = fp->layers[numlayers]->neurons;
  fp->layers[numlayers]->predecessors = predecessors;
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  size_t i;
  for (i = 0; i < num_out_neurons; i++) {
    double *weight_i = weights[i];
    double bias_i = bias[i];
    out_neurons[i]->expr = create_dense_expr(weight_i, bias_i, num_in_neurons);
  }
  update_state_using_previous_layers_parallel(man, fp, numlayers,
                                              use_area_heuristic);
  if (activation == RELU) {
    handle_final_relu_layer(pr, fp->out, out_neurons, num_out_neurons,
                            has_activation);
  } else if (activation == SIGMOID) {
    handle_final_sigmoid_layer(pr, fp->out, out_neurons, num_out_neurons,
                               has_activation);
  } else if (activation == TANH) {
    handle_final_tanh_layer(pr, fp->out, out_neurons, num_out_neurons,
                            has_activation);
  }

  else {
    for (i = 0; i < num_out_neurons; i++) {
      out->output_inf[i] = out_neurons[i]->lb;
      out->output_sup[i] = out_neurons[i]->ub;
    }
  }
  // printf("finish\n");
  // fppoly_fprint(stdout,man,fp,NULL);
  // fflush(stdout);
  // layer_fprint(stdout,fp->layers[1],NULL);
  // layer_fprint(stdout,fp->layers[fp->numlayers-1],NULL);
  return;
}

void ffn_handle_last_relu_layer(elina_manager_t *man,
                                elina_abstract0_t *element, double **weights,
                                double *bias, size_t num_out_neurons,
                                size_t num_in_neurons, size_t *predecessors,
                                bool has_relu, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_relu, RELU, true,
                        use_area_heuristic);
}

void ffn_handle_last_sigmoid_layer(elina_manager_t *man,
                                   elina_abstract0_t *element, double **weights,
                                   double *bias, size_t num_out_neurons,
                                   size_t num_in_neurons, size_t *predecessors,
                                   bool has_sigmoid, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_sigmoid, SIGMOID,
                        true, use_area_heuristic);
}

void ffn_handle_last_tanh_layer(elina_manager_t *man,
                                elina_abstract0_t *element, double **weights,
                                double *bias, size_t num_out_neurons,
                                size_t num_in_neurons, size_t *predecessors,
                                bool has_tanh, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_tanh, TANH, true,
                        use_area_heuristic);
}

void ffn_handle_last_parabola_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool has_parabola, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_parabola, PARABOLA,
                        true, use_area_heuristic);
}

void ffn_handle_last_log_layer(elina_manager_t *man, elina_abstract0_t *element,
                               double **weights, double *bias,
                               size_t num_out_neurons, size_t num_in_neurons,
                               size_t *predecessors, bool has_log,
                               bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_log, LOG, true,
                        use_area_heuristic);
}

void ffn_handle_last_relu_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool has_relu, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_relu, RELU, false,
                        use_area_heuristic);
}

void ffn_handle_last_sigmoid_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool has_sigmoid, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_sigmoid, SIGMOID,
                        false, use_area_heuristic);
}

void ffn_handle_last_tanh_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool has_tanh, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_tanh, TANH, false,
                        use_area_heuristic);
}

void ffn_handle_last_parabola_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool has_parabola, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_parabola, PARABOLA,
                        false, use_area_heuristic);
}

void ffn_handle_last_log_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool has_log, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_log, LOG, false,
                        use_area_heuristic);
}

double get_lb_using_predecessor_layer(fppoly_internal_t *pr, fppoly_t *fp,
                                      expr_t **lexpr_ptr, size_t k,
                                      bool use_area_heuristic) {
  expr_t *tmp_l;
  neuron_t **aux_neurons = fp->layers[k]->neurons;
  expr_t *lexpr = *lexpr_ptr;
  double res = INFINITY;
  if (fp->layers[k]->type == FFN || fp->layers[k]->type == CONV) {

    if (fp->layers[k]->activation == RELU) {
      tmp_l = lexpr;
      lexpr =
          lexpr_replace_relu_bounds(pr, lexpr, aux_neurons, use_area_heuristic);
      free_expr(tmp_l);
    } else if (fp->layers[k]->activation == SIGMOID) {
      tmp_l = lexpr;
      // printf("start\n");
      // fflush(stdout);
      lexpr = lexpr_replace_sigmoid_bounds(pr, lexpr, aux_neurons);
      // printf("finish\n");
      // fflush(stdout);
      free_expr(tmp_l);
    } else if (fp->layers[k]->activation == TANH) {
      tmp_l = lexpr;
      lexpr = lexpr_replace_tanh_bounds(pr, lexpr, aux_neurons);
      free_expr(tmp_l);
    }

    else if (fp->layers[k]->activation == PARABOLA) {
      tmp_l = lexpr;
      lexpr = lexpr_replace_parabola_bounds(pr, lexpr, aux_neurons);
      free_expr(tmp_l);
    } else if (fp->layers[k]->activation == LOG) {
      tmp_l = lexpr;
      lexpr = lexpr_replace_log_bounds(pr, lexpr, aux_neurons);
      free_expr(tmp_l);
    }
    tmp_l = lexpr;
    res = compute_lb_from_expr(pr, lexpr, fp, k);
    *lexpr_ptr = expr_from_previous_layer(pr, lexpr, fp->layers[k]);
    free_expr(tmp_l);
  } else {
    expr_t *tmp_l = lexpr;
    *lexpr_ptr = lexpr_replace_maxpool_or_lstm_bounds(pr, lexpr, aux_neurons);
    free_expr(tmp_l);
  }
  return res;
}

double get_ub_using_predecessor_layer(fppoly_internal_t *pr, fppoly_t *fp,
                                      expr_t **uexpr_ptr, size_t k,
                                      bool use_area_heuristic) {
  expr_t *tmp_u;
  neuron_t **aux_neurons = fp->layers[k]->neurons;
  expr_t *uexpr = *uexpr_ptr;
  double res = INFINITY;
  if (fp->layers[k]->type == FFN || fp->layers[k]->type == CONV) {

    if (fp->layers[k]->activation == RELU) {
      tmp_u = uexpr;
      uexpr =
          uexpr_replace_relu_bounds(pr, uexpr, aux_neurons, use_area_heuristic);
      free_expr(tmp_u);
    } else if (fp->layers[k]->activation == SIGMOID) {
      tmp_u = uexpr;

      uexpr = uexpr_replace_sigmoid_bounds(pr, uexpr, aux_neurons);
      free_expr(tmp_u);
    } else if (fp->layers[k]->activation == TANH) {
      tmp_u = uexpr;
      uexpr = uexpr_replace_tanh_bounds(pr, uexpr, aux_neurons);
      free_expr(tmp_u);
    }

    else if (fp->layers[k]->activation == PARABOLA) {
      tmp_u = uexpr;
      uexpr = uexpr_replace_parabola_bounds(pr, uexpr, aux_neurons);
      free_expr(tmp_u);
    } else if (fp->layers[k]->activation == LOG) {
      tmp_u = uexpr;
      uexpr = uexpr_replace_log_bounds(pr, uexpr, aux_neurons);
      free_expr(tmp_u);
    }
    tmp_u = uexpr;
    res = compute_ub_from_expr(pr, uexpr, fp, k);
    *uexpr_ptr = expr_from_previous_layer(pr, uexpr, fp->layers[k]);
    free_expr(tmp_u);
  } else {
    expr_t *tmp_u = uexpr;
    *uexpr_ptr = uexpr_replace_maxpool_or_lstm_bounds(pr, uexpr, aux_neurons);
    free_expr(tmp_u);
  }
  return res;
}

double get_lb_using_previous_layers(elina_manager_t *man, fppoly_t *fp,
                                    expr_t *expr, size_t layerno,
                                    bool use_area_heuristic) {
  size_t i;
  int k;
  // size_t numlayers = fp->numlayers;
  expr_t *lexpr = copy_expr(expr);
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  if (fp->numlayers == layerno) {

    k = layerno - 1;
  } else {
    k = fp->layers[layerno]->predecessors[0] - 1;
  }
  double res = INFINITY;
  while (k >= 0) {

    if (fp->layers[k]->type == RESIDUAL) {
      if (fp->layers[k]->activation == RELU) {
        neuron_t **aux_neurons = fp->layers[k]->neurons;
        expr_t *tmp_l = lexpr;
        lexpr = lexpr_replace_relu_bounds(pr, lexpr, aux_neurons,
                                          use_area_heuristic);
        free_expr(tmp_l);
      }
      expr_t *lexpr_copy = copy_expr(lexpr);
      lexpr_copy->inf_cst = 0;
      lexpr_copy->sup_cst = 0;
      size_t predecessor1 = fp->layers[k]->predecessors[0] - 1;
      size_t predecessor2 = fp->layers[k]->predecessors[1] - 1;

      char *predecessor_map = (char *)calloc(k, sizeof(char));
      // Assume no nested residual layers
      int iter = fp->layers[predecessor1]->predecessors[0] - 1;
      while (iter >= 0) {
        predecessor_map[iter] = 1;
        iter = fp->layers[iter]->predecessors[0] - 1;
      }
      iter = fp->layers[predecessor2]->predecessors[0] - 1;
      int common_predecessor = 0;
      while (iter >= 0) {
        if (predecessor_map[iter] == 1) {
          common_predecessor = iter;
          break;
        }
        iter = fp->layers[iter]->predecessors[0] - 1;
      }

      iter = predecessor1;
      while (iter != common_predecessor) {
        get_lb_using_predecessor_layer(pr, fp, &lexpr, iter,
                                       use_area_heuristic);
        iter = fp->layers[iter]->predecessors[0] - 1;
      }
      iter = predecessor2;
      while (iter != common_predecessor) {
        get_lb_using_predecessor_layer(pr, fp, &lexpr_copy, iter,
                                       use_area_heuristic);
        iter = fp->layers[iter]->predecessors[0] - 1;
      }
      free(predecessor_map);
      add_expr(pr, lexpr, lexpr_copy);

      free_expr(lexpr_copy);

      // Assume at least one non-residual layer between two residual layers
      k = common_predecessor;

      continue;
    } else {

      res = fmin(res, get_lb_using_predecessor_layer(pr, fp, &lexpr, k,
                                                     use_area_heuristic));
      k = fp->layers[k]->predecessors[0] - 1;
    }
  }

  res = fmin(res, compute_lb_from_expr(pr, lexpr, fp, -1));
  free_expr(lexpr);
  return res;
}

double get_ub_using_previous_layers(elina_manager_t *man, fppoly_t *fp,
                                    expr_t *expr, size_t layerno,
                                    bool use_area_heuristic) {
  size_t i;
  int k;
  // size_t numlayers = fp->numlayers;
  expr_t *uexpr = copy_expr(expr);
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  if (fp->numlayers == layerno) {
    k = layerno - 1;
  } else {
    k = fp->layers[layerno]->predecessors[0] - 1;
  }
  double res = INFINITY;
  while (k >= 0) {
    if (fp->layers[k]->type == RESIDUAL) {
      if (fp->layers[k]->activation == RELU) {
        neuron_t **aux_neurons = fp->layers[k]->neurons;
        expr_t *tmp_u = uexpr;
        uexpr = uexpr_replace_relu_bounds(pr, uexpr, aux_neurons,
                                          use_area_heuristic);
        free_expr(tmp_u);
      }
      expr_t *uexpr_copy = copy_expr(uexpr);
      uexpr_copy->inf_cst = 0;
      uexpr_copy->sup_cst = 0;
      size_t predecessor1 = fp->layers[k]->predecessors[0] - 1;
      size_t predecessor2 = fp->layers[k]->predecessors[1] - 1;

      char *predecessor_map = (char *)calloc(k, sizeof(char));
      // Assume no nested residual layers
      int iter = fp->layers[predecessor1]->predecessors[0] - 1;
      while (iter >= 0) {
        predecessor_map[iter] = 1;
        iter = fp->layers[iter]->predecessors[0] - 1;
      }
      iter = fp->layers[predecessor2]->predecessors[0] - 1;
      int common_predecessor = 0;
      while (iter >= 0) {
        if (predecessor_map[iter] == 1) {
          common_predecessor = iter;
          break;
        }
        iter = fp->layers[iter]->predecessors[0] - 1;
      }

      iter = predecessor1;
      while (iter != common_predecessor) {
        get_ub_using_predecessor_layer(pr, fp, &uexpr, iter,
                                       use_area_heuristic);
        iter = fp->layers[iter]->predecessors[0] - 1;
      }
      iter = predecessor2;
      while (iter != common_predecessor) {
        get_ub_using_predecessor_layer(pr, fp, &uexpr_copy, iter,
                                       use_area_heuristic);
        iter = fp->layers[iter]->predecessors[0] - 1;
      }
      free(predecessor_map);
      add_expr(pr, uexpr, uexpr_copy);

      free_expr(uexpr_copy);

      // Assume at least one non-residual layer between two residual layers
      k = common_predecessor;

      continue;
    } else {

      res = fmin(res, get_ub_using_predecessor_layer(pr, fp, &uexpr, k,
                                                     use_area_heuristic));
      k = fp->layers[k]->predecessors[0] - 1;
    }
  }

  res = fmin(res, compute_ub_from_expr(pr, uexpr, fp, -1));
  free_expr(uexpr);
  return res;
}

void coeff_to_interval(elina_coeff_t *coeff, double *inf, double *sup) {
  double d;
  if (coeff->discr == ELINA_COEFF_SCALAR) {
    elina_scalar_t *scalar = coeff->val.scalar;
    d = scalar->val.dbl;
    *inf = -d;
    *sup = d;
  } else {
    elina_interval_t *interval = coeff->val.interval;
    d = interval->inf->val.dbl;
    *inf = -d;
    d = interval->sup->val.dbl;
    *sup = d;
  }
}

expr_t *elina_linexpr0_to_expr(elina_linexpr0_t *linexpr0) {
  size_t size = linexpr0->size;
  size_t i;
  expr_t *res = (expr_t *)malloc(sizeof(expr_t));
  res->inf_coeff = (double *)malloc(size * sizeof(double));
  res->sup_coeff = (double *)malloc(size * sizeof(double));
  res->size = size;
  if (linexpr0->discr == ELINA_LINEXPR_SPARSE) {
    res->type = SPARSE;
    res->dim = (size_t *)malloc(size * sizeof(size_t));
  } else {
    res->type = DENSE;
    res->dim = NULL;
  }
  size_t k;
  for (i = 0; i < size; i++) {
    elina_coeff_t *coeff;
    if (res->type == SPARSE) {
      k = linexpr0->p.linterm[i].dim;
      res->dim[i] = k;
      coeff = &linexpr0->p.linterm[i].coeff;
      coeff_to_interval(coeff, &res->inf_coeff[i], &res->sup_coeff[i]);
    } else {
      k = i;
      coeff = &linexpr0->p.coeff[k];
      coeff_to_interval(coeff, &res->inf_coeff[k], &res->sup_coeff[k]);
    }
  }
  elina_coeff_t *cst = &linexpr0->cst;
  coeff_to_interval(cst, &res->inf_cst, &res->sup_cst);
  return res;
}

elina_interval_t *get_bounds_for_linexpr0(elina_manager_t *man,
                                          elina_abstract0_t *element,
                                          elina_linexpr0_t *linexpr0,
                                          size_t layerno) {

  elina_interval_t *res = elina_interval_alloc();
  fppoly_t *fp = fppoly_of_abstract0(element);
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  // printf("start %p %lu\n",fp->layers[layerno],layerno);
  // fflush(stdout);
  expr_t *tmp = elina_linexpr0_to_expr(linexpr0);
  // printf("coming here\n");
  // fflush(stdout);
  expr_t *expr = expr_from_previous_layer(pr, tmp, fp->layers[layerno]);
  // printf("end\n");
  // fflush(stdout);
  expr_t *expr2 = copy_expr(expr);

  double lb = get_lb_using_previous_layers(man, fp, expr, layerno, true);

  double ub = get_ub_using_previous_layers(man, fp, expr2, layerno, true);

  elina_interval_set_double(res, -lb, ub);
  free_expr(expr);
  free_expr(expr2);
  free_expr(tmp);
  return res;
}

void create_lstm_layer(elina_manager_t *man, elina_abstract0_t *abs, size_t h,
                       size_t *predecessors) {
  fppoly_t *fp = fppoly_of_abstract0(abs);
  size_t numlayers = fp->numlayers;
  fppoly_add_new_layer(fp, h, LSTM, NONE);
  fp->lstm_index = numlayers;
}

void handle_lstm_layer(elina_manager_t *man, elina_abstract0_t *abs,
                       double **weights, double *bias, size_t d, size_t h,
                       size_t *predecessors, bool use_area_heuristic) {
  fppoly_t *fp = fppoly_of_abstract0(abs);
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  size_t lstm_index = fp->lstm_index;
  layer_t *layer = fp->layers[lstm_index];
  neuron_t **out_neurons = fp->layers[lstm_index]->neurons;
  fp->layers[lstm_index]->predecessors = predecessors;
  size_t i;
  neuron_t *neuron = neuron_alloc();
  bool first_time_step = (layer->h_t_inf == NULL && layer->h_t_sup == NULL);
  size_t k = h + d;
  if (first_time_step) {
    layer->h_t_inf = (double *)malloc(h * sizeof(double));
    layer->h_t_sup = (double *)malloc(h * sizeof(double));
    layer->c_t_inf = (double *)malloc(h * sizeof(double));
    layer->c_t_sup = (double *)malloc(h * sizeof(double));
  }

  // TODO: Fix, debug: 	for(i=0; i< h; i++){
  for (i = 0; i < 1; i++) {
    // printf("i = %d\n",(int)i);
    expr_t *f_t_lexpr, *i_t_lexpr, *o_t_lexpr, *c_t_lexpr;
    if (first_time_step) {
      i_t_lexpr = create_dense_expr(weights[i], bias[i], d);
      c_t_lexpr = create_dense_expr(weights[h + i], bias[h + i], d);
      f_t_lexpr = create_dense_expr(weights[2 * h + i], bias[2 * h + i], d);
      o_t_lexpr = create_dense_expr(weights[3 * h + i], bias[3 * h + i], d);
    } else {
      expr_t *tmp1 = create_dense_expr(weights[i], bias[i], d + h);
      expr_t *tmp2 = create_dense_expr(weights[h + i], bias[h + i], d + h);
      expr_t *tmp3 =
          create_dense_expr(weights[2 * h + i], bias[2 * h + i], d + h);
      expr_t *tmp4 =
          create_dense_expr(weights[3 * h + i], bias[3 * h + i], d + h);
      i_t_lexpr = concretize_dense_sub_expr(pr, tmp1, layer->h_t_inf,
                                            layer->h_t_sup, d, d + h);
      c_t_lexpr = concretize_dense_sub_expr(pr, tmp2, layer->h_t_inf,
                                            layer->h_t_sup, d, d + h);
      f_t_lexpr = concretize_dense_sub_expr(pr, tmp3, layer->h_t_inf,
                                            layer->h_t_sup, d, d + h);
      o_t_lexpr = concretize_dense_sub_expr(pr, tmp4, layer->h_t_inf,
                                            layer->h_t_sup, d, d + h);
      free_expr(tmp1);
      free_expr(tmp2);
      free_expr(tmp3);
      free_expr(tmp4);
    }

    // expr_print(f_t_lexpr);

    // printf("computing forget...\n");
    expr_t *f_t_uexpr = copy_expr(f_t_lexpr);
    expr_t *tmp_f_t_lexpr = copy_expr(f_t_lexpr);
    expr_t *tmp_f_t_uexpr = copy_expr(f_t_uexpr);
    double lb_f_t = get_lb_using_previous_layers(
        man, fp, tmp_f_t_lexpr, lstm_index, use_area_heuristic);
    double ub_f_t = get_ub_using_previous_layers(
        man, fp, tmp_f_t_uexpr, lstm_index, use_area_heuristic);
    /* free_expr(tmp_f_t_lexpr); */
    /* free_expr(tmp_f_t_uexpr); */

    neuron->lb = lb_f_t;
    neuron->ub = ub_f_t;
    // printf("forget gate before sigmoid: lb = %lf, ub = %lf\n",neuron->lb,
    // neuron->ub); expr_print(f_t_lexpr); expr_print(f_t_uexpr);
    lb_f_t = apply_sigmoid_lexpr(pr, &f_t_lexpr, neuron);
    ub_f_t = apply_sigmoid_uexpr(pr, &f_t_uexpr, neuron);
    // printf("forget gate after sigmoid: lb_f_t = %lf, ub_f_t =
    // %lf\n",lb_f_t,ub_f_t); expr_print(f_t_lexpr); expr_print(f_t_uexpr);
    // printf("forget gate done\n\n");

    // printf("computing input...\n");
    expr_t *i_t_uexpr = copy_expr(i_t_lexpr);
    expr_t *tmp_i_t_lexpr = copy_expr(i_t_lexpr);
    expr_t *tmp_i_t_uexpr = copy_expr(i_t_uexpr);
    double lb_i_t = get_lb_using_previous_layers(
        man, fp, tmp_i_t_lexpr, lstm_index, use_area_heuristic);
    double ub_i_t = get_ub_using_previous_layers(
        man, fp, tmp_i_t_uexpr, lstm_index, use_area_heuristic);
    /* free_expr(tmp_i_t_lexpr); */
    /* free_expr(tmp_i_t_uexpr); */
    neuron->lb = lb_i_t;
    neuron->ub = ub_i_t;
    // printf("input gate before sigmoid: lb = %lf, ub = %lf\n",neuron->lb,
    // neuron->ub); expr_print(i_t_uexpr);
    lb_i_t = apply_sigmoid_lexpr(pr, &i_t_lexpr, neuron);
    ub_i_t = apply_sigmoid_uexpr(pr, &i_t_uexpr, neuron);
    // expr_print(i_t_uexpr);
    // printf("input gate after sigmoid: lb_i_t = %lf, ub_i_t =
    // %lf\n",lb_i_t,ub_i_t); printf("input gate done\n\n");

    // printf("computing output...\n");
    expr_t *o_t_uexpr = copy_expr(o_t_lexpr);
    expr_t *tmp_o_t_lexpr = copy_expr(o_t_lexpr);
    expr_t *tmp_o_t_uexpr = copy_expr(o_t_uexpr);
    double lb_o_t = get_lb_using_previous_layers(
        man, fp, tmp_o_t_lexpr, lstm_index, use_area_heuristic);
    double ub_o_t = get_ub_using_previous_layers(
        man, fp, tmp_o_t_uexpr, lstm_index, use_area_heuristic);
    /* free_expr(tmp_o_t_lexpr); */
    /* free_expr(tmp_o_t_uexpr); */

    neuron->lb = lb_o_t;
    neuron->ub = ub_o_t;
    // printf("output gate before sigmoid: lb = %lf, ub = %lf\n",neuron->lb,
    // neuron->ub);
    lb_o_t = apply_sigmoid_lexpr(pr, &o_t_lexpr, neuron);
    ub_o_t = apply_sigmoid_uexpr(pr, &o_t_uexpr, neuron);
    // printf("output gate after sigmoid: lb = %lf, ub = %lf\n",lb_o_t,ub_o_t);
    out_neurons[i]->lb = lb_o_t;
    out_neurons[i]->ub = ub_o_t;
    out_neurons[i]->lexpr = o_t_lexpr;
    out_neurons[i]->uexpr = o_t_uexpr;
    // printf("output gate done\n\n");

    // printf("computing control state...\n");
    // printf("control expression:\n");
    // expr_print(c_t_lexpr);
    // printf("...\n");
    expr_t *c_t_uexpr = copy_expr(c_t_lexpr);
    expr_t *tmp_c_t_lexpr = copy_expr(c_t_lexpr);
    expr_t *tmp_c_t_uexpr = copy_expr(c_t_uexpr);
    double lb_c_t = get_lb_using_previous_layers(
        man, fp, tmp_c_t_lexpr, lstm_index, use_area_heuristic);
    double ub_c_t = get_ub_using_previous_layers(
        man, fp, tmp_c_t_uexpr, lstm_index, use_area_heuristic);
    neuron->lb = lb_c_t;
    neuron->ub = ub_c_t;
    // expr_print(c_t_lexpr);
    // expr_print(c_t_uexpr);
    // printf("control before tanh: lb = %lf, ub =
    // %lf\n",neuron->lb,neuron->ub);
    lb_c_t = apply_tanh_lexpr(pr, &c_t_lexpr, neuron);
    ub_c_t = apply_tanh_uexpr(pr, &c_t_uexpr, neuron);
    // printf("control after tanh: lb = %lf, ub = %lf\n",lb_c_t,ub_c_t);
    // printf("control expression:\n");
    // expr_print(c_t_lexpr);
    // expr_print(c_t_uexpr);

    // printf("=======================\n");

    // printf("multiplying control by input:\n");
    expr_t *tmp_l, *tmp_u;
    double width1 = ub_i_t + lb_i_t;
    double width2 = ub_c_t + lb_c_t;
    tmp_l = c_t_lexpr;
    tmp_u = c_t_uexpr;
    // printf("control: [%lf %lf], input: [%lf
    // %lf]\n",lb_c_t,ub_c_t,lb_i_t,ub_i_t); printf("control before multiplying
    // by input:\n"); expr_print(c_t_lexpr); expr_print(c_t_uexpr);
    if (width1 < width2) {
      // printf("concretize input\n");
      c_t_lexpr = multiply_expr(pr, c_t_lexpr, lb_i_t, ub_i_t);
      c_t_uexpr = multiply_expr(pr, c_t_uexpr, lb_i_t, ub_i_t);
    } else {
      // printf("concretize control\n");
      if (lb_c_t < 0) {
        c_t_lexpr = multiply_expr(pr, i_t_lexpr, lb_c_t, ub_c_t);
        c_t_uexpr = multiply_expr(pr, i_t_uexpr, lb_c_t, ub_c_t);
      } else if (ub_c_t < 0) {
        c_t_lexpr = multiply_expr(pr, i_t_uexpr, lb_c_t, ub_c_t);
        c_t_uexpr = multiply_expr(pr, i_t_lexpr, lb_c_t, ub_c_t);
      } else {
        c_t_lexpr = multiply_expr(pr, i_t_lexpr, 0, 0);
        c_t_uexpr = multiply_expr(pr, i_t_uexpr, 0, 0);
        double tmp1, tmp2;
        elina_double_interval_mul_expr_coeff(pr, &tmp1, &tmp2, lb_i_t, ub_i_t,
                                             lb_c_t, ub_c_t);
        c_t_lexpr->inf_cst += tmp1;
        c_t_lexpr->sup_cst += tmp2;
        c_t_uexpr->inf_cst += tmp1;
        c_t_uexpr->sup_cst += tmp2;
      }
    }

    // printf("control after multiplying by input:\n");
    // expr_print(c_t_lexpr);
    // expr_print(c_t_uexpr);

    free_expr(tmp_l);
    free_expr(tmp_u);

    // printf("here\n\n\n");
    // printf("====================================\n");

    if (!first_time_step) {
      if (layer->c_t_inf[i] < 0) {
        tmp_l =
            multiply_expr(pr, f_t_lexpr, layer->c_t_inf[i], layer->c_t_sup[i]);
        tmp_u =
            multiply_expr(pr, f_t_uexpr, layer->c_t_inf[i], layer->c_t_sup[i]);
      } else if (layer->c_t_sup[i] < 0) {
        tmp_l =
            multiply_expr(pr, f_t_uexpr, layer->c_t_inf[i], layer->c_t_sup[i]);
        tmp_u =
            multiply_expr(pr, f_t_lexpr, layer->c_t_inf[i], layer->c_t_sup[i]);
      } else {
        tmp_l = multiply_expr(pr, f_t_lexpr, 0, 0);
        tmp_u = multiply_expr(pr, f_t_uexpr, 0, 0);
        double tmp1, tmp2;
        elina_double_interval_mul_expr_coeff(pr, &tmp1, &tmp2, lb_f_t, ub_f_t,
                                             layer->c_t_inf[i],
                                             layer->c_t_sup[i]);
        tmp_l->inf_cst += tmp1;
        tmp_l->sup_cst += tmp2;
        tmp_u->inf_cst += tmp1;
        tmp_u->sup_cst += tmp2;
      }
      add_expr(pr, c_t_lexpr, tmp_l);
      add_expr(pr, c_t_uexpr, tmp_u);
      free_expr(tmp_l);
      free_expr(tmp_u);
    }
    layer->c_t_inf[i] = get_lb_using_previous_layers(
        man, fp, c_t_lexpr, lstm_index, use_area_heuristic);
    layer->c_t_sup[i] = get_ub_using_previous_layers(
        man, fp, c_t_uexpr, lstm_index, use_area_heuristic);

    neuron->lb = layer->c_t_inf[i];
    neuron->ub = layer->c_t_sup[i];

    // printf("c_t ---> lb = %lf, ub = %lf\n", neuron->lb, neuron->ub);

    lb_c_t = apply_tanh_lexpr(pr, &c_t_lexpr, neuron);
    ub_c_t = apply_tanh_uexpr(pr, &c_t_uexpr, neuron);

    width1 = ub_o_t + lb_o_t;
    width2 = ub_c_t + lb_c_t;

    expr_t *h_t_lexpr, *h_t_uexpr;
    if (width1 < width2) {
      h_t_lexpr = multiply_expr(pr, c_t_lexpr, lb_o_t, ub_o_t);
      h_t_uexpr = multiply_expr(pr, c_t_uexpr, lb_o_t, ub_o_t);
    } else {
      h_t_lexpr = multiply_expr(pr, o_t_lexpr, lb_c_t, ub_c_t);
      h_t_uexpr = multiply_expr(pr, o_t_uexpr, lb_c_t, ub_c_t);
    }

    layer->h_t_inf[i] = get_lb_using_previous_layers(
        man, fp, h_t_lexpr, lstm_index, use_area_heuristic);
    layer->h_t_sup[i] = get_ub_using_previous_layers(
        man, fp, h_t_uexpr, lstm_index, use_area_heuristic);

    free_expr(f_t_lexpr);
    free_expr(f_t_uexpr);
    free_expr(i_t_lexpr);
    free_expr(i_t_uexpr);
    free_expr(c_t_lexpr);
    free_expr(c_t_uexpr);
    free_expr(h_t_lexpr);
    free_expr(h_t_uexpr);
  }
  free_neuron(neuron);
  // update_state_using_previous_layers_parallel(man,fp,numlayers);
  return;
}

bool is_greater(elina_manager_t *man, elina_abstract0_t *element, elina_dim_t y,
                elina_dim_t x, bool use_area_heuristic) {
  fppoly_t *fp = fppoly_of_abstract0(element);
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  if (1) {

    expr_t *sub = (expr_t *)malloc(sizeof(expr_t));
    // sub->size = size;
    sub->inf_cst = 0;
    sub->sup_cst = 0;
    sub->inf_coeff = (double *)malloc(2 * sizeof(double));
    sub->sup_coeff = (double *)malloc(2 * sizeof(double));
    sub->dim = (size_t *)malloc(2 * sizeof(size_t));
    sub->size = 2;
    sub->type = SPARSE;
    sub->inf_coeff[0] = -1;
    sub->sup_coeff[0] = 1;
    sub->dim[0] = y;
    sub->inf_coeff[1] = 1;
    sub->sup_coeff[1] = -1;
    sub->dim[1] = x;

    // layer_fprint(stdout,fp->layers[3],NULL);
    double lb = get_lb_using_previous_layers(man, fp, sub, fp->numlayers,
                                             use_area_heuristic);

    // free_expr(sub);

    if (lb < 0) {
      return true;
    } else {
      return false;
    }
  }
  output_abstract_t *out = fp->out;
  expr_t *exprA = out->lexpr[y];
  expr_t *exprB = out->uexpr[x];
  if (exprA == NULL) {
    return false;
  } else {
    if (exprB == NULL) {
      if (out->output_inf[y] < 0) {
        return true;
      } else {
        return false;
      }

    } else {
      // printf("before access %zu %zu\n", exprA->size,exprB->size);
      // fflush(stdout);
      size_t sizeA = exprA->size;
      size_t sizeB = exprB->size;
      // printf("after access\n");
      // fflush(stdout);
      size_t i, k;
      expr_t *sub = (expr_t *)malloc(sizeof(expr_t));
      //
      // sub->size = size;
      sub->inf_cst = exprA->inf_cst + exprB->sup_cst;
      sub->sup_cst = exprA->sup_cst + exprB->inf_cst;
      // printf("getting here\n");
      // expr_print(exprA);
      // expr_print(exprB);
      // fflush(stdout);
      if (exprA->type == DENSE) {
        sub->inf_coeff = (double *)malloc(sizeA * sizeof(double));
        sub->sup_coeff = (double *)malloc(sizeA * sizeof(double));
        sub->dim = NULL;
        sub->size = sizeA;
        sub->type = DENSE;
        if (exprB->type == DENSE) {
          for (i = 0; i < sizeA; i++) {
            sub->inf_coeff[i] = exprA->inf_coeff[i] + exprB->sup_coeff[i];
            sub->sup_coeff[i] = exprA->sup_coeff[i] + exprB->inf_coeff[i];
          }
        } else {
          k = 0;
          for (i = 0; i < sizeA; i++) {
            if (k < sizeB && exprB->dim[k] == i) {
              sub->inf_coeff[i] = exprA->inf_coeff[i] + exprB->sup_coeff[k];
              sub->sup_coeff[i] = exprA->sup_coeff[i] + exprB->inf_coeff[k];
              k++;
            } else {
              sub->inf_coeff[i] = exprA->inf_coeff[i];
              sub->sup_coeff[i] = exprA->sup_coeff[i];
            }
          }
        }

      } else {
        if (exprB->type == DENSE) {
          sub->inf_coeff = (double *)malloc(sizeB * sizeof(double));
          sub->sup_coeff = (double *)malloc(sizeB * sizeof(double));
          sub->dim = NULL;
          sub->size = sizeB;
          sub->type = DENSE;
          i = 0;
          for (k = 0; k < sizeB; k++) {
            if (i < sizeA && exprA->dim[i] == k) {
              sub->inf_coeff[k] = exprA->inf_coeff[i] + exprB->sup_coeff[k];
              sub->sup_coeff[k] = exprA->sup_coeff[i] + exprB->inf_coeff[k];
              i++;
            } else {
              sub->inf_coeff[i] = exprB->sup_coeff[k];
              sub->sup_coeff[i] = exprB->inf_coeff[k];
            }
          }
        } else {
          sub->inf_coeff = (double *)malloc((sizeA + sizeB) * sizeof(double));
          sub->sup_coeff = (double *)malloc((sizeA + sizeB) * sizeof(double));
          sub->dim = NULL;

          sub->type = SPARSE;
          size_t l = 0;
          i = 0;
          k = 0;
          sub->dim = (size_t *)malloc((sizeA + sizeB) * sizeof(size_t));
          while (i < sizeA && k < sizeB) {
            if (exprA->dim[i] < exprB->dim[k]) {
              sub->inf_coeff[l] = exprA->inf_coeff[i];
              sub->sup_coeff[l] = exprA->sup_coeff[i];
              sub->dim[l] = exprA->dim[i];
              i++;

            } else if (exprB->dim[k] < exprA->dim[i]) {
              sub->inf_coeff[l] = exprB->sup_coeff[k];
              sub->sup_coeff[l] = exprB->inf_coeff[k];
              sub->dim[l] = exprB->dim[k];
              k++;
            } else {
              sub->inf_coeff[l] = exprA->inf_coeff[i] + exprB->sup_coeff[k];
              sub->sup_coeff[l] = exprA->sup_coeff[i] + exprB->inf_coeff[k];
              sub->dim[l] = exprA->dim[i];
              i++;
              k++;
            }
            l++;
          }
          while (i < sizeA) {
            sub->inf_coeff[l] = exprA->inf_coeff[i];
            sub->sup_coeff[l] = exprA->sup_coeff[i];
            sub->dim[l] = exprA->dim[i];
            i++;
            l++;
          }
          while (k < sizeB) {
            sub->inf_coeff[l] = exprB->inf_coeff[k];
            sub->sup_coeff[l] = exprB->sup_coeff[k];
            sub->dim[l] = exprB->dim[k];
            k++;
            l++;
          }
          sub->size = l;
          sub->inf_coeff =
              (double *)realloc(sub->inf_coeff, l * sizeof(double));
          sub->sup_coeff =
              (double *)realloc(sub->sup_coeff, l * sizeof(double));
          sub->dim = (size_t *)realloc(sub->dim, l * sizeof(size_t));
        }
      }

      // expr_print(sub);
      // fflush(stdout);
      double lb = compute_lb_from_expr(pr, sub, fp, -1);
      // printf("y: %zu x: %zu lb: %g\n",y,x,lb);
      // fflush(stdout);
      free_expr(sub);
      // double lb = -out->output_inf[y] - out->output_sup[x];
      if (lb < 0) {
        return true;
      } else {
        return false;
      }
    }
  }
}

long int max(long int a, long int b){
	return a> b? a : b;

}

void conv_handle_first_layer(elina_manager_t *man, elina_abstract0_t *abs,
                             double *filter_weights, double *filter_bias,
                             size_t *input_size, size_t *filter_size,
                             size_t num_filters, size_t *strides,
                             size_t *output_size, size_t pad_top,
                             size_t pad_left, bool has_bias,
                             size_t *predecessors) {

  size_t i, j;
  size_t num_pixels = input_size[0] * input_size[1] * input_size[2];

  size_t size = output_size[0] * output_size[1] * output_size[2];

  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  fppoly_t *res = fppoly_of_abstract0(abs);
  fppoly_alloc_first_layer(res, size, CONV, RELU);

  neuron_t **neurons = res->layers[0]->neurons;
  res->layers[0]->predecessors = predecessors;
  size_t out_x, out_y, out_z;
  size_t inp_x, inp_y, inp_z;
  size_t x_shift, y_shift;

  for (out_x = 0; out_x < output_size[0]; out_x++) {
    for (out_y = 0; out_y < output_size[1]; out_y++) {
      for (out_z = 0; out_z < output_size[2]; out_z++) {
        size_t mat_x = out_x * output_size[1] * output_size[2] +
                       out_y * output_size[2] + out_z;
        size_t num_coeff = input_size[2] * filter_size[0] * filter_size[1];
        size_t actual_coeff = 0;
        double *coeff = (double *)malloc(num_coeff * sizeof(double));
        size_t *dim = (size_t *)malloc(num_coeff * sizeof(double));
        i = 0;
        for (inp_z = 0; inp_z < input_size[2]; inp_z++) {
          for (x_shift = 0; x_shift < filter_size[0]; x_shift++) {
            for (y_shift = 0; y_shift < filter_size[1]; y_shift++) {
              long int x_val = out_x * strides[0] + x_shift - pad_top;
              long int y_val = out_y * strides[1] + y_shift - pad_left;
              if (y_val < 0 || y_val >= (long int)input_size[1]) {
                continue;
              }

              if (x_val < 0 || x_val >= (long int)input_size[0]) {
                continue;
              }
              size_t mat_y = x_val * input_size[1] * input_size[2] +
                             y_val * input_size[2] + inp_z;
              if (mat_y >= num_pixels) {
                continue;
              }
              size_t filter_index =
                  x_shift * filter_size[1] * input_size[2] * output_size[2] +
                  y_shift * input_size[2] * output_size[2] +
                  inp_z * output_size[2] + out_z;
              coeff[i] = filter_weights[filter_index];
              dim[i] = mat_y;
              i++;
              actual_coeff++;
            }
          }
        }
        double cst = has_bias ? filter_bias[out_z] : 0;
        neurons[mat_x]->expr =
            create_sparse_expr(coeff, cst, dim, actual_coeff);
        sort_sparse_expr(neurons[mat_x]->expr);

        neurons[mat_x]->lb =
            compute_lb_from_expr(pr, neurons[mat_x]->expr, res, -1);
        neurons[mat_x]->ub =
            compute_ub_from_expr(pr, neurons[mat_x]->expr, res, -1);
        free(coeff);
        free(dim);
      }
    }
  }

  // printf("return here\n");
  // fppoly_fprint(stdout,man,res,NULL);
  // fflush(stdout);
  return;
}

void conv_handle_intermediate_layer(
    elina_manager_t *man, elina_abstract0_t *element, double *filter_weights,
    double *filter_bias, size_t *input_size, size_t *filter_size,
    size_t num_filters, size_t *strides, size_t *output_size, size_t pad_top,
    size_t pad_left, bool has_bias, size_t *predecessors,
    activation_type_t activation, bool use_area_heuristic) {
  // printf("conv intermediate starts here\n");
  // fflush(stdout);
  fppoly_t *fp = fppoly_of_abstract0(element);
  size_t numlayers = fp->numlayers;
  size_t i, j;
  size_t num_pixels = input_size[0] * input_size[1] * input_size[2];

  output_size[2] = num_filters;
  size_t num_out_neurons = output_size[0] * output_size[1] * output_size[2];
  // printf("num_out_neurons: %zu %zu\n",num_out_neurons,num_pixels);
  // fflush(stdout);
  fppoly_add_new_layer(fp, num_out_neurons, CONV, activation);
  neuron_t **out_neurons = fp->layers[numlayers]->neurons;
  fp->layers[numlayers]->predecessors = predecessors;
  size_t out_x, out_y, out_z;
  size_t inp_x, inp_y, inp_z;
  size_t x_shift, y_shift;

  for (out_x = 0; out_x < output_size[0]; out_x++) {
    for (out_y = 0; out_y < output_size[1]; out_y++) {
      for (out_z = 0; out_z < output_size[2]; out_z++) {
        size_t mat_x = out_x * output_size[1] * output_size[2] +
                       out_y * output_size[2] + out_z;
        size_t num_coeff = input_size[2] * filter_size[0] * filter_size[1];
        size_t actual_coeff = 0;
        double *coeff = (double *)malloc(num_coeff * sizeof(double));
        size_t *dim = (size_t *)malloc(num_coeff * sizeof(double));
        i = 0;
        for (inp_z = 0; inp_z < input_size[2]; inp_z++) {
          for (x_shift = 0; x_shift < filter_size[0]; x_shift++) {
            for (y_shift = 0; y_shift < filter_size[1]; y_shift++) {
              long int x_val = out_x * strides[0] + x_shift - pad_top;
              long int y_val = out_y * strides[1] + y_shift - pad_left;
              if (y_val < 0 || y_val >= (long int)input_size[1]) {
                continue;
              }

              if (x_val < 0 || x_val >= (long int)input_size[0]) {
                continue;
              }
              size_t mat_y = x_val * input_size[1] * input_size[2] +
                             y_val * input_size[2] + inp_z;
              if (mat_y >= num_pixels) {
                continue;
              }
              size_t filter_index =
                  x_shift * filter_size[1] * input_size[2] * output_size[2] +
                  y_shift * input_size[2] * output_size[2] +
                  inp_z * output_size[2] + out_z;
              coeff[i] = filter_weights[filter_index];
              dim[i] = mat_y;
              actual_coeff++;
              i++;
            }
          }
        }
        double cst = has_bias ? filter_bias[out_z] : 0;
        out_neurons[mat_x]->expr =
            create_sparse_expr(coeff, cst, dim, actual_coeff);
        sort_sparse_expr(out_neurons[mat_x]->expr);
        free(coeff);
        free(dim);
      }
    }
  }

  update_state_using_previous_layers_parallel(man, fp, numlayers,
                                              use_area_heuristic);

  // printf("return here2\n");
  // fppoly_fprint(stdout,man,fp,NULL);
  // fflush(stdout);
  return;
}

void conv_handle_intermediate_relu_layer(
    elina_manager_t *man, elina_abstract0_t *element, double *filter_weights,
    double *filter_bias, size_t *input_size, size_t *filter_size,
    size_t num_filters, size_t *strides, size_t *output_size, size_t pad_top,
    size_t pad_left, bool has_bias, size_t *predecessors,
    bool use_area_heuristic) {

  conv_handle_intermediate_layer(man, element, filter_weights, filter_bias,
                                 input_size, filter_size, num_filters, strides,
                                 output_size, pad_top, pad_left, has_bias,
                                 predecessors, RELU, use_area_heuristic);
}

void conv_handle_intermediate_affine_layer(
    elina_manager_t *man, elina_abstract0_t *element, double *filter_weights,
    double *filter_bias, size_t *input_size, size_t *filter_size,
    size_t num_filters, size_t *strides, size_t *output_size, size_t pad_top,
    size_t pad_left, bool has_bias, size_t *predecessors,
    bool use_area_heuristic) {

  conv_handle_intermediate_layer(man, element, filter_weights, filter_bias,
                                 input_size, filter_size, num_filters, strides,
                                 output_size, pad_top, pad_left, has_bias,
                                 predecessors, NONE, use_area_heuristic);
}

size_t handle_maxpool_layer(elina_manager_t *man, elina_abstract0_t *element,
                            size_t *pool_size, size_t *input_size,
                            size_t *predecessors) {
  // assert(dimensionality==3);
  // printf("maxpool start\n");
  // fflush(stdout);
  // printf("maxpool start\n");
  // fflush(stdout);
  assert(pool_size[0] == 2 && pool_size[1] == 2 && pool_size[2] == 1);
  // assert(stride[0]==2 && stride[1]==2 && stride[2]==1);

  size_t i, j, k;
  size_t *output_size = (size_t *)malloc(3 * sizeof(size_t));
  for (i = 0; i < 3; i++) {
    output_size[i] = input_size[i] / pool_size[i];
  }

  size_t num_input_neurons = input_size[0] * input_size[1] * input_size[2];
  size_t num_out_neurons = output_size[0] * output_size[1] * output_size[2];

  size_t o12 = output_size[1] * output_size[2];
  size_t i12 = input_size[1] * input_size[2];
  size_t p01 = pool_size[0] * pool_size[1];

  fppoly_t *fp = fppoly_of_abstract0(element);
  size_t numlayers = fp->numlayers;
  fppoly_add_new_layer(fp, num_out_neurons, MAXPOOL, NONE);
  size_t out_pos;
  double *inf = (double *)calloc(p01, sizeof(double));
  double *sup = (double *)calloc(p01, sizeof(double));
  size_t *pool_map = (size_t *)calloc(p01, sizeof(size_t));
  neuron_t **out_neurons = fp->layers[numlayers]->neurons;
  fp->layers[numlayers]->predecessors = predecessors;
  size_t count = 0;
  for (out_pos = 0; out_pos < num_out_neurons; out_pos++) {
    size_t out_x = out_pos / o12;
    size_t out_y = (out_pos - out_x * o12) / output_size[2];
    size_t out_z = out_pos - out_x * o12 - out_y * output_size[2];
    size_t inp_x = out_x * pool_size[0];
    size_t inp_y = out_y * pool_size[1];
    size_t inp_z = out_z;
    size_t inp_pos = inp_x * i12 + inp_y * input_size[2] + inp_z;
    size_t pool_start_dim = out_pos * pool_size[0] * pool_size[1];
    // printf("inpXYZ: %zu, %zu, %zu %zu %zu\n", inp_x, inp_y, inp_z, out_pos,
    // num_out_neurons); printf("outXYZ: %zu, %zu, %zu\n", out_x, out_y, out_z);
    // fflush(stdout);
    size_t x_shift, y_shift, l = 0;
    double sum_u = 0.0;
    double sum_l = 0.0;
    double max_u = -INFINITY;
    double max_l = -INFINITY;

    size_t max_l_var = 0.0;
    size_t max_u_var = 0.0;
    size_t min_width_var = 0.0;
    double min_width = INFINITY;
    for (x_shift = 0; x_shift < pool_size[0]; x_shift++) {
      for (y_shift = 0; y_shift < pool_size[1]; y_shift++) {
        size_t pool_cur_dim = inp_pos + x_shift * i12 + y_shift * input_size[2];
        // printf("pool_cur_dim %zu %zu
        // %zu\n",pool_cur_dim,fp->layers[numlayers-1]->dims,numlayers);
        // fflush(stdout);
        pool_map[l] = pool_cur_dim;
        // use the ReLU bounds from the previous layer
        double lb = -fp->layers[numlayers - 1]->neurons[pool_cur_dim]->lb;
        double ub = fp->layers[numlayers - 1]->neurons[pool_cur_dim]->ub;
        if (ub <= 0) {
          inf[l] = 0.0;
          sup[l] = 0.0;
        } else if (lb > 0) {
          inf[l] = lb;
          sup[l] = ub;
        } else {
          inf[l] = 0;
          sup[l] = ub;
        }
        // printf("inf: %g %g\n",inf[l],sup[l]);
        // fflush(stdout);
        sum_u = sum_u + sup[l];
        sum_l = sum_l + inf[l];
        if (sup[l] > max_u) {
          max_u = sup[l];
          max_u_var = pool_map[l];
        }
        if (inf[l] > max_l) {
          max_l = inf[l];
          max_l_var = pool_map[l];
        }
        if ((ub - lb) < min_width) {
          min_width = ub - lb;
          min_width_var = pool_map[l];
        }
        l++;
      }
    }

    bool flag = false;
    size_t var = 0;
    for (j = 0; j < p01; j++) {
      bool is_greater = true;
      for (k = 0; k < p01; k++) {
        if (k == j)
          continue;
        if ((inf[k] == sup[k]) && (inf[j] >= sup[k])) {
          continue;
        } else if ((inf[j] == inf[k]) && (sup[j] == sup[k]) &&
                   (inf[j] == sup[j])) {
          continue;
        } else if (inf[j] <= sup[k]) {
          is_greater = false;
          break;
        }
      }
      if (is_greater) {
        flag = true;
        var = pool_map[j];
        break;
      }
    }
    // printf("max_l: %gmax_u: %g\n",max_l,max_u);
    // fflush(stdout);
    if (flag) {
      // if(0){
      // x_new = x_var
      count++;
      // printf("out_pos: %zu\n",out_pos);
      // fflush(stdout);
      double coeff[1];
      size_t dim[1];
      coeff[0] = 1;
      dim[0] = var;
      out_neurons[out_pos]->lexpr = create_sparse_expr(coeff, 0, dim, 1);
      out_neurons[out_pos]->uexpr = create_sparse_expr(coeff, 0, dim, 1);
      // out_neurons[out_pos]->expr = create_sparse_expr(coeff,0,dim,1);
    } else {
      // max_l	<= x_new <= max_u
      double lcoeff[1];
      size_t ldim[1];
      lcoeff[0] = 1;
      ldim[0] = max_l_var;
      // lcoeff[0] = 0;
      // ldim[0] = 0;
      // printf("max_l: %gmax_u: %g\n",max_l,max_u);
      // fflush(stdout);
      // expr_t * expr = (expr_t *)malloc(sizeof(expr_t));
      // expr->inf_coeff= expr->sup_coeff = expr->dim = NULL;
      // expr->size = 0;
      // expr->inf_cst = -max_l;

      // out_neurons[out_pos]->lexpr =
      // NULL;//create_sparse_expr(NULL,max_l,NULL,0);
      out_neurons[out_pos]->lexpr = create_sparse_expr(lcoeff, 0, ldim, 1);
      // double *ucoeff = (double *)malloc(p01*sizeof(double));
      // size_t *udim = (size_t *)malloc(p01*sizeof(size_t));
      // for(j=0; j < p01; j++){
      //	ucoeff[j] = 1.0;
      //	udim[j] = pool_map[j];
      //}
      double ucoeff[1];
      size_t udim[1];
      ucoeff[0] = 0;
      udim[0] = 0;
      out_neurons[out_pos]->uexpr = create_sparse_expr(ucoeff, max_u, udim, 1);
      // out_neurons[out_pos]->uexpr =
      // create_sparse_expr(ucoeff,max_l-sum_l,udim,p01);
      // sort_sparse_expr(out_neurons[out_pos]->uexpr);
      // free(ucoeff);
      // free(udim);
      // out_neurons[out_pos]->lexpr = create_cst_expr(-max_l,max_l);
      // out_neurons[out_pos]->uexpr = create_cst_expr(-max_u,max_u);
    }
    out_neurons[out_pos]->lb = -max_l;
    out_neurons[out_pos]->ub = max_u;
  }
  // update_state_using_previous_layers_parallel(man,fp,numlayers);
  free(inf);
  free(sup);
  free(pool_map);
  free(output_size);
  // printf("count: %zu\n",count);
  // fflush(stdout);
  // printf("return here2\n");
  // fppoly_fprint(stdout,man,fp,NULL);
  // fflush(stdout);
  return num_out_neurons;
}

void handle_residual_layer(elina_manager_t *man, elina_abstract0_t *element,
                           size_t num_neurons, size_t *predecessors,
                           activation_type_t activation,
                           bool use_area_heuristic) {
  fppoly_t *fp = fppoly_of_abstract0(element);
  size_t numlayers = fp->numlayers;
  fppoly_add_new_layer(fp, num_neurons, RESIDUAL, activation);
  fp->layers[numlayers]->predecessors = predecessors;
  size_t i;
  neuron_t **neurons = fp->layers[numlayers]->neurons;
  for (i = 0; i < num_neurons; i++) {
    double *coeff = (double *)malloc(sizeof(double));
    coeff[0] = 1;
    size_t *dim = (size_t *)malloc(sizeof(size_t));
    dim[0] = i;
    neurons[i]->expr = create_sparse_expr(coeff, 0, dim, 1);
  }
  update_state_using_previous_layers_parallel(man, fp, numlayers,
                                              use_area_heuristic);
}

void handle_residual_relu_layer(elina_manager_t *man,
                                elina_abstract0_t *element, size_t num_neurons,
                                size_t *predecessors, bool use_area_heuristic) {
  handle_residual_layer(man, element, num_neurons, predecessors, RELU,
                        use_area_heuristic);
}

void handle_residual_affine_layer(elina_manager_t *man,
                                  elina_abstract0_t *element,
                                  size_t num_neurons, size_t *predecessors,
                                  bool use_area_heuristic) {
  handle_residual_layer(man, element, num_neurons, predecessors, NONE,
                        use_area_heuristic);
}

void free_neuron(neuron_t *neuron){
  if (neuron->expr) {
    free_expr(neuron->expr);
  }
        if(neuron->lexpr){
		free_expr(neuron->lexpr);
	}
        if (neuron->uexpr) {
          free_expr(neuron->uexpr);
        }
        free(neuron);
}

void free_non_lstm_layer_expr(elina_manager_t *man, elina_abstract0_t *abs, size_t layerno){
    fppoly_t *fp = fppoly_of_abstract0(abs);
    if(layerno >= fp->numlayers){
        fprintf(stdout,"the layer does not exist\n");
        return;
    }
    layer_t * layer = fp->layers[layerno];
    size_t dims = layer->dims;
    size_t i;
    for(i=0; i < dims; i++){
        neuron_t *neuron = layer->neurons[i];
        if (neuron->expr) {
          free_expr(neuron->expr);
        }
        if(neuron->lexpr){
            free_expr(neuron->lexpr);
        }
        if (neuron->uexpr) {
          free_expr(neuron->uexpr);
        }
    }
}

void layer_free(layer_t * layer){
	size_t dims = layer->dims;
	size_t i;
	for(i=0; i < dims; i++){
		free_neuron(layer->neurons[i]);
	}
	free(layer->neurons);
	layer->neurons = NULL;
	if(layer->h_t_inf!=NULL){
		free(layer->h_t_inf);
		layer->h_t_inf = NULL;
	}

	if(layer->h_t_sup!=NULL){
		free(layer->h_t_sup);
		layer->h_t_sup = NULL;
	}
	
	if(layer->c_t_inf!=NULL){
		free(layer->c_t_inf);
		layer->c_t_inf = NULL;
	}

	if(layer->c_t_sup!=NULL){
		free(layer->c_t_sup);
		layer->c_t_sup = NULL;
	}

        free(layer);
	layer = NULL;
}



void fppoly_free(elina_manager_t *man, fppoly_t *fp){
	size_t i;
	size_t output_size = fp->layers[fp->numlayers-1]->dims;
	for(i=0; i < fp->numlayers; i++){
		layer_free(fp->layers[i]);
	}

        for (i = 0; i < output_size; i++) {
          if (fp->out->lexpr[i]) {
            free_expr(fp->out->lexpr[i]);
          }
          if (fp->out->uexpr[i]) {
            free_expr(fp->out->uexpr[i]);
          }
        }

        free(fp->layers);
	fp->layers = NULL;
	free(fp->input_inf);
	fp->input_inf = NULL;
        if(fp->input_lexpr!=NULL && fp->input_uexpr!=NULL){
		for(i=0; i < fp->num_pixels; i++){
			free(fp->input_lexpr[i]);
			free(fp->input_uexpr[i]);
		}
	
		free(fp->input_lexpr);
		fp->input_lexpr = NULL;
		free(fp->input_uexpr);
		fp->input_uexpr = NULL;
        }
	free(fp->input_sup);
	fp->input_sup = NULL;
        free(fp->out->output_inf);
        fp->out->output_inf = NULL;
        free(fp->out->output_sup);
        fp->out->output_sup = NULL;
        free(fp->out->lexpr);
        fp->out->lexpr = NULL;
        free(fp->out->uexpr);
        fp->out->uexpr = NULL;
        free(fp->out);
        fp->out = NULL;
        free(fp);
	fp = NULL;
}






void fppoly_fprint(FILE* stream, elina_manager_t* man, fppoly_t* fp, char** name_of_dim){
	size_t i;
	for(i = 0; i < fp->numlayers; i++){
		fprintf(stream,"layer: %zu\n", i);
		layer_fprint(stream, fp->layers[i], name_of_dim);
	}
	size_t output_size = fp->layers[fp->numlayers-1]->dims;
        if (fp->out != NULL) {
          fprintf(stream, "OUTPUT bounds: \n");
          for (i = 0; i < output_size; i++) {
            fprintf(stream, "%zu: [%g,%g] \n", i, -fp->out->output_inf[i],
                    fp->out->output_sup[i]);
          }
        }
}


elina_interval_t * box_for_neuron(elina_manager_t* man, elina_abstract0_t * abs, size_t layerno, size_t neuron_no){
	fppoly_t *fp = fppoly_of_abstract0(abs);
	if(layerno >= fp->numlayers){
		fprintf(stdout,"the layer does not exist\n");
		return NULL;
	}
	layer_t * layer = fp->layers[layerno];
	size_t dims = layer->dims;
	if(neuron_no >= dims){
		fprintf(stdout,"the neuron does not exist\n");
		return NULL;
	}
	neuron_t * neuron = layer->neurons[neuron_no];
	elina_interval_t * res = elina_interval_alloc();
	elina_interval_set_double(res,-neuron->lb,neuron->ub);
	return res;
}

elina_interval_t ** box_for_layer(elina_manager_t* man, elina_abstract0_t * abs, size_t layerno){
	fppoly_t *fp = fppoly_of_abstract0(abs);
	if(layerno >= fp->numlayers){
		fprintf(stdout,"the layer does not exist\n");
		return NULL;
	}
	layer_t * layer = fp->layers[layerno];
	size_t dims = layer->dims;
	elina_interval_t ** itv_arr = (elina_interval_t **)malloc(dims*sizeof(elina_interval_t *));
	size_t i;
	for(i=0; i< dims; i++){
		itv_arr[i] = box_for_neuron(man, abs, layerno, i);
		
	}
	return itv_arr;
}


size_t get_num_neurons_in_layer(elina_manager_t* man, elina_abstract0_t * abs, size_t layerno){
	fppoly_t *fp = fppoly_of_abstract0(abs);
	if(layerno >= fp->numlayers){
		fprintf(stdout,"the layer does not exist\n");
		return 0;
	}
	layer_t * layer = fp->layers[layerno];
	size_t dims = layer->dims;
	
	return dims;
}

elina_linexpr0_t * get_expr_for_output_neuron(elina_manager_t *man, elina_abstract0_t *abs, size_t i, bool is_lower){
	fppoly_internal_t *pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
	fppoly_t *fp = fppoly_of_abstract0(abs);
	
	size_t output_size = fp->layers[fp->numlayers-1]->dims;
	if(i >= output_size){
		return NULL;
	}
	size_t num_pixels = fp->num_pixels;
	expr_t * expr = NULL;
	if(is_lower){
          expr = fp->out->lexpr[i];
        }
	else{
          expr = fp->out->uexpr[i];
        }
	elina_linexpr0_t * res = NULL;
	size_t j,k;
	if((fp->input_lexpr!=NULL) && (fp->input_uexpr!=NULL)){
		if(is_lower){
			expr =  replace_input_poly_cons_in_lexpr(pr, expr, fp);
		}
		else{
			expr =  replace_input_poly_cons_in_uexpr(pr, expr, fp);
		}
	}
	size_t expr_size = expr->size;
	if(expr->type==SPARSE){
		sort_sparse_expr(expr);
		res = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,expr_size);
	}
	else{
		res = elina_linexpr0_alloc(ELINA_LINEXPR_DENSE,expr_size);
	}
	elina_linexpr0_set_cst_interval_double(res,-expr->inf_cst,expr->sup_cst);
	
	for(j=0;j < expr_size; j++){
		if(expr->type==DENSE){
			k = j;
		}
		else{
			k = expr->dim[j];
		}
		elina_linexpr0_set_coeff_interval_double(res,k,-expr->inf_coeff[j],expr->sup_coeff[j]);
	}
	if((fp->input_lexpr!=NULL) && (fp->input_uexpr!=NULL)){
		free_expr(expr);
	}
	return res;
}

elina_linexpr0_t * get_lexpr_for_output_neuron(elina_manager_t *man,elina_abstract0_t *abs, size_t i){
	return get_expr_for_output_neuron(man,abs,i, true);
}

elina_linexpr0_t * get_uexpr_for_output_neuron(elina_manager_t *man,elina_abstract0_t *abs, size_t i){
	return get_expr_for_output_neuron(man,abs,i, false);
}

void update_bounds_for_neuron(elina_manager_t *man, elina_abstract0_t *abs, size_t layerno, size_t neuron_no, double lb, double ub){
	fppoly_t *fp = fppoly_of_abstract0(abs);
	if(layerno >= fp->numlayers){
		fprintf(stdout,"the layer does not exist\n");
		return;
	}
	layer_t * layer = fp->layers[layerno];
	neuron_t * neuron = layer->neurons[neuron_no];
	neuron->lb = -lb;
	neuron->ub = ub;
}

/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
 *  This software is distributed under GNU Lesser General Public License
 * Version 3.0. For more information, see the ELINA project website at:
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

#include "fppoly_gpu.h"

#include <cuda.h>
#include <iostream>

__device__ __host__ void elina_double_interval_mul(double *a_inf, double *a_sup,
                                                   double b_inf, double b_sup,
                                                   double c_inf, double c_sup) {
  if (c_inf <= 0) {
    /* interval c is positive */
    if (b_inf <= 0) {
      /*interval b is positive*/
      if ((b_inf == 0) || (c_inf == 0)) {
        *a_inf = 0.0;
      } else {
        *a_inf = b_inf * -c_inf;
      }

      if ((b_sup == 0) || (c_sup == 0)) {
        *a_sup = 0.0;
      } else {
        *a_sup = b_sup * c_sup;
      }
    } else if (b_sup <= 0) {
      /* interval b is negative */
      if ((c_sup == 0) || (b_inf == 0)) {
        *a_inf = 0.0;
      } else {
        *a_inf = c_sup * b_inf;
      }

      if ((c_inf == 0) || (b_sup == 0)) {
        *a_sup = 0.0;
      } else {
        *a_sup = -c_inf * b_sup;
      }
    } else {
      /* there is 0 in between for b */
      if ((c_sup == 0) || (b_inf == 0)) {
        *a_inf = 0.0;
      } else {
        *a_inf = b_inf * c_sup;
      }

      if ((c_sup == 0) || (b_sup == 0)) {
        *a_sup = 0.0;
      } else {
        *a_sup = b_sup * c_sup;
      }
    }
  } else if (c_sup <= 0) {
    /* interval c is negative */
    if (b_inf <= 0) {
      /*interval b is positive*/
      if ((b_sup == 0) || (c_inf == 0)) {
        *a_inf = 0.0;
      } else {
        *a_inf = b_sup * c_inf;
      }

      if ((b_inf == 0) || (c_sup == 0)) {
        *a_sup = 0.0;
      } else {
        *a_sup = -b_inf * c_sup;
      }
    } else if (b_sup <= 0) {
      /* interval b is negative */
      if ((b_sup == 0) || (c_sup == 0)) {
        *a_inf = 0.0;
      } else {
        *a_inf = b_sup * -c_sup;
      }

      if ((b_inf == 0) || (c_inf == 0)) {
        *a_sup = 0.0;
      } else {
        *a_sup = b_inf * c_inf;
      }
    } else {
      /* there is 0 in between for b */
      if ((c_inf == 0) || (b_sup == 0)) {
        *a_inf = 0.0;
      } else {
        *a_inf = b_sup * c_inf;
      }

      if ((c_inf == 0) || (b_inf == 0)) {
        *a_sup = 0.0;
      } else {
        *a_sup = b_inf * c_inf;
      }
    }
  } else if (b_inf <= 0) {
    /* interval b is positive */
    if (c_inf <= 0) {
      /*interval c is positive */
      if ((b_inf == 0) || (c_inf == 0)) {
        *a_inf = 0.0;
      } else {
        *a_inf = -b_inf * c_inf;
      }

      if ((b_sup == 0) || (c_sup == 0)) {
        *a_sup = 0.0;
      } else {
        *a_sup = b_sup * c_sup;
      }
    } else if (c_sup <= 0) {
      /* interval c is negative */
      if ((b_sup == 0) || (c_inf == 0)) {
        *a_inf = 0.0;
      } else {
        *a_inf = b_sup * c_inf;
      }

      if ((b_inf == 0) || (c_sup == 0)) {
        *a_sup = 0.0;
      } else {
        *a_sup = -b_inf * c_sup;
      }
    } else {
      /* there is 0 in between for c */
      if ((b_sup == 0) || (c_inf == 0)) {
        *a_inf = 0.0;
      } else {
        *a_inf = b_sup * c_inf;
      }

      if ((b_sup == 0) || (c_sup == 0)) {
        *a_sup = 0.0;
      } else {
        *a_sup = b_sup * c_sup;
      }
    }
  } else if (b_sup <= 0) {
    /* interval b is negative */
    if (c_inf <= 0) {
      /* interval c is positive */
      if ((b_inf == 0) || (c_sup == 0)) {
        *a_inf = 0.0;
      } else {
        *a_inf = b_inf * c_sup;
      }

      if ((b_sup == 0) || (c_inf == 0)) {
        *a_sup = 0.0;
      } else {
        *a_sup = b_sup * -c_inf;
      }
    } else if (c_sup <= 0) {
      /* interval c is negative */
      if ((b_sup == 0) || (c_sup == 0)) {
        *a_inf = 0.0;
      } else {
        *a_inf = -b_sup * c_sup;
      }

      if ((b_inf == 0) || (c_inf == 0)) {
        *a_sup = 0.0;
      } else {
        *a_sup = b_inf * c_inf;
      }
    } else {
      /* there is 0 in between for c */
      if ((b_inf == 0) || (c_sup == 0)) {
        *a_inf = 0.0;
      } else {
        *a_inf = b_inf * c_sup;
      }

      if ((b_inf == 0) || (c_inf == 0)) {
        *a_sup = 0.0;
      } else {
        *a_sup = b_inf * c_inf;
      }
    }
  } else {
    /* there is 0 in between for both b and c */
    double tmp_inf1 = b_sup * c_inf;
    double tmp_sup1 = b_inf * c_inf;
    double tmp_inf2 = b_inf * c_sup;
    double tmp_sup2 = b_sup * c_sup;
    *a_inf = fmax(tmp_inf1, tmp_inf2);
    *a_sup = fmax(tmp_sup1, tmp_sup2);
  }
}

__device__ __host__ void elina_double_interval_div(double *a_inf, double *a_sup,
                                                   double b_inf, double b_sup,
                                                   double c_inf, double c_sup) {
  if (c_inf < 0) {
    /* c is positive */
    if (b_inf <= 0) {
      /* b is positive */
      *a_inf = b_inf / c_sup;
      *a_sup = b_sup / -c_inf;
    } else if (b_sup <= 0) {
      /* b is negative */
      *a_inf = -b_inf / c_inf;
      *a_sup = b_sup / c_sup;
    } else {
      /* 0 is in the middle of b: one divides b by c->inf */
      *a_inf = b_inf / -c_inf;
      *a_sup = b_sup / -c_inf;
    }
  } else if (c_sup < 0) {
    /* c is negative */
    if (b_inf <= 0) {
      /* b is positive */
      *a_sup = b_inf / c_inf;
      *a_inf = -b_sup / c_sup;
    } else if (b_sup <= 0) {
      /* b is negative */
      *a_inf = b_sup / c_inf;
      *a_sup = -b_inf / c_sup;
    } else {
      /* 0 is in the middle of b: one cross-divide b by c->sup */
      *a_inf = b_sup / c_sup;
      *a_sup = b_inf / c_sup;
    }
  } else if ((b_inf == 0) && (b_sup == 0)) {
    /* b is [0,0] */
    *a_inf = b_inf;
    *a_sup = b_sup;
  } else {
    *a_inf = INFINITY;
    *a_sup = INFINITY;
  }
}

fppoly_t *fppoly_of_abstract0(elina_abstract0_t *a) {
  return (fppoly_t *)a->value;
}

elina_abstract0_t *abstract0_of_fppoly(elina_manager_t *man, fppoly_t *fp) {
  elina_abstract0_t *r = (elina_abstract0_t *)malloc(sizeof(elina_abstract0_t));
  assert(r);
  r->value = fp;
  r->man = elina_manager_copy(man);

  return r;
}

static inline void fppoly_internal_free(fppoly_internal_t *pr) {
  if (pr) {
    pr->funid = ELINA_FUNID_UNKNOWN;
    free(pr);
    pr = nullptr;
  }
}

static inline fppoly_internal_t *fppoly_internal_alloc() {
  fppoly_internal_t *pr =
      (fppoly_internal_t *)malloc(sizeof(fppoly_internal_t));
  pr->funid = ELINA_FUNID_UNKNOWN;
  pr->man = nullptr;
  pr->funopt = nullptr;
  pr->min_denormal = ldexpl(1.0, -1074);
  pr->ulp = ldexpl(1.0, -52);

  return pr;
}

/* back pointer to our internal structure from the manager */
fppoly_internal_t *fppoly_init_from_manager(elina_manager_t *man,
                                            elina_funid_t funid) {
  fppoly_internal_t *pr = (fppoly_internal_t *)man->internal;
  pr->funid = funid;

  if (!(pr->man)) {
    pr->man = man;
  }

  return pr;
}

elina_manager_t *fppoly_manager_alloc() {
  std::cout << "This is the GPU version of fppoly!" << std::endl;

  void **funptr;
  // fesetround(FE_UPWARD);
  fppoly_internal_t *pr = fppoly_internal_alloc();

  elina_manager_t *man = elina_manager_alloc(
      "fppoly",                              /* Library name */
      "1.0",                                 /* version */
      pr,                                    /* internal structure */
      (void (*)(void *))fppoly_internal_free /* free function for internal */
  );

  funptr = man->funptr;
  funptr[ELINA_FUNID_FREE] = (void *)&fppoly_free;
  /* 3.Printing */
  funptr[ELINA_FUNID_FPRINT] = (void *)&fppoly_fprint;

  return man;
}

void expr_fprint(FILE *stream, expr_t *expr) {
  if ((expr->inf_coeff == nullptr) || (expr->sup_coeff == nullptr)) {
    fprintf(stdout, "+ [%g, %g]\n", -expr->inf_cst, expr->sup_cst);

    return;
  }

  size_t size = expr->size;

  for (size_t i = 0; i < size; i++) {
    if (i == 0) {
      if (expr->type == DENSE) {
        fprintf(stream, "[%g, %g]x0 ", -expr->inf_coeff[0], expr->sup_coeff[0]);
      } else {
        fprintf(stream, "[%g, %g]x%zu ", -expr->inf_coeff[0],
                expr->sup_coeff[0], expr->dim[0]);
      }
    } else {
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

void expr_print(expr_t *expr) { expr_fprint(stdout, expr); }

__device__ __host__ expr_t *alloc_expr() {
  expr_t *expr = (expr_t *)malloc(sizeof(expr_t));

  expr->inf_coeff = nullptr;
  expr->sup_coeff = nullptr;

  expr->dim = nullptr;

  return expr;
}

__device__ __host__ expr_t *create_dense_expr(double *coeff, double cst,
                                              size_t size) {
  expr_t *expr = (expr_t *)malloc(sizeof(expr_t));

  expr->inf_coeff = (double *)malloc(size * sizeof(double));
  expr->sup_coeff = (double *)malloc(size * sizeof(double));

  expr->dim = nullptr;
  expr->size = size;
  expr->inf_cst = -cst;
  expr->sup_cst = cst;
  expr->type = DENSE;

  for (size_t i = 0; i < size; i++) {
    expr->inf_coeff[i] = -coeff[i];
    expr->sup_coeff[i] = coeff[i];
  }

  return expr;
}

__device__ __host__ expr_t *create_cst_expr(double l, double u) {
  expr_t *expr = (expr_t *)malloc(sizeof(expr_t));

  expr->inf_coeff = nullptr;
  expr->sup_coeff = nullptr;
  expr->dim = nullptr;
  expr->type = SPARSE;
  expr->size = 0;
  expr->inf_cst = l;
  expr->sup_cst = u;

  return expr;
}

__device__ __host__ expr_t *create_sparse_expr(double *coeff, double cst,
                                               size_t *dim, size_t size) {
  expr_t *expr = (expr_t *)malloc(sizeof(expr_t));

  if (size > 0) {
    expr->inf_coeff = (double *)malloc(size * sizeof(double));
    expr->sup_coeff = (double *)malloc(size * sizeof(double));
    expr->dim = (size_t *)malloc(size * sizeof(size_t));
  } else {
    expr->inf_coeff = nullptr;
    expr->sup_coeff = nullptr;
    expr->dim = nullptr;
  }

  expr->size = size;
  expr->inf_cst = -cst;
  expr->sup_cst = cst;
  expr->type = SPARSE;

  for (size_t i = 0; i < size; i++) {
    expr->inf_coeff[i] = -coeff[i];
    expr->sup_coeff[i] = coeff[i];
    expr->dim[i] = dim[i];
  }

  return expr;
}

__device__ __host__ void free_expr(expr_t *expr) {
  if (expr->inf_coeff) {
    free(expr->inf_coeff);
    expr->inf_coeff = nullptr;
  }

  if (expr->sup_coeff) {
    free(expr->sup_coeff);
    expr->sup_coeff = nullptr;
  }

  if ((expr->type == SPARSE) && expr->dim) {
    free(expr->dim);
  }

  expr->dim = nullptr;
  free(expr);
  expr = nullptr;
}

__device__ __host__ expr_t *copy_cst_expr(expr_t *src) {
  expr_t *dst = (expr_t *)malloc(sizeof(expr_t));

  dst->inf_coeff = nullptr;
  dst->sup_coeff = nullptr;
  dst->inf_cst = src->inf_cst;
  dst->sup_cst = src->sup_cst;
  dst->type = src->type;
  dst->dim = nullptr;
  dst->size = src->size;

  return dst;
}

__device__ __host__ expr_t *copy_expr(expr_t *src) {
  expr_t *dst = (expr_t *)malloc(sizeof(expr_t));

  dst->inf_coeff = (double *)malloc(src->size * sizeof(double));
  dst->sup_coeff = (double *)malloc(src->size * sizeof(double));

  dst->inf_cst = src->inf_cst;
  dst->sup_cst = src->sup_cst;
  dst->type = src->type;

  for (size_t i = 0; i < src->size; i++) {
    dst->inf_coeff[i] = src->inf_coeff[i];
    dst->sup_coeff[i] = src->sup_coeff[i];
  }

  if (src->type == SPARSE) {
    dst->dim = (size_t *)malloc(src->size * sizeof(size_t));

    for (size_t i = 0; i < src->size; i++) {
      dst->dim[i] = src->dim[i];
    }
  }

  dst->size = src->size;

  return dst;
}

__device__ __host__ void merge_sparse_expr(expr_t *expr, size_t l, size_t m,
                                           size_t r) {
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

  /* Copy the remaining elements of L[], if there are any */
  while (i < n1) {
    expr->dim[k] = L[i];
    expr->inf_coeff[k] = L2[i];
    expr->sup_coeff[k] = L3[i];
    i++;
    k++;
  }

  /* Copy the remaining elements of R[], if there are any */
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
__device__ __host__ void merge_sort_sparse_expr(expr_t *expr, size_t l,
                                                size_t r) {
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

__device__ __host__ void sort_sparse_expr(expr_t *expr) {
  merge_sort_sparse_expr(expr, 0, expr->size - 1);
}

__device__ __host__ neuron_t *neuron_alloc() {
  neuron_t *res = (neuron_t *)malloc(sizeof(neuron_t));
  res->expr = nullptr;
  res->lb = -INFINITY;
  res->ub = INFINITY;
  res->maxpool_lexpr = nullptr;
  res->maxpool_uexpr = nullptr;

  return res;
}

__device__ __host__ layer_t *create_layer(size_t size, layertype_t type,
                                          activation_type_t activation) {
  layer_t *layer = (layer_t *)malloc(sizeof(layer_t));
  layer->dims = size;
  layer->type = type;
  layer->activation = activation;
  layer->neurons = (neuron_t **)malloc(size * sizeof(neuron_t *));

  for (size_t i = 0; i < size; i++) {
    layer->neurons[i] = neuron_alloc();
  }

  return layer;
}

void fppoly_from_network_input_box(fppoly_t *res, size_t intdim, size_t realdim,
                                   double *inf_array, double *sup_array) {
  res->layers = nullptr;
  res->numlayers = 0;
  size_t num_pixels = intdim + realdim;
  res->input_inf = (double *)malloc(num_pixels * sizeof(double));
  res->input_sup = (double *)malloc(num_pixels * sizeof(double));
  res->input_lexpr = nullptr;
  res->input_uexpr = nullptr;

  for (size_t i = 0; i < num_pixels; i++) {
    res->input_inf[i] = -inf_array[i];
    res->input_sup[i] = sup_array[i];
  }

  res->num_pixels = num_pixels;
  res->out = nullptr;
}

elina_abstract0_t *fppoly_from_network_input(elina_manager_t *man,
                                             size_t intdim, size_t realdim,
                                             double *inf_array,
                                             double *sup_array) {
  fppoly_t *res = (fppoly_t *)malloc(sizeof(fppoly_t));
  fppoly_from_network_input_box(res, intdim, realdim, inf_array, sup_array);

  return abstract0_of_fppoly(man, res);
}

/*
elina_abstract0_t* fppoly_from_network_input_poly(elina_manager_t *man, size_t
intdim, size_t realdim, double *inf_array, double *sup_array, double *
lexpr_weights, double * lexpr_cst, size_t * lexpr_dim, double * uexpr_weights,
                                                  double * uexpr_cst, size_t *
uexpr_dim, size_t expr_size)
{
    fppoly_t * res = (fppoly_t *)malloc(sizeof(fppoly_t));
    fppoly_from_network_input_box(res, intdim, realdim, inf_array, sup_array);
    size_t num_pixels = intdim + realdim;
    res->input_lexpr = (expr_t **)malloc(num_pixels*sizeof(expr_t *));
    res->input_uexpr = (expr_t **)malloc(num_pixels*sizeof(expr_t *));

    size_t i;
        double * tmp_weights = (double*)malloc(expr_size*sizeof(double));
    size_t * tmp_dim = (size_t*)malloc(expr_size*sizeof(size_t));
    for(i = 0; i < num_pixels; i++){

        size_t j;
        for(j=0; j < expr_size; j++){
            tmp_weights[j] = lexpr_weights[i*expr_size+j];
            tmp_dim[j] = lexpr_dim[i*expr_size+j];
        }
        res->input_lexpr[i] = create_sparse_expr(tmp_weights, lexpr_cst[i],
tmp_dim, expr_size); sort_sparse_expr(res->input_lexpr[i]);
    //printf("w: %p %g %g %g cst: %g dim: %p %zu %zu
%zu\n",lexpr_weights[i],lexpr_weights[i][0],lexpr_weights[i][1],
lexpr_weights[i][2],lexpr_cst[i],lexpr_dim[i],lexpr_dim[i][0],lexpr_dim[i][1],
lexpr_dim[i][2]);
        //expr_print(res->input_lexpr[i]);
        //fflush(stdout);
        for(j=0; j < expr_size; j++){
            tmp_weights[j] = uexpr_weights[i*expr_size+j];
            tmp_dim[j] = uexpr_dim[i*expr_size+j];
        }
        res->input_uexpr[i] = create_sparse_expr(tmp_weights, uexpr_cst[i],
tmp_dim, expr_size); sort_sparse_expr(res->input_uexpr[i]);
    //    expr_print(res->input_uexpr[i]);
    //    fflush(stdout);
    }
    free(tmp_weights);
    free(tmp_dim);
    return abstract0_of_fppoly(man,res);

}
*/

void fppoly_alloc_first_layer(fppoly_t *fp, size_t size, size_t num_pixels,
                              layertype_t type, activation_type_t activation) {
  layer_t *layer = create_layer(size, type, activation);
  fp->layers = (layer_t **)malloc(20 * sizeof(layer_t *));
  fp->layers[0] = layer;
  fp->numlayers = 1;
}

void fppoly_add_new_layer(fppoly_t *fp, size_t size, layertype_t type,
                          activation_type_t activation) {
  size_t numlayers = fp->numlayers;
  fp->layers[numlayers] = create_layer(size, type, activation);
  fp->numlayers++;
}

__device__ __host__ void
elina_double_interval_add_expr_coeff(fppoly_internal_t *pr, double *res_inf,
                                     double *res_sup, double inf, double sup,
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

__device__ __host__ void
elina_double_interval_add_cst_coeff(fppoly_internal_t *pr, double *res_inf,
                                    double *res_sup, double inf, double sup,
                                    double inf_expr, double sup_expr) {
  elina_double_interval_add_expr_coeff(pr, res_inf, res_sup, inf, sup, inf_expr,
                                       sup_expr);
  *res_inf += pr->min_denormal;
  *res_sup += pr->min_denormal;
}

__device__ __host__ void
elina_double_interval_mul_expr_coeff(fppoly_internal_t *pr, double *res_inf,
                                     double *res_sup, double inf, double sup,
                                     double inf_expr, double sup_expr) {
  elina_double_interval_mul(res_inf, res_sup, inf, sup, inf_expr, sup_expr);
  double maxA = fmax(fabs(inf_expr), fabs(sup_expr));
  double tmp1, tmp2;
  elina_double_interval_mul(&tmp1, &tmp2, inf, sup, maxA * pr->ulp,
                            maxA * pr->ulp);
  *res_inf += tmp1;
  *res_sup += tmp2;
}

__device__ __host__ void
elina_double_interval_mul_cst_coeff(fppoly_internal_t *pr, double *res_inf,
                                    double *res_sup, double inf, double sup,
                                    double inf_expr, double sup_expr) {
  elina_double_interval_mul_expr_coeff(pr, res_inf, res_sup, inf, sup, inf_expr,
                                       sup_expr);
  *res_inf += pr->min_denormal;
  *res_sup += pr->min_denormal;
}

__device__ __host__ expr_t *multiply_expr(fppoly_internal_t *pr, expr_t *expr,
                                          double mul_inf, double mul_sup) {
  expr_t *res = alloc_expr();

  if (expr->size > 0) {
    res->inf_coeff = (double *)malloc(expr->size * sizeof(double));
    res->sup_coeff = (double *)malloc(expr->size * sizeof(double));
  } else {
    res->inf_coeff = nullptr;
    res->sup_coeff = nullptr;
  }

  res->type = expr->type;

  for (size_t i = 0; i < expr->size; i++) {
    // res->coeff[i] = mul_coeff*expr->coeff[i];
    elina_double_interval_mul_expr_coeff(
        pr, &res->inf_coeff[i], &res->sup_coeff[i], mul_inf, mul_sup,
        expr->inf_coeff[i], expr->sup_coeff[i]);
  }

  if (expr->type == SPARSE) {
    if (expr->size > 0) {
      res->dim = (size_t *)malloc(expr->size * sizeof(size_t));

      for (size_t i = 0; i < expr->size; i++) {
        res->dim[i] = expr->dim[i];
      }
    } else {
      res->dim = nullptr;
    }
  }

  res->size = expr->size;

  elina_double_interval_mul_cst_coeff(pr, &res->inf_cst, &res->sup_cst, mul_inf,
                                      mul_sup, expr->inf_cst, expr->sup_cst);

  // res->cst = mul_coeff*expr->cst;
  return res;
}

__device__ __host__ expr_t *multiply_cst_expr(fppoly_internal_t *pr,
                                              expr_t *expr, double mul_inf,
                                              double mul_sup) {
  expr_t *res = alloc_expr();
  res->inf_coeff = nullptr;
  res->sup_coeff = nullptr;
  res->dim = nullptr;
  res->type = expr->type;
  res->size = expr->size;
  elina_double_interval_mul_cst_coeff(pr, &res->inf_cst, &res->sup_cst, mul_inf,
                                      mul_sup, expr->inf_cst, expr->sup_cst);
  // res->cst = mul_coeff*expr->cst;

  return res;
}

__device__ __host__ void add_cst_expr(fppoly_internal_t *pr, expr_t *exprA,
                                      expr_t *exprB) {
  double maxA = fmax(fabs(exprA->inf_cst), fabs(exprA->sup_cst));
  double maxB = fmax(fabs(exprB->inf_cst), fabs(exprB->sup_cst));
  exprA->inf_cst = exprA->inf_cst + exprB->inf_cst + (maxA + maxB) * pr->ulp +
                   pr->min_denormal;
  exprA->sup_cst = exprA->sup_cst + exprB->sup_cst + (maxA + maxB) * pr->ulp +
                   pr->min_denormal;
}

// A = A + B
__device__ __host__ void add_expr(fppoly_internal_t *pr, expr_t *exprA,
                                  expr_t *exprB) {
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
          if (k < sizeB && (exprB->dim[k] == i)) {
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
          if ((i < sizeA) && (exprA->dim[i] == k)) {
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
        exprA->dim = nullptr;
      } else {
        i = 0;
        k = 0;
        size_t l = 0;
        // printf("before sort\n");
        // expr_print(exprA);
        // fflush(stdout);
        if (exprA->size > 0) {
          sort_sparse_expr(exprA);
        }
        // printf("after sort\n");
        // expr_print(exprA);
        // fflush(stdout);
        sort_sparse_expr(exprB);
        new_inf_coeff = (double *)malloc((sizeA + sizeB) * sizeof(double));
        new_sup_coeff = (double *)malloc((sizeA + sizeB) * sizeof(double));
        size_t *new_dim = (size_t *)malloc((sizeA + sizeB) * sizeof(size_t));

        while ((i < sizeA) && (k < sizeB)) {
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

        double *tmp_inf_coeff = (double *)malloc(l * sizeof(double));
        double *tmp_sup_coeff = (double *)malloc(l * sizeof(double));

        for (size_t n = 0; n < l; n++) {
          tmp_inf_coeff[n] = new_inf_coeff[n];
          tmp_sup_coeff[n] = new_sup_coeff[n];
        }

        free(new_inf_coeff);
        free(new_sup_coeff);

        new_inf_coeff = tmp_inf_coeff;
        new_sup_coeff = tmp_sup_coeff;

        size_t *tmp_dim = (size_t *)malloc(l * sizeof(size_t));

        for (size_t n = 0; n < l; n++) {
          tmp_dim[n] = new_dim[n];
        }

        free(new_dim);
        new_dim = tmp_dim;

        free(exprA->dim);
        exprA->dim = new_dim;
        exprA->size = l;
      }

      if (exprA->inf_coeff) {
        free(exprA->inf_coeff);
        exprA->inf_coeff = nullptr;
      }

      if (exprA->sup_coeff) {
        free(exprA->sup_coeff);
        exprA->sup_coeff = nullptr;
      }

      exprA->inf_coeff = new_inf_coeff;
      exprA->sup_coeff = new_sup_coeff;
    }
  }
}

__device__ __host__ expr_t *
replace_input_poly_cons_in_lexpr(fppoly_internal_t *pr, expr_t *expr,
                                 fppoly_t *fp) {
  size_t dims = expr->size;
  double tmp1, tmp2;
  expr_t *res;

  size_t k;

  if (expr->type == DENSE) {
    k = 0;
  } else {
    k = expr->dim[0];
  }

  expr_t *mul_expr = nullptr;

  if (expr->sup_coeff[0] < 0) {
    mul_expr = fp->input_uexpr[k];
  } else if (expr->inf_coeff[0] < 0) {
    mul_expr = fp->input_lexpr[k];
  }

  if (mul_expr != nullptr) {
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

  for (size_t i = 1; i < dims; i++) {
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }

    expr_t *mul_expr = nullptr;
    expr_t *sum_expr = nullptr;

    if (expr->sup_coeff[i] < 0) {
      mul_expr = fp->input_uexpr[k];
    } else if (expr->inf_coeff[i] < 0) {
      mul_expr = fp->input_lexpr[k];
    }

    if (mul_expr != nullptr) {
      if (mul_expr->size == 0) {
        sum_expr = multiply_cst_expr(pr, mul_expr, expr->inf_coeff[i],
                                     expr->sup_coeff[i]);
        add_cst_expr(pr, res, sum_expr);
      } else if ((expr->inf_coeff[i] != 0) && (expr->sup_coeff[i] != 0)) {
        sum_expr =
            multiply_expr(pr, mul_expr, expr->inf_coeff[i], expr->sup_coeff[i]);
        add_expr(pr, res, sum_expr);
      }
      // free_expr(mul_expr);
      if (sum_expr != nullptr) {
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

__device__ __host__ expr_t *
replace_input_poly_cons_in_uexpr(fppoly_internal_t *pr, expr_t *expr,
                                 fppoly_t *fp) {
  size_t dims = expr->size;
  double tmp1, tmp2;
  expr_t *res;

  size_t k;

  if (expr->type == DENSE) {
    k = 0;
  } else {
    k = expr->dim[0];
  }

  expr_t *mul_expr = nullptr;

  if (expr->sup_coeff[0] < 0) {
    mul_expr = fp->input_lexpr[k];
  } else if (expr->inf_coeff[0] < 0) {
    mul_expr = fp->input_uexpr[k];
  }

  if (mul_expr != nullptr) {
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
  for (size_t i = 1; i < dims; i++) {
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }

    expr_t *mul_expr = nullptr;
    expr_t *sum_expr = nullptr;

    if (expr->sup_coeff[i] < 0) {
      mul_expr = fp->input_lexpr[k];
    } else if (expr->inf_coeff[i] < 0) {
      mul_expr = fp->input_uexpr[k];
    }

    if (mul_expr != nullptr) {
      if (mul_expr->size == 0) {
        sum_expr = multiply_cst_expr(pr, mul_expr, expr->inf_coeff[i],
                                     expr->sup_coeff[i]);
        add_cst_expr(pr, res, sum_expr);
      } else if ((expr->inf_coeff[i] != 0) && (expr->sup_coeff[i] != 0)) {
        sum_expr =
            multiply_expr(pr, mul_expr, expr->inf_coeff[i], expr->sup_coeff[i]);
        add_expr(pr, res, sum_expr);
      }
      // free_expr(mul_expr);
      if (sum_expr != nullptr) {
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

__device__ __host__ double compute_lb_from_expr(fppoly_internal_t *pr,
                                                expr_t *expr, fppoly_t *fp) {
  // printf("start\n");
  // fflush(stdout);
  if ((fp->input_lexpr != nullptr) && (fp->input_uexpr != nullptr)) {
    expr = replace_input_poly_cons_in_lexpr(pr, expr, fp);
  }
  // expr_print(expr);
  // fflush(stdout);
  size_t dims = expr->size;
  double res_inf = expr->inf_cst;

  if ((expr->inf_coeff == nullptr) || (expr->sup_coeff == nullptr)) {
    return 0;
  }

  double tmp1, tmp2;
  size_t k;

  for (size_t i = 0; i < dims; i++) {
    // if(expr->inf_coeff[i]<0){
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }

    elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                              expr->sup_coeff[i], fp->input_inf[k],
                              fp->input_sup[k]);
    // printf("tmp1: %g\n",tmp1);
    res_inf = res_inf + tmp1;
  }
  //    printf("inf: %g\n",-res_inf);
  //    fflush(stdout);
  if ((fp->input_lexpr != nullptr) && (fp->input_uexpr != nullptr)) {
    free_expr(expr);
  }
  // printf("finish\n");
  // fflush(stdout);
  return res_inf;
}

__device__ __host__ double compute_ub_from_expr(fppoly_internal_t *pr,
                                                expr_t *expr, fppoly_t *fp) {
  if ((fp->input_lexpr != nullptr) && (fp->input_uexpr != nullptr)) {
    expr = replace_input_poly_cons_in_uexpr(pr, expr, fp);
  }

  size_t dims = expr->size;
  double res_sup = expr->sup_cst;
  if ((expr->inf_coeff == nullptr) || (expr->sup_coeff == nullptr)) {
    return 0;
  }

  double tmp1, tmp2;
  size_t k;

  for (size_t i = 0; i < dims; i++) {
    // if(expr->inf_coeff[i]<0){
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }

    elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                              expr->sup_coeff[i], fp->input_inf[k],
                              fp->input_sup[k]);
    res_sup = res_sup + tmp2;
  }
  // printf("sup: %g\n",res_sup);
  // fflush(stdout);
  if ((fp->input_lexpr != nullptr) && (fp->input_uexpr != nullptr)) {
    free_expr(expr);
  }

  return res_sup;
}

void ffn_handle_first_layer(elina_manager_t *man, elina_abstract0_t *abs,
                            double **weights, double *bias, size_t size,
                            size_t num_pixels, activation_type_t activation) {
  // printf("start \n");
  // fflush(stdout);
  fppoly_t *res = fppoly_of_abstract0(abs);
  fppoly_alloc_first_layer(res, size, num_pixels, FFN, activation);
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  // for(size_t i = 0; i < num_pixels; i++)
  //{
  //    elina_interval_print(itv[i]);
  //    printf("\n");
  //}
  // fflush(stdout);
  neuron_t **neurons = res->layers[0]->neurons;

  for (size_t i = 0; i < size; i++) {
    neuron_t *neuron = neurons[i];
    double *weight_i = weights[i];
    double bias_i = bias[i];
    neuron->expr = create_dense_expr(weight_i, bias_i, num_pixels);
    neuron->lb = compute_lb_from_expr(pr, neuron->expr, res);
    neuron->ub = compute_ub_from_expr(pr, neuron->expr, res);
  }

  // printf("return here\n");
  // fppoly_fprint(stdout,man,res,nullptr);
  // fflush(stdout);
}

void ffn_handle_first_relu_layer(elina_manager_t *man, elina_abstract0_t *abs,
                                 double **weights, double *bias, size_t size,
                                 size_t num_pixels) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels, RELU);
}

void ffn_handle_first_sigmoid_layer(elina_manager_t *man,
                                    elina_abstract0_t *abs, double **weights,
                                    double *bias, size_t size,
                                    size_t num_pixels) {
  // ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels, SIGMOID);
}

void ffn_handle_first_tanh_layer(elina_manager_t *man, elina_abstract0_t *abs,
                                 double **weights, double *bias, size_t size,
                                 size_t num_pixels) {
  // ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels, TANH);
}

__device__ __host__ expr_t *lexpr_replace_relu_bounds(fppoly_internal_t *pr,
                                                      expr_t *expr,
                                                      neuron_t **neurons) {
  size_t num_neurons = expr->size;
  expr_t *res = alloc_expr();
  res->inf_coeff = (double *)malloc(num_neurons * sizeof(double));
  res->sup_coeff = (double *)malloc(num_neurons * sizeof(double));
  res->inf_cst = expr->inf_cst;
  res->sup_cst = expr->sup_cst;
  res->type = expr->type;
  res->size = num_neurons;

  size_t k;

  for (size_t i = 0; i < num_neurons; i++) {
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

      if ((area1 < area2) && (area1 < area3)) {
        // if(1){
        // res->coeff[i] = lambda*expr->coeff[i];
        elina_double_interval_mul_expr_coeff(
            pr, &res->inf_coeff[i], &res->sup_coeff[i], lambda_inf, lambda_sup,
            expr->inf_coeff[i], expr->sup_coeff[i]);
      } else if ((area2 < area1) && (area2 < area3)) {
        res->inf_coeff[i] = 0.0;
        res->sup_coeff[i] = 0.0;
      } else {
        res->inf_coeff[i] = expr->inf_coeff[i];
        res->sup_coeff[i] = expr->sup_coeff[i];
      }
    } else {
      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;
      double tmp1, tmp2;
      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i], 0, ub);
      res->inf_cst = res->inf_cst + tmp1;
      res->sup_cst = res->sup_cst - tmp1;
    }
  }

  if (expr->type == SPARSE) {
    res->dim = (size_t *)malloc(num_neurons * sizeof(size_t));

    for (size_t i = 0; i < num_neurons; i++) {
      res->dim[i] = expr->dim[i];
    }
  }

  return res;
}

__device__ __host__ expr_t *uexpr_replace_relu_bounds(fppoly_internal_t *pr,
                                                      expr_t *expr,
                                                      neuron_t **neurons) {
  size_t num_neurons = expr->size;
  expr_t *res = alloc_expr();
  res->inf_coeff = (double *)malloc(num_neurons * sizeof(double));
  res->sup_coeff = (double *)malloc(num_neurons * sizeof(double));
  res->inf_cst = expr->inf_cst;
  res->sup_cst = expr->sup_cst;
  res->type = expr->type;
  res->size = num_neurons;

  size_t k;

  for (size_t i = 0; i < num_neurons; i++) {
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

      if ((area1 < area2) && (area1 < area3)) {
        // if(1){
        // res->coeff[i] = lambda*expr->coeff[i];
        elina_double_interval_mul_expr_coeff(
            pr, &res->inf_coeff[i], &res->sup_coeff[i], lambda_inf, lambda_sup,
            expr->inf_coeff[i], expr->sup_coeff[i]);
      } else if ((area2 < area1) && (area2 < area3)) {
        res->inf_coeff[i] = 0.0;
        res->sup_coeff[i] = 0.0;
      } else {
        res->inf_coeff[i] = expr->inf_coeff[i];
        res->sup_coeff[i] = expr->sup_coeff[i];
      }
    } else {
      res->inf_coeff[i] = 0.0;
      res->sup_coeff[i] = 0.0;
      double tmp1, tmp2;
      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i], 0, ub);
      res->inf_cst = res->inf_cst - tmp2;
      res->sup_cst = res->sup_cst + tmp2;
    }
  }

  if (expr->type == SPARSE) {
    res->dim = (size_t *)malloc(num_neurons * sizeof(size_t));

    for (size_t i = 0; i < num_neurons; i++) {
      res->dim[i] = expr->dim[i];
    }
  }

  return res;
}

/*
void compute_chord_slope(double *slope_inf, double *slope_sup, double f_sup_l,
double f_sup_u, double f_inf_l, double f_inf_u, double inf_l, double inf_u,
double sup_l, double sup_u){ double num_l =  f_sup_l + f_inf_u; double num_u =
f_sup_u + f_inf_l;

    double den_l = sup_l + inf_u;
    double den_u = sup_u + inf_l;

        elina_double_interval_div(slope_inf, slope_sup, num_l, num_u, den_l,
den_u);
}
*/

/*
void compute_derivative(double *slope_inf, double *slope_sup, double s_curve_l,
double s_curve_u, double sq_l, double sq_u, bool is_sigmoid){ double
sq_den_sup_l, sq_den_sup_u; elina_double_interval_mul(&sq_den_sup_l,
&sq_den_sup_u, sq_l, sq_u, sq_l, sq_u); if(is_sigmoid){
            elina_double_interval_div(slope_inf, slope_sup, s_curve_l,
s_curve_u, sq_den_sup_l, sq_den_sup_u);
    }
    else{
            *slope_inf = -1 + sq_den_sup_u;
            *slope_sup = 1 + sq_den_sup_l;
    }

}
*/

/*
expr_t * lexpr_replace_s_curve_bounds(fppoly_internal_t * pr, expr_t * expr,
neuron_t ** neurons, bool is_sigmoid){ size_t num_neurons = expr->size; size_t
i,k; expr_t * res = alloc_expr(); res->inf_coeff = (double
*)malloc(num_neurons*sizeof(double)); res->sup_coeff = (double
*)malloc(num_neurons*sizeof(double)); res->inf_cst = expr->inf_cst; res->sup_cst
= expr->sup_cst; res->type = expr->type; res->size = num_neurons;

    for(i = 0; i < num_neurons; i++){
        if(expr->type==DENSE){
            k = i;
        }
        else{
            k = expr->dim[i];
        }
        double lb = neurons[k]->lb;
        double ub = neurons[k]->ub;
        //fesetround(FE_DOWNWARD);
        double e_sup_l = is_sigmoid ? -exp(ub) : -tanh(ub);
        double e_inf_l = is_sigmoid ? -exp(-lb) : -tanh(-lb);

        //fesetround(FE_UPWARD);
        double e_sup_u = is_sigmoid ? exp(ub) : tanh(ub);
        double e_inf_u = is_sigmoid ? exp(-lb) : tanh(-lb);

        double f_sup_l, f_sup_u;
        double f_inf_l, f_inf_u;
        double den_sup_l, den_sup_u;
        double den_inf_l, den_inf_u;


        if(is_sigmoid){
            den_sup_l = -1 + e_sup_l;
            den_sup_u = 1 + e_sup_u;
            den_inf_l = -1 + e_inf_l;
            den_inf_u = 1 + e_inf_u;
            elina_double_interval_div(&f_sup_l, &f_sup_u, e_sup_l, e_sup_u,
den_sup_l, den_sup_u); elina_double_interval_div(&f_inf_l, &f_inf_u, e_inf_l,
e_inf_u, den_inf_l, den_inf_u);
        }
        else{
            f_inf_l = e_inf_l;
            f_inf_u = e_inf_u;
            f_sup_l = e_sup_l;
            f_sup_u = e_sup_u;
            den_inf_l = e_inf_l;
            den_inf_u = e_inf_u;
            den_sup_l = e_sup_l;
            den_sup_u = e_sup_u;
        }

        if((-lb==ub)|| (-f_inf_l==f_sup_u)){

            res->inf_coeff[i] = 0.0;
            res->sup_coeff[i] = 0.0;
            double tmp1, tmp2;
            elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],f_inf_l,f_sup_u);
            res->inf_cst = res->inf_cst + tmp1;
            res->sup_cst = res->sup_cst - tmp1;

        }

        else if(expr->sup_coeff[i]<0 || expr->inf_coeff[i] < 0){
            double slope_inf, slope_sup;
            double intercept_inf, intercept_sup;
            double add_inf, add_sup;
                    double mul_inf, mul_sup;
                    double x_l, x_u;
                    double f_x_l, f_x_u;
                    bool boxify = false;
                    if(expr->sup_coeff[i] < 0){
                        if(ub<0){

                           compute_chord_slope(&slope_inf, &slope_sup, f_sup_l,
f_sup_u, f_inf_l, f_inf_u, lb, -lb, -ub, ub);
                       //if(slope_inf>0){
                    //boxify=true;
                       //}
                       //else{
                    //compute_derivative( &slope_inf, &slope_sup, e_inf_l,
e_inf_u, den_inf_l, den_inf_u, is_sigmoid);
                    //printf("slope2: %.30f
%.30f\n",slope_inf*(ub+lb),slope_sup*(ub+lb));
                       // }    //fflush(stdout);

                    x_l = ub;
                    x_u = -ub;
                    f_x_l = f_sup_l;
                    f_x_u = f_sup_u;
                    //elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,-lb,
lb, slope_inf,slope_sup);
                    //elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_inf_l,
f_inf_u, intercept_inf, intercept_sup);
                    elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,x_l,
x_u, slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_x_l,
f_x_u, intercept_inf, intercept_sup); double tmp1, tmp2, tmp3, tmp4;
                           elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,lb,
-lb, slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&tmp3,&tmp4,add_inf,
add_sup, tmp1, tmp2); if(tmp4<f_inf_u){

                                boxify = true;
                                }

                   //}
                //printf("slope: %.30f %.30f %.30f
%.30f\n",slope_inf*(ub+lb),slope_sup*(ub+lb),- tmp3,tmp4);
                //    fflush(stdout);
                        //elina_double_interval_add_cst_coeff(pr,&tmp3,&tmp4,add_inf,
add_sup, tmp1, tmp2);

                        }
                         else if(lb<=0){

                            compute_derivative( &slope_inf, &slope_sup, e_sup_l,
e_sup_u, den_sup_l, den_sup_u, is_sigmoid); x_l = ub; x_u = -ub; f_x_l =
f_sup_l; f_x_u = f_sup_u;
                            elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,x_l,
x_u, slope_inf,slope_sup);
                            elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_x_l,
f_x_u, intercept_inf, intercept_sup); double tmp1, tmp2, tmp3, tmp4;
                            elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,lb,
-lb, slope_inf,slope_sup);
                            elina_double_interval_add_cst_coeff(pr,&tmp3,&tmp4,add_inf,
add_sup, tmp1, tmp2); if(-tmp3<f_inf_u){

                                boxify = true;
                            }
                      }
                    else{

                            if(lb<=ub){
                                //double slope_inf1, slope_sup1;
                                //double slope_inf2, slope_sup2;
                                compute_derivative( &slope_inf, &slope_sup,
e_sup_l, e_sup_u, den_sup_l, den_sup_u, is_sigmoid);

                            }
                            else{
                                compute_derivative( &slope_inf, &slope_sup,
e_inf_l, e_inf_u, den_inf_l, den_inf_u, is_sigmoid);
                            }
                   // if(slope_inf1>=slope_sup2){
                     //   slope_inf = slope_inf2;
                       // slope_sup = slope_sup2;
                    //}
                    //else if(slope_inf2>=slope_sup1){
                      //  slope_inf = slope_inf1;
                        //slope_sup = slope_sup1;
                    //}
                    //else{
                      //  boxify = true;
                    //}
                    x_l = ub;
                    x_u = -ub;
                    f_x_l = f_sup_l;
                    f_x_u = f_sup_u;
                    elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,x_l,
x_u, slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_x_l,
f_x_u, intercept_inf, intercept_sup); double tmp1, tmp2, tmp3, tmp4;
                    elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,lb, -lb,
slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&tmp3,&tmp4,add_inf,
add_sup, tmp1, tmp2); if(-tmp3<f_inf_u){

                        boxify = true;
                    }
                    }




                }
                else{
                    if(ub < 0){

                    compute_derivative( &slope_inf, &slope_sup, e_inf_l,
e_inf_u, den_inf_l, den_inf_u, is_sigmoid); x_l = -lb; x_u = lb; f_x_l =
f_inf_l; f_x_u = f_inf_u;
                    elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,x_l,
x_u, slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_x_l,
f_x_u, intercept_inf, intercept_sup); double tmp1, tmp2, tmp3, tmp4;
                    elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,-ub, ub,
slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&tmp3,&tmp4,add_inf,
add_sup, tmp1, tmp2); if(tmp4>-f_sup_l){

                        boxify = true;
                    }
                        }
                else if(lb<=0){
                    compute_chord_slope(&slope_inf, &slope_sup, f_sup_l,
f_sup_u, f_inf_l, f_inf_u, lb, -lb, -ub, ub);
                //if(slope_inf>0){
                    //boxify=true;
                //}
               // compute_derivative( &slope_inf, &slope_sup, e_sup_l, e_sup_u,
den_sup_l, den_sup_u, is_sigmoid);
                //}
                //else{
                    x_l = -lb;
                    x_u = lb;
                    f_x_l = f_inf_l;
                    f_x_u = f_inf_u;
                    //elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,x_l,
x_u, slope_inf,slope_sup);
                    elina_double_interval_mul(&intercept_inf,&intercept_sup,x_l,
x_u, slope_inf,slope_sup);
                    //add_inf = f_x_l + intercept_inf;
                    //add_sup = f_x_u + intercept_sup;

                           elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_x_l,
f_x_u, intercept_inf, intercept_sup); double tmp1, tmp2, tmp3, tmp4;
                           elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,-ub,
ub, slope_inf,slope_sup);
                           elina_double_interval_add_cst_coeff(pr,&tmp3,&tmp4,add_inf,
add_sup, tmp1, tmp2);

                   if(-tmp3>f_sup_u){

                            boxify = true;
                           }
                   //}
                }
                else{

                    //double slope_inf1, slope_sup1;
                    //double slope_inf2, slope_sup2;
                    if(lb<=ub){
                    compute_derivative( &slope_inf, &slope_sup, e_sup_l,
e_sup_u, den_sup_l, den_sup_u, is_sigmoid);

                    }
                    else{
                    compute_derivative( &slope_inf, &slope_sup, e_inf_l,
e_inf_u, den_inf_l, den_inf_u, is_sigmoid);
                    }
                    //if(slope_inf1>=slope_sup2){
                      //  slope_inf = slope_inf2;
                        //slope_sup = slope_sup2;
                   // }
                    //else if(slope_inf2>=slope_sup1){
                      //  slope_inf = slope_inf1;
                       // slope_sup = slope_sup1;
                    //}
                    //else{
                      //  boxify = true;
                    //}
                    x_l = -lb;
                    x_u = lb;
                    f_x_l = f_inf_l;
                    f_x_u = f_inf_u;
                    elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,x_l,
x_u, slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_x_l,
f_x_u, intercept_inf, intercept_sup); double tmp1, tmp2, tmp3, tmp4;
                    elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,-ub, ub,
slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&tmp3,&tmp4,add_inf,
add_sup, tmp1, tmp2); if(tmp4>-f_sup_l){

                        boxify = true;
                    }
                }


            }
            if(boxify){
                //printf("boxify lexpr\n");
            //fflush(stdout);
                res->inf_coeff[i] = 0.0;
                res->sup_coeff[i] = 0.0;
                double tmp1, tmp2;
                elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],f_inf_l,f_sup_u);
                res->inf_cst = res->inf_cst + tmp1;
                res->sup_cst = res->sup_cst - tmp1;

            }
            else{
                elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],slope_inf,slope_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
                    elina_double_interval_mul_cst_coeff(pr, &mul_inf, &mul_sup,
add_inf, add_sup, expr->inf_coeff[i], expr->sup_coeff[i] );
                    elina_double_interval_add_cst_coeff(pr,&res->inf_cst,&res->sup_cst,mul_inf,
mul_sup, res->inf_cst, res->sup_cst);
            }

        }


        else{

            res->inf_coeff[i] = 0.0;
            res->sup_coeff[i] = 0.0;
            double tmp1, tmp2;
            elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],f_inf_l,f_sup_u);
            res->inf_cst = res->inf_cst + tmp1;
            res->sup_cst = res->sup_cst - tmp1;
        }
    }
    if(expr->type==SPARSE){
        res->dim = (size_t*)malloc(num_neurons*sizeof(size_t));
        for(i=0; i < num_neurons; i++){
            res->dim[i] = expr->dim[i];
        }
    }
    return res;
}
*/

/*
expr_t * uexpr_replace_s_curve_bounds(fppoly_internal_t *pr, expr_t * expr,
neuron_t ** neurons, bool is_sigmoid){ size_t num_neurons = expr->size; size_t
i,k; expr_t * res = alloc_expr(); res->inf_coeff = (double
*)malloc(num_neurons*sizeof(double)); res->sup_coeff = (double
*)malloc(num_neurons*sizeof(double)); res->inf_cst = expr->inf_cst; res->sup_cst
= expr->sup_cst; res->type = expr->type; res->size = num_neurons;

    for(i = 0; i < num_neurons; i++){
        if(expr->type==DENSE){
            k = i;
        }
        else{
            k = expr->dim[i];
        }
        double lb = neurons[k]->lb;
        double ub = neurons[k]->ub;
        //fesetround(FE_DOWNWARD);
        double e_sup_l = is_sigmoid ? -exp(ub) : -tanh(ub);
        double e_inf_l = is_sigmoid ? -exp(-lb) : -tanh(-lb);

        //fesetround(FE_UPWARD);
        double e_sup_u = is_sigmoid ? exp(ub) : tanh(ub);
        double e_inf_u = is_sigmoid ? exp(-lb) : tanh(-lb);

        double f_sup_l, f_sup_u;
        double f_inf_l, f_inf_u;
        double den_sup_l, den_sup_u;
        double den_inf_l, den_inf_u;

        if(is_sigmoid){
            den_sup_l = -1 + e_sup_l;
            den_sup_u = 1 + e_sup_u;
            den_inf_l = -1 + e_inf_l;
            den_inf_u = 1 + e_inf_u;
            elina_double_interval_div(&f_sup_l, &f_sup_u, e_sup_l, e_sup_u,
den_sup_l, den_sup_u); elina_double_interval_div(&f_inf_l, &f_inf_u, e_inf_l,
e_inf_u, den_inf_l, den_inf_u);
        }
        else{
            f_inf_l = e_inf_l;
            f_inf_u = e_inf_u;
            f_sup_l = e_sup_l;
            f_sup_u = e_sup_u;
            den_inf_l = e_inf_l;
            den_inf_u = e_inf_u;
            den_sup_l = e_sup_l;
            den_sup_u = e_sup_u;
        }

        if((-lb==ub) || (-f_inf_l==f_sup_u)){
            //printf("boxify %.30f %.30f %.30f
%.30f\n",-f_inf_l,f_inf_u,-f_sup_l,f_sup_u); res->inf_coeff[i] = 0.0;
            res->sup_coeff[i] = 0.0;
            double tmp1, tmp2;
            elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],f_inf_l,f_sup_u);
            res->inf_cst = res->inf_cst - tmp2;
            res->sup_cst = res->sup_cst + tmp2;
        }
       // else if(lb+ub<0.000000000001){
    //    printf("boxify %g %g %g %g\n",lb,ub,f_inf_u,f_sup_u);
    //    fflush(stdout);
          //  res->inf_coeff[i] = 0.0;
            //res->sup_coeff[i] = 0.0;
            //double tmp1, tmp2;
            //elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],f_inf_l,f_sup_u);
            //res->inf_cst = res->inf_cst + tmp1;
            //res->sup_cst = res->sup_cst - tmp1;
       // }
        else if(expr->sup_coeff[i]<0 || expr->inf_coeff[i] < 0){
            double slope_inf, slope_sup;
            double intercept_inf, intercept_sup;
            double add_inf, add_sup;
            double mul_inf, mul_sup;
            double x_l, x_u;
            double f_x_l, f_x_u;
            bool boxify = false;
            if(expr->sup_coeff[i] < 0){
                if(ub<0){

                    compute_derivative( &slope_inf, &slope_sup, e_inf_l,
e_inf_u, den_inf_l, den_inf_u, is_sigmoid); x_l = -lb; x_u = lb; f_x_l =
f_inf_l; f_x_u = f_inf_u;
                    //elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,-lb,
lb, slope_inf,slope_sup);
                    //elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_inf_l,
f_inf_u, intercept_inf, intercept_sup);
                    elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,x_l,
x_u, slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_x_l,
f_x_u, intercept_inf, intercept_sup); double tmp1, tmp2, tmp3, tmp4;
                    elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,-ub, ub,
slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&tmp3,&tmp4,add_inf,
add_sup, tmp1, tmp2); if(tmp4>-f_sup_l){ boxify = true;
                    }
                }
                else if(lb<=0){

                    compute_chord_slope(&slope_inf, &slope_sup, f_sup_l,
f_sup_u, f_inf_l, f_inf_u, lb, -lb, -ub, ub);
                //if(slope_inf>0){
                //boxify = true;
                //}
                //else{
                        //compute_derivative( &slope_inf, &slope_sup, e_sup_l,
e_sup_u, den_sup_l, den_sup_u, is_sigmoid);
               //}
                    x_l = -lb;
                    x_u = lb;
                    f_x_l = f_inf_l;
                    f_x_u = f_inf_u;
                    elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,x_l,
x_u, slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_x_l,
f_x_u, intercept_inf, intercept_sup); double tmp1, tmp2, tmp3, tmp4;
                                elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,-ub,
ub, slope_inf,slope_sup);
                                elina_double_interval_add_cst_coeff(pr,&tmp3,&tmp4,add_inf,
add_sup, tmp1, tmp2); if(-tmp3>f_sup_u){ boxify = true;
                                }
                    //}
                }
                else{
                    //double slope_inf1, slope_sup1;
                    //double slope_inf2, slope_sup2;
                    if(lb<=ub){
                    compute_derivative( &slope_inf, &slope_sup, e_sup_l,
e_sup_u, den_sup_l, den_sup_u, is_sigmoid);

                    }
                    else{
                    compute_derivative( &slope_inf, &slope_sup, e_inf_l,
e_inf_u, den_inf_l, den_inf_u, is_sigmoid);
                    }
                    //if(slope_inf1>=slope_sup2){
                      //  slope_inf = slope_inf2;
                        //slope_sup = slope_sup2;
                    //}
                    //else if(slope_inf2>=slope_sup1){
                      //  slope_inf = slope_inf1;
                        //slope_sup = slope_sup1;
                    //}
                    //else{
                      //  boxify = true;
                    //}
                    x_l = -lb;
                    x_u = lb;
                    f_x_l = f_inf_l;
                    f_x_u = f_inf_u;
                    elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,x_l,
x_u, slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_x_l,
f_x_u, intercept_inf, intercept_sup); double tmp1, tmp2, tmp3, tmp4;
                    elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,-ub, ub,
slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&tmp3,&tmp4,add_inf,
add_sup, tmp1, tmp2); if(tmp4>-f_sup_l){ boxify = true;
                    }
                }

            }
            else{
                if(ub < 0){

                    compute_chord_slope(&slope_inf, &slope_sup, f_sup_l,
f_sup_u, f_inf_l, f_inf_u, lb, -lb, -ub, ub);
                    //if(slope_inf>0){
                //boxify = true;
                //}
                //else{
                    //compute_derivative( &slope_inf, &slope_sup, e_inf_l,
e_inf_u, den_inf_l, den_inf_u, is_sigmoid);
                //}
                    x_l = ub;
                    x_u = -ub;
                    f_x_l = f_sup_l;
                    f_x_u = f_sup_u;
                    elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,x_l,
x_u, slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_x_l,
f_x_u, intercept_inf, intercept_sup); double tmp1, tmp2, tmp3, tmp4;
                                elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,lb,
-lb, slope_inf,slope_sup);
                                elina_double_interval_add_cst_coeff(pr,&tmp3,&tmp4,add_inf,
add_sup, tmp1, tmp2);

                                if(tmp4<f_inf_u){
                                boxify = true;
                                }
                   // }


                }
                else if(lb<=0){

                    compute_derivative( &slope_inf, &slope_sup, e_sup_l,
e_sup_u, den_sup_l, den_sup_u, is_sigmoid);

                    x_l = ub;
                    x_u = -ub;
                    f_x_l = f_sup_l;
                    f_x_u = f_sup_u;
                    elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,x_l,
x_u, slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_x_l,
f_x_u, intercept_inf, intercept_sup); double tmp1, tmp2, tmp3, tmp4;
                    elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,lb, -lb,
slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&tmp3,&tmp4,add_inf,
add_sup, tmp1, tmp2);

                    if(-tmp3<f_inf_u){
                        boxify = true;
                    }


                }
                else{
                    //double slope_inf1, slope_sup1;
                    //double slope_inf2, slope_sup2;
                    if(lb<=ub){
                    compute_derivative( &slope_inf, &slope_sup, e_sup_l,
e_sup_u, den_sup_l, den_sup_u, is_sigmoid);

                    }
                    else{
                    compute_derivative( &slope_inf, &slope_sup, e_inf_l,
e_inf_u, den_inf_l, den_inf_u, is_sigmoid);
                    }
                    //if(slope_inf1>=slope_sup2){
                      //  slope_inf = slope_inf2;
                        //slope_sup = slope_sup2;
                    //}
                    //else if(slope_inf2>=slope_sup1){
                      //  slope_inf = slope_inf1;
                        //slope_sup = slope_sup1;
                    //}
                    //else{
                      //  boxify = true;
                    //}

                    x_l = ub;
                    x_u = -ub;
                    f_x_l = f_sup_l;
                    f_x_u = f_sup_u;
                    elina_double_interval_mul_cst_coeff(pr,&intercept_inf,&intercept_sup,x_l,
x_u, slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&add_inf,&add_sup,f_x_l,
f_x_u, intercept_inf, intercept_sup); double tmp1, tmp2, tmp3, tmp4;
                    elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,lb, -lb,
slope_inf,slope_sup);
                    elina_double_interval_add_cst_coeff(pr,&tmp3,&tmp4,add_inf,
add_sup, tmp1, tmp2); if(-tmp3<f_inf_u){ boxify = true;
                    }

                }


            }
       // printf("slope inf: %g slope sup: %g\n",slope_inf,slope_sup);
       // fflush(stdout);
            if(boxify){
                //printf("boxify uexpr\n");
                //fflush(stdout);
                res->inf_coeff[i] = 0.0;
                res->sup_coeff[i] = 0.0;
                double tmp1, tmp2;
                elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],f_inf_l,f_sup_u);
                res->inf_cst = res->inf_cst - tmp2;
                res->sup_cst = res->sup_cst + tmp2;

            }
            else{
            elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],slope_inf,slope_sup,expr->inf_coeff[i],expr->sup_coeff[i]);

            elina_double_interval_mul_cst_coeff(pr, &mul_inf, &mul_sup, add_inf,
add_sup, expr->inf_coeff[i], expr->sup_coeff[i] );
            elina_double_interval_add_cst_coeff(pr,&res->inf_cst,&res->sup_cst,mul_inf,
mul_sup, res->inf_cst, res->sup_cst);
            }
           // printf("cst: %g %g\n",res->inf_cst,res->sup_cst);
        //fflush(stdout);
        }


        else{

            res->inf_coeff[i] = 0.0;
            res->sup_coeff[i] = 0.0;
            double tmp1, tmp2;
            elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],f_inf_l,f_sup_u);
            res->inf_cst = res->inf_cst - tmp2;
            res->sup_cst = res->sup_cst + tmp2;
        }
    }
    if(expr->type==SPARSE){
        res->dim = (size_t*)malloc(num_neurons*sizeof(size_t));
        for(i=0; i < num_neurons; i++){
            res->dim[i] = expr->dim[i];
        }
    }
    return res;
}
*/

/*
expr_t * uexpr_replace_sigmoid_bounds(fppoly_internal_t *pr, expr_t * expr,
neuron_t ** neurons){ return uexpr_replace_s_curve_bounds(pr, expr, neurons,
true);
}

expr_t * uexpr_replace_tanh_bounds(fppoly_internal_t *pr, expr_t * expr,
neuron_t ** neurons){ return uexpr_replace_s_curve_bounds(pr, expr, neurons,
false);
}

expr_t * lexpr_replace_sigmoid_bounds(fppoly_internal_t *pr, expr_t * expr,
neuron_t ** neurons){ return lexpr_replace_s_curve_bounds(pr, expr, neurons,
true);
}

expr_t * lexpr_replace_tanh_bounds(fppoly_internal_t *pr, expr_t * expr,
neuron_t ** neurons){ return lexpr_replace_s_curve_bounds(pr, expr, neurons,
false);
}
*/

__device__ __host__ expr_t *lexpr_replace_maxpool_bounds(fppoly_internal_t *pr,
                                                         expr_t *expr,
                                                         neuron_t **neurons) {
  // printf("begin\n");
  // fflush(stdout);
  size_t num_neurons = expr->size;
  expr_t *res;

  size_t k;

  if (expr->type == DENSE) {
    k = 0;
  } else {
    k = expr->dim[0];
  }

  neuron_t *neuron_k = neurons[k];

  if (expr->sup_coeff[0] < 0) {
    // expr_print(neuron_k->maxpool_uexpr);
    if (neuron_k->maxpool_uexpr == nullptr) {
      res = (expr_t *)malloc(sizeof(expr_t));
      res->inf_coeff = res->sup_coeff = nullptr;
      res->dim = nullptr;
      res->size = 0;
      res->type = SPARSE;
      elina_double_interval_mul_cst_coeff(
          pr, &res->inf_cst, &res->sup_cst, neuron_k->lb, neuron_k->ub,
          expr->inf_coeff[0], expr->sup_coeff[0]);
    } else {
      res = multiply_expr(pr, neuron_k->maxpool_uexpr, expr->inf_coeff[0],
                          expr->sup_coeff[0]);
    }
    // printf("multiply end %zu \n",k);
    // expr_print(res);
    // fflush(stdout);
  } else if (expr->inf_coeff[0] < 0) {
    // expr_print(neuron_k->maxpool_lexpr);
    if (neuron_k->maxpool_lexpr == nullptr) {
      res = (expr_t *)malloc(sizeof(expr_t));
      res->inf_coeff = nullptr;
      res->sup_coeff = nullptr;
      res->dim = nullptr;
      res->size = 0;
      res->type = SPARSE;
      elina_double_interval_mul_cst_coeff(
          pr, &res->inf_cst, &res->sup_cst, neuron_k->lb, neuron_k->ub,
          expr->inf_coeff[0], expr->sup_coeff[0]);
    } else {
      res = multiply_expr(pr, neuron_k->maxpool_lexpr, expr->inf_coeff[0],
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

  for (size_t i = 1; i < num_neurons; i++) {
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }

    neuron_t *neuron_k = neurons[k];

    if (expr->sup_coeff[i] < 0) {
      // expr_print(neuron_k->maxpool_uexpr);
      // printf("add start %zu %zu\n",k,i);

      // expr_print(res);

      if (neuron_k->maxpool_uexpr == nullptr) {
        mul_expr = (expr_t *)malloc(sizeof(expr_t));
        mul_expr->inf_coeff = nullptr;
        mul_expr->sup_coeff = nullptr;
        mul_expr->dim = nullptr;
        mul_expr->size = 0;
        mul_expr->type = SPARSE;
        // printf("lb: %g %g\n");
        elina_double_interval_mul_cst_coeff(
            pr, &mul_expr->inf_cst, &mul_expr->sup_cst, neuron_k->lb,
            neuron_k->ub, expr->inf_coeff[i], expr->sup_coeff[i]);
        res->inf_cst += mul_expr->inf_cst;
        res->sup_cst += mul_expr->sup_cst;
      } else {
        mul_expr = multiply_expr(pr, neuron_k->maxpool_uexpr,
                                 expr->inf_coeff[i], expr->sup_coeff[i]);
        add_expr(pr, res, mul_expr);
      }
      // expr_print(mul_expr);
      // fflush(stdout);
      // printf("add finish\n");
      // expr_print(res);
      // fflush(stdout);
      free_expr(mul_expr);
    } else if (expr->inf_coeff[i] < 0) {
      // expr_print(neuron_k->maxpool_lexpr);
      // printf("add start %zu %zu\n",k,i);

      // expr_print(res);

      if (neuron_k->maxpool_lexpr == nullptr) {
        mul_expr = (expr_t *)malloc(sizeof(expr_t));
        mul_expr->inf_coeff = nullptr;
        mul_expr->sup_coeff = nullptr;
        mul_expr->dim = nullptr;
        mul_expr->size = 0;
        mul_expr->type = SPARSE;
        elina_double_interval_mul_cst_coeff(
            pr, &mul_expr->inf_cst, &mul_expr->sup_cst, neuron_k->lb,
            neuron_k->ub, expr->inf_coeff[i], expr->sup_coeff[i]);
        res->inf_cst += mul_expr->inf_cst;
        res->sup_cst += mul_expr->sup_cst;
      } else {
        mul_expr = multiply_expr(pr, neuron_k->maxpool_lexpr,
                                 expr->inf_coeff[i], expr->sup_coeff[i]);
        // printf("add start1 %zu %zu\n",k,i);
        // expr_print(res);
        // expr_print(mul_expr);
        // fflush(stdout);
        add_expr(pr, res, mul_expr);
      }
      // expr_print(mul_expr);
      //    fflush(stdout);
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

__device__ __host__ expr_t *uexpr_replace_maxpool_bounds(fppoly_internal_t *pr,
                                                         expr_t *expr,
                                                         neuron_t **neurons) {
  size_t num_neurons = expr->size;
  expr_t *res;

  size_t k;

  if (expr->type == DENSE) {
    k = 0;
  } else {
    k = expr->dim[0];
  }

  neuron_t *neuron_k = neurons[k];

  if (expr->sup_coeff[0] < 0) {
    res = multiply_expr(pr, neuron_k->maxpool_lexpr, expr->inf_coeff[0],
                        expr->sup_coeff[0]);
  } else if (expr->inf_coeff[0] < 0) {
    res = multiply_expr(pr, neuron_k->maxpool_uexpr, expr->inf_coeff[0],
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

  for (size_t i = 1; i < num_neurons; i++) {
    if (expr->type == DENSE) {
      k = i;
    } else {
      k = expr->dim[i];
    }

    neuron_t *neuron_k = neurons[k];

    if (expr->sup_coeff[i] < 0) {
      expr_t *mul_expr = multiply_expr(pr, neuron_k->maxpool_lexpr,
                                       expr->inf_coeff[i], expr->sup_coeff[i]);
      add_expr(pr, res, mul_expr);
      free_expr(mul_expr);
    } else if (expr->inf_coeff[i] < 0) {
      expr_t *mul_expr = multiply_expr(pr, neuron_k->maxpool_uexpr,
                                       expr->inf_coeff[i], expr->sup_coeff[i]);
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

__device__ __host__ expr_t *expr_from_previous_layer(fppoly_internal_t *pr,
                                                     expr_t *expr,
                                                     layer_t *prev_layer) {
  if (expr->size == 0) {
    return copy_cst_expr(expr);
  }

  if ((expr->inf_coeff == nullptr) || (expr->sup_coeff == nullptr)) {
    return alloc_expr();
  }

  // printf("coming here %zu\n",expr->size);
  //    fflush(stdout);
  neuron_t **prev_neurons = prev_layer->neurons;
  // size_t out_num_neurons = prev_layer->dims;
  size_t in_num_neurons = expr->size;
  size_t k;
  expr_t *res;

  if (expr->type == DENSE) {
    k = 0;
  } else {
    k = expr->dim[0];
  }

  // printf("start\n");
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

  for (size_t i = 1; i < in_num_neurons; i++) {
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
    } else if ((expr->inf_coeff[i] != 0) || (expr->sup_coeff[i] != 0)) {
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

void update_state_using_previous_layers(elina_manager_t *man, fppoly_t *fp,
                                        size_t layerno) {
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  int k;

  neuron_t **out_neurons = fp->layers[layerno]->neurons;
  size_t num_out_neurons = fp->layers[layerno]->dims;

  for (size_t i = 0; i < num_out_neurons; i++) {
    // printf("i %zu \n",i);
    // fflush(stdout);
    expr_t *lexpr = copy_expr(out_neurons[i]->expr);
    expr_t *uexpr = copy_expr(out_neurons[i]->expr);
    // printf("end\n");
    //    fflush(stdout);
    for (k = layerno - 1; k >= 0; k--) {
      neuron_t **aux_neurons = fp->layers[k]->neurons;
      // size_t dims = fp->layers[k]->dims;
      expr_t *tmp_l;
      expr_t *tmp_u;
      if (fp->layers[k]->type != MAXPOOL) {
        tmp_l = lexpr;
        tmp_u = uexpr;
        // if(i==0 && layerno==1){
        //    printf("before relu %zu\n",lexpr->size);
        //    expr_print(lexpr);
        //    expr_print(uexpr);
        //    fflush(stdout);
        //}
        if (fp->layers[k]->activation == RELU) {
          lexpr = lexpr_replace_relu_bounds(pr, lexpr, aux_neurons);
          // printf("doesnt return\n");

          // fflush(stdout);
          uexpr = uexpr_replace_relu_bounds(pr, uexpr, aux_neurons);
        }
        /*
                else if(fp->layers[k]->activation==SIGMOID){
                        //printf("start\n");
                        //fflush(stdout);
                    lexpr = lexpr_replace_sigmoid_bounds(pr,lexpr,aux_neurons);
                        //printf("finish\n");
                        //fflush(stdout);
                    uexpr = uexpr_replace_sigmoid_bounds(pr,uexpr,aux_neurons);
                }
                else if(fp->layers[k]->activation==TANH){
                    lexpr = lexpr_replace_tanh_bounds(pr,lexpr,aux_neurons);
                    uexpr = uexpr_replace_tanh_bounds(pr,uexpr,aux_neurons);
                }
        */
        free_expr(tmp_l);
        free_expr(tmp_u);

        // if(i==0 && layerno==1){
        //    printf("after relu %zu\n",lexpr->size);
        //    expr_print(lexpr);
        //   expr_print(uexpr);
        //   fflush(stdout);
        // }
        // if(k>=1){
        tmp_l = lexpr;
        tmp_u = uexpr;
        // if(fp->layers[k]->type==MAXPOOL){
        // printf("before replacing %zu %zu\n",lexpr->size,uexpr->size);
        // expr_print(lexpr);
        // expr_print(uexpr);
        // fflush(stdout);
        //}
        lexpr = expr_from_previous_layer(pr, lexpr, fp->layers[k]);
        // printf("doesnt return\n");
        // expr_print(uexpr);
        // fflush(stdout);
        uexpr = expr_from_previous_layer(pr, uexpr, fp->layers[k]);
        // if(fp->layers[k]->type==MAXPOOL){
        // printf("replacing expression from previous layer\n");
        // expr_print(lexpr);
        // expr_print(uexpr);
        ////fflush(stdout);
        //}
        free_expr(tmp_l);
        free_expr(tmp_u);
      } else {
        expr_t *tmp_l = lexpr;
        expr_t *tmp_u = uexpr;
        // printf("before maxpool %zu %zu\n",lexpr->size,i);
        // expr_print(lexpr);
        // expr_print(uexpr);
        // fflush(stdout);
        lexpr = lexpr_replace_maxpool_bounds(pr, lexpr, aux_neurons);

        // fflush(stdout);
        uexpr = uexpr_replace_maxpool_bounds(pr, uexpr, aux_neurons);
        // printf("after maxpool %zu\n",lexpr->size);
        // expr_print(lexpr);
        // fflush(stdout);
        free_expr(tmp_l);

        free_expr(tmp_u);
      }
      //}
      // else{
      //    tmp_l = lexpr;
      //    tmp_u = uexpr;
      //    lexpr = expr_from_previous_layer(lexpr, fp->layers[0]);

      //    uexpr = expr_from_previous_layer(uexpr, fp->layers[0]);

      //    free_expr(tmp_l);
      //    free_expr(tmp_u);
      //}
      // printf("coming here\n");
      // fflush(stdout);
    }
    // expr_print(lexpr);
    // printf("uexpr: %zu %zu\n",uexpr->size,i);
    // expr_print(uexpr);
    // fflush(stdout);

    out_neurons[i]->lb = compute_lb_from_expr(pr, lexpr, fp);
    //- bias_i;
    out_neurons[i]->ub = compute_ub_from_expr(pr, uexpr, fp); //+ bias_i;
    // printf("lb: %g ub: %g\n",out_neurons[i]->lb,out_neurons[i]->ub);
    if (fp->out != nullptr) {
      fp->out->lexpr[i] = lexpr;
      fp->out->uexpr[i] = uexpr;
    } else {
      free_expr(lexpr);
      free_expr(uexpr);
    }
  }
}

void ffn_handle_intermediate_layer(elina_manager_t *man,
                                   elina_abstract0_t *element, double **weights,
                                   double *bias, size_t num_out_neurons,
                                   size_t num_in_neurons,
                                   activation_type_t activation) {
  // printf("ReLU start here %zu %zu\n",num_in_neurons,num_out_neurons);
  // fflush(stdout);
  fppoly_t *fp = fppoly_of_abstract0(element);
  size_t numlayers = fp->numlayers;
  fppoly_add_new_layer(fp, num_out_neurons, FFN, activation);
  neuron_t **out_neurons = fp->layers[numlayers]->neurons;

  for (size_t i = 0; i < num_out_neurons; i++) {
    double *weight_i = weights[i];
    double bias_i = bias[i];
    out_neurons[i]->expr = create_dense_expr(weight_i, bias_i, num_in_neurons);
  }

  update_state_using_previous_layers(man, fp, numlayers);

  // printf("return here2\n");
  // fppoly_fprint(stdout,man,fp,nullptr);
  // fflush(stdout);
}

void ffn_handle_intermediate_relu_layer(elina_manager_t *man,
                                        elina_abstract0_t *element,
                                        double **weights, double *bias,
                                        size_t num_out_neurons,
                                        size_t num_in_neurons) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, RELU);
}

void ffn_handle_intermediate_sigmoid_layer(elina_manager_t *man,
                                           elina_abstract0_t *element,
                                           double **weights, double *bias,
                                           size_t num_out_neurons,
                                           size_t num_in_neurons) {
  // ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
  // num_in_neurons, SIGMOID);
}

void ffn_handle_intermediate_tanh_layer(elina_manager_t *man,
                                        elina_abstract0_t *element,
                                        double **weights, double *bias,
                                        size_t num_out_neurons,
                                        size_t num_in_neurons) {
  // ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
  // num_in_neurons, TANH);
}

__device__ __host__ double
apply_relu_lexpr(fppoly_internal_t *pr, expr_t **lexpr_p, neuron_t *neuron) {
  expr_t *lexpr = *lexpr_p;
  size_t size = lexpr->size;
  double lb = neuron->lb;
  double ub = neuron->ub;
  double width = lb + ub;

  if (ub < 0) {
    free_expr(*lexpr_p);
    *lexpr_p = nullptr;

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

    for (size_t i = 0; i < size; i++) {
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
    *lexpr_p = nullptr;

    return 0;
  }
}

__device__ __host__ double
apply_relu_uexpr(fppoly_internal_t *pr, expr_t **uexpr_p, neuron_t *neuron) {
  expr_t *uexpr = *uexpr_p;
  size_t size = uexpr->size;
  double lb = neuron->lb;
  double ub = neuron->ub;
  double width = lb + ub;

  if (ub < 0) {
    free_expr(*uexpr_p);
    *uexpr_p = nullptr;

    return 0;
  }

  if (lb < 0) {
    return ub;
  }

  double lambda_inf = -ub / width;
  double lambda_sup = ub / width;

  for (size_t i = 0; i < size; i++) {
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
  if (has_relu) {
    for (size_t i = 0; i < size; i++) {
      out->output_inf[i] = apply_relu_lexpr(pr, &out->lexpr[i], neurons[i]);
      out->output_sup[i] = apply_relu_uexpr(pr, &out->uexpr[i], neurons[i]);
    }
  } else {
    for (size_t i = 0; i < size; i++) {
      out->output_inf[i] = neurons[i]->lb;
      out->output_sup[i] = neurons[i]->ub;
    }
  }
}

void ffn_handle_last_layer(elina_manager_t *man, elina_abstract0_t *element,
                           double **weights, double *bias,
                           size_t num_out_neurons, size_t num_in_neurons,
                           bool has_activation, activation_type_t activation) {
  // printf("last\n");
  // fflush(stdout);
  fppoly_t *fp = fppoly_of_abstract0(element);
  size_t numlayers = fp->numlayers;

  if (has_activation) {
    fppoly_add_new_layer(fp, num_out_neurons, FFN, activation);
  } else {
    fppoly_add_new_layer(fp, num_out_neurons, FFN, NONE);
  }

  output_abstract_t *out =
      (output_abstract_t *)malloc(sizeof(output_abstract_t));
  out->output_inf = (double *)malloc(num_out_neurons * sizeof(double));
  out->output_sup = (double *)malloc(num_out_neurons * sizeof(double));
  out->lexpr = (expr_t **)malloc(num_out_neurons * sizeof(expr_t *));
  out->uexpr = (expr_t **)malloc(num_out_neurons * sizeof(expr_t *));
  fp->out = out;
  neuron_t **out_neurons = fp->layers[numlayers]->neurons;
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  for (size_t i = 0; i < num_out_neurons; i++) {
    double *weight_i = weights[i];
    double bias_i = bias[i];
    out_neurons[i]->expr = create_dense_expr(weight_i, bias_i, num_in_neurons);
  }

  update_state_using_previous_layers(man, fp, numlayers);

  if (activation == RELU) {
    handle_final_relu_layer(pr, fp->out, out_neurons, num_out_neurons,
                            has_activation);
  } else {
    for (size_t i = 0; i < num_out_neurons; i++) {
      out->output_inf[i] = out_neurons[i]->lb;
      out->output_sup[i] = out_neurons[i]->ub;
    }
  }
  // fppoly_fprint(stdout,man,fp,nullptr);
}

void ffn_handle_last_relu_layer(elina_manager_t *man,
                                elina_abstract0_t *element, double **weights,
                                double *bias, size_t num_out_neurons,
                                size_t num_in_neurons, bool has_relu) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, has_relu, RELU);
}

void ffn_handle_last_sigmoid_layer(elina_manager_t *man,
                                   elina_abstract0_t *element, double **weights,
                                   double *bias, size_t num_out_neurons,
                                   size_t num_in_neurons, bool has_sigmoid) {
  // ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
  // num_in_neurons, has_sigmoid, SIGMOID);
}

void ffn_handle_last_tanh_layer(elina_manager_t *man,
                                elina_abstract0_t *element, double **weights,
                                double *bias, size_t num_out_neurons,
                                size_t num_in_neurons, bool has_tanh) {
  // ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
  // num_in_neurons, has_tanh, TANH);
}

double get_lb_using_previous_layers(elina_manager_t *man, fppoly_t *fp,
                                    expr_t *expr) {
  size_t numlayers = fp->numlayers;
  expr_t *lexpr = expr;
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  for (int k = numlayers - 1; k >= 0; k--) {
    //    printf("k: %zu\n",k);
    //    fflush(stdout);
    neuron_t **aux_neurons = fp->layers[k]->neurons;
    expr_t *tmp_l;

    if (fp->layers[k]->type != MAXPOOL) {
      if (fp->layers[k]->activation == RELU) {
        tmp_l = lexpr;
        lexpr = lexpr_replace_relu_bounds(pr, lexpr, aux_neurons);
        free_expr(tmp_l);
      }
      /*
          else if(fp->layers[k]->activation==SIGMOID){
              tmp_l = lexpr;
      //printf("start\n");
      //fflush(stdout);
              lexpr = lexpr_replace_sigmoid_bounds(pr,lexpr,aux_neurons);
      //printf("finish\n");
      //fflush(stdout);
              free_expr(tmp_l);
          }
          else if(fp->layers[k]->activation==TANH){
               tmp_l = lexpr;
              lexpr = lexpr_replace_tanh_bounds(pr,lexpr,aux_neurons);
              free_expr(tmp_l);
          }
      */
      tmp_l = lexpr;
      lexpr = expr_from_previous_layer(pr, lexpr, fp->layers[k]);
      free_expr(tmp_l);
    } else {
      expr_t *tmp_l = lexpr;
      lexpr = lexpr_replace_maxpool_bounds(pr, lexpr, aux_neurons);
      free_expr(tmp_l);
    }
  }

  double res = compute_lb_from_expr(pr, lexpr, fp);

  return res;
}

bool is_greater(elina_manager_t *man, elina_abstract0_t *element, elina_dim_t y,
                elina_dim_t x) {
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

    double lb = get_lb_using_previous_layers(man, fp, sub);

    // free_expr(sub);

    if (lb < 0) {
      return true;
    } else {
      return false;
    }
  }
  /*
  output_abstract_t * out = fp->out;
  expr_t * exprA = out->lexpr[y];
  expr_t * exprB = out->uexpr[x];
  if(exprA==nullptr){
      return false;
  }
  else{
      if(exprB==nullptr){
          if(out->output_inf[y]<0){
              return true;
          }
          else{
              return false;
          }

      }
      else{
          //printf("before access %zu %zu\n", exprA->size,exprB->size);
          //fflush(stdout);
          size_t sizeA = exprA->size;
          size_t sizeB = exprB->size;
          //printf("after access\n");
          //fflush(stdout);
          size_t i,k;
          expr_t * sub = (expr_t *)malloc(sizeof(expr_t));
          //
          //sub->size = size;
          sub->inf_cst = exprA->inf_cst + exprB->sup_cst;
          sub->sup_cst = exprA->sup_cst + exprB->inf_cst;
          //printf("getting here\n");
          //expr_print(exprA);
          //expr_print(exprB);
          //fflush(stdout);
          if(exprA->type==DENSE){
              sub->inf_coeff = (double*)malloc(sizeA*sizeof(double));
              sub->sup_coeff = (double*)malloc(sizeA*sizeof(double));
              sub->dim = nullptr;
              sub->size = sizeA;
              sub->type = DENSE;
              if(exprB->type==DENSE){
                      for(i=0; i < sizeA; i++){
                          sub->inf_coeff[i] = exprA->inf_coeff[i] +
  exprB->sup_coeff[i]; sub->sup_coeff[i] = exprA->sup_coeff[i] +
  exprB->inf_coeff[i];
                      }
              }
              else{
                  k = 0;
                  for(i=0; i < sizeA; i++){
                      if(k < sizeB && exprB->dim[k]==i){
                          sub->inf_coeff[i] = exprA->inf_coeff[i] +
  exprB->sup_coeff[k]; sub->sup_coeff[i] = exprA->sup_coeff[i] +
  exprB->inf_coeff[k]; k++;
                      }
                      else{
                          sub->inf_coeff[i] = exprA->inf_coeff[i];
                          sub->sup_coeff[i] = exprA->sup_coeff[i];
                      }
                  }
              }

          }
          else{
              if(exprB->type==DENSE){
                  sub->inf_coeff = (double*)malloc(sizeB*sizeof(double));
                  sub->sup_coeff = (double*)malloc(sizeB*sizeof(double));
                  sub->dim = nullptr;
                  sub->size = sizeB;
                  sub->type = DENSE;
                  i = 0;
                  for(k=0; k < sizeB; k++){
                      if(i < sizeA && exprA->dim[i]==k){
                          sub->inf_coeff[k] = exprA->inf_coeff[i] +
  exprB->sup_coeff[k]; sub->sup_coeff[k] = exprA->sup_coeff[i] +
  exprB->inf_coeff[k]; i++;
                      }
                      else{
                          sub->inf_coeff[i] = exprB->sup_coeff[k];
                          sub->sup_coeff[i] = exprB->inf_coeff[k];
                      }
                  }
              }
              else{
                  sub->inf_coeff =
  (double*)malloc((sizeA+sizeB)*sizeof(double)); sub->sup_coeff =
  (double*)malloc((sizeA+sizeB)*sizeof(double)); sub->dim = nullptr;

                  sub->type = SPARSE;
                  size_t l = 0;
                  i=0;
                  k=0;
                  sub->dim = (size_t *)malloc((sizeA+sizeB)*sizeof(size_t));
                  while(i < sizeA && k < sizeB){
                      if(exprA->dim[i] < exprB->dim[k]){
                          sub->inf_coeff[l] = exprA->inf_coeff[i];
                          sub->sup_coeff[l] = exprA->sup_coeff[i];
                          sub->dim[l] = exprA->dim[i];
                          i++;

                      }
                      else if(exprB->dim[k] < exprA->dim[i]){
                          sub->inf_coeff[l] = exprB->sup_coeff[k];
                          sub->sup_coeff[l] = exprB->inf_coeff[k];
                          sub->dim[l] = exprB->dim[k];
                          k++;
                      }
                      else{
                          sub->inf_coeff[l] = exprA->inf_coeff[i] +
  exprB->sup_coeff[k]; sub->sup_coeff[l] = exprA->sup_coeff[i] +
  exprB->inf_coeff[k]; sub->dim[l] = exprA->dim[i]; i++; k++;
                      }
                      l++;
                  }
                  while(i < sizeA){
                      sub->inf_coeff[l] = exprA->inf_coeff[i];
                      sub->sup_coeff[l] = exprA->sup_coeff[i];
                      sub->dim[l] = exprA->dim[i];
                      i++;
                      l++;
                  }
                  while(k < sizeB){
                      sub->inf_coeff[l] = exprB->inf_coeff[k];
                      sub->sup_coeff[l] = exprB->sup_coeff[k];
                      sub->dim[l] = exprB->dim[k];
                      k++;
                      l++;
                  }
                  sub->size = l;
                  sub->inf_coeff =
  (double*)realloc(sub->inf_coeff,l*sizeof(double)); sub->sup_coeff =
  (double*)realloc(sub->sup_coeff,l*sizeof(double)); sub->dim = (size_t
  *)realloc(sub->dim,l*sizeof(size_t));
              }
          }

          //expr_print(sub);
          //fflush(stdout);
          double lb = compute_lb_from_expr(pr,sub,fp);
          //printf("y: %zu x: %zu lb: %g\n",y,x,lb);
          //fflush(stdout);
          free_expr(sub);
          //double lb = -out->output_inf[y] - out->output_sup[x];
          if(lb<0){
              return true;
          }
          else{
              return false;
          }
      }
  }
*/
}

/*
void conv_handle_first_layer(elina_manager_t *man, elina_abstract0_t *abs,
double *filter_weights, double *filter_bias, size_t *input_size, size_t
*filter_size, size_t num_filters, size_t *strides, bool is_valid_padding, bool
has_bias){

    size_t i;
    size_t num_pixels = input_size[0]*input_size[1]*input_size[2];

    size_t output_size[3];
    if(is_valid_padding){
        output_size[0] = ceil((double)(input_size[0] - filter_size[0] + 1) /
(double)strides[0]); output_size[1] = ceil((double)(input_size[1] -
filter_size[1] + 1) / (double)strides[1]);
    }
    else{
        output_size[0] = ceil((double)input_size[0] / (double)strides[0]);
        output_size[1] = ceil((double)input_size[1] / (double)strides[1]);
    }

    output_size[2] = num_filters;
    size_t size = output_size[0]*output_size[1]*output_size[2];

    fppoly_internal_t * pr =
fppoly_init_from_manager(man,ELINA_FUNID_ASSIGN_LINEXPR_ARRAY); fppoly_t *res =
fppoly_of_abstract0(abs); fppoly_alloc_first_layer(res,size, num_pixels, CONV,
RELU);

    neuron_t ** neurons = res->layers[0]->neurons;
    size_t out_x, out_y, out_z;
        //size_t inp_x, inp_y;
        size_t inp_z;
    size_t x_shift, y_shift;

    long int pad_along_height=0, pad_along_width=0;
    long int pad_top=0,  pad_left=0,  tmp=0;
    if(!is_valid_padding){
        if (input_size[0] % strides[0] == 0){
            long int tmp = filter_size[0] - strides[0];
            pad_along_height = max(tmp, long(0));
        }
        else{
            tmp = filter_size[0] - (input_size[0] % strides[0]);
            pad_along_height = max(tmp, long(0));
        }
        if (input_size[1] % strides[1] == 0){
            tmp = filter_size[1] - strides[1];
            pad_along_width = max(tmp, long(0));
        }
        else{
            tmp = filter_size[1] - (input_size[1] % strides[1]);
            pad_along_width = max(tmp, long(0));
        }
        pad_top = pad_along_height / 2;

        pad_left = pad_along_width / 2;

    }

    for(out_x=0; out_x < output_size[0]; out_x++) {
        for(out_y = 0; out_y < output_size[1]; out_y++) {
         for(out_z=0; out_z < output_size[2]; out_z++) {
             size_t mat_x = out_x*output_size[1]*output_size[2] +
out_y*output_size[2] + out_z; size_t num_coeff =
input_size[2]*filter_size[0]*filter_size[1]; size_t actual_coeff = 0; double
*coeff = (double *)malloc(num_coeff*sizeof(double)); size_t *dim = (size_t
*)malloc(num_coeff*sizeof(double)); i=0; for(inp_z=0; inp_z <input_size[2];
inp_z++) { for(x_shift = 0; x_shift < filter_size[0]; x_shift++) { for(y_shift
=0; y_shift < filter_size[1]; y_shift++) { long int x_val =
out_x*strides[0]+x_shift-pad_top; long int y_val =
out_y*strides[1]+y_shift-pad_left; if(y_val<0 || y_val >= (long
int)input_size[1]){ continue;
                       }

                       if(x_val<0 || x_val >= (long int)input_size[0]){
                             continue;
                       }
                     size_t mat_y = x_val*input_size[1]*input_size[2] +
y_val*input_size[2] + inp_z; if(mat_y>=num_pixels){ continue;
                           }
                     size_t filter_index =
x_shift*filter_size[1]*input_size[2]*output_size[2] +
y_shift*input_size[2]*output_size[2] + inp_z*output_size[2] + out_z; coeff[i] =
filter_weights[filter_index]; dim[i] = mat_y; i++; actual_coeff++;
                 }
            }
            }
            double cst = has_bias? filter_bias[out_z] : 0;
                neurons[mat_x]->expr =
create_sparse_expr(coeff,cst,dim,actual_coeff);
            sort_sparse_expr(neurons[mat_x]->expr);


           neurons[mat_x]->lb = compute_lb_from_expr(pr,
neurons[mat_x]->expr,res); neurons[mat_x]->ub = compute_ub_from_expr(pr,
neurons[mat_x]->expr,res); free(coeff); free(dim);
            }
         }
    }


    //printf("return here\n");
    //fppoly_fprint(stdout,man,res,nullptr);
    //fflush(stdout);
    return;
}


void conv_handle_intermediate_relu_layer(elina_manager_t* man,
elina_abstract0_t* element, double *filter_weights, double * filter_bias, size_t
* input_size, size_t *filter_size, size_t num_filters, size_t *strides, bool
is_valid_padding, bool has_bias){
    //printf("conv intermediate starts here\n");
    //fflush(stdout);
    fppoly_t *fp = fppoly_of_abstract0(element);
    size_t numlayers = fp->numlayers;
    size_t i;
    size_t num_pixels = input_size[0]*input_size[1]*input_size[2];
    size_t output_size[3];
    if(is_valid_padding){
        output_size[0] = ceil((double)(input_size[0] - filter_size[0] + 1) /
(double)strides[0]); output_size[1] = ceil((double)(input_size[1] -
filter_size[1] + 1) / (double)strides[1]);
    }
    else{
        output_size[0] = ceil((double)input_size[0] / (double)strides[0]);
        output_size[1] = ceil((double)input_size[1] / (double)strides[1]);
    }

    output_size[2] = num_filters;
    size_t num_out_neurons = output_size[0]*output_size[1]*output_size[2];
    //printf("num_out_neurons: %zu %zu\n",num_out_neurons,num_pixels);
    //fflush(stdout);
    fppoly_add_new_layer(fp,num_out_neurons, CONV, RELU);
    neuron_t ** out_neurons = fp->layers[numlayers]->neurons;
    size_t out_x, out_y, out_z;
        //size_t inp_x, inp_y;
        size_t inp_z;
    size_t x_shift, y_shift;

    long int pad_along_height=0, pad_along_width=0;
    long int pad_top=0,  pad_left=0,  tmp=0;
    if(!is_valid_padding){
        if (input_size[0] % strides[0] == 0){
            long int tmp = filter_size[0] - strides[0];
            pad_along_height = max(tmp, long(0));
        }
        else{
            tmp = filter_size[0] - (input_size[0] % strides[0]);
            pad_along_height = max(tmp, long(0));
        }
        if (input_size[1] % strides[1] == 0){
            tmp = filter_size[1] - strides[1];
            pad_along_width = max(tmp, long(0));
        }
        else{
            tmp = filter_size[1] - (input_size[1] % strides[1]);
            pad_along_width = max(tmp, long(0));
        }
        pad_top = pad_along_height / 2;

        pad_left = pad_along_width / 2;

    }

    for(out_x=0; out_x < output_size[0]; out_x++) {
        for(out_y = 0; out_y < output_size[1]; out_y++) {
         for(out_z=0; out_z < output_size[2]; out_z++) {
             size_t mat_x = out_x*output_size[1]*output_size[2] +
out_y*output_size[2] + out_z; size_t num_coeff =
input_size[2]*filter_size[0]*filter_size[1]; size_t actual_coeff = 0; double
*coeff = (double *)malloc(num_coeff*sizeof(double)); size_t *dim = (size_t
*)malloc(num_coeff*sizeof(double)); i=0; for(inp_z=0; inp_z <input_size[2];
inp_z++) { for(x_shift = 0; x_shift < filter_size[0]; x_shift++) { for(y_shift
=0; y_shift < filter_size[1]; y_shift++) { long int x_val =
out_x*strides[0]+x_shift-pad_top; long int y_val =
out_y*strides[1]+y_shift-pad_left; if(y_val<0 || y_val >= (long
int)input_size[1]){ continue;
                       }

                       if(x_val<0 || x_val >= (long int)input_size[0]){
                             continue;
                       }
                     size_t mat_y = x_val*input_size[1]*input_size[2] +
y_val*input_size[2] + inp_z; if(mat_y>=num_pixels){ continue;
                           }
                     size_t filter_index =
x_shift*filter_size[1]*input_size[2]*output_size[2] +
y_shift*input_size[2]*output_size[2] + inp_z*output_size[2] + out_z; coeff[i] =
filter_weights[filter_index]; dim[i] = mat_y; actual_coeff++; i++;
                 }
            }
            }
           double cst = has_bias? filter_bias[out_z] : 0;
               out_neurons[mat_x]->expr =
create_sparse_expr(coeff,cst,dim,actual_coeff);
           sort_sparse_expr(out_neurons[mat_x]->expr);
           free(coeff);
           free(dim);
            }
         }
    }
    //printf("gets till here %zu\n",numlayers);
    //fflush(stdout);
    update_state_using_previous_layers(man,fp,numlayers);

    //printf("return here2\n");
    //fppoly_fprint(stdout,man,fp,nullptr);
    //fflush(stdout);
    return;
}
*/

/*
size_t handle_maxpool_layer(elina_manager_t *man, elina_abstract0_t *element,
               size_t *pool_size, size_t *input_size){
    //assert(dimensionality==3);
    //printf("maxpool start\n");
    //fflush(stdout);
    //printf("maxpool start\n");
    //fflush(stdout);
    assert(pool_size[0]==2 && pool_size[1]==2 && pool_size[2]==1);
    //assert(stride[0]==2 && stride[1]==2 && stride[2]==1);

    size_t i,j,k;
    size_t * output_size = (size_t *)malloc(3*sizeof(size_t));
    for(i=0; i < 3; i++){
        output_size[i] = input_size[i]/pool_size[i];
    }


    //size_t num_input_neurons = input_size[0]*input_size[1]*input_size[2];
    size_t num_out_neurons = output_size[0]*output_size[1]*output_size[2];

        size_t o12 = output_size[1]*output_size[2];
       size_t i12 = input_size[1]*input_size[2];
        size_t p01 = pool_size[0]*pool_size[1];

    fppoly_t * fp = fppoly_of_abstract0(element);
    size_t numlayers = fp->numlayers;
    fppoly_add_new_layer(fp,num_out_neurons, MAXPOOL, NONE);
    size_t out_pos;
    double * inf = (double *) calloc(p01,sizeof(double));
    double * sup = (double *) calloc(p01,sizeof(double));
    size_t * pool_map = (size_t *)calloc(p01,sizeof(double));
    neuron_t ** out_neurons = fp->layers[numlayers]->neurons;
    size_t count = 0;
    for(out_pos=0; out_pos<num_out_neurons; out_pos++){
        size_t out_x = out_pos / o12;
        size_t out_y = (out_pos-out_x*o12) / output_size[2];
        size_t out_z = out_pos-out_x*o12 - out_y*output_size[2];
        size_t inp_x = out_x*pool_size[0];
        size_t inp_y = out_y*pool_size[1];
        size_t inp_z = out_z;
        size_t inp_pos = inp_x*i12 + inp_y*input_size[2] + inp_z;
        //size_t pool_start_dim = out_pos*pool_size[0]*pool_size[1];
        //printf("inpXYZ: %zu, %zu, %zu %zu %zu\n", inp_x, inp_y, inp_z,
out_pos, num_out_neurons);
            //printf("outXYZ: %zu, %zu, %zu\n", out_x, out_y, out_z);
        //fflush(stdout);
        size_t x_shift, y_shift, l = 0;
        double sum_u = 0.0;
        double sum_l = 0.0;
        double max_u = -INFINITY;
        double max_l = -INFINITY;

        size_t max_l_var = 0;
        size_t max_u_var = 0;
        size_t min_width_var = 0;
        double min_width = INFINITY;
        for(x_shift = 0; x_shift < pool_size[0]; x_shift++){
            for(y_shift = 0; y_shift < pool_size[1]; y_shift++){
                size_t pool_cur_dim = inp_pos + x_shift*i12 +
y_shift*input_size[2];
                //printf("pool_cur_dim %zu %zu
%zu\n",pool_cur_dim,fp->layers[numlayers - 1]->dims,numlayers);
                //fflush(stdout);
                pool_map[l] = pool_cur_dim;
                // use the ReLU bounds from the previous layer
                double lb = -fp->layers[numlayers -
1]->neurons[pool_cur_dim]->lb; double ub = fp->layers[numlayers -
1]->neurons[pool_cur_dim]->ub; if(ub<=0){ inf[l] = 0.0; sup[l] = 0.0;
                }
                else if(lb>0){
                    inf[l] = lb;
                    sup[l] = ub;
                }
                else{
                    inf[l] = 0;
                    sup[l] = ub;
                }
                //printf("inf: %g %g\n",inf[l],sup[l]);
                //fflush(stdout);
                sum_u = sum_u + sup[l];
                sum_l = sum_l + inf[l];
                if(sup[l]>max_u){
                    max_u = sup[l];
                    max_u_var = pool_map[l];
                }
                if(inf[l] > max_l){
                    max_l = inf[l];
                    max_l_var = pool_map[l];
                }
                if((ub-lb) < min_width){
                    min_width = ub - lb;
                    min_width_var = pool_map[l];
                }
                l++;
            }
        }

        bool flag = false;
        size_t var = 0;
        for(j=0; j < p01; j++){
            bool is_greater = true;
            for(k = 0;  k < p01; k++){
                if(k==j)continue;
                if((inf[k]==sup[k]) && (inf[j]>=sup[k])){
                    continue;
                }
                else if((inf[j]==inf[k]) && (sup[j]==sup[k]) &&
(inf[j]==sup[j])){ continue;
                }
                else if(inf[j]<=sup[k]){
                    is_greater = false;
                    break;
                }
            }
            if(is_greater){
                flag = true;
                var =pool_map[j];
                break;
            }
        }
        //printf("max_l: %gmax_u: %g\n",max_l,max_u);
        //fflush(stdout);
        if(flag){
        //if(0){
            //x_new = x_var
            count++;
            //printf("out_pos: %zu\n",out_pos);
            //fflush(stdout);
            double coeff[1];
            size_t dim[1];
            coeff[0] = 1;
            dim[0] = var;
            out_neurons[out_pos]->maxpool_lexpr =
create_sparse_expr(coeff,0,dim,1); out_neurons[out_pos]->maxpool_uexpr =
create_sparse_expr(coeff,0,dim,1);
            //out_neurons[out_pos]->expr = create_sparse_expr(coeff,0,dim,1);
        }
        else{
            //max_l    <= x_new <= max_u
            double lcoeff[1];
            size_t ldim[1];
            lcoeff[0] = 1;
            ldim[0] = max_l_var;
            //lcoeff[0] = 0;
            //ldim[0] = 0;
            //printf("max_l: %gmax_u: %g\n",max_l,max_u);
            //fflush(stdout);
            //expr_t * expr = (expr_t *)malloc(sizeof(expr_t));
            //expr->inf_coeff= expr->sup_coeff = expr->dim = nullptr;
            //expr->size = 0;
            //expr->inf_cst = -max_l;

            //out_neurons[out_pos]->maxpool_lexpr =
nullptr;//create_sparse_expr(nullptr,max_l,nullptr,0);
            out_neurons[out_pos]->maxpool_lexpr =
create_sparse_expr(lcoeff,0,ldim,1);
            //double *ucoeff = (double *)malloc(p01*sizeof(double));
            //size_t *udim = (size_t *)malloc(p01*sizeof(size_t));
            //for(j=0; j < p01; j++){
            //    ucoeff[j] = 1.0;
            //    udim[j] = pool_map[j];
            //}
            double ucoeff[1];
            size_t udim[1];
            ucoeff[0] = 0;
            udim[0] = 0;
            out_neurons[out_pos]->maxpool_uexpr =
create_sparse_expr(ucoeff,max_u,udim,1);
            //out_neurons[out_pos]->maxpool_uexpr =
create_sparse_expr(ucoeff,max_l-sum_l,udim,p01);
            //sort_sparse_expr(out_neurons[out_pos]->maxpool_uexpr);
            //free(ucoeff);
            //free(udim);
            //out_neurons[out_pos]->maxpool_lexpr =
create_cst_expr(-max_l,max_l);
            //out_neurons[out_pos]->maxpool_uexpr =
create_cst_expr(-max_u,max_u);
        }
        out_neurons[out_pos]->lb = -max_l;
        out_neurons[out_pos]->ub = max_u;
    }
    //update_state_using_previous_layers(man,fp,numlayers);
        free(inf);
    free(sup);
    free(pool_map);
    free(output_size);
        //printf("count: %zu\n",count);
    //fflush(stdout);
    //printf("return here2\n");
    //fppoly_fprint(stdout,man,fp,nullptr);
    //fflush(stdout);
    return num_out_neurons;
}
*/

__device__ __host__ void free_neuron(neuron_t *neuron) {
  if (neuron->expr) {
    free_expr(neuron->expr);
  }

  if (neuron->maxpool_lexpr) {
    free_expr(neuron->maxpool_lexpr);
  }

  if (neuron->maxpool_uexpr) {
    free_expr(neuron->maxpool_uexpr);
  }

  free(neuron);
}

__device__ __host__ void layer_free(layer_t *layer) {
  for (size_t i = 0; i < layer->dims; i++) {
    free_neuron(layer->neurons[i]);
  }

  free(layer->neurons);
  layer->neurons = nullptr;

  free(layer);
  layer = nullptr;
}

void fppoly_free(elina_manager_t *man, fppoly_t *fp) {
  size_t output_size = fp->layers[fp->numlayers - 1]->dims;

  for (size_t i = 0; i < fp->numlayers; i++) {
    layer_free(fp->layers[i]);
  }

  for (size_t i = 0; i < output_size; i++) {
    if (fp->out->lexpr[i]) {
      free_expr(fp->out->lexpr[i]);
    }

    if (fp->out->uexpr[i]) {
      free_expr(fp->out->uexpr[i]);
    }
  }

  free(fp->layers);
  fp->layers = nullptr;
  free(fp->input_inf);
  fp->input_inf = nullptr;

  if ((fp->input_lexpr != nullptr) && (fp->input_uexpr != nullptr)) {
    for (size_t i = 0; i < fp->num_pixels; i++) {
      free(fp->input_lexpr[i]);
      free(fp->input_uexpr[i]);
    }

    free(fp->input_lexpr);
    fp->input_lexpr = nullptr;
    free(fp->input_uexpr);
    fp->input_uexpr = nullptr;
  }

  free(fp->input_sup);
  fp->input_sup = nullptr;
  free(fp->out->output_inf);
  fp->out->output_inf = nullptr;
  free(fp->out->output_sup);
  fp->out->output_sup = nullptr;
  free(fp->out->lexpr);
  fp->out->lexpr = nullptr;
  free(fp->out->uexpr);
  fp->out->uexpr = nullptr;
  free(fp->out);
  fp->out = nullptr;
  free(fp);
  fp = nullptr;
}

void neuron_fprint(FILE *stream, neuron_t *neuron, char **name_of_dim) {
  // expr_fprint(stream,neuron->expr);
  fprintf(stream, "[%g, %g]\n", -neuron->lb, neuron->ub);
}

void layer_fprint(FILE *stream, layer_t *layer, char **name_of_dim) {
  size_t dims = layer->dims;

  for (size_t i = 0; i < dims; i++) {
    fprintf(stream, "neuron: %zu ", i);
    neuron_fprint(stream, layer->neurons[i], name_of_dim);
  }
}

void fppoly_fprint(FILE *stream, elina_manager_t *man, fppoly_t *fp,
                   char **name_of_dim) {
  for (size_t i = 0; i < fp->numlayers; i++) {
    fprintf(stream, "layer: %zu\n", i);
    layer_fprint(stream, fp->layers[i], name_of_dim);
  }

  size_t output_size = fp->layers[fp->numlayers - 1]->dims;

  if (fp->out != nullptr) {
    fprintf(stream, "OUTPUT bounds: \n");

    for (size_t i = 0; i < output_size; i++) {
      fprintf(stream, "%zu: [%g,%g] \n", i, -fp->out->output_inf[i],
              fp->out->output_sup[i]);
    }
  }
}

/*
elina_interval_t * box_for_neuron(elina_abstract0_t * abs, size_t layerno,
size_t neuron_no){ fppoly_t *fp = fppoly_of_abstract0(abs); if(layerno >=
fp->numlayers){ fprintf(stdout,"the layer does not exist\n"); return nullptr;
    }
    layer_t * layer = fp->layers[layerno];
    size_t dims = layer->dims;
    if(neuron_no >= dims){
        fprintf(stdout,"the neuron does not exist\n");
        return nullptr;
    }
    neuron_t * neuron = layer->neurons[neuron_no];
    elina_interval_t * res = elina_interval_alloc();
    elina_interval_set_double(res,-neuron->lb,neuron->ub);
    return res;
}

elina_interval_t ** box_for_layer(elina_abstract0_t * abs, size_t layerno){
    fppoly_t *fp = fppoly_of_abstract0(abs);
    if(layerno >= fp->numlayers){
        fprintf(stdout,"the layer does not exist\n");
        return nullptr;
    }
    layer_t * layer = fp->layers[layerno];
    size_t dims = layer->dims;
    elina_interval_t ** itv_arr = (elina_interval_t
**)malloc(dims*sizeof(elina_interval_t *)); size_t i; for(i=0; i< dims; i++){
        itv_arr[i] = box_for_neuron(abs, layerno, i);
    }
    return itv_arr;
}
*/

elina_linexpr0_t *get_expr_for_output_neuron(elina_manager_t *man,
                                             elina_abstract0_t *abs, size_t i,
                                             bool is_lower) {
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  fppoly_t *fp = fppoly_of_abstract0(abs);

  size_t output_size = fp->layers[fp->numlayers - 1]->dims;

  if (i >= output_size) {
    return nullptr;
  }

  expr_t *expr = nullptr;

  if (is_lower) {
    expr = fp->out->lexpr[i];
  } else {
    expr = fp->out->uexpr[i];
  }

  elina_linexpr0_t *res = nullptr;

  if ((fp->input_lexpr != nullptr) && (fp->input_uexpr != nullptr)) {
    if (is_lower) {
      expr = replace_input_poly_cons_in_lexpr(pr, expr, fp);
    } else {
      expr = replace_input_poly_cons_in_uexpr(pr, expr, fp);
    }
  }

  size_t expr_size = expr->size;

  if (expr->type == SPARSE) {
    sort_sparse_expr(expr);
    res = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, expr_size);
  } else {
    res = elina_linexpr0_alloc(ELINA_LINEXPR_DENSE, expr_size);
  }

  elina_linexpr0_set_cst_interval_double(res, -expr->inf_cst, expr->sup_cst);

  size_t k;

  for (size_t j = 0; j < expr_size; j++) {
    if (expr->type == DENSE) {
      k = j;
    } else {
      k = expr->dim[j];
    }

    elina_linexpr0_set_coeff_interval_double(res, k, -expr->inf_coeff[j],
                                             expr->sup_coeff[j]);
  }

  if ((fp->input_lexpr != nullptr) && (fp->input_uexpr != nullptr)) {
    free_expr(expr);
  }

  return res;
}

elina_linexpr0_t *get_lexpr_for_output_neuron(elina_manager_t *man,
                                              elina_abstract0_t *abs,
                                              size_t i) {
  return get_expr_for_output_neuron(man, abs, i, true);
}

elina_linexpr0_t *get_uexpr_for_output_neuron(elina_manager_t *man,
                                              elina_abstract0_t *abs,
                                              size_t i) {
  return get_expr_for_output_neuron(man, abs, i, false);
}

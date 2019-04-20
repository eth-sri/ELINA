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

#include <chrono>
#include <cuda.h>
#include <iostream>

bool results[90];
bool results_calculated;
size_t output_counter;

#define gpuErrchk(ans)                                                         \
  { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line,
                      bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
            line);
    if (abort)
      exit(code);
  }
}

__constant__ const double min_denormal = 4.940656458412465441766e-324;
__constant__ const double ulp = 2.220446049250313080848e-16;

bool initialized = false;

__device__ void
elina_double_interval_mul(double *const a_inf, double *const a_sup,
                          const double b_inf, const double b_sup,
                          const double c_inf, const double c_sup) {
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

__device__ void
elina_double_interval_div(double *const a_inf, double *const a_sup,
                          const double b_inf, const double b_sup,
                          const double c_inf, const double c_sup) {
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
  results_calculated = false;
  output_counter = 1;
  size_t set_limit = 3 * size_t(2u << 30u);
  size_t limit;

  if (!initialized) {
    // SET THIS ONLY ONCE BEFORE FIRST KERNEL CALL, ELSE CUDA ERROR (INVALID
    // ARGUMENT ERROR)
    cudaDeviceSetLimit(cudaLimitMallocHeapSize, set_limit);
    initialized = true;
  }

  cudaDeviceGetLimit(&limit, cudaLimitMallocHeapSize);

  std::cout << "The Heap Limit is " << limit << std::endl;

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

__device__ void expr_print(const expr_t *const expr) {
  if ((expr->inf_coeff == nullptr) || (expr->sup_coeff == nullptr)) {
    printf("+ [%g, %g]\n", -expr->inf_cst, expr->sup_cst);

    return;
  }

  size_t size = expr->size;

  for (size_t i = 0; i < size; i++) {
    if (i == 0) {
      if (expr->type == DENSE) {
        printf("[%g, %g]x0 ", -expr->inf_coeff[0], expr->sup_coeff[0]);
      } else {
        printf("[%g, %g]x%zu ", -expr->inf_coeff[0], expr->sup_coeff[0],
               expr->dim[0]);
      }
    } else {
      if (expr->type == DENSE) {
        printf("+ [%g, %g]x%zu ", -expr->inf_coeff[i], expr->sup_coeff[i], i);
      } else {
        printf("+ [%g, %g]x%zu ", -expr->inf_coeff[i], expr->sup_coeff[i],
               expr->dim[i]);
      }
    }
  }

  printf("+ [%g, %g]\n", -expr->inf_cst, expr->sup_cst);
}

__device__ expr_t *alloc_expr() {
  expr_t *expr = (expr_t *)malloc(sizeof(expr_t));

  expr->inf_coeff = nullptr;
  expr->sup_coeff = nullptr;

  expr->dim = nullptr;

  return expr;
}

__device__ expr_t *alloc_dense_expr_and_arrays(const size_t size) {
  expr_t *expr = (expr_t *)malloc(sizeof(expr_t));

  expr->inf_coeff = (double *)malloc(size * sizeof(double));
  expr->sup_coeff = (double *)malloc(size * sizeof(double));

  expr->size = size;
  expr->type = DENSE;

  expr->inf_cst = 0;
  expr->sup_cst = 0;

  expr->dim = nullptr;

  return expr;
}

__device__ __host__ expr_t *
create_dense_expr(const double *coeff, const double cst, const size_t size) {
  expr_t *expr = (expr_t *)malloc(sizeof(expr_t));

  expr->inf_coeff = (double *)malloc(size * sizeof(double));
  expr->sup_coeff = (double *)malloc(size * sizeof(double));

  expr->size = size;
  expr->type = DENSE;

  expr->inf_cst = -cst;
  expr->sup_cst = cst;

  for (size_t i = 0; i < size; i++) {
    expr->inf_coeff[i] = -coeff[i];
    expr->sup_coeff[i] = coeff[i];
  }

  expr->dim = nullptr;

  return expr;
}

__device__ expr_t *create_cst_expr(const double l, const double u) {
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

__device__ __host__ expr_t *create_sparse_expr(const double *coeff,
                                               const double cst,
                                               const size_t *dim,
                                               const size_t size) {
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

__device__ expr_t *copy_cst_expr(const expr_t *const src) {
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

__device__ expr_t *copy_expr(const expr_t *const src) {
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

// merges an array of size_t and mirrors the same merge operation to two arrays
// of double
__device__ __host__ void merge(size_t *source, size_t *target,
                               double *source_mirror, double *target_mirror,
                               size_t llength, size_t ulength) {
  size_t length = llength + ulength;

  size_t *lsource_ptr = source;
  size_t *usource_ptr = source + llength;
  size_t *target_ptr = target;

  double *lsource_mirror_ptr = source_mirror;
  double *usource_mirror_ptr = source_mirror + llength;
  double *target_mirror_ptr = target_mirror;

  while (lsource_ptr < (source + llength) && usource_ptr < (source + length)) {
    if (*lsource_ptr < *usource_ptr) {
      *target_ptr = *lsource_ptr;
      lsource_ptr++;
      target_ptr++;

      *target_mirror_ptr = *lsource_mirror_ptr;
      lsource_mirror_ptr++;
      target_mirror_ptr++;
    } else {
      *target_ptr = *usource_ptr;

      usource_ptr++;
      target_ptr++;

      *target_mirror_ptr = *usource_mirror_ptr;

      usource_mirror_ptr++;
      target_mirror_ptr++;
    }
  }

  while (lsource_ptr < (source + llength)) {
    *target_ptr = *lsource_ptr;

    lsource_ptr++;
    target_ptr++;

    *target_mirror_ptr = *lsource_mirror_ptr;

    lsource_mirror_ptr++;
    target_mirror_ptr++;
  }

  while (usource_ptr < (source + length)) {
    *target_ptr = *usource_ptr;

    usource_ptr++;
    target_ptr++;

    *target_mirror_ptr = *usource_mirror_ptr;

    usource_mirror_ptr++;
    target_mirror_ptr++;
  }
}

template <typename T> __device__ __host__ void swap(T &a, T &b) {
  T tmp = b;
  b = a;
  a = tmp;
}

// sorts an array of size_t and mirrors the same permutation to two arrays of
// double
__device__ __host__ void sort_sparse_expr(size_t *array, double *mirror_array,
                                          const size_t length) {
  size_t *a = array;
  size_t *b = nullptr;
  double *mirror_a = mirror_array;
  double *mirror_b = nullptr;

  bool sorted = true;

  size_t startl = 0;
  size_t startu;
  size_t endu;

  bool first = true;
  size_t i = 0;

  // check if the array is sorted before going in the main loop
  for (; i < length - 1; i++) {
    if (a[i] > a[i + 1]) {
      startu = i + 1;
      first = false;
      sorted = false;

      // only allocate memory if the array is not sorted
      b = (size_t *)malloc(length * sizeof(size_t));
      mirror_b = (double *)malloc(length * sizeof(double));

      i++;

      break;
    }
  }

  if (sorted) {
    return;
  }

  while (true) {
    for (; i < length - 1; i++) {
      if (a[i] > a[i + 1]) {
        if (first) {
          // mark end of first chunk when encountering a descending number
          startu = i + 1;
          first = false;
          sorted = false;
        } else {
          // mark end of second chunk when encountering a descending number and
          // merge the two chunks
          endu = i + 1;
          merge(a + startl, b + startl, mirror_a + startl, mirror_b + startl,
                startu - startl, endu - startu);
          first = true;
          startl = endu;
        }
      }
    }

    if (sorted == true) {
      break;
    }

    if (first) {
      // for odd numbers of chunks, need to copy the last chunk over
      memcpy(b + startl, a + startl, (length - startl) * sizeof(size_t));
      memcpy(mirror_b + startl, mirror_a + startl,
             (length - startl) * sizeof(double));
    } else {
      // for even numbers of chunks, need to mark past-the-end pointer as endu
      // and do a merge
      endu = length;
      merge(a + startl, b + startl, mirror_a + startl, mirror_b + startl,
            startu - startl, endu - startu);
      first = true;
    }

    swap(a, b);
    swap(mirror_a, mirror_b);

    i = 0;
    startl = 0;
    sorted = true;
  }

  if (array == a) {
    // if a points to the original array, just free the temporaries
    free(b);
    free(mirror_b);
  } else {
    // if a points to the temporaries created in merge_sort, swap, then memcpy
    // from b to a, then free the temporaries
    swap(a, b);
    swap(mirror_a, mirror_b);

    memcpy(a, b, length * sizeof(size_t));
    memcpy(mirror_a, mirror_b, length * sizeof(double));

    free(b);
    free(mirror_b);
  }
}

layer_t *create_layer(const size_t num_out_neurons, const size_t num_in_neurons,
                      const layertype_t type,
                      const activation_type_t activation) {
  layer_t *layer = (layer_t *)malloc(sizeof(layer_t));

  layer->num_out_neurons = num_out_neurons;
  layer->num_in_neurons = num_in_neurons;

  layer->type = type;
  layer->activation = activation;

  cudaMalloc((void **)&layer->lb_array, num_out_neurons * sizeof(double));
  cudaMalloc((void **)&layer->ub_array, num_out_neurons * sizeof(double));

  cudaMalloc((void **)&layer->expr_array, num_out_neurons * sizeof(expr_t *));
  layer->maxpool_lexpr_array = nullptr;
  layer->maxpool_uexpr_array = nullptr;

  return layer;
}

void fppoly_from_network_input_box(fppoly_t *const res, const size_t intdim,
                                   const size_t realdim,
                                   const double *inf_array,
                                   const double *sup_array) {
  res->layers = nullptr;
  res->numlayers = 0;

  size_t num_pixels = intdim + realdim;

  double *tmp_input_inf = (double *)malloc(num_pixels * sizeof(double));
  double *tmp_input_sup = (double *)malloc(num_pixels * sizeof(double));

  for (size_t i = 0; i < num_pixels; i++) {
    tmp_input_inf[i] = -inf_array[i];
    tmp_input_sup[i] = sup_array[i];
  }

  cudaMalloc((void **)&(res->input_inf), num_pixels * sizeof(double));
  cudaMalloc((void **)&(res->input_sup), num_pixels * sizeof(double));

  res->input_lexpr = nullptr;
  res->input_uexpr = nullptr;

  cudaMemcpy(res->input_inf, tmp_input_inf, num_pixels * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(res->input_sup, tmp_input_sup, num_pixels * sizeof(double),
             cudaMemcpyHostToDevice);

  free(tmp_input_inf);
  free(tmp_input_sup);

  res->num_pixels = num_pixels;
  res->out = nullptr;
}

elina_abstract0_t *fppoly_from_network_input(elina_manager_t *man,
                                             const size_t intdim,
                                             const size_t realdim,
                                             const double *inf_array,
                                             const double *sup_array) {
  fppoly_t *res = (fppoly_t *)malloc(sizeof(fppoly_t));
  fppoly_from_network_input_box(res, intdim, realdim, inf_array, sup_array);

  return abstract0_of_fppoly(man, res);
}

void fppoly_add_new_layer(fppoly_t *const fp, const size_t num_out_neurons,
                          const size_t num_in_neurons, const layertype_t type,
                          const activation_type_t activation) {
  const size_t numlayers = fp->numlayers;
  fp->layers[numlayers] =
      create_layer(num_out_neurons, num_in_neurons, type, activation);
  fp->numlayers++;
}

__device__ void elina_double_interval_add_expr_coeff(
    double *const res_inf, double *const res_sup, const double inf,
    const double sup, const double inf_expr, const double sup_expr) {
  *res_inf = inf + inf_expr;
  *res_sup = sup + sup_expr;
  const double maxA = fmax(fabs(inf_expr), fabs(sup_expr));
  double tmp1, tmp2;
  elina_double_interval_mul(&tmp1, &tmp2, inf, sup, maxA * ulp, maxA * ulp);
  *res_inf += tmp1;
  *res_sup += tmp2;
}

__device__ void elina_double_interval_add_cst_coeff(
    double *const res_inf, double *const res_sup, const double inf,
    const double sup, const double inf_expr, const double sup_expr) {
  elina_double_interval_add_expr_coeff(res_inf, res_sup, inf, sup, inf_expr,
                                       sup_expr);
  *res_inf += min_denormal;
  *res_sup += min_denormal;
}

__device__ void elina_double_interval_mul_expr_coeff(
    double *const res_inf, double *const res_sup, const double inf,
    const double sup, const double inf_expr, const double sup_expr) {
  elina_double_interval_mul(res_inf, res_sup, inf, sup, inf_expr, sup_expr);
  const double maxA = fmax(fabs(inf_expr), fabs(sup_expr));
  double tmp1, tmp2;
  elina_double_interval_mul(&tmp1, &tmp2, inf, sup, maxA * ulp, maxA * ulp);
  *res_inf += tmp1;
  *res_sup += tmp2;
}

__device__ void elina_double_interval_mul_cst_coeff(
    double *const res_inf, double *const res_sup, const double inf,
    const double sup, const double inf_expr, const double sup_expr) {
  elina_double_interval_mul_expr_coeff(res_inf, res_sup, inf, sup, inf_expr,
                                       sup_expr);
  *res_inf += min_denormal;
  *res_sup += min_denormal;
}

__device__ expr_t *multiply_expr(expr_t *const expr, const double mul_inf,
                                 const double mul_sup) {
  expr_t *bla = alloc_expr();

  bla->inf_coeff = (double *)malloc(expr->size * sizeof(double));
  bla->sup_coeff = (double *)malloc(expr->size * sizeof(double));

  bla->type = expr->type;

  for (size_t j = 0; j < expr->size; j++) {
    elina_double_interval_mul_expr_coeff(&bla->inf_coeff[j], &bla->sup_coeff[j],
                                         mul_inf, mul_sup, expr->inf_coeff[j],
                                         expr->sup_coeff[j]);
  }

  bla->size = expr->size;

  elina_double_interval_mul_cst_coeff(&bla->inf_cst, &bla->sup_cst, mul_inf,
                                      mul_sup, expr->inf_cst, expr->sup_cst);

  return bla;
}

__device__ expr_t *multiply_cst_expr(expr_t *const expr, const double mul_inf,
                                     const double mul_sup) {
  expr_t *res = alloc_expr();
  res->inf_coeff = nullptr;
  res->sup_coeff = nullptr;
  res->dim = nullptr;
  res->type = expr->type;
  res->size = expr->size;
  elina_double_interval_mul_cst_coeff(&res->inf_cst, &res->sup_cst, mul_inf,
                                      mul_sup, expr->inf_cst, expr->sup_cst);
  // res->cst = mul_coeff*expr->cst;

  return res;
}

__device__ void add_cst_expr(expr_t *const exprA, expr_t *const exprB) {
  double maxA = fmax(fabs(exprA->inf_cst), fabs(exprA->sup_cst));
  double maxB = fmax(fabs(exprB->inf_cst), fabs(exprB->sup_cst));
  exprA->inf_cst =
      exprA->inf_cst + exprB->inf_cst + (maxA + maxB) * ulp + min_denormal;
  exprA->sup_cst =
      exprA->sup_cst + exprB->sup_cst + (maxA + maxB) * ulp + min_denormal;
}

// A = A + B
__device__ void add_expr(expr_t *const exprA, expr_t *const exprB) {
  assert(exprA->size == exprB->size);

  double maxA = fmax(fabs(exprA->inf_cst), fabs(exprA->sup_cst));
  double maxB = fmax(fabs(exprB->inf_cst), fabs(exprB->sup_cst));

  exprA->inf_cst += exprB->inf_cst + (maxA + maxB) * ulp + min_denormal;
  exprA->sup_cst += exprB->sup_cst + (maxA + maxB) * ulp + min_denormal;

  for (size_t j = 0; j < exprB->size; j++) {
    maxA = fmax(fabs(exprA->inf_coeff[j]), fabs(exprA->sup_coeff[j]));
    maxB = fmax(fabs(exprB->inf_coeff[j]), fabs(exprB->sup_coeff[j]));

    exprA->inf_coeff[j] =
        exprA->inf_coeff[j] + exprB->inf_coeff[j] + (maxA + maxB) * ulp;
    exprA->sup_coeff[j] =
        exprA->sup_coeff[j] + exprB->sup_coeff[j] + (maxA + maxB) * ulp;
  }
}

__global__ void compute_lb_from_expr(double *lb_array, expr_t **expr_array,
                                     double *input_inf, double *input_sup,
                                     const size_t size) {

  size_t n = blockIdx.x;

  if (n < size) {
    expr_t *expr = expr_array[n];

    double res_inf = expr->inf_cst;

    if ((expr->inf_coeff == nullptr) || (expr->sup_coeff == nullptr)) {
      lb_array[n] = 0;

      return;
    }

    double tmp1, tmp2;
    size_t k;

    for (size_t i = 0; i < expr->size; i++) {
      if (expr->type == DENSE) {
        k = i;
      } else {
        k = expr->dim[i];
      }

      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i], input_inf[k], input_sup[k]);
      res_inf = res_inf + tmp1;
    }

    lb_array[n] = res_inf;
  }
}

__global__ void compute_ub_from_expr(double *ub_array, expr_t **expr_array,
                                     double *input_inf, double *input_sup,
                                     const size_t size) {
  size_t n = blockIdx.x;

  if (n < size) {
    expr_t *expr = expr_array[n];

    double res_sup = expr->sup_cst;

    if ((expr->inf_coeff == nullptr) || (expr->sup_coeff == nullptr)) {
      ub_array[n] = 0;

      return;
    }

    double tmp1, tmp2;
    size_t k;

    for (size_t i = 0; i < expr->size; i++) {
      if (expr->type == DENSE) {
        k = i;
      } else {
        k = expr->dim[i];
      }

      elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i],
                                expr->sup_coeff[i], input_inf[k], input_sup[k]);
      res_sup = res_sup + tmp2;
    }

    ub_array[n] = res_sup;
  }
}

__global__ void device_layer_create_dense_expr(expr_t **expr_array,
                                               const double *weights,
                                               const double *bias,
                                               const size_t num_out_neurons,
                                               const size_t num_in_neurons) {
  size_t i = blockIdx.x;

  if (i < num_out_neurons) {
    const double *weight_i = weights + i * num_in_neurons;
    const double bias_i = bias[i];
    expr_array[i] = create_dense_expr(weight_i, bias_i, num_in_neurons);
  }
}

void layer_create_dense_exprs(expr_t **expr_array, const double **weights,
                              const double *bias, const size_t num_out_neurons,
                              const size_t num_in_neurons) {
  double *tmp_weights;
  cudaMalloc((void **)&tmp_weights,
             num_out_neurons * num_in_neurons * sizeof(double));

  double *tmp_bias;
  cudaMalloc((void **)&tmp_bias, num_out_neurons * sizeof(double));

  for (size_t i = 0; i < num_out_neurons; i++) {
    cudaMemcpy(tmp_weights + i * num_in_neurons, weights[i],
               num_in_neurons * sizeof(double), cudaMemcpyHostToDevice);
  }

  cudaMemcpy(tmp_bias, bias, num_out_neurons * sizeof(double),
             cudaMemcpyHostToDevice);

  device_layer_create_dense_expr<<<num_out_neurons, 1>>>(
      expr_array, tmp_weights, tmp_bias, num_out_neurons, num_in_neurons);

  cudaFree(tmp_weights);
  cudaFree(tmp_bias);
}

__global__ void layer_copy_exprs(expr_t **expr_array, expr_t **lexpr_array,
                                 expr_t **uexpr_array, const size_t size) {
  size_t i = blockIdx.x;

  if (i < size) {
    lexpr_array[i] = copy_expr(expr_array[i]);
    uexpr_array[i] = copy_expr(expr_array[i]);
  }
}

__global__ void free_expr_array(expr_t **expr_array, const size_t size) {
  size_t i = blockIdx.x;

  if (i < size) {
    free_expr(expr_array[i]);
  }
}

void layer_compute_bounds_from_exprs(expr_t **expr_array, double *lb_array,
                                     double *ub_array, double *input_inf,
                                     double *input_sup, expr_t **input_lexpr,
                                     expr_t **input_uexpr, const size_t size) {
  // allocate
  expr_t **lexpr_array;
  expr_t **uexpr_array;

  cudaMalloc((void **)&lexpr_array, size * sizeof(expr_t *));
  cudaMalloc((void **)&uexpr_array, size * sizeof(expr_t *));

  layer_copy_exprs<<<size, 1>>>(expr_array, lexpr_array, uexpr_array, size);

  compute_lb_from_expr<<<size, 1>>>(lb_array, lexpr_array, input_inf, input_sup,
                                    size);
  compute_ub_from_expr<<<size, 1>>>(ub_array, uexpr_array, input_inf, input_sup,
                                    size);

  // free
  free_expr_array<<<size, 1>>>(lexpr_array, size);
  free_expr_array<<<size, 1>>>(uexpr_array, size);

  cudaFree(lexpr_array);
  cudaFree(uexpr_array);
}

void ffn_handle_first_layer(elina_manager_t *man, elina_abstract0_t *abs,
                            const double **weights, const double *bias,
                            const size_t size, const size_t num_pixels,
                            const activation_type_t activation) {
  fppoly_t *res = fppoly_of_abstract0(abs);
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  res->layers = (layer_t **)malloc(20 * sizeof(layer_t *));
  fppoly_add_new_layer(res, size, num_pixels, FFN, activation);

  expr_t **expr_array = res->layers[0]->expr_array;

  layer_create_dense_exprs(expr_array, weights, bias, size, num_pixels);
  layer_compute_bounds_from_exprs(
      expr_array, res->layers[0]->lb_array, res->layers[0]->ub_array,
      res->input_inf, res->input_sup, res->input_lexpr, res->input_uexpr, size);
}

void ffn_handle_first_relu_layer(elina_manager_t *man, elina_abstract0_t *abs,
                                 const double **weights, const double *bias,
                                 const size_t size, const size_t num_pixels) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels, RELU);
}

void ffn_handle_first_sigmoid_layer(elina_manager_t *man,
                                    elina_abstract0_t *abs,
                                    const double **weights, const double *bias,
                                    const size_t size,
                                    const size_t num_pixels) {
  // ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels, SIGMOID);
}

void ffn_handle_first_tanh_layer(elina_manager_t *man, elina_abstract0_t *abs,
                                 const double **weights, const double *bias,
                                 const size_t size, const size_t num_pixels) {
  // ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels, TANH);
}

__global__ void
lexpr_replace_relu_bounds(expr_t **expr_array, double *lb_array,
                          double *ub_array,
                          const size_t num_out_neurons_last_layer,
                          const size_t num_out_neurons_current_layer) {
  size_t n = blockIdx.x;
  size_t i = blockIdx.y * blockDim.x + threadIdx.x;

  if (n < num_out_neurons_last_layer) {
    expr_t *expr = expr_array[n];

    if (i < num_out_neurons_current_layer) {
      const double lb = lb_array[i];
      const double ub = ub_array[i];
      const double width = ub + lb;
      const double lambda_inf = -ub / width;
      const double lambda_sup = ub / width;

      const double old_inf_coeff = expr->inf_coeff[i];
      const double old_sup_coeff = expr->sup_coeff[i];

      if ((old_sup_coeff == 0) && (old_inf_coeff == 0)) {
        expr->inf_coeff[i] = 0.0;
        expr->sup_coeff[i] = 0.0;

        return;
      } else if (ub <= 0) {
        expr->inf_coeff[i] = 0.0;
        expr->sup_coeff[i] = 0.0;

        return;
      } else if (lb < 0) {
        expr->inf_coeff[i] = old_inf_coeff;
        expr->sup_coeff[i] = old_sup_coeff;
      } else if (old_sup_coeff < 0) {
        const double mu_inf = lambda_inf * lb;
        const double mu_sup = lambda_sup * lb;
        elina_double_interval_mul_expr_coeff(
            &expr->inf_coeff[i], &expr->sup_coeff[i], lambda_inf, lambda_sup,
            old_inf_coeff, old_sup_coeff);
        double tmp1, tmp2;
        elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, mu_inf, mu_sup,
                                            old_inf_coeff, old_sup_coeff);
        expr->inf_cst = expr->inf_cst + tmp1 + min_denormal;
        expr->sup_cst = expr->sup_cst + tmp2 + min_denormal;
      } else if (old_inf_coeff < 0) {
        const double area1 = lb * ub;
        const double area2 = 0.5 * ub * width;
        const double area3 = 0.5 * lb * width;

        if ((area1 < area2) && (area1 < area3)) {
          elina_double_interval_mul_expr_coeff(
              &expr->inf_coeff[i], &expr->sup_coeff[i], lambda_inf, lambda_sup,
              old_inf_coeff, old_sup_coeff);
        } else if ((area2 < area1) && (area2 < area3)) {
          expr->inf_coeff[i] = 0.0;
          expr->sup_coeff[i] = 0.0;
        } else {
          expr->inf_coeff[i] = old_inf_coeff;
          expr->sup_coeff[i] = old_sup_coeff;
        }
      } else {
        expr->inf_coeff[i] = 0.0;
        expr->sup_coeff[i] = 0.0;
        double tmp1, tmp2;
        elina_double_interval_mul(&tmp1, &tmp2, old_inf_coeff, old_sup_coeff, 0,
                                  ub);
        expr->inf_cst = expr->inf_cst + tmp1;
        expr->sup_cst = expr->sup_cst - tmp1;
      }
    }
  }
}

__global__ void
uexpr_replace_relu_bounds(expr_t **expr_array, double *lb_array,
                          double *ub_array,
                          const size_t num_out_neurons_last_layer,
                          const size_t num_out_neurons_current_layer) {
  size_t n = blockIdx.x;
  size_t i = blockIdx.y * blockDim.x + threadIdx.x;

  if (n < num_out_neurons_last_layer) {
    expr_t *expr = expr_array[n];

    if (i < num_out_neurons_current_layer) {
      const double lb = lb_array[i];
      const double ub = ub_array[i];
      const double width = ub + lb;
      const double lambda_inf = -ub / width;
      const double lambda_sup = ub / width;

      const double old_inf_coeff = expr->inf_coeff[i];
      const double old_sup_coeff = expr->sup_coeff[i];

      if ((old_sup_coeff == 0) && (old_inf_coeff == 0)) {
        expr->inf_coeff[i] = 0.0;
        expr->sup_coeff[i] = 0.0;

        return;
      } else if (ub <= 0) {
        expr->inf_coeff[i] = 0.0;
        expr->sup_coeff[i] = 0.0;

        return;
      } else if (lb < 0) {
        expr->inf_coeff[i] = old_inf_coeff;
        expr->sup_coeff[i] = old_sup_coeff;
      } else if (old_inf_coeff < 0) {
        const double mu_inf = lambda_inf * lb;
        const double mu_sup = lambda_sup * lb;
        elina_double_interval_mul_expr_coeff(
            &expr->inf_coeff[i], &expr->sup_coeff[i], lambda_inf, lambda_sup,
            old_inf_coeff, old_sup_coeff);
        double tmp1, tmp2;
        elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, mu_inf, mu_sup,
                                            old_inf_coeff, old_sup_coeff);
        expr->inf_cst = expr->inf_cst + tmp1 + min_denormal;
        expr->sup_cst = expr->sup_cst + tmp2 + min_denormal;
      } else if (old_sup_coeff < 0) {
        const double area1 = lb * ub;
        const double area2 = 0.5 * ub * width;
        const double area3 = 0.5 * lb * width;

        if ((area1 < area2) && (area1 < area3)) {
          elina_double_interval_mul_expr_coeff(
              &expr->inf_coeff[i], &expr->sup_coeff[i], lambda_inf, lambda_sup,
              old_inf_coeff, old_sup_coeff);
        } else if ((area2 < area1) && (area2 < area3)) {
          expr->inf_coeff[i] = 0.0;
          expr->sup_coeff[i] = 0.0;
        } else {
          expr->inf_coeff[i] = old_inf_coeff;
          expr->sup_coeff[i] = old_sup_coeff;
        }
      } else {
        expr->inf_coeff[i] = 0.0;
        expr->sup_coeff[i] = 0.0;
        double tmp1, tmp2;
        elina_double_interval_mul(&tmp1, &tmp2, old_inf_coeff, old_sup_coeff, 0,
                                  ub);
        expr->inf_cst = expr->inf_cst - tmp2;
        expr->sup_cst = expr->sup_cst + tmp2;
      }
    }
  }
}

__global__ void
expr_from_previous_layer(expr_t **expr_array, expr_t **res_array,
                         expr_t **aux_expr_array,
                         const size_t num_out_neurons_last_layer,
                         const size_t num_out_neurons_current_layer,
                         const size_t num_in_neurons_current_layer) {
  size_t n = blockIdx.x;
  size_t j = blockIdx.y * blockDim.x + threadIdx.x;

  if (n < num_out_neurons_last_layer) {
    expr_t *expr = expr_array[n];
    expr_t *res = res_array[n];

    size_t i;

    if (j < num_in_neurons_current_layer) {
      if (j == 0) {
        i = 0;

        expr_t *aux_expr = aux_expr_array[i];

        elina_double_interval_mul_cst_coeff(
            &res->inf_cst, &res->sup_cst, expr->inf_coeff[0],
            expr->sup_coeff[0], aux_expr->inf_cst, aux_expr->sup_cst);

        double tmp1, tmp2;
        double maxRes, maxMul;

        for (i = 1; i < num_out_neurons_current_layer; i++) {
          if ((expr->inf_coeff[i] != 0) || (expr->sup_coeff[i] != 0)) {
            aux_expr = aux_expr_array[i];

            elina_double_interval_mul_cst_coeff(
                &tmp1, &tmp2, expr->inf_coeff[i], expr->sup_coeff[i],
                aux_expr->inf_cst, aux_expr->sup_cst);

            maxRes = fmax(fabs(res->inf_cst), fabs(res->sup_cst));
            maxMul = fmax(fabs(tmp1), fabs(tmp2));

            res->inf_cst += tmp1 + (maxRes + maxMul) * ulp + min_denormal;
            res->sup_cst += tmp2 + (maxRes + maxMul) * ulp + min_denormal;
          }
        }

        res->inf_cst = res->inf_cst + expr->inf_cst;
        res->sup_cst = res->sup_cst + expr->sup_cst;
      }

      i = 0;

      expr_t *aux_expr = aux_expr_array[i];

      elina_double_interval_mul_expr_coeff(
          &res->inf_coeff[j], &res->sup_coeff[j], expr->inf_coeff[0],
          expr->sup_coeff[0], aux_expr->inf_coeff[j], aux_expr->sup_coeff[j]);

      double tmp1, tmp2;
      double maxRes, maxMul;

      for (i = 1; i < num_out_neurons_current_layer; i++) {
        if ((expr->inf_coeff[i] != 0) || (expr->sup_coeff[i] != 0)) {
          aux_expr = aux_expr_array[i];

          elina_double_interval_mul_expr_coeff(
              &tmp1, &tmp2, expr->inf_coeff[i], expr->sup_coeff[i],
              aux_expr->inf_coeff[j], aux_expr->sup_coeff[j]);

          maxRes = fmax(fabs(res->inf_coeff[j]), fabs(res->sup_coeff[j]));
          maxMul = fmax(fabs(tmp1), fabs(tmp2));

          res->inf_coeff[j] =
              res->inf_coeff[j] + tmp1 + (maxRes + maxMul) * ulp;
          res->sup_coeff[j] =
              res->sup_coeff[j] + tmp2 + (maxRes + maxMul) * ulp;
        }
      }
    }
  }
}

__global__ void layer_allocate_exprs(expr_t **expr_array,
                                     const size_t array_size,
                                     const size_t expr_size) {
  size_t i = blockIdx.x;

  if (i < array_size) {
    expr_array[i] = alloc_dense_expr_and_arrays(expr_size);
  }
}

void update_state_using_previous_layers(elina_manager_t *man, fppoly_t *fp,
                                        const size_t layerno) {
  auto start = std::chrono::system_clock::now();

  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  const size_t num_out_neurons_last_layer =
      fp->layers[layerno]->num_out_neurons;
  std::cout << "num_out_neurons_last " << num_out_neurons_last_layer
            << std::endl;

  expr_t **expr_array = fp->layers[layerno]->expr_array;

  double *lb_array = fp->layers[layerno]->lb_array;
  double *ub_array = fp->layers[layerno]->ub_array;

  expr_t **lexpr_array;
  expr_t **uexpr_array;

  cudaMalloc((void **)&lexpr_array,
             num_out_neurons_last_layer * sizeof(expr_t *));
  cudaMalloc((void **)&uexpr_array,
             num_out_neurons_last_layer * sizeof(expr_t *));

  expr_t **lexpr_array_tmp;
  expr_t **uexpr_array_tmp;

  cudaMalloc((void **)&lexpr_array_tmp,
             num_out_neurons_last_layer * sizeof(expr_t *));
  cudaMalloc((void **)&uexpr_array_tmp,
             num_out_neurons_last_layer * sizeof(expr_t *));

  layer_copy_exprs<<<num_out_neurons_last_layer, 1>>>(
      expr_array, lexpr_array, uexpr_array, num_out_neurons_last_layer);

  for (int k = layerno - 1; k >= 0; k--) {
    const size_t num_out_neurons_current_layer = fp->layers[k]->num_out_neurons;
    const size_t num_in_neurons_current_layer = fp->layers[k]->num_in_neurons;
    std::cout << "num_out_neurons_current " << num_out_neurons_current_layer
              << " num_in_neurons_current " << num_in_neurons_current_layer
              << std::endl;

    const size_t num_threads = 512;

    const dim3 num_blocks_relu(num_out_neurons_last_layer,
                               num_out_neurons_current_layer / 512 + 1, 1);
    const dim3 num_blocks_linear(num_out_neurons_last_layer,
                                 num_in_neurons_current_layer / 512 + 1, 1);

    std::cout << "num_threads" << num_threads << " num_blocks_relu "
              << num_blocks_relu.y << " num_blocks_linear "
              << num_blocks_linear.y << std::endl;

    expr_t **aux_expr_array = fp->layers[k]->expr_array;

    double *aux_lb_array = fp->layers[k]->lb_array;
    double *aux_ub_array = fp->layers[k]->ub_array;

    if (fp->layers[k]->activation == RELU) {
      lexpr_replace_relu_bounds<<<num_blocks_relu, num_threads>>>(
          lexpr_array, aux_lb_array, aux_ub_array, num_out_neurons_last_layer,
          num_out_neurons_current_layer);
      uexpr_replace_relu_bounds<<<num_blocks_relu, num_threads>>>(
          uexpr_array, aux_lb_array, aux_ub_array, num_out_neurons_last_layer,
          num_out_neurons_current_layer);
    }

    layer_allocate_exprs<<<num_out_neurons_last_layer, 1>>>(
        lexpr_array_tmp, num_out_neurons_last_layer,
        num_in_neurons_current_layer);
    layer_allocate_exprs<<<num_out_neurons_last_layer, 1>>>(
        uexpr_array_tmp, num_out_neurons_last_layer,
        num_in_neurons_current_layer);

    expr_from_previous_layer<<<num_blocks_linear, num_threads>>>(
        lexpr_array, lexpr_array_tmp, aux_expr_array,
        num_out_neurons_last_layer, num_out_neurons_current_layer,
        num_in_neurons_current_layer);
    expr_from_previous_layer<<<num_blocks_linear, num_threads>>>(
        uexpr_array, uexpr_array_tmp, aux_expr_array,
        num_out_neurons_last_layer, num_out_neurons_current_layer,
        num_in_neurons_current_layer);

    std::swap(lexpr_array, lexpr_array_tmp);
    std::swap(uexpr_array, uexpr_array_tmp);

    free_expr_array<<<num_out_neurons_last_layer, 1>>>(
        lexpr_array_tmp, num_out_neurons_current_layer);
    free_expr_array<<<num_out_neurons_last_layer, 1>>>(
        uexpr_array_tmp, num_out_neurons_current_layer);
  }

  compute_lb_from_expr<<<num_out_neurons_last_layer, 1>>>(
      lb_array, lexpr_array, fp->input_inf, fp->input_sup,
      num_out_neurons_last_layer);
  compute_ub_from_expr<<<num_out_neurons_last_layer, 1>>>(
      ub_array, uexpr_array, fp->input_inf, fp->input_sup,
      num_out_neurons_last_layer);

  if (fp->out != nullptr) {
    fp->out->lexpr = lexpr_array;
    fp->out->uexpr = uexpr_array;
  } else {
    free_expr_array<<<num_out_neurons_last_layer, 1>>>(
        lexpr_array, num_out_neurons_last_layer);
    free_expr_array<<<num_out_neurons_last_layer, 1>>>(
        uexpr_array, num_out_neurons_last_layer);

    cudaFree(lexpr_array);
    cudaFree(uexpr_array);
  }

  cudaFree(lexpr_array_tmp);
  cudaFree(uexpr_array_tmp);

  cudaDeviceSynchronize();

  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl
            << std::endl;
}

void ffn_handle_intermediate_layer(elina_manager_t *man,
                                   elina_abstract0_t *element,
                                   const double **weights, const double *bias,
                                   const size_t num_out_neurons,
                                   const size_t num_in_neurons,
                                   const activation_type_t activation) {
  fppoly_t *fp = fppoly_of_abstract0(element);
  fppoly_add_new_layer(fp, num_out_neurons, num_in_neurons, FFN, activation);
  expr_t **expr_array = fp->layers[fp->numlayers - 1]->expr_array;

  layer_create_dense_exprs(expr_array, weights, bias, num_out_neurons,
                           num_in_neurons);

  update_state_using_previous_layers(man, fp, fp->numlayers - 1);
}

void ffn_handle_intermediate_relu_layer(elina_manager_t *man,
                                        elina_abstract0_t *element,
                                        const double **weights,
                                        const double *bias,
                                        const size_t num_out_neurons,
                                        const size_t num_in_neurons) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, RELU);
}

void ffn_handle_intermediate_sigmoid_layer(elina_manager_t *man,
                                           elina_abstract0_t *element,
                                           const double **weights,
                                           const double *bias,
                                           const size_t num_out_neurons,
                                           const size_t num_in_neurons) {
  // ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
  // num_in_neurons, SIGMOID);
}

void ffn_handle_intermediate_tanh_layer(elina_manager_t *man,
                                        elina_abstract0_t *element,
                                        const double **weights,
                                        const double *bias,
                                        const size_t num_out_neurons,
                                        const size_t num_in_neurons) {
  // ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
  // num_in_neurons, TANH);
}

__global__ void apply_relu_lexpr(double *lb_array, expr_t **lexpr_array,
                                 double *aux_lb_array, double *aux_ub_array,
                                 const size_t num_out_neurons) {
  size_t n = blockIdx.x;

  if (n < num_out_neurons) {
    expr_t *lexpr = lexpr_array[n];

    const size_t size = lexpr->size;
    const double lb = aux_lb_array[n];
    const double ub = aux_ub_array[n];
    const double width = lb + ub;

    if (ub < 0) {
      free_expr(lexpr_array[n]);
      lexpr_array[n] = nullptr;

      lb_array[n] = 0;

      return;
    }

    if (lb < 0) {
      lb_array[n] = lb;

      return;
    }

    const double area1 = lb * ub;
    const double area2 = 0.5 * ub * width;

    if (area1 < area2) {
      const double lambda_inf = -ub / width;
      const double lambda_sup = ub / width;

      for (size_t i = 0; i < size; i++) {
        // lexpr->coeff[i] = lexpr->coeff[i]*lambda;
        elina_double_interval_mul_expr_coeff(
            &lexpr->inf_coeff[i], &lexpr->sup_coeff[i], lambda_inf, lambda_sup,
            lexpr->inf_coeff[i], lexpr->sup_coeff[i]);
      }
      // lexpr->cst = lexpr->cst*lambda;
      elina_double_interval_mul_cst_coeff(&lexpr->inf_cst, &lexpr->sup_cst,
                                          lambda_inf, lambda_sup,
                                          lexpr->inf_cst, lexpr->sup_cst);
      // double res, res1;

      lb_array[n] = -(lambda_inf * lb);
    } else {
      free_expr(lexpr_array[n]);
      lexpr_array[n] = nullptr;

      lb_array[n] = 0;
    }
  }
}

__global__ void apply_relu_uexpr(double *ub_array, expr_t **uexpr_array,
                                 double *aux_lb_array, double *aux_ub_array,
                                 const size_t num_out_neurons) {
  size_t n = blockIdx.x;

  if (n < num_out_neurons) {
    expr_t *uexpr = uexpr_array[n];

    const size_t size = uexpr->size;
    const double lb = aux_lb_array[n];
    const double ub = aux_ub_array[n];
    const double width = lb + ub;

    if (ub < 0) {
      free_expr(uexpr_array[n]);
      uexpr_array[n] = nullptr;

      ub_array[n] = 0;

      return;
    }

    if (lb < 0) {
      ub_array[n] = ub;

      return;
    }

    const double lambda_inf = -ub / width;
    const double lambda_sup = ub / width;

    for (size_t i = 0; i < size; i++) {
      // uexpr->coeff[i] = uexpr->coeff[i]*lambda;
      elina_double_interval_mul_expr_coeff(
          &uexpr->inf_coeff[i], &uexpr->sup_coeff[i], lambda_inf, lambda_sup,
          uexpr->inf_coeff[i], uexpr->sup_coeff[i]);
    }

    elina_double_interval_mul_cst_coeff(&uexpr->inf_cst, &uexpr->sup_cst,
                                        lambda_inf, lambda_sup, uexpr->inf_cst,
                                        uexpr->sup_cst);
    const double mu_inf = lambda_inf * lb;
    const double mu_sup = lambda_sup * lb;
    // uexpr->cst = uexpr->cst*lambda;
    // uexpr->cst = uexpr->cst + lambda*lb;
    uexpr->inf_cst += mu_inf;
    uexpr->sup_cst += mu_sup;

    ub_array[n] = ub;
  }
}

__global__ void assign_output_inf_sup(double *lb_array, double *ub_array,
                                      double *aux_lb_array,
                                      double *aux_ub_array,
                                      const size_t num_out_neurons) {
  size_t i = blockIdx.x;

  if (i < num_out_neurons) {
    lb_array[i] = aux_lb_array[i];
    ub_array[i] = aux_ub_array[i];
    printf("out inf number %i is: %g\n", i, lb_array[i]);
    printf("out sup number %i is: %g\n", i, ub_array[i]);
  }
}

void handle_final_relu_layer(fppoly_internal_t *pr, output_abstract_t *out,
                             double *aux_lb_array, double *aux_ub_array,
                             const size_t size, const bool has_relu) {
  if (has_relu) {
    apply_relu_lexpr<<<size, 1>>>(out->output_inf, out->lexpr, aux_lb_array,
                                  aux_ub_array, size);
    apply_relu_uexpr<<<size, 1>>>(out->output_sup, out->uexpr, aux_lb_array,
                                  aux_ub_array, size);
  } else {
    assign_output_inf_sup<<<size, 1>>>(out->output_inf, out->output_sup,
                                       aux_lb_array, aux_ub_array, size);
  }
}

void handle_final_non_relu_layer(fppoly_internal_t *pr, output_abstract_t *out,
                                 double *aux_lb_array, double *aux_ub_array,
                                 const size_t size) {
  assign_output_inf_sup<<<size, 1>>>(out->output_inf, out->output_sup,
                                     aux_lb_array, aux_ub_array, size);
}

output_abstract_t *allocate_output_abstract(const size_t num_out_neurons) {
  output_abstract_t *out =
      (output_abstract_t *)malloc(sizeof(output_abstract_t));

  cudaMalloc((void **)&(out->output_inf), num_out_neurons * sizeof(double));
  cudaMalloc((void **)&(out->output_sup), num_out_neurons * sizeof(double));

  out->lexpr = nullptr;
  out->uexpr = nullptr;

  return out;
}

void ffn_handle_last_layer(elina_manager_t *man, elina_abstract0_t *element,
                           const double **weights, const double *bias,
                           const size_t num_out_neurons,
                           const size_t num_in_neurons,
                           const bool has_activation,
                           const activation_type_t activation) {
  fppoly_t *fp = fppoly_of_abstract0(element);

  if (has_activation) {
    fppoly_add_new_layer(fp, num_out_neurons, num_in_neurons, FFN, activation);
  } else {
    fppoly_add_new_layer(fp, num_out_neurons, num_in_neurons, FFN, NONE);
  }

  fp->out = allocate_output_abstract(num_out_neurons);

  expr_t **expr_array = fp->layers[fp->numlayers - 1]->expr_array;

  double *lb_array = fp->layers[fp->numlayers - 1]->lb_array;
  double *ub_array = fp->layers[fp->numlayers - 1]->ub_array;

  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  layer_create_dense_exprs(expr_array, weights, bias, num_out_neurons,
                           num_in_neurons);

  update_state_using_previous_layers(man, fp, fp->numlayers - 1);

  if (activation == RELU) {
    handle_final_relu_layer(pr, fp->out, lb_array, ub_array, num_out_neurons,
                            has_activation);
  } else {
    handle_final_non_relu_layer(pr, fp->out, lb_array, ub_array,
                                num_out_neurons);
  }
}

void ffn_handle_last_relu_layer(elina_manager_t *man,
                                elina_abstract0_t *element,
                                const double **weights, const double *bias,
                                const size_t num_out_neurons,
                                const size_t num_in_neurons,
                                const bool has_relu) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, has_relu, RELU);
}

void ffn_handle_last_sigmoid_layer(elina_manager_t *man,
                                   elina_abstract0_t *element,
                                   const double **weights, const double *bias,
                                   const size_t num_out_neurons,
                                   const size_t num_in_neurons,
                                   const bool has_sigmoid) {
  // ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
  // num_in_neurons, has_sigmoid, SIGMOID);
}

void ffn_handle_last_tanh_layer(elina_manager_t *man,
                                elina_abstract0_t *element,
                                const double **weights, const double *bias,
                                const size_t num_out_neurons,
                                const size_t num_in_neurons,
                                const bool has_tanh) {
  // ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
  // num_in_neurons, has_tanh, TANH);
}

__global__ void create_sub_expr(expr_t **sub, const size_t index,
                                const elina_dim_t y, const elina_dim_t x) {

  expr_t *expr = (expr_t *)malloc(sizeof(expr_t));

  expr->inf_cst = 0;
  expr->sup_cst = 0;

  expr->inf_coeff = (double *)malloc(10 * sizeof(double));
  expr->sup_coeff = (double *)malloc(10 * sizeof(double));
  expr->dim = nullptr;

  expr->size = 10;
  expr->type = DENSE;

  for (size_t i = 0; i < 10; i++) {
    expr->inf_coeff[i] = 0.;
    expr->sup_coeff[i] = 0.;
  }

  expr->inf_coeff[y] = -1.;
  expr->sup_coeff[y] = 1.;

  expr->inf_coeff[x] = 1.;
  expr->sup_coeff[x] = -1.;

  sub[index] = expr;
}

void get_lb_using_previous_layers(elina_manager_t *man,
                                  const fppoly_t *const fp) {
  const size_t numlayers = fp->numlayers;
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  const size_t num_out_neurons_last_layer = 90;

  double *lb_dev;
  cudaMalloc((void **)&lb_dev, num_out_neurons_last_layer * sizeof(double));

  expr_t **lexpr_array;
  cudaMalloc((void **)&lexpr_array,
             num_out_neurons_last_layer * sizeof(expr_t *));

  size_t index = 0;

  for (elina_dim_t y = 0; y < 10; y++) {
    for (elina_dim_t x = 0; x < 10; x++) {
      if (y != x) {
        create_sub_expr<<<1, 1>>>(lexpr_array, index, y, x);
        index++;
      }
    }
  }

  expr_t **lexpr_array_tmp;
  cudaMalloc((void **)&lexpr_array_tmp,
             num_out_neurons_last_layer * sizeof(expr_t *));

  for (int k = numlayers - 1; k >= 0; k--) {
    const size_t num_out_neurons_current_layer = fp->layers[k]->num_out_neurons;
    const size_t num_in_neurons_current_layer = fp->layers[k]->num_in_neurons;

    const size_t num_threads = 512;

    const dim3 num_blocks_relu(num_out_neurons_last_layer,
                               num_out_neurons_current_layer / 512 + 1, 1);
    const dim3 num_blocks_linear(num_out_neurons_last_layer,
                                 num_in_neurons_current_layer / 512 + 1, 1);

    expr_t **aux_expr_array = fp->layers[k]->expr_array;

    double *aux_lb_array = fp->layers[k]->lb_array;
    double *aux_ub_array = fp->layers[k]->ub_array;

    if (fp->layers[k]->activation == RELU) {
      lexpr_replace_relu_bounds<<<num_blocks_relu, num_threads>>>(
          lexpr_array, aux_lb_array, aux_ub_array, num_out_neurons_last_layer,
          num_out_neurons_current_layer);
    }

    layer_allocate_exprs<<<num_out_neurons_last_layer, 1>>>(
        lexpr_array_tmp, num_out_neurons_last_layer,
        num_in_neurons_current_layer);

    expr_from_previous_layer<<<num_blocks_linear, num_threads>>>(
        lexpr_array, lexpr_array_tmp, aux_expr_array,
        num_out_neurons_last_layer, num_out_neurons_current_layer,
        num_in_neurons_current_layer);

    std::swap(lexpr_array, lexpr_array_tmp);

    free_expr_array<<<num_out_neurons_last_layer, 1>>>(
        lexpr_array_tmp, num_out_neurons_current_layer);
  }

  compute_lb_from_expr<<<num_out_neurons_last_layer, 1>>>(
      lb_dev, lexpr_array, fp->input_inf, fp->input_sup,
      num_out_neurons_last_layer);

  double lb[num_out_neurons_last_layer];
  cudaMemcpy(&lb, lb_dev, num_out_neurons_last_layer * sizeof(double),
             cudaMemcpyDeviceToHost);

  free_expr_array<<<num_out_neurons_last_layer, 1>>>(
      lexpr_array, num_out_neurons_last_layer);
  cudaFree(lexpr_array);
  cudaFree(lexpr_array_tmp);
  cudaFree(lb_dev);

  for (size_t i = 0; i < num_out_neurons_last_layer; i++) {
    if (lb[i] < 0) {
      results[i] = true;
    } else {
      results[i] = false;
    }
  }
}

bool is_greater(elina_manager_t *man, elina_abstract0_t *element,
                const elina_dim_t y, const elina_dim_t x) {
  const fppoly_t *fp = fppoly_of_abstract0(element);
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  if (!results_calculated) {
    get_lb_using_previous_layers(man, fp);
    results_calculated = true;

    return results[0];
  } else {
    bool result = results[output_counter];
    output_counter++;

    return result;
  }
}

__global__ void
create_dense_expr_device(expr_t **expr_array, const size_t index,
                         const double *inf_coeff, const double *sup_coeff,
                         const double inf_cst, const double sup_cst,
                         const size_t size, const exprtype_t type) {
  expr_t *expr = (expr_t *)malloc(sizeof(expr_t));

  expr->inf_cst = inf_cst;
  expr->sup_cst = sup_cst;

  expr->inf_coeff = (double *)malloc(size * sizeof(double));
  expr->sup_coeff = (double *)malloc(size * sizeof(double));

  for (size_t i = 0; i < size; i++) {
    expr->inf_coeff[i] = inf_coeff[i];
    expr->sup_coeff[i] = sup_coeff[i];
  }

  expr->dim = nullptr;
  expr->size = size;
  expr->type = type;

  expr_array[index] = expr;
}

void copy_dense_expr_host_to_device(expr_t **expr_array, size_t index,
                                    const expr_t *const src) {
  double *inf_coeff_tmp;
  double *sup_coeff_tmp;

  cudaMalloc((void **)&inf_coeff_tmp, src->size * sizeof(double));
  cudaMalloc((void **)&sup_coeff_tmp, src->size * sizeof(double));

  cudaMemcpy(inf_coeff_tmp, src->inf_coeff, src->size * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(sup_coeff_tmp, src->sup_coeff, src->size * sizeof(double),
             cudaMemcpyHostToDevice);

  create_dense_expr_device<<<1, 1>>>(expr_array, index, inf_coeff_tmp,
                                     sup_coeff_tmp, src->inf_cst, src->sup_cst,
                                     src->size, src->type);

  cudaFree(inf_coeff_tmp);
  cudaFree(sup_coeff_tmp);
}

void device_layer_create_sparse_exprs(
    expr_t **expr_array, const double *filter_weights,
    const double *filter_bias, const size_t *input_size,
    const size_t *output_size, const size_t *filter_size, const size_t *strides,
    const bool has_bias, const long int pad_top, const long int pad_left,
    const size_t num_pixels) {
  for (size_t out_x = 0; out_x < output_size[0]; out_x++) {
    for (size_t out_y = 0; out_y < output_size[1]; out_y++) {
      for (size_t out_z = 0; out_z < output_size[2]; out_z++) {
        const size_t mat_x = out_x * output_size[1] * output_size[2] +
                             out_y * output_size[2] + out_z;
        const size_t num_coeff =
            input_size[2] * filter_size[0] * filter_size[1];
        size_t actual_coeff = 0;
        double *coeff = (double *)malloc(num_coeff * sizeof(double));
        size_t *dim = (size_t *)malloc(num_coeff * sizeof(size_t));
        size_t i = 0;

        for (size_t inp_z = 0; inp_z < input_size[2]; inp_z++) {
          for (size_t x_shift = 0; x_shift < filter_size[0]; x_shift++) {
            for (size_t y_shift = 0; y_shift < filter_size[1]; y_shift++) {
              const long int x_val = out_x * strides[0] + x_shift - pad_top;
              const long int y_val = out_y * strides[1] + y_shift - pad_left;

              if ((y_val < 0) || (y_val >= (long int)input_size[1])) {
                continue;
              }

              if ((x_val < 0) || (x_val >= (long int)input_size[0])) {
                continue;
              }

              const size_t mat_y = x_val * input_size[1] * input_size[2] +
                                   y_val * input_size[2] + inp_z;

              if (mat_y >= num_pixels) {
                continue;
              }

              const size_t filter_index =
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

        sort_sparse_expr(dim, coeff, actual_coeff);

        const double cst = has_bias ? filter_bias[out_z] : 0;

        double *dense_coeff = (double *)malloc(num_pixels * sizeof(double));

        size_t k = 0;
        for (size_t i = 0; i < num_pixels; i++) {
          if ((dim[k] == i) && (k < actual_coeff)) {
            dense_coeff[i] = coeff[k];
            k++;
          } else {
            dense_coeff[i] = 0.;
          }
        }

        expr_t *res = create_dense_expr(dense_coeff, cst, num_pixels);

        copy_dense_expr_host_to_device(expr_array, mat_x, res);

        free_expr(res);
        free(coeff);
        free(dense_coeff);
        free(dim);
      }
    }
  }
}

void layer_create_sparse_exprs(fppoly_t *const fp, const double *filter_weights,
                               const double *filter_bias,
                               const size_t *input_size,
                               const size_t *filter_size,
                               const size_t num_filters, const size_t *strides,
                               const bool is_valid_padding,
                               const bool has_bias) {
  const size_t num_pixels = input_size[0] * input_size[1] * input_size[2];

  size_t output_size[3];

  if (is_valid_padding) {
    output_size[0] =
        ceil((double)(input_size[0] - filter_size[0] + 1) / (double)strides[0]);
    output_size[1] =
        ceil((double)(input_size[1] - filter_size[1] + 1) / (double)strides[1]);
  } else {
    output_size[0] = ceil((double)input_size[0] / (double)strides[0]);
    output_size[1] = ceil((double)input_size[1] / (double)strides[1]);
  }

  output_size[2] = num_filters;

  const size_t num_out_neurons =
      output_size[0] * output_size[1] * output_size[2];
  fppoly_add_new_layer(fp, num_out_neurons, num_pixels, CONV, RELU);
  expr_t **expr_array = fp->layers[fp->numlayers - 1]->expr_array;

  long int pad_along_height = 0;
  long int pad_along_width = 0;
  long int pad_top = 0;
  long int pad_left = 0;

  if (!is_valid_padding) {
    if (input_size[0] % strides[0] == 0) {
      const long int tmp = filter_size[0] - strides[0];
      pad_along_height = max(tmp, long(0));
    } else {
      const long int tmp = filter_size[0] - (input_size[0] % strides[0]);
      pad_along_height = max(tmp, long(0));
    }

    if (input_size[1] % strides[1] == 0) {
      const long int tmp = filter_size[1] - strides[1];
      pad_along_width = max(tmp, long(0));
    } else {
      const long int tmp = filter_size[1] - (input_size[1] % strides[1]);
      pad_along_width = max(tmp, long(0));
    }

    pad_top = pad_along_height / 2;
    pad_left = pad_along_width / 2;
  }

  const size_t size =
      filter_size[0] * filter_size[1] * input_size[2] * output_size[2];

  double *filter_weights_tmp = (double *)malloc(size * sizeof(double));
  double *filter_bias_tmp = (double *)malloc(output_size[2] * sizeof(double));

  size_t *input_size_tmp = (size_t *)malloc(3 * sizeof(size_t));
  size_t *output_size_tmp = (size_t *)malloc(3 * sizeof(size_t));
  size_t *filter_size_tmp = (size_t *)malloc(2 * sizeof(size_t));
  size_t *strides_tmp = (size_t *)malloc(2 * sizeof(size_t));

  cudaMemcpy(filter_weights_tmp, filter_weights, size * sizeof(double),
             cudaMemcpyHostToHost);
  cudaMemcpy(filter_bias_tmp, filter_bias, output_size[2] * sizeof(double),
             cudaMemcpyHostToHost);

  cudaMemcpy(input_size_tmp, input_size, 3 * sizeof(size_t),
             cudaMemcpyHostToHost);
  cudaMemcpy(output_size_tmp, output_size, 3 * sizeof(size_t),
             cudaMemcpyHostToHost);
  cudaMemcpy(filter_size_tmp, filter_size, 2 * sizeof(size_t),
             cudaMemcpyHostToHost);
  cudaMemcpy(strides_tmp, strides, 2 * sizeof(size_t), cudaMemcpyHostToHost);

  device_layer_create_sparse_exprs(
      expr_array, filter_weights_tmp, filter_bias_tmp, input_size_tmp,
      output_size_tmp, filter_size_tmp, strides_tmp, has_bias, pad_top,
      pad_left, num_pixels);

  free(filter_weights_tmp);
  free(filter_bias_tmp);

  free(input_size_tmp);
  free(output_size_tmp);
  free(filter_size_tmp);
  free(strides_tmp);
}

void conv_handle_first_layer(elina_manager_t *man, elina_abstract0_t *element,
                             const double *filter_weights,
                             const double *filter_bias,
                             const size_t *input_size,
                             const size_t *filter_size,
                             const size_t num_filters, const size_t *strides,
                             const bool is_valid_padding, const bool has_bias) {
  fppoly_t *const fp = fppoly_of_abstract0(element);
  fp->layers = (layer_t **)malloc(20 * sizeof(layer_t *));

  layer_create_sparse_exprs(fp, filter_weights, filter_bias, input_size,
                            filter_size, num_filters, strides, is_valid_padding,
                            has_bias);

  const size_t num_out_neurons = fp->layers[fp->numlayers - 1]->num_out_neurons;
  expr_t **expr_array = fp->layers[fp->numlayers - 1]->expr_array;

  layer_compute_bounds_from_exprs(expr_array, fp->layers[0]->lb_array,
                                  fp->layers[0]->ub_array, fp->input_inf,
                                  fp->input_sup, fp->input_lexpr,
                                  fp->input_uexpr, num_out_neurons);
}

void conv_handle_intermediate_relu_layer(
    elina_manager_t *man, elina_abstract0_t *element,
    const double *filter_weights, const double *filter_bias,
    const size_t *input_size, const size_t *filter_size,
    const size_t num_filters, const size_t *strides,
    const bool is_valid_padding, const bool has_bias) {
  fppoly_t *const fp = fppoly_of_abstract0(element);

  layer_create_sparse_exprs(fp, filter_weights, filter_bias, input_size,
                            filter_size, num_filters, strides, is_valid_padding,
                            has_bias);

  update_state_using_previous_layers(man, fp, fp->numlayers - 1);
}

void free_layer(layer_t *layer) {
  free_expr_array<<<layer->num_out_neurons, 1>>>(layer->expr_array,
                                                 layer->num_out_neurons);

  if (layer->maxpool_lexpr_array != nullptr) {
    free_expr_array<<<layer->num_out_neurons, 1>>>(layer->maxpool_lexpr_array,
                                                   layer->num_out_neurons);
  }

  if (layer->maxpool_uexpr_array != nullptr) {
    free_expr_array<<<layer->num_out_neurons, 1>>>(layer->maxpool_uexpr_array,
                                                   layer->num_out_neurons);
  }

  cudaFree(layer->expr_array);
  cudaFree(layer->maxpool_lexpr_array);
  cudaFree(layer->maxpool_uexpr_array);
  layer->expr_array = nullptr;
  layer->maxpool_lexpr_array = nullptr;
  layer->maxpool_uexpr_array = nullptr;

  cudaFree(layer->lb_array);
  cudaFree(layer->ub_array);
  layer->lb_array = nullptr;
  layer->ub_array = nullptr;

  free(layer);
  layer = nullptr;
}

void fppoly_free(elina_manager_t *man, fppoly_t *fp) {
  size_t output_size = fp->layers[fp->numlayers - 1]->num_out_neurons;

  for (size_t i = 0; i < fp->numlayers; i++) {
    free_layer(fp->layers[i]);
  }

  free(fp->layers);
  fp->layers = nullptr;

  free_expr_array<<<output_size, 1>>>(fp->out->lexpr, output_size);
  free_expr_array<<<output_size, 1>>>(fp->out->uexpr, output_size);

  if ((fp->input_lexpr != nullptr) && (fp->input_uexpr != nullptr)) {
    free_expr_array<<<fp->num_pixels, 1>>>(fp->input_lexpr, fp->num_pixels);
    free_expr_array<<<fp->num_pixels, 1>>>(fp->input_uexpr, fp->num_pixels);

    cudaFree(fp->input_lexpr);
    fp->input_lexpr = nullptr;
    cudaFree(fp->input_uexpr);
    fp->input_uexpr = nullptr;
  }

  cudaFree(fp->input_inf);
  fp->input_inf = nullptr;
  cudaFree(fp->input_sup);
  fp->input_sup = nullptr;

  cudaFree(fp->out->output_inf);
  fp->out->output_inf = nullptr;
  cudaFree(fp->out->output_sup);
  fp->out->output_sup = nullptr;

  cudaFree(fp->out->lexpr);
  fp->out->lexpr = nullptr;
  cudaFree(fp->out->uexpr);
  fp->out->uexpr = nullptr;

  free(fp->out);
  fp->out = nullptr;

  free(fp);
  fp = nullptr;
}

void layer_print(const layer_t *layer) {
  // neurons_print<<<1, 1>>>(layer->neurons, layer->num_out_neurons);
}

__global__ void output_print(const double *output_inf, const double *output_sup,
                             const size_t output_size) {
  for (size_t i = 0; i < output_size; i++) {
    printf("%zu: %g \n", i, output_inf[i]);
    printf("%zu: %g \n", i, output_sup[i]);
  }
}

void fppoly_fprint(FILE *const stream, elina_manager_t *man,
                   const fppoly_t *const fp, const char **name_of_dim) {
  for (size_t i = 0; i < fp->numlayers; i++) {
    printf("layer: %zu\n", i);
    layer_print(fp->layers[i]);
  }

  const size_t output_size = fp->layers[fp->numlayers - 1]->num_out_neurons;

  if (fp->out != nullptr) {
    printf("OUTPUT bounds: \n");

    output_print<<<1, 1>>>(fp->out->output_inf, fp->out->output_sup,
                           output_size);
  }
}

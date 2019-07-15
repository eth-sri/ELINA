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

const size_t num_threads = 128;

bool results[90];
bool results_calculated;
size_t output_counter;

constexpr int DIGS = DECIMAL_DIG;

#ifdef single
__constant__ const float_type min_denormal = 1.40129846e-45;
__constant__ const float_type ulp = 1.1920929e-07;
#else
__constant__ const float_type min_denormal = 4.940656458412465441766e-324;
__constant__ const float_type ulp = 2.220446049250313080848e-16;
#endif

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

long int offset_y = 0;
long int offset_x = 0;

long int length_y = 0;
long int length_x = 0;

long int shift_y = 0;
long int shift_x = 0;

__device__ void
elina_double_interval_mul(float_type *const a_inf, float_type *const a_sup,
                          const float_type b_inf, const float_type b_sup,
                          const float_type c_inf, const float_type c_sup) {
  if (c_inf >= 0) {
    /* interval c is positive */
    if (b_inf >= 0) {
      /*interval b is positive*/
      *a_inf = b_inf * c_inf;
      *a_sup = b_sup * c_sup;
    } else if (b_sup <= 0) {
      /* interval b is negative */
      *a_inf = c_sup * b_inf;
      *a_sup = c_inf * b_sup;
    } else {
      /* there is 0 in between for b */
      *a_inf = b_inf * c_sup;
      *a_sup = b_sup * c_sup;
    }
  } else if (c_sup <= 0) {
    /* interval c is negative */
    if (b_inf >= 0) {
      /*interval b is positive*/
      *a_inf = b_sup * c_inf;
      *a_sup = b_inf * c_sup;
    } else if (b_sup <= 0) {
      /* interval b is negative */
      *a_inf = b_sup * c_sup;
      *a_sup = b_inf * c_inf;
    } else {
      /* there is 0 in between for b */
      *a_inf = b_sup * c_inf;
      *a_sup = b_inf * c_inf;
    }
  }
  /* there is 0 in between for c */
  else if (b_inf >= 0) {
    /*interval b is positive*/
    *a_inf = b_sup * c_inf;
    *a_sup = b_sup * c_sup;
  } else if (b_sup <= 0) {
    /* interval b is negative */
    *a_inf = b_inf * c_sup;
    *a_sup = b_inf * c_inf;
  } else {
    /* there is 0 in between for both b and c */
    float_type tmp_inf1 = b_sup * c_inf;
    float_type tmp_sup1 = b_inf * c_inf;
    float_type tmp_inf2 = b_inf * c_sup;
    float_type tmp_sup2 = b_sup * c_sup;
    *a_inf = min(tmp_inf1, tmp_inf2);
    *a_sup = max(tmp_sup1, tmp_sup2);
  }
}

__device__ void
elina_double_interval_mul2(float_type *const a_inf, float_type *const a_sup,
                           const float_type b_inf, const float_type b_sup,
                           const float_type c_inf, const float_type c_sup) {
  float_type inf_inf = b_inf * c_inf;
  float_type inf_sup = b_inf * c_sup;
  float_type sup_inf = b_sup * c_inf;
  float_type sup_sup = b_sup * c_sup;

  if (inf_inf < inf_sup) {
    if (sup_inf < sup_sup) {
      if (inf_inf < sup_inf) {
        *a_inf = inf_inf;
      } else {
        *a_inf = sup_inf;
      }

      if (inf_sup < sup_sup) {
        *a_sup = sup_sup;
      } else {
        *a_sup = inf_sup;
      }
    } else {
      if (inf_inf < sup_sup) {
        *a_inf = inf_inf;
      } else {
        *a_inf = sup_sup;
      }

      if (inf_sup < sup_inf) {
        *a_sup = sup_inf;
      } else {
        *a_sup = inf_sup;
      }
    }
  } else {
    if (sup_inf < sup_sup) {
      if (inf_sup < sup_inf) {
        *a_inf = inf_sup;
      } else {
        *a_inf = sup_inf;
      }

      if (inf_inf < sup_sup) {
        *a_sup = sup_sup;
      } else {
        *a_sup = inf_inf;
      }
    } else {
      if (inf_sup < sup_sup) {
        *a_inf = inf_sup;
      } else {
        *a_inf = sup_sup;
      }

      if (inf_inf < sup_inf) {
        *a_sup = sup_inf;
      } else {
        *a_sup = inf_inf;
      }
    }
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
  std::cout << std::endl << "This is the GPU version of fppoly!" << std::endl;
  results_calculated = false;
  output_counter = 1;

  void **funptr;
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

void fppoly_from_network_input_box(fppoly_t *const res, const size_t intdim,
                                   const size_t realdim,
                                   const double *inf_array,
                                   const double *sup_array) {
  res->layers = nullptr;
  res->numlayers = 0;

  size_t num_pixels = intdim + realdim;

  float_type *tmp_input_inf =
      (float_type *)malloc(num_pixels * sizeof(float_type));
  float_type *tmp_input_sup =
      (float_type *)malloc(num_pixels * sizeof(float_type));

  for (size_t i = 0; i < num_pixels; i++) {
    tmp_input_inf[i] = inf_array[i];
    tmp_input_sup[i] = sup_array[i];
  }

  cudaMalloc((void **)&(res->input_inf), num_pixels * sizeof(float_type));
  cudaMalloc((void **)&(res->input_sup), num_pixels * sizeof(float_type));

  cudaMemcpy(res->input_inf, tmp_input_inf, num_pixels * sizeof(float_type),
             cudaMemcpyHostToDevice);
  cudaMemcpy(res->input_sup, tmp_input_sup, num_pixels * sizeof(float_type),
             cudaMemcpyHostToDevice);

  free(tmp_input_inf);
  free(tmp_input_sup);

  res->num_pixels = num_pixels;
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

void ffn_add_layer(fppoly_t *const fp, const size_t num_out_neurons,
                   const size_t num_in_neurons, const layertype_t type,
                   const activation_type_t activation) {
  layer_t *layer = (layer_t *)malloc(sizeof(layer_t));

  layer->num_out_neurons = num_out_neurons;
  layer->num_in_neurons = num_in_neurons;

  layer->type = type;
  layer->activation = activation;

  cudaMalloc((void **)&layer->lb_array, num_out_neurons * sizeof(float_type));
  cudaMalloc((void **)&layer->ub_array, num_out_neurons * sizeof(float_type));

  cudaMalloc((void **)&layer->coeffs,
             num_out_neurons * num_in_neurons * sizeof(float_type));
  cudaMalloc((void **)&layer->csts, num_out_neurons * sizeof(float_type));

  layer->filter_weights = nullptr;
  layer->filter_bias = nullptr;

  layer->input_size = nullptr;
  layer->output_size = nullptr;
  layer->filter_size = nullptr;
  layer->strides = nullptr;
  layer->pad = nullptr;

  fp->layers[fp->numlayers] = layer;

  fp->numlayers++;
}

__device__ void elina_double_interval_mul_expr_coeff(
    float_type *const res_inf, float_type *const res_sup, const float_type inf,
    const float_type sup, const float_type inf_expr,
    const float_type sup_expr) {
  elina_double_interval_mul2(res_inf, res_sup, inf, sup, inf_expr, sup_expr);

  const float_type maxA = max(fabs(inf_expr), fabs(sup_expr));
  float_type tmp1, tmp2;

  elina_double_interval_mul2(&tmp1, &tmp2, inf, sup, -maxA * ulp, maxA * ulp);

  *res_inf += tmp1;
  *res_sup += tmp2;
}

__device__ void elina_double_interval_mul_cst_coeff(float_type *const res_inf,
                                                    float_type *const res_sup,
                                                    const float_type inf,
                                                    const float_type sup,
                                                    const float_type inf_expr,
                                                    const float_type sup_expr) {
  elina_double_interval_mul_expr_coeff(res_inf, res_sup, inf, sup, inf_expr,
                                       sup_expr);

  *res_inf -= min_denormal;
  *res_sup += min_denormal;
}

__global__ void compute_lb_from_expr(float_type *__restrict__ lb_array,
                                     const float_type *__restrict__ inf_coeff,
                                     const float_type *__restrict__ sup_coeff,
                                     const float_type *__restrict__ inf_cst,
                                     const float_type *__restrict__ input_inf,
                                     const float_type *__restrict__ input_sup,
                                     const size_t expr_size) {
  const size_t n = blockIdx.x;

  float_type res_inf = inf_cst[n];

  float_type tmp1, tmp2;

  for (size_t i = 0; i < expr_size; i++) {
    elina_double_interval_mul2(&tmp1, &tmp2, inf_coeff[n * expr_size + i],
                               sup_coeff[n * expr_size + i], input_inf[i],
                               input_sup[i]);

    res_inf = res_inf + tmp1;
  }

  lb_array[n] = res_inf;
}

__global__ void compute_ub_from_expr(float_type *__restrict__ ub_array,
                                     const float_type *__restrict__ inf_coeff,
                                     const float_type *__restrict__ sup_coeff,
                                     const float_type *__restrict__ sup_cst,
                                     const float_type *__restrict__ input_inf,
                                     const float_type *__restrict__ input_sup,
                                     const size_t expr_size) {
  const size_t n = blockIdx.x;

  float_type res_sup = sup_cst[n];

  float_type tmp1, tmp2;

  for (size_t i = 0; i < expr_size; i++) {
    elina_double_interval_mul2(&tmp1, &tmp2, inf_coeff[n * expr_size + i],
                               sup_coeff[n * expr_size + i], input_inf[i],
                               input_sup[i]);

    res_sup = res_sup + tmp2;
  }

  ub_array[n] = res_sup;
}

__global__ void compute_lb_from_expr_conv_sparse(
    float_type *__restrict__ lb_array, const float_type *__restrict__ inf_coeff,
    const float_type *__restrict__ sup_coeff,
    const float_type *__restrict__ inf_cst,
    const float_type *__restrict__ input_inf,
    const float_type *__restrict__ input_sup, const size_t num_chunks,
    const size_t chunk_counter, const size_t output_size_x,
    const size_t output_size_y, const size_t output_size_z, long int offset_x,
    long int offset_y, long int length_x, long int length_y, long int shift_x,
    long int shift_y) {
  const size_t last_x = blockIdx.x;
  const size_t last_y = blockIdx.y;
  const size_t last_z = blockIdx.z;

  const size_t local_n =
      last_x * gridDim.y * gridDim.z + last_y * gridDim.z + last_z;
  const size_t global_n = last_x * gridDim.y * num_chunks * gridDim.z +
                          last_y * num_chunks * gridDim.z +
                          chunk_counter * gridDim.z + last_z;

  const long int min_out_x = offset_x + last_x * shift_x;
  const long int min_out_y = offset_y + last_y * shift_y;

  float_type res_inf = inf_cst[local_n];

  float_type tmp1, tmp2;

  for (long int out_x = 0; out_x < length_x; out_x++) {
    if ((out_x + min_out_x < 0) || (out_x + min_out_x >= output_size_x)) {
      continue;
    }

    for (long int out_y = 0; out_y < length_y; out_y++) {
      if ((out_y + min_out_y < 0) || (out_y + min_out_y >= output_size_y)) {
        continue;
      }

      for (size_t out_z = 0; out_z < output_size_z; out_z++) {
        size_t i = (out_x + min_out_x) * output_size_y * output_size_z +
                   (out_y + min_out_y) * output_size_z + out_z;
        size_t j =
            out_x * length_y * output_size_z + out_y * output_size_z + out_z;

        size_t mat_out = local_n * length_x * length_y * output_size_z + j;

        elina_double_interval_mul2(&tmp1, &tmp2, inf_coeff[mat_out],
                                   sup_coeff[mat_out], input_inf[i],
                                   input_sup[i]);

        res_inf = res_inf + tmp1;
      }
    }
  }

  lb_array[global_n] = res_inf;
}

__global__ void compute_ub_from_expr_conv_sparse(
    float_type *__restrict__ ub_array, const float_type *__restrict__ inf_coeff,
    const float_type *__restrict__ sup_coeff,
    const float_type *__restrict__ sup_cst,
    const float_type *__restrict__ input_inf,
    const float_type *__restrict__ input_sup, const size_t num_chunks,
    const size_t chunk_counter, const size_t output_size_x,
    const size_t output_size_y, const size_t output_size_z, long int offset_x,
    long int offset_y, long int length_x, long int length_y, long int shift_x,
    long int shift_y) {
  const size_t last_x = blockIdx.x;
  const size_t last_y = blockIdx.y;
  const size_t last_z = blockIdx.z;

  const size_t local_n =
      last_x * gridDim.y * gridDim.z + last_y * gridDim.z + last_z;
  const size_t global_n = last_x * gridDim.y * num_chunks * gridDim.z +
                          last_y * num_chunks * gridDim.z +
                          chunk_counter * gridDim.z + last_z;
  const long int min_out_x = offset_x + last_x * shift_x;
  const long int min_out_y = offset_y + last_y * shift_y;

  float_type res_sup = sup_cst[local_n];

  float_type tmp1, tmp2;

  for (long int out_x = 0; out_x < length_x; out_x++) {
    if ((out_x + min_out_x < 0) || (out_x + min_out_x >= output_size_x)) {
      continue;
    }

    for (long int out_y = 0; out_y < length_y; out_y++) {
      if ((out_y + min_out_y < 0) || (out_y + min_out_y >= output_size_y)) {
        continue;
      }

      for (size_t out_z = 0; out_z < output_size_z; out_z++) {
        size_t i = (out_x + min_out_x) * output_size_y * output_size_z +
                   (out_y + min_out_y) * output_size_z + out_z;
        size_t j =
            out_x * length_y * output_size_z + out_y * output_size_z + out_z;

        size_t mat_out = local_n * length_x * length_y * output_size_z + j;

        elina_double_interval_mul2(&tmp1, &tmp2, inf_coeff[mat_out],
                                   sup_coeff[mat_out], input_inf[i],
                                   input_sup[i]);

        res_sup = res_sup + tmp2;
      }
    }
  }

  ub_array[global_n] = res_sup;
}

template <typename FloatType>
__global__ void device_layer_create_dense_expr(
    float_type *__restrict__ coeffs, float_type *__restrict__ csts,
    const FloatType *__restrict__ weights, const FloatType *__restrict__ bias,
    const size_t num_out_neurons, const size_t num_in_neurons) {
  const size_t i = blockIdx.x;

  const FloatType *weight_i = weights + i * num_in_neurons;
  const FloatType bias_i = bias[i];

  csts[i] = bias_i;

  for (size_t j = 0; j < num_in_neurons; j++) {
    coeffs[i * num_in_neurons + j] = weight_i[j];
  }
}

void layer_create_dense_exprs(float_type *coeffs, float_type *csts,
                              const double **weights, const double *bias,
                              const size_t num_out_neurons,
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
      coeffs, csts, tmp_weights, tmp_bias, num_out_neurons, num_in_neurons);

  cudaFree(tmp_weights);
  cudaFree(tmp_bias);
}

__global__ void copy_expr_array(float_type *__restrict__ target_inf_coeff,
                                float_type *__restrict__ target_sup_coeff,
                                float_type *__restrict__ target_inf_cst,
                                float_type *__restrict__ target_sup_cst,
                                const float_type *__restrict__ source_coeffs,
                                const float_type *__restrict__ source_csts,
                                const size_t num_exprs,
                                const size_t expr_size) {
  const size_t i = blockIdx.x;

  for (size_t j = 0; j < expr_size; j++) {
    target_inf_coeff[i * expr_size + j] = source_coeffs[i * expr_size + j];
    target_sup_coeff[i * expr_size + j] = source_coeffs[i * expr_size + j];
  }

  target_inf_cst[i] = source_csts[i];
  target_sup_cst[i] = source_csts[i];
}

__global__ void layer_compute_bounds_from_exprs_conv(
    const float_type *__restrict__ coeffs, const float_type *__restrict__ csts,
    float_type *__restrict__ lb_array, float_type *__restrict__ ub_array,
    const float_type *__restrict__ input_inf,
    const float_type *__restrict__ input_sup, const size_t output_size_x,
    const size_t output_size_y, const size_t output_size_z,
    const size_t input_size_x, const size_t input_size_y,
    const size_t input_size_z, const size_t filter_size_x,
    const size_t filter_size_y, const size_t stride_x, const size_t stride_y,
    const long int pad_x, const long int pad_y) {
  size_t out_x = blockIdx.x;
  size_t out_y = blockIdx.y;
  size_t out_z = blockIdx.z;

  const size_t mat_out =
      out_x * output_size_y * output_size_z + out_y * output_size_z + out_z;

  float_type res_inf = csts[out_z];
  float_type res_sup = csts[out_z];

  float_type tmp1, tmp2;

  for (size_t x_shift = 0; x_shift < filter_size_x; x_shift++) {
    for (size_t y_shift = 0; y_shift < filter_size_y; y_shift++) {
      for (size_t inp_z = 0; inp_z < input_size_z; inp_z++) {
        const long int x_val = out_x * stride_x + x_shift - pad_x;
        const long int y_val = out_y * stride_y + y_shift - pad_y;

        if ((y_val < 0) || (y_val >= (long int)input_size_y)) {
          continue;
        }

        if ((x_val < 0) || (x_val >= (long int)input_size_x)) {
          continue;
        }

        const size_t mat_in =
            x_val * input_size_y * input_size_z + y_val * input_size_z + inp_z;

        const size_t filter_index =
            out_z * filter_size_x * filter_size_y * input_size_z +
            x_shift * filter_size_y * input_size_z + y_shift * input_size_z +
            inp_z;

        elina_double_interval_mul2(&tmp1, &tmp2, coeffs[filter_index],
                                   coeffs[filter_index], input_inf[mat_in],
                                   input_sup[mat_in]);

        res_inf = res_inf + tmp1;
        res_sup = res_sup + tmp2;
      }
    }
  }

  lb_array[mat_out] = res_inf;
  ub_array[mat_out] = res_sup;
}

__global__ void layer_compute_bounds_from_exprs(
    const float_type *__restrict__ coeffs, const float_type *__restrict__ csts,
    float_type *lb_array, float_type *ub_array,
    const float_type *__restrict__ input_inf,
    const float_type *__restrict__ input_sup, const size_t num_out_neurons,
    const size_t num_in_neurons) {
  const size_t n = blockIdx.x;

  float_type res_inf = csts[n];
  float_type res_sup = csts[n];

  float_type tmp1, tmp2;

  for (size_t i = 0; i < num_in_neurons; i++) {
    elina_double_interval_mul2(&tmp1, &tmp2, coeffs[n * num_in_neurons + i],
                               coeffs[n * num_in_neurons + i], input_inf[i],
                               input_sup[i]);

    res_inf = res_inf + tmp1;
    res_sup = res_sup + tmp2;
  }

  lb_array[n] = res_inf;
  ub_array[n] = res_sup;
}

void ffn_handle_first_layer(elina_manager_t *man, elina_abstract0_t *abs,
                            const double **weights, const double *bias,
                            const size_t size, const size_t num_pixels,
                            const activation_type_t activation) {
  fppoly_t *res = fppoly_of_abstract0(abs);
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  res->layers = (layer_t **)malloc(2000 * sizeof(layer_t *));
  ffn_add_layer(res, size, num_pixels, FFN, activation);

  float_type *coeffs = res->layers[0]->coeffs;
  float_type *csts = res->layers[0]->csts;

  layer_create_dense_exprs(coeffs, csts, weights, bias, size, num_pixels);
  layer_compute_bounds_from_exprs<<<res->layers[0]->num_out_neurons, 1>>>(
      coeffs, csts, res->layers[0]->lb_array, res->layers[0]->ub_array,
      res->input_inf, res->input_sup, res->layers[0]->num_out_neurons,
      res->layers[0]->num_in_neurons);
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

// TODO: Assess if non-determinism introduced by atomics is problematic for
// analyzer.
__global__ void lexpr_replace_relu_bounds(
    float_type *__restrict__ inf_coeff, float_type *__restrict__ sup_coeff,
    float_type *__restrict__ inf_cst, float_type *__restrict__ sup_cst,
    const float_type *__restrict__ lb_array,
    const float_type *__restrict__ ub_array,
    const size_t num_out_neurons_current_layer) {
  const size_t n = blockIdx.x;
  const size_t i = blockIdx.y * blockDim.x + threadIdx.x;

  if (i < num_out_neurons_current_layer) {
    const size_t a = n * num_out_neurons_current_layer + i;

    const float_type lb = lb_array[i];
    const float_type ub = ub_array[i];
    const float_type width = ub - lb;
    const float_type lambda_inf = ub / width;
    const float_type lambda_sup = ub / width;

    const float_type old_inf_coeff = inf_coeff[a];
    const float_type old_sup_coeff = sup_coeff[a];

    if ((old_sup_coeff == 0) && (old_inf_coeff == 0)) {
      inf_coeff[a] = 0.0;
      sup_coeff[a] = 0.0;

      return;
    } else if (ub <= 0) {
      inf_coeff[a] = 0.0;
      sup_coeff[a] = 0.0;

      return;
    } else if (lb > 0) {
      inf_coeff[a] = old_inf_coeff;
      sup_coeff[a] = old_sup_coeff;
    } else if (old_sup_coeff < 0) {
      const float_type mu_inf = -lambda_inf * lb;
      const float_type mu_sup = -lambda_sup * lb;
      elina_double_interval_mul_expr_coeff(&inf_coeff[a], &sup_coeff[a],
                                           lambda_inf, lambda_sup,
                                           old_inf_coeff, old_sup_coeff);
      float_type tmp1, tmp2;
      elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, mu_inf, mu_sup,
                                          old_inf_coeff, old_sup_coeff);

      atomicAdd(&inf_cst[n], tmp1 - min_denormal);
      atomicAdd(&sup_cst[n], tmp2 + min_denormal);
    } else if (old_inf_coeff > 0) {
      const float_type area1 = -lb * ub;
      const float_type area2 = 0.5 * ub * width;
      const float_type area3 = -0.5 * lb * width;

      if ((area1 < area2) && (area1 < area3)) {
        elina_double_interval_mul_expr_coeff(&inf_coeff[a], &sup_coeff[a],
                                             lambda_inf, lambda_sup,
                                             old_inf_coeff, old_sup_coeff);
      } else if ((area2 < area1) && (area2 < area3)) {
        inf_coeff[a] = 0.0;
        sup_coeff[a] = 0.0;
      } else {
        inf_coeff[a] = old_inf_coeff;
        sup_coeff[a] = old_sup_coeff;
      }
    } else {
      inf_coeff[a] = 0.0;
      sup_coeff[a] = 0.0;
      float_type tmp1, tmp2;
      elina_double_interval_mul2(&tmp1, &tmp2, old_inf_coeff, old_sup_coeff, 0,
                                 ub);

      atomicAdd(&inf_cst[n], tmp1);
      atomicAdd(&sup_cst[n], tmp1);
    }
  }
}

__global__ void uexpr_replace_relu_bounds(
    float_type *__restrict__ inf_coeff, float_type *__restrict__ sup_coeff,
    float_type *__restrict__ inf_cst, float_type *__restrict__ sup_cst,
    const float_type *__restrict__ lb_array,
    const float_type *__restrict__ ub_array,
    const size_t num_out_neurons_current_layer) {
  const size_t n = blockIdx.x;
  const size_t i = blockIdx.y * blockDim.x + threadIdx.x;

  if (i < num_out_neurons_current_layer) {
    const size_t a = n * num_out_neurons_current_layer + i;

    const float_type lb = lb_array[i];
    const float_type ub = ub_array[i];
    const float_type width = ub - lb;
    const float_type lambda_inf = ub / width;
    const float_type lambda_sup = ub / width;

    const float_type old_inf_coeff = inf_coeff[a];
    const float_type old_sup_coeff = sup_coeff[a];

    if ((old_sup_coeff == 0) && (old_inf_coeff == 0)) {
      inf_coeff[a] = 0.0;
      sup_coeff[a] = 0.0;

      return;
    } else if (ub <= 0) {
      inf_coeff[a] = 0.0;
      sup_coeff[a] = 0.0;

      return;
    } else if (lb > 0) {
      inf_coeff[a] = old_inf_coeff;
      sup_coeff[a] = old_sup_coeff;
    } else if (old_inf_coeff > 0) {
      const float_type mu_inf = -lambda_inf * lb;
      const float_type mu_sup = -lambda_sup * lb;
      elina_double_interval_mul_expr_coeff(&inf_coeff[a], &sup_coeff[a],
                                           lambda_inf, lambda_sup,
                                           old_inf_coeff, old_sup_coeff);
      float_type tmp1, tmp2;
      elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, mu_inf, mu_sup,
                                          old_inf_coeff, old_sup_coeff);

      atomicAdd(&inf_cst[n], tmp1 - min_denormal);
      atomicAdd(&sup_cst[n], tmp2 + min_denormal);
    } else if (old_sup_coeff < 0) {
      const float_type area1 = -lb * ub;
      const float_type area2 = 0.5 * ub * width;
      const float_type area3 = -0.5 * lb * width;

      if ((area1 < area2) && (area1 < area3)) {
        elina_double_interval_mul_expr_coeff(&inf_coeff[a], &sup_coeff[a],
                                             lambda_inf, lambda_sup,
                                             old_inf_coeff, old_sup_coeff);
      } else if ((area2 < area1) && (area2 < area3)) {
        inf_coeff[a] = 0.0;
        sup_coeff[a] = 0.0;
      } else {
        inf_coeff[a] = old_inf_coeff;
        sup_coeff[a] = old_sup_coeff;
      }
    } else {
      inf_coeff[a] = 0.0;
      sup_coeff[a] = 0.0;
      float_type tmp1, tmp2;
      elina_double_interval_mul2(&tmp1, &tmp2, old_inf_coeff, old_sup_coeff, 0,
                                 ub);

      atomicAdd(&inf_cst[n], tmp2);
      atomicAdd(&sup_cst[n], tmp2);
    }
  }
}

// TODO: Assess if non-determinism introduced by atomics is problematic for
// analyzer.
__global__ void lexpr_replace_relu_bounds_conv_sparse(
    float_type *__restrict__ inf_coeff, float_type *__restrict__ sup_coeff,
    float_type *__restrict__ inf_cst, float_type *__restrict__ sup_cst,
    const float_type *__restrict__ lb_array,
    const float_type *__restrict__ ub_array, const size_t output_size_x,
    const size_t output_size_y, const size_t output_size_z, long int offset_x,
    long int offset_y, long int length_x, long int length_y, long int shift_x,
    long int shift_y) {
  const size_t last_x = blockIdx.x;
  const size_t last_y = blockIdx.y;
  const size_t last_z = blockIdx.z;

  const size_t n = last_x * gridDim.y * gridDim.z + last_y * gridDim.z + last_z;

  const long int min_out_x = offset_x + last_x * shift_x;
  const long int min_out_y = offset_y + last_y * shift_y;

  for (long int out_x = 0; out_x < length_x; out_x++) {
    if ((out_x + min_out_x < 0) || (out_x + min_out_x >= output_size_x)) {
      continue;
    }

    for (long int out_y = 0; out_y < length_y; out_y++) {
      if ((out_y + min_out_y < 0) || (out_y + min_out_y >= output_size_y)) {
        continue;
      }

      for (size_t out_z = 0; out_z < output_size_z; out_z++) {
        size_t i = (out_x + min_out_x) * output_size_y * output_size_z +
                   (out_y + min_out_y) * output_size_z + out_z;
        size_t j =
            out_x * length_y * output_size_z + out_y * output_size_z + out_z;

        size_t a = n * length_x * length_y * output_size_z + j;

        const float_type lb = lb_array[i];
        const float_type ub = ub_array[i];
        const float_type width = ub - lb;
        const float_type lambda_inf = ub / width;
        const float_type lambda_sup = ub / width;

        const float_type old_inf_coeff = inf_coeff[a];
        const float_type old_sup_coeff = sup_coeff[a];

        if ((old_sup_coeff == 0) && (old_inf_coeff == 0)) {
          inf_coeff[a] = 0.0;
          sup_coeff[a] = 0.0;

          continue;
        } else if (ub <= 0) {
          inf_coeff[a] = 0.0;
          sup_coeff[a] = 0.0;

          continue;
        } else if (lb > 0) {
          inf_coeff[a] = old_inf_coeff;
          sup_coeff[a] = old_sup_coeff;
        } else if (old_sup_coeff < 0) {
          const float_type mu_inf = -lambda_inf * lb;
          const float_type mu_sup = -lambda_sup * lb;
          elina_double_interval_mul_expr_coeff(&inf_coeff[a], &sup_coeff[a],
                                               lambda_inf, lambda_sup,
                                               old_inf_coeff, old_sup_coeff);
          float_type tmp1, tmp2;
          elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, mu_inf, mu_sup,
                                              old_inf_coeff, old_sup_coeff);

          inf_cst[n] += tmp1 - min_denormal;
          sup_cst[n] += tmp2 + min_denormal;
        } else if (old_inf_coeff > 0) {
          const float_type area1 = -lb * ub;
          const float_type area2 = 0.5 * ub * width;
          const float_type area3 = -0.5 * lb * width;

          if ((area1 < area2) && (area1 < area3)) {
            elina_double_interval_mul_expr_coeff(&inf_coeff[a], &sup_coeff[a],
                                                 lambda_inf, lambda_sup,
                                                 old_inf_coeff, old_sup_coeff);
          } else if ((area2 < area1) && (area2 < area3)) {
            inf_coeff[a] = 0.0;
            sup_coeff[a] = 0.0;
          } else {
            inf_coeff[a] = old_inf_coeff;
            sup_coeff[a] = old_sup_coeff;
          }
        } else {
          inf_coeff[a] = 0.0;
          sup_coeff[a] = 0.0;
          float_type tmp1, tmp2;
          elina_double_interval_mul2(&tmp1, &tmp2, old_inf_coeff, old_sup_coeff,
                                     0, ub);

          inf_cst[n] += tmp1;
          sup_cst[n] += tmp1;
        }
      }
    }
  }
}

__global__ void uexpr_replace_relu_bounds_conv_sparse(
    float_type *__restrict__ inf_coeff, float_type *__restrict__ sup_coeff,
    float_type *__restrict__ inf_cst, float_type *__restrict__ sup_cst,
    const float_type *__restrict__ lb_array,
    const float_type *__restrict__ ub_array, const size_t output_size_x,
    const size_t output_size_y, const size_t output_size_z, long int offset_x,
    long int offset_y, long int length_x, long int length_y, long int shift_x,
    long int shift_y) {
  const size_t last_x = blockIdx.x;
  const size_t last_y = blockIdx.y;
  const size_t last_z = blockIdx.z;

  const size_t n = last_x * gridDim.y * gridDim.z + last_y * gridDim.z + last_z;

  const long int min_out_x = offset_x + last_x * shift_x;
  const long int min_out_y = offset_y + last_y * shift_y;

  for (long int out_x = 0; out_x < length_x; out_x++) {
    if ((out_x + min_out_x < 0) || (out_x + min_out_x >= output_size_x)) {
      continue;
    }

    for (long int out_y = 0; out_y < length_y; out_y++) {
      if ((out_y + min_out_y < 0) || (out_y + min_out_y >= output_size_y)) {
        continue;
      }

      for (size_t out_z = 0; out_z < output_size_z; out_z++) {
        size_t i = (out_x + min_out_x) * output_size_y * output_size_z +
                   (out_y + min_out_y) * output_size_z + out_z;
        size_t j =
            out_x * length_y * output_size_z + out_y * output_size_z + out_z;

        size_t a = n * length_x * length_y * output_size_z + j;

        const float_type lb = lb_array[i];
        const float_type ub = ub_array[i];
        const float_type width = ub - lb;
        const float_type lambda_inf = ub / width;
        const float_type lambda_sup = ub / width;

        const float_type old_inf_coeff = inf_coeff[a];
        const float_type old_sup_coeff = sup_coeff[a];

        if ((old_sup_coeff == 0) && (old_inf_coeff == 0)) {
          inf_coeff[a] = 0.0;
          sup_coeff[a] = 0.0;

          continue;
        } else if (ub <= 0) {
          inf_coeff[a] = 0.0;
          sup_coeff[a] = 0.0;

          continue;
        } else if (lb > 0) {
          inf_coeff[a] = old_inf_coeff;
          sup_coeff[a] = old_sup_coeff;
        } else if (old_inf_coeff > 0) {
          const float_type mu_inf = -lambda_inf * lb;
          const float_type mu_sup = -lambda_sup * lb;
          elina_double_interval_mul_expr_coeff(&inf_coeff[a], &sup_coeff[a],
                                               lambda_inf, lambda_sup,
                                               old_inf_coeff, old_sup_coeff);
          float_type tmp1, tmp2;
          elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, mu_inf, mu_sup,
                                              old_inf_coeff, old_sup_coeff);

          inf_cst[n] += tmp1 - min_denormal;
          sup_cst[n] += tmp2 + min_denormal;
        } else if (old_sup_coeff < 0) {
          const float_type area1 = -lb * ub;
          const float_type area2 = 0.5 * ub * width;
          const float_type area3 = -0.5 * lb * width;

          if ((area1 < area2) && (area1 < area3)) {
            elina_double_interval_mul_expr_coeff(&inf_coeff[a], &sup_coeff[a],
                                                 lambda_inf, lambda_sup,
                                                 old_inf_coeff, old_sup_coeff);
          } else if ((area2 < area1) && (area2 < area3)) {
            inf_coeff[a] = 0.0;
            sup_coeff[a] = 0.0;
          } else {
            inf_coeff[a] = old_inf_coeff;
            sup_coeff[a] = old_sup_coeff;
          }
        } else {
          inf_coeff[a] = 0.0;
          sup_coeff[a] = 0.0;
          float_type tmp1, tmp2;
          elina_double_interval_mul2(&tmp1, &tmp2, old_inf_coeff, old_sup_coeff,
                                     0, ub);

          inf_cst[n] += tmp2;
          sup_cst[n] += tmp2;
        }
      }
    }
  }
}

__global__ void
coeffs_from_previous_layer(const float_type *__restrict__ expr_inf_coeff,
                           const float_type *__restrict__ expr_sup_coeff,
                           float_type *__restrict__ res_inf_coeff,
                           float_type *__restrict__ res_sup_coeff,
                           const float_type *__restrict__ aux_coeffs,
                           const size_t num_out_neurons_current_layer,
                           const size_t num_in_neurons_current_layer) {
  const size_t n = blockIdx.x;
  const size_t j = blockIdx.y * blockDim.x + threadIdx.x;

  if (j < num_in_neurons_current_layer) {
    const size_t b = n * num_in_neurons_current_layer + j;

    float_type inf_coeff = 0;
    float_type sup_coeff = 0;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    for (size_t i = 0; i < num_out_neurons_current_layer; i++) {
      size_t a = n * num_out_neurons_current_layer + i;
      size_t c = i * num_in_neurons_current_layer + j;

      const float_type prev_inf_coeff = expr_inf_coeff[a];
      const float_type prev_sup_coeff = expr_sup_coeff[a];

      if ((prev_inf_coeff != 0) || (prev_sup_coeff != 0)) {
        elina_double_interval_mul_expr_coeff(&tmp1, &tmp2, prev_inf_coeff,
                                             prev_sup_coeff, aux_coeffs[c],
                                             aux_coeffs[c]);

        maxRes = max(fabs(inf_coeff), fabs(sup_coeff));
        maxMul = max(fabs(tmp1), fabs(tmp2));

        inf_coeff = inf_coeff + tmp1 - (maxRes + maxMul) * ulp;
        sup_coeff = sup_coeff + tmp2 + (maxRes + maxMul) * ulp;
      }
    }

    res_inf_coeff[b] = inf_coeff;
    res_sup_coeff[b] = sup_coeff;
  }
}

__global__ void coeffs_from_previous_layer_conv(
    const float_type *__restrict__ expr_inf_coeff,
    const float_type *__restrict__ expr_sup_coeff,
    float_type *__restrict__ res_inf_coeff,
    float_type *__restrict__ res_sup_coeff,
    const float_type *__restrict__ aux_coeffs, const size_t output_size_x,
    const size_t output_size_y, const size_t output_size_z,
    const size_t input_size_x, const size_t input_size_y,
    const size_t input_size_z, const size_t filter_size_x,
    const size_t filter_size_y, const size_t stride_x, const size_t stride_y,
    const size_t pad_x, const size_t pad_y) {
  const size_t n = blockIdx.x;

  const size_t inp_z = threadIdx.x;
  const size_t y_shift = threadIdx.y;
  const size_t x_shift = threadIdx.z;

  float_type tmp1, tmp2;
  float_type maxRes, maxMul;

  if (x_shift * filter_size_y * input_size_z + y_shift * input_size_z + inp_z <
      filter_size_x * filter_size_y * input_size_z) {
    for (size_t out_x = 0; out_x < output_size_x; out_x++) {
      for (size_t out_y = 0; out_y < output_size_y; out_y++) {
        const long int x_val = out_x * stride_x + x_shift - pad_x;
        const long int y_val = out_y * stride_y + y_shift - pad_y;

        if (!((y_val < 0) || (y_val >= (long int)input_size_y))) {
          if (!((x_val < 0) || (x_val >= (long int)input_size_x))) {
            const size_t mat_in = x_val * input_size_y * input_size_z +
                                  y_val * input_size_z + inp_z;
            const size_t b =
                n * input_size_x * input_size_y * input_size_z + mat_in;

            float_type inf_coeff = res_inf_coeff[b];
            float_type sup_coeff = res_sup_coeff[b];

            for (size_t out_z = 0; out_z < output_size_z; out_z++) {
              const size_t mat_out = out_x * output_size_y * output_size_z +
                                     out_y * output_size_z + out_z;

              const size_t a =
                  n * output_size_x * output_size_y * output_size_z + mat_out;

              const float_type prev_inf_coeff = expr_inf_coeff[a];
              const float_type prev_sup_coeff = expr_sup_coeff[a];

              if ((prev_inf_coeff != 0) || (prev_sup_coeff != 0)) {
                const size_t filter_index =
                    out_z * filter_size_x * filter_size_y * input_size_z +
                    x_shift * filter_size_y * input_size_z +
                    y_shift * input_size_z + inp_z;

                const float_type aux_coeff = aux_coeffs[filter_index];
                elina_double_interval_mul_expr_coeff(
                    &tmp1, &tmp2, prev_inf_coeff, prev_sup_coeff, aux_coeff,
                    aux_coeff);

                maxRes = max(fabs(inf_coeff), fabs(sup_coeff));
                maxMul = max(fabs(tmp1), fabs(tmp2));

                inf_coeff = inf_coeff + tmp1 - (maxRes + maxMul) * ulp;
                sup_coeff = sup_coeff + tmp2 + (maxRes + maxMul) * ulp;
              }
            }

            res_inf_coeff[b] = inf_coeff;
            res_sup_coeff[b] = sup_coeff;
          }
        }

        __syncthreads();
      }
    }
  }
}

__global__ void coeffs_from_previous_layer_conv_filter_serial(
    const float_type *__restrict__ expr_inf_coeff,
    const float_type *__restrict__ expr_sup_coeff,
    float_type *__restrict__ res_inf_coeff,
    float_type *__restrict__ res_sup_coeff,
    const float_type *__restrict__ aux_coeffs, const size_t output_size_x,
    const size_t output_size_y, const size_t output_size_z,
    const size_t input_size_x, const size_t input_size_y,
    const size_t input_size_z, const size_t filter_size_x,
    const size_t filter_size_y, const size_t stride_x, const size_t stride_y,
    const size_t pad_x, const size_t pad_y) {
  const size_t n = blockIdx.x;

  const size_t inp_z = threadIdx.x;

  float_type tmp1, tmp2;
  float_type maxRes, maxMul;

  if (inp_z < input_size_z) {
    for (size_t out_x = 0; out_x < output_size_x; out_x++) {
      for (size_t out_y = 0; out_y < output_size_y; out_y++) {
        for (size_t x_shift = 0; x_shift < filter_size_x; x_shift++) {
          for (size_t y_shift = 0; y_shift < filter_size_y; y_shift++) {
            const long int x_val = out_x * stride_x + x_shift - pad_x;
            const long int y_val = out_y * stride_y + y_shift - pad_y;

            if (!((y_val < 0) || (y_val >= (long int)input_size_y))) {
              if (!((x_val < 0) || (x_val >= (long int)input_size_x))) {
                const size_t mat_in = x_val * input_size_y * input_size_z +
                                      y_val * input_size_z + inp_z;
                const size_t b =
                    n * input_size_x * input_size_y * input_size_z + mat_in;

                float_type inf_coeff = res_inf_coeff[b];
                float_type sup_coeff = res_sup_coeff[b];

                for (size_t out_z = 0; out_z < output_size_z; out_z++) {
                  const size_t mat_out = out_x * output_size_y * output_size_z +
                                         out_y * output_size_z + out_z;

                  const size_t a =
                      n * output_size_x * output_size_y * output_size_z +
                      mat_out;

                  const float_type prev_inf_coeff = expr_inf_coeff[a];
                  const float_type prev_sup_coeff = expr_sup_coeff[a];

                  if ((prev_inf_coeff != 0) || (prev_sup_coeff != 0)) {
                    const size_t filter_index =
                        out_z * filter_size_x * filter_size_y * input_size_z +
                        x_shift * filter_size_y * input_size_z +
                        y_shift * input_size_z + inp_z;

                    const float_type aux_coeff = aux_coeffs[filter_index];
                    elina_double_interval_mul_expr_coeff(
                        &tmp1, &tmp2, prev_inf_coeff, prev_sup_coeff, aux_coeff,
                        aux_coeff);

                    maxRes = max(fabs(inf_coeff), fabs(sup_coeff));
                    maxMul = max(fabs(tmp1), fabs(tmp2));

                    inf_coeff = inf_coeff + tmp1 - (maxRes + maxMul) * ulp;
                    sup_coeff = sup_coeff + tmp2 + (maxRes + maxMul) * ulp;
                  }
                }

                res_inf_coeff[b] = inf_coeff;
                res_sup_coeff[b] = sup_coeff;
              }
            }
          }
        }
      }
    }
  }
}

__global__ void coeffs_from_previous_layer_conv_sparse(
    const float_type *__restrict__ expr_inf_coeff,
    const float_type *__restrict__ expr_sup_coeff,
    float_type *__restrict__ res_inf_coeff,
    float_type *__restrict__ res_sup_coeff,
    const float_type *__restrict__ aux_coeffs, const size_t output_size_x,
    const size_t output_size_y, const size_t output_size_z,
    const size_t input_size_x, const size_t input_size_y,
    const size_t input_size_z, long int offset_x, long int offset_y,
    long int length_x, long int length_y, long int shift_x, long int shift_y,
    const size_t filter_size_x, const size_t filter_size_y,
    const size_t stride_x, const size_t stride_y, const size_t pad_x,
    const size_t pad_y) {
  const size_t last_x = blockIdx.x;
  const size_t last_y = blockIdx.y;
  const size_t last_z = blockIdx.z;

  const size_t n = last_x * gridDim.y * gridDim.z + last_y * gridDim.z + last_z;

  const long int min_out_x = offset_x + last_x * shift_x;
  const long int min_out_y = offset_y + last_y * shift_y;

  const size_t inp_z = threadIdx.x;
  const size_t y_shift = threadIdx.y;
  const size_t x_shift = threadIdx.z;

  float_type tmp1, tmp2;
  float_type maxRes, maxMul;

  if (x_shift * filter_size_y * input_size_z + y_shift * input_size_z + inp_z <
      filter_size_x * filter_size_y * input_size_z) {
    for (long int out_x = 0; out_x < length_x; out_x++) {
      if ((out_x + min_out_x < 0) || (out_x + min_out_x >= output_size_x)) {
        continue;
      }

      for (long int out_y = 0; out_y < length_y; out_y++) {
        if ((out_y + min_out_y < 0) || (out_y + min_out_y >= output_size_y)) {
          continue;
        }

        const long int x_val =
            out_x * stride_x + min_out_x * stride_x + x_shift - pad_x;
        const long int y_val =
            out_y * stride_y + min_out_y * stride_y + y_shift - pad_y;

        if (!((x_val < 0) || (x_val >= (long int)input_size_x))) {
          if (!((y_val < 0) || (y_val >= (long int)input_size_y))) {
            const size_t mat_in =
                (out_x * stride_x + x_shift) *
                    ((length_y - 1) * stride_y + filter_size_y) * input_size_z +
                (out_y * stride_y + y_shift) * input_size_z + inp_z;
            const size_t b = n * ((length_x - 1) * stride_x + filter_size_x) *
                                 ((length_y - 1) * stride_y + filter_size_y) *
                                 input_size_z +
                             mat_in;

            float_type inf_coeff = res_inf_coeff[b];
            float_type sup_coeff = res_sup_coeff[b];

            for (size_t out_z = 0; out_z < output_size_z; out_z++) {
              const size_t mat_out = out_x * length_y * output_size_z +
                                     out_y * output_size_z + out_z;
              const size_t a =
                  n * length_x * length_y * output_size_z + mat_out;

              const float_type prev_inf_coeff = expr_inf_coeff[a];
              const float_type prev_sup_coeff = expr_sup_coeff[a];

              if ((prev_inf_coeff != 0) || (prev_sup_coeff != 0)) {
                const size_t filter_index =
                    out_z * filter_size_x * filter_size_y * input_size_z +
                    x_shift * filter_size_y * input_size_z +
                    y_shift * input_size_z + inp_z;

                const float_type aux_coeff = aux_coeffs[filter_index];
                elina_double_interval_mul_expr_coeff(
                    &tmp1, &tmp2, prev_inf_coeff, prev_sup_coeff, aux_coeff,
                    aux_coeff);

                maxRes = max(fabs(inf_coeff), fabs(sup_coeff));
                maxMul = max(fabs(tmp1), fabs(tmp2));

                inf_coeff = inf_coeff + tmp1 - (maxRes + maxMul) * ulp;
                sup_coeff = sup_coeff + tmp2 + (maxRes + maxMul) * ulp;
              }
            }

            res_inf_coeff[b] = inf_coeff;
            res_sup_coeff[b] = sup_coeff;
          }
        }

        __syncthreads();
      }
    }
  }
}

__global__ void coeffs_from_previous_layer_conv_sparse_filter_serial(
    const float_type *__restrict__ expr_inf_coeff,
    const float_type *__restrict__ expr_sup_coeff,
    float_type *__restrict__ res_inf_coeff,
    float_type *__restrict__ res_sup_coeff,
    const float_type *__restrict__ aux_coeffs, const size_t output_size_x,
    const size_t output_size_y, const size_t output_size_z,
    const size_t input_size_x, const size_t input_size_y,
    const size_t input_size_z, long int offset_x, long int offset_y,
    long int length_x, long int length_y, long int shift_x, long int shift_y,
    const size_t filter_size_x, const size_t filter_size_y,
    const size_t stride_x, const size_t stride_y, const size_t pad_x,
    const size_t pad_y) {
  const size_t last_x = blockIdx.x;
  const size_t last_y = blockIdx.y;
  const size_t last_z = blockIdx.z;

  const size_t n = last_x * gridDim.y * gridDim.z + last_y * gridDim.z + last_z;

  const long int min_out_x = offset_x + last_x * shift_x;
  const long int min_out_y = offset_y + last_y * shift_y;

  const size_t inp_z = threadIdx.x;

  float_type tmp1, tmp2;
  float_type maxRes, maxMul;

  if (inp_z < input_size_z) {
    for (long int out_x = 0; out_x < length_x; out_x++) {
      if ((out_x + min_out_x < 0) || (out_x + min_out_x >= output_size_x)) {
        continue;
      }

      for (long int out_y = 0; out_y < length_y; out_y++) {
        if ((out_y + min_out_y < 0) || (out_y + min_out_y >= output_size_y)) {
          continue;
        }

        for (size_t x_shift = 0; x_shift < filter_size_x; x_shift++) {
          for (size_t y_shift = 0; y_shift < filter_size_y; y_shift++) {
            const long int x_val =
                out_x * stride_x + min_out_x * stride_x + x_shift - pad_x;
            const long int y_val =
                out_y * stride_y + min_out_y * stride_y + y_shift - pad_y;

            if (!((x_val < 0) || (x_val >= (long int)input_size_x))) {
              if (!((y_val < 0) || (y_val >= (long int)input_size_y))) {
                const size_t mat_in =
                    (out_x * stride_x + x_shift) *
                        ((length_y - 1) * stride_y + filter_size_y) *
                        input_size_z +
                    (out_y * stride_y + y_shift) * input_size_z + inp_z;
                const size_t b =
                    n * ((length_x - 1) * stride_x + filter_size_x) *
                        ((length_y - 1) * stride_y + filter_size_y) *
                        input_size_z +
                    mat_in;

                float_type inf_coeff = res_inf_coeff[b];
                float_type sup_coeff = res_sup_coeff[b];

                for (size_t out_z = 0; out_z < output_size_z; out_z++) {
                  const size_t mat_out = out_x * length_y * output_size_z +
                                         out_y * output_size_z + out_z;
                  const size_t a =
                      n * length_x * length_y * output_size_z + mat_out;

                  const float_type prev_inf_coeff = expr_inf_coeff[a];
                  const float_type prev_sup_coeff = expr_sup_coeff[a];

                  if ((prev_inf_coeff != 0) || (prev_sup_coeff != 0)) {
                    const size_t filter_index =
                        out_z * filter_size_x * filter_size_y * input_size_z +
                        x_shift * filter_size_y * input_size_z +
                        y_shift * input_size_z + inp_z;

                    const float_type aux_coeff = aux_coeffs[filter_index];
                    elina_double_interval_mul_expr_coeff(
                        &tmp1, &tmp2, prev_inf_coeff, prev_sup_coeff, aux_coeff,
                        aux_coeff);

                    maxRes = max(fabs(inf_coeff), fabs(sup_coeff));
                    maxMul = max(fabs(tmp1), fabs(tmp2));

                    inf_coeff = inf_coeff + tmp1 - (maxRes + maxMul) * ulp;
                    sup_coeff = sup_coeff + tmp2 + (maxRes + maxMul) * ulp;
                  }
                }

                res_inf_coeff[b] = inf_coeff;
                res_sup_coeff[b] = sup_coeff;
              }
            }
          }
        }
      }
    }
  }
}

__global__ void
csts_from_previous_layer(const float_type *__restrict__ expr_inf_coeff,
                         const float_type *__restrict__ expr_sup_coeff,
                         const float_type *__restrict__ expr_inf_cst,
                         const float_type *__restrict__ expr_sup_cst,
                         float_type *__restrict__ res_inf_cst,
                         float_type *__restrict__ res_sup_cst,
                         const float_type *__restrict__ aux_csts,
                         const size_t num_out_neurons_current_layer) {
  const size_t n = blockIdx.x;

  float_type inf_cst = expr_inf_cst[n];
  float_type sup_cst = expr_sup_cst[n];

  float_type tmp1, tmp2;
  float_type maxRes, maxMul;

  for (size_t i = 0; i < num_out_neurons_current_layer; i++) {
    size_t a = n * num_out_neurons_current_layer + i;

    const float_type prev_inf_coeff = expr_inf_coeff[a];
    const float_type prev_sup_coeff = expr_sup_coeff[a];

    if ((prev_inf_coeff != 0) || (prev_sup_coeff != 0)) {
      elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, prev_inf_coeff,
                                          prev_sup_coeff, aux_csts[i],
                                          aux_csts[i]);

      maxRes = max(fabs(inf_cst), fabs(sup_cst));
      maxMul = max(fabs(tmp1), fabs(tmp2));

      inf_cst += tmp1 - (maxRes + maxMul) * ulp - min_denormal;
      sup_cst += tmp2 + (maxRes + maxMul) * ulp + min_denormal;
    }
  }

  res_inf_cst[n] = inf_cst;
  res_sup_cst[n] = sup_cst;
}

__global__ void
csts_from_previous_layer_conv(const float_type *__restrict__ expr_inf_coeff,
                              const float_type *__restrict__ expr_sup_coeff,
                              const float_type *__restrict__ expr_inf_cst,
                              const float_type *__restrict__ expr_sup_cst,
                              float_type *__restrict__ res_inf_cst,
                              float_type *__restrict__ res_sup_cst,
                              const float_type *__restrict__ aux_csts,
                              const size_t current_layer_out_size_x,
                              const size_t current_layer_out_size_y,
                              const size_t current_layer_out_size_z) {
  const size_t n = blockIdx.x;

  float_type inf_cst = expr_inf_cst[n];
  float_type sup_cst = expr_sup_cst[n];

  float_type tmp1, tmp2;
  float_type maxRes, maxMul;

  for (size_t i = 0; i < current_layer_out_size_x * current_layer_out_size_y;
       i++) {
    for (size_t j = 0; j < current_layer_out_size_z; j++) {
      size_t a = n * current_layer_out_size_x * current_layer_out_size_y *
                     current_layer_out_size_z +
                 i * current_layer_out_size_z + j;

      const float_type prev_inf_coeff = expr_inf_coeff[a];
      const float_type prev_sup_coeff = expr_sup_coeff[a];

      if ((prev_inf_coeff != 0) || (prev_sup_coeff != 0)) {
        elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, prev_inf_coeff,
                                            prev_sup_coeff, aux_csts[j],
                                            aux_csts[j]);

        maxRes = max(fabs(inf_cst), fabs(sup_cst));
        maxMul = max(fabs(tmp1), fabs(tmp2));

        inf_cst += tmp1 - (maxRes + maxMul) * ulp - min_denormal;
        sup_cst += tmp2 + (maxRes + maxMul) * ulp + min_denormal;
      }
    }
  }

  res_inf_cst[n] = inf_cst;
  res_sup_cst[n] = sup_cst;
}

__global__ void csts_from_previous_layer_conv_sparse(
    const float_type *__restrict__ expr_inf_coeff,
    const float_type *__restrict__ expr_sup_coeff,
    const float_type *__restrict__ expr_inf_cst,
    const float_type *__restrict__ expr_sup_cst,
    float_type *__restrict__ res_inf_cst, float_type *__restrict__ res_sup_cst,
    const float_type *__restrict__ aux_csts, const size_t output_size_x,
    const size_t output_size_y, const size_t output_size_z, long int offset_x,
    long int offset_y, long int length_x, long int length_y, long int shift_x,
    long int shift_y) {
  const size_t last_x = blockIdx.x;
  const size_t last_y = blockIdx.y;
  const size_t last_z = blockIdx.z;

  const size_t n = last_x * gridDim.y * gridDim.z + last_y * gridDim.z + last_z;

  const long int min_out_x = offset_x + last_x * shift_x;
  const long int min_out_y = offset_y + last_y * shift_y;

  float_type inf_cst = expr_inf_cst[n];
  float_type sup_cst = expr_sup_cst[n];

  float_type tmp1, tmp2;
  float_type maxRes, maxMul;

  for (long int out_x = 0; out_x < length_x; out_x++) {
    if ((out_x + min_out_x < 0) || (out_x + min_out_x >= output_size_x)) {
      continue;
    }

    for (long int out_y = 0; out_y < length_y; out_y++) {
      if ((out_y + min_out_y < 0) || (out_y + min_out_y >= output_size_y)) {
        continue;
      }

      for (size_t out_z = 0; out_z < output_size_z; out_z++) {
        size_t mat_out = n * length_x * length_y * output_size_z +
                         out_x * length_y * output_size_z +
                         out_y * output_size_z + out_z;

        const float_type prev_inf_coeff = expr_inf_coeff[mat_out];
        const float_type prev_sup_coeff = expr_sup_coeff[mat_out];

        if ((prev_inf_coeff != 0) || (prev_sup_coeff != 0)) {
          float_type aux_cst = aux_csts[out_z];

          elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, prev_inf_coeff,
                                              prev_sup_coeff, aux_cst, aux_cst);

          maxRes = max(fabs(inf_cst), fabs(sup_cst));
          maxMul = max(fabs(tmp1), fabs(tmp2));

          inf_cst += tmp1 - (maxRes + maxMul) * ulp - min_denormal;
          sup_cst += tmp2 + (maxRes + maxMul) * ulp + min_denormal;
        }
      }
    }
  }

  res_inf_cst[n] = inf_cst;
  res_sup_cst[n] = sup_cst;
}

__global__ void device_layer_create_sparse_exprs(
    float_type *dense_coeff, float_type *bias, const float_type *filter_weights,
    const float_type *filter_bias, const size_t chunk_counter,
    const size_t input_size_x, const size_t input_size_y,
    const size_t input_size_z, const size_t filter_size_x,
    const size_t filter_size_y, const size_t stride_x, const size_t stride_y,
    const long int pad_x, const long int pad_y, const size_t num_pixels) {
  size_t out_x = blockIdx.x;
  size_t out_y = blockIdx.y;
  size_t out_z = blockIdx.z;

  const size_t local_mat_x =
      out_x * gridDim.y * gridDim.z + out_y * gridDim.z + out_z;

  for (size_t x_shift = 0; x_shift < filter_size_x; x_shift++) {
    for (size_t y_shift = 0; y_shift < filter_size_y; y_shift++) {
      for (size_t inp_z = 0; inp_z < input_size_z; inp_z++) {
        const long int x_val = out_x * stride_x + x_shift - pad_x;
        const long int y_val = out_y * stride_y + y_shift - pad_y;

        if ((y_val < 0) || (y_val >= (long int)input_size_y)) {
          continue;
        }

        if ((x_val < 0) || (x_val >= (long int)input_size_x)) {
          continue;
        }

        const size_t mat_y = x_shift * filter_size_y * input_size_z +
                             y_shift * input_size_z + inp_z;
        const size_t filter_index = (chunk_counter * gridDim.z + out_z) *
                                        filter_size_x * filter_size_y *
                                        input_size_z +
                                    x_shift * filter_size_y * input_size_z +
                                    y_shift * input_size_z + inp_z;

        dense_coeff[local_mat_x * filter_size_x * filter_size_y * input_size_z +
                    mat_y] = filter_weights[filter_index];
      }
    }
  }

  bias[local_mat_x] = filter_bias[chunk_counter * gridDim.z + out_z];
}

void update_state_using_previous_layers(elina_manager_t *man, fppoly_t *fp,
                                        const size_t layerno) {
  auto start = std::chrono::system_clock::now();

  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  const size_t num_out_neurons_last_layer =
      fp->layers[layerno]->num_out_neurons;
  const size_t num_in_neurons_last_layer = fp->layers[layerno]->num_in_neurons;

  const size_t num_in_neurons_first_layer = fp->layers[0]->num_in_neurons;

  std::cout << "num_out_neurons_last " << num_out_neurons_last_layer
            << std::endl;

  float_type *coeffs = fp->layers[layerno]->coeffs;
  float_type *csts = fp->layers[layerno]->csts;

  float_type *lb_array = fp->layers[layerno]->lb_array;
  float_type *ub_array = fp->layers[layerno]->ub_array;

  float_type *linf_coeff;
  float_type *lsup_coeff;
  float_type *linf_cst;
  float_type *lsup_cst;

  cudaMalloc((void **)&linf_coeff, num_out_neurons_last_layer *
                                       num_in_neurons_last_layer *
                                       sizeof(float_type));
  cudaMalloc((void **)&lsup_coeff, num_out_neurons_last_layer *
                                       num_in_neurons_last_layer *
                                       sizeof(float_type));

  cudaMalloc((void **)&linf_cst,
             num_out_neurons_last_layer * sizeof(float_type));
  cudaMalloc((void **)&lsup_cst,
             num_out_neurons_last_layer * sizeof(float_type));

  float_type *uinf_coeff;
  float_type *usup_coeff;
  float_type *uinf_cst;
  float_type *usup_cst;

  cudaMalloc((void **)&uinf_coeff, num_out_neurons_last_layer *
                                       num_in_neurons_last_layer *
                                       sizeof(float_type));
  cudaMalloc((void **)&usup_coeff, num_out_neurons_last_layer *
                                       num_in_neurons_last_layer *
                                       sizeof(float_type));

  cudaMalloc((void **)&uinf_cst,
             num_out_neurons_last_layer * sizeof(float_type));
  cudaMalloc((void **)&usup_cst,
             num_out_neurons_last_layer * sizeof(float_type));

  copy_expr_array<<<num_out_neurons_last_layer, 1>>>(
      linf_coeff, lsup_coeff, linf_cst, lsup_cst, coeffs, csts,
      num_out_neurons_last_layer, num_in_neurons_last_layer);
  copy_expr_array<<<num_out_neurons_last_layer, 1>>>(
      uinf_coeff, usup_coeff, uinf_cst, usup_cst, coeffs, csts,
      num_out_neurons_last_layer, num_in_neurons_last_layer);

  float_type *linf_coeff_tmp;
  float_type *lsup_coeff_tmp;
  float_type *linf_cst_tmp;
  float_type *lsup_cst_tmp;

  float_type *uinf_coeff_tmp;
  float_type *usup_coeff_tmp;
  float_type *uinf_cst_tmp;
  float_type *usup_cst_tmp;

  cudaMalloc((void **)&linf_cst_tmp,
             num_out_neurons_last_layer * sizeof(float_type));
  cudaMalloc((void **)&lsup_cst_tmp,
             num_out_neurons_last_layer * sizeof(float_type));
  cudaMalloc((void **)&uinf_cst_tmp,
             num_out_neurons_last_layer * sizeof(float_type));
  cudaMalloc((void **)&usup_cst_tmp,
             num_out_neurons_last_layer * sizeof(float_type));

  for (int k = layerno - 1; k >= 0; k--) {
    const size_t num_out_neurons_current_layer = fp->layers[k]->num_out_neurons;
    const size_t num_in_neurons_current_layer = fp->layers[k]->num_in_neurons;
    std::cout << "num_out_neurons_current " << num_out_neurons_current_layer
              << " num_in_neurons_current " << num_in_neurons_current_layer
              << std::endl;

    const dim3 num_blocks_relu(num_out_neurons_last_layer,
                               num_out_neurons_current_layer / num_threads + 1,
                               1);
    const dim3 num_blocks_linear(num_out_neurons_last_layer,
                                 num_in_neurons_current_layer / num_threads + 1,
                                 1);

    std::cout << "num_threads" << num_threads << " num_blocks_relu "
              << num_blocks_relu.y << " num_blocks_linear "
              << num_blocks_linear.y << std::endl;

    float_type *aux_coeffs;
    float_type *aux_csts;

    if (fp->layers[k]->type == CONV) {
      aux_coeffs = fp->layers[k]->filter_weights;
      aux_csts = fp->layers[k]->filter_bias;
    } else {
      aux_coeffs = fp->layers[k]->coeffs;
      aux_csts = fp->layers[k]->csts;
    }

    float_type *aux_lb_array = fp->layers[k]->lb_array;
    float_type *aux_ub_array = fp->layers[k]->ub_array;

    if (fp->layers[k]->activation == RELU) {
      lexpr_replace_relu_bounds<<<num_blocks_relu, num_threads>>>(
          linf_coeff, lsup_coeff, linf_cst, lsup_cst, aux_lb_array,
          aux_ub_array, num_out_neurons_current_layer);
      uexpr_replace_relu_bounds<<<num_blocks_relu, num_threads>>>(
          uinf_coeff, usup_coeff, uinf_cst, usup_cst, aux_lb_array,
          aux_ub_array, num_out_neurons_current_layer);
    }

    cudaMalloc((void **)&linf_coeff_tmp, num_out_neurons_last_layer *
                                             num_in_neurons_current_layer *
                                             sizeof(float_type));
    cudaMalloc((void **)&lsup_coeff_tmp, num_out_neurons_last_layer *
                                             num_in_neurons_current_layer *
                                             sizeof(float_type));
    cudaMalloc((void **)&uinf_coeff_tmp, num_out_neurons_last_layer *
                                             num_in_neurons_current_layer *
                                             sizeof(float_type));
    cudaMalloc((void **)&usup_coeff_tmp, num_out_neurons_last_layer *
                                             num_in_neurons_current_layer *
                                             sizeof(float_type));

    cudaMemset(linf_coeff_tmp, 0,
               num_out_neurons_last_layer * num_in_neurons_current_layer *
                   sizeof(float_type));
    cudaMemset(lsup_coeff_tmp, 0,
               num_out_neurons_last_layer * num_in_neurons_current_layer *
                   sizeof(float_type));
    cudaMemset(uinf_coeff_tmp, 0,
               num_out_neurons_last_layer * num_in_neurons_current_layer *
                   sizeof(float_type));
    cudaMemset(usup_coeff_tmp, 0,
               num_out_neurons_last_layer * num_in_neurons_current_layer *
                   sizeof(float_type));

    if (fp->layers[k]->type == CONV) {
      if (fp->layers[k]->input_size[2] >= 128) {
        coeffs_from_previous_layer_conv_filter_serial<<<
            num_out_neurons_last_layer, fp->layers[k]->input_size[2]>>>(
            linf_coeff, lsup_coeff, linf_coeff_tmp, lsup_coeff_tmp, aux_coeffs,
            fp->layers[k]->output_size[0], fp->layers[k]->output_size[1],
            fp->layers[k]->output_size[2], fp->layers[k]->input_size[0],
            fp->layers[k]->input_size[1], fp->layers[k]->input_size[2],
            fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1],
            fp->layers[k]->strides[0], fp->layers[k]->strides[1],
            fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
        coeffs_from_previous_layer_conv_filter_serial<<<
            num_out_neurons_last_layer, fp->layers[k]->input_size[2]>>>(
            uinf_coeff, usup_coeff, uinf_coeff_tmp, usup_coeff_tmp, aux_coeffs,
            fp->layers[k]->output_size[0], fp->layers[k]->output_size[1],
            fp->layers[k]->output_size[2], fp->layers[k]->input_size[0],
            fp->layers[k]->input_size[1], fp->layers[k]->input_size[2],
            fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1],
            fp->layers[k]->strides[0], fp->layers[k]->strides[1],
            fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
      } else {
        coeffs_from_previous_layer_conv<<<
            num_out_neurons_last_layer,
            dim3(fp->layers[k]->input_size[2], fp->layers[k]->filter_size[1],
                 fp->layers[k]->filter_size[0])>>>(
            linf_coeff, lsup_coeff, linf_coeff_tmp, lsup_coeff_tmp, aux_coeffs,
            fp->layers[k]->output_size[0], fp->layers[k]->output_size[1],
            fp->layers[k]->output_size[2], fp->layers[k]->input_size[0],
            fp->layers[k]->input_size[1], fp->layers[k]->input_size[2],
            fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1],
            fp->layers[k]->strides[0], fp->layers[k]->strides[1],
            fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
        coeffs_from_previous_layer_conv<<<
            num_out_neurons_last_layer,
            dim3(fp->layers[k]->input_size[2], fp->layers[k]->filter_size[1],
                 fp->layers[k]->filter_size[0])>>>(
            uinf_coeff, usup_coeff, uinf_coeff_tmp, usup_coeff_tmp, aux_coeffs,
            fp->layers[k]->output_size[0], fp->layers[k]->output_size[1],
            fp->layers[k]->output_size[2], fp->layers[k]->input_size[0],
            fp->layers[k]->input_size[1], fp->layers[k]->input_size[2],
            fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1],
            fp->layers[k]->strides[0], fp->layers[k]->strides[1],
            fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
      }

      csts_from_previous_layer_conv<<<num_out_neurons_last_layer, 1>>>(
          linf_coeff, lsup_coeff, linf_cst, lsup_cst, linf_cst_tmp,
          lsup_cst_tmp, aux_csts, fp->layers[k]->output_size[0],
          fp->layers[k]->output_size[1], fp->layers[k]->output_size[2]);
      csts_from_previous_layer_conv<<<num_out_neurons_last_layer, 1>>>(
          uinf_coeff, usup_coeff, uinf_cst, usup_cst, uinf_cst_tmp,
          usup_cst_tmp, aux_csts, fp->layers[k]->output_size[0],
          fp->layers[k]->output_size[1], fp->layers[k]->output_size[2]);

    } else {
      coeffs_from_previous_layer<<<num_blocks_linear, num_threads>>>(
          linf_coeff, lsup_coeff, linf_coeff_tmp, lsup_coeff_tmp, aux_coeffs,
          num_out_neurons_current_layer, num_in_neurons_current_layer);
      coeffs_from_previous_layer<<<num_blocks_linear, num_threads>>>(
          uinf_coeff, usup_coeff, uinf_coeff_tmp, usup_coeff_tmp, aux_coeffs,
          num_out_neurons_current_layer, num_in_neurons_current_layer);

      csts_from_previous_layer<<<num_out_neurons_last_layer, 1>>>(
          linf_coeff, lsup_coeff, linf_cst, lsup_cst, linf_cst_tmp,
          lsup_cst_tmp, aux_csts, num_out_neurons_current_layer);
      csts_from_previous_layer<<<num_out_neurons_last_layer, 1>>>(
          uinf_coeff, usup_coeff, uinf_cst, usup_cst, uinf_cst_tmp,
          usup_cst_tmp, aux_csts, num_out_neurons_current_layer);
    }

    std::swap(linf_coeff, linf_coeff_tmp);
    std::swap(lsup_coeff, lsup_coeff_tmp);
    std::swap(linf_cst, linf_cst_tmp);
    std::swap(lsup_cst, lsup_cst_tmp);

    std::swap(uinf_coeff, uinf_coeff_tmp);
    std::swap(usup_coeff, usup_coeff_tmp);
    std::swap(uinf_cst, uinf_cst_tmp);
    std::swap(usup_cst, usup_cst_tmp);

    cudaFree(linf_coeff_tmp);
    cudaFree(lsup_coeff_tmp);
    cudaFree(uinf_coeff_tmp);
    cudaFree(usup_coeff_tmp);
  }

  compute_lb_from_expr<<<num_out_neurons_last_layer, 1>>>(
      lb_array, linf_coeff, lsup_coeff, linf_cst, fp->input_inf, fp->input_sup,
      num_in_neurons_first_layer);
  compute_ub_from_expr<<<num_out_neurons_last_layer, 1>>>(
      ub_array, uinf_coeff, usup_coeff, usup_cst, fp->input_inf, fp->input_sup,
      num_in_neurons_first_layer);

  cudaFree(linf_coeff);
  cudaFree(lsup_coeff);
  cudaFree(linf_cst);
  cudaFree(lsup_cst);

  cudaFree(uinf_coeff);
  cudaFree(usup_coeff);
  cudaFree(uinf_cst);
  cudaFree(usup_cst);

  cudaFree(linf_cst_tmp);
  cudaFree(lsup_cst_tmp);

  cudaFree(uinf_cst_tmp);
  cudaFree(usup_cst_tmp);

  cudaDeviceSynchronize();

  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl
            << std::endl;
}

size_t predict_size(fppoly_t *fp, const size_t layerno) {
  size_t free;
  size_t total;

  cudaMemGetInfo(&free, &total);

  std::cout << "USABLE " << free - (2 << 27) << std::endl;

  const size_t num_out_neurons_last_layer =
      fp->layers[layerno]->num_out_neurons;

  offset_x = -fp->layers[layerno]->pad[0];
  offset_y = -fp->layers[layerno]->pad[1];

  std::cout << "offset_x " << offset_x << " offset_y " << offset_y << std::endl;

  length_x = fp->layers[layerno]->filter_size[0];
  length_y = fp->layers[layerno]->filter_size[1];

  std::cout << "length_x " << length_x << " length_y " << length_y << std::endl;

  shift_x = fp->layers[layerno]->strides[0];
  shift_y = fp->layers[layerno]->strides[1];

  std::cout << "shift_x " << shift_x << " shift_y " << shift_y << std::endl;

  size_t current_size = num_out_neurons_last_layer * length_x * length_y *
                        fp->layers[layerno]->input_size[2] * sizeof(float_type);
  size_t last_size;

  size_t maximum_size = 0;

  std::cout << "Starting size: " << current_size << " bytes" << std::endl;

  for (int k = layerno - 1; k >= 0; k--) {
    const size_t num_out_neurons_current_layer = fp->layers[k]->num_out_neurons;
    const size_t num_in_neurons_current_layer = fp->layers[k]->num_in_neurons;
    std::cout << "num_out_neurons_current " << num_out_neurons_current_layer
              << " num_in_neurons_current " << num_in_neurons_current_layer
              << std::endl;

    offset_x = fp->layers[k]->strides[0] * offset_x - fp->layers[k]->pad[0];
    offset_y = fp->layers[k]->strides[1] * offset_y - fp->layers[k]->pad[1];

    std::cout << "offset_x " << offset_x << " offset_y " << offset_y
              << std::endl;

    length_x = (length_x - 1) * fp->layers[k]->strides[0] +
               fp->layers[k]->filter_size[0];
    length_y = (length_y - 1) * fp->layers[k]->strides[1] +
               fp->layers[k]->filter_size[1];

    std::cout << "length_x " << length_x << " length_y " << length_y
              << std::endl;

    shift_x = fp->layers[k]->strides[0] * shift_x;
    shift_y = fp->layers[k]->strides[1] * shift_y;

    std::cout << "shift_x " << shift_x << " shift_y " << shift_y << std::endl;

    size_t missing_length = length_x * length_y * fp->layers[k]->input_size[2];

    std::cout << "POST SIZE " << missing_length << std::endl;

    std::cout << std::endl;

    last_size = current_size;
    current_size = num_out_neurons_last_layer * length_x * length_y *
                   fp->layers[k]->input_size[2] * sizeof(float_type);

    std::cout << "Size: " << current_size << " bytes" << std::endl;

    if (4 * (last_size + current_size) > maximum_size) {
      maximum_size = 4 * (last_size + current_size);
    }
  }

  size_t num_chunks = 1;

  while (maximum_size > free) {
    maximum_size /= 2;
    num_chunks *= 2;
  }

  if (num_chunks > 1) {
    std::cout << "Number of chunks " << num_chunks << std::endl;
  }

  return num_chunks;
}

void update_state_using_previous_layers_conv_chunk(elina_manager_t *man,
                                                   fppoly_t *fp,
                                                   const size_t layerno,
                                                   const size_t num_chunks,
                                                   const size_t chunk_counter) {
  auto start = std::chrono::system_clock::now();

  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  const size_t x_y_size_last_layer =
      fp->layers[layerno]->output_size[0] * fp->layers[layerno]->output_size[1];
  const size_t num_filters_last_layer =
      fp->layers[layerno]->output_size[2] / num_chunks;

  std::cout << "num_out_neurons_last "
            << x_y_size_last_layer * num_filters_last_layer << std::endl;

  offset_x = -fp->layers[layerno]->pad[0];
  offset_y = -fp->layers[layerno]->pad[1];

  length_x = fp->layers[layerno]->filter_size[0];
  length_y = fp->layers[layerno]->filter_size[1];

  shift_x = fp->layers[layerno]->strides[0];
  shift_y = fp->layers[layerno]->strides[1];

  float_type *coeffs;
  float_type *csts;

  cudaMalloc((void **)&coeffs, x_y_size_last_layer * num_filters_last_layer *
                                   length_x * length_y *
                                   fp->layers[layerno]->input_size[2] *
                                   sizeof(float_type));
  cudaMalloc((void **)&csts,
             x_y_size_last_layer * num_filters_last_layer * sizeof(float_type));

  cudaMemset(coeffs, 0,
             x_y_size_last_layer * num_filters_last_layer * length_x *
                 length_y * fp->layers[layerno]->input_size[2] *
                 sizeof(float_type));
  cudaMemset(csts, 0,
             x_y_size_last_layer * num_filters_last_layer * sizeof(float_type));

  device_layer_create_sparse_exprs<<<dim3(fp->layers[layerno]->output_size[0],
                                          fp->layers[layerno]->output_size[1],
                                          fp->layers[layerno]->output_size[2] /
                                              num_chunks),
                                     1>>>(
      coeffs, csts, fp->layers[layerno]->filter_weights,
      fp->layers[layerno]->filter_bias, chunk_counter,
      fp->layers[layerno]->input_size[0], fp->layers[layerno]->input_size[1],
      fp->layers[layerno]->input_size[2], fp->layers[layerno]->filter_size[0],
      fp->layers[layerno]->filter_size[1], fp->layers[layerno]->strides[0],
      fp->layers[layerno]->strides[1], fp->layers[layerno]->pad[0],
      fp->layers[layerno]->pad[1], fp->layers[layerno]->num_in_neurons);

  float_type *lb_array = fp->layers[layerno]->lb_array;
  float_type *ub_array = fp->layers[layerno]->ub_array;

  float_type *linf_coeff;
  float_type *lsup_coeff;
  float_type *linf_cst;
  float_type *lsup_cst;

  cudaMalloc((void **)&linf_coeff,
             x_y_size_last_layer * num_filters_last_layer * length_x *
                 length_y * fp->layers[layerno]->input_size[2] *
                 sizeof(float_type));
  cudaMalloc((void **)&lsup_coeff,
             x_y_size_last_layer * num_filters_last_layer * length_x *
                 length_y * fp->layers[layerno]->input_size[2] *
                 sizeof(float_type));

  cudaMalloc((void **)&linf_cst,
             x_y_size_last_layer * num_filters_last_layer * sizeof(float_type));
  cudaMalloc((void **)&lsup_cst,
             x_y_size_last_layer * num_filters_last_layer * sizeof(float_type));

  float_type *uinf_coeff;
  float_type *usup_coeff;
  float_type *uinf_cst;
  float_type *usup_cst;

  cudaMalloc((void **)&uinf_coeff,
             x_y_size_last_layer * num_filters_last_layer * length_x *
                 length_y * fp->layers[layerno]->input_size[2] *
                 sizeof(float_type));
  cudaMalloc((void **)&usup_coeff,
             x_y_size_last_layer * num_filters_last_layer * length_x *
                 length_y * fp->layers[layerno]->input_size[2] *
                 sizeof(float_type));

  cudaMalloc((void **)&uinf_cst,
             x_y_size_last_layer * num_filters_last_layer * sizeof(float_type));
  cudaMalloc((void **)&usup_cst,
             x_y_size_last_layer * num_filters_last_layer * sizeof(float_type));

  copy_expr_array<<<x_y_size_last_layer * num_filters_last_layer, 1>>>(
      linf_coeff, lsup_coeff, linf_cst, lsup_cst, coeffs, csts,
      x_y_size_last_layer * num_filters_last_layer,
      length_x * length_y * fp->layers[layerno]->input_size[2]);
  copy_expr_array<<<x_y_size_last_layer * num_filters_last_layer, 1>>>(
      uinf_coeff, usup_coeff, uinf_cst, usup_cst, coeffs, csts,
      x_y_size_last_layer * num_filters_last_layer,
      length_x * length_y * fp->layers[layerno]->input_size[2]);

  cudaFree(coeffs);
  cudaFree(csts);

  float_type *linf_coeff_tmp;
  float_type *lsup_coeff_tmp;
  float_type *linf_cst_tmp;
  float_type *lsup_cst_tmp;

  float_type *uinf_coeff_tmp;
  float_type *usup_coeff_tmp;
  float_type *uinf_cst_tmp;
  float_type *usup_cst_tmp;

  cudaMalloc((void **)&linf_cst_tmp,
             x_y_size_last_layer * num_filters_last_layer * sizeof(float_type));
  cudaMalloc((void **)&lsup_cst_tmp,
             x_y_size_last_layer * num_filters_last_layer * sizeof(float_type));
  cudaMalloc((void **)&uinf_cst_tmp,
             x_y_size_last_layer * num_filters_last_layer * sizeof(float_type));
  cudaMalloc((void **)&usup_cst_tmp,
             x_y_size_last_layer * num_filters_last_layer * sizeof(float_type));

  for (int k = layerno - 1; k >= 0; k--) {
    const size_t num_out_neurons_current_layer = fp->layers[k]->num_out_neurons;
    const size_t num_in_neurons_current_layer = fp->layers[k]->num_in_neurons;
    std::cout << "num_out_neurons_current " << num_out_neurons_current_layer
              << " num_in_neurons_current " << num_in_neurons_current_layer
              << std::endl;

    float_type *aux_coeffs;
    float_type *aux_csts;

    if (fp->layers[k]->type == CONV) {
      aux_coeffs = fp->layers[k]->filter_weights;
      aux_csts = fp->layers[k]->filter_bias;
    } else {
      aux_coeffs = fp->layers[k]->coeffs;
      aux_csts = fp->layers[k]->csts;
    }

    float_type *aux_lb_array = fp->layers[k]->lb_array;
    float_type *aux_ub_array = fp->layers[k]->ub_array;

    if (fp->layers[k]->activation == RELU) {
      lexpr_replace_relu_bounds_conv_sparse<<<
          dim3(fp->layers[layerno]->output_size[0],
               fp->layers[layerno]->output_size[1],
               fp->layers[layerno]->output_size[2] / num_chunks),
          1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, aux_lb_array,
               aux_ub_array, fp->layers[k]->output_size[0],
               fp->layers[k]->output_size[1], fp->layers[k]->output_size[2],
               offset_x, offset_y, length_x, length_y, shift_x, shift_y);
      uexpr_replace_relu_bounds_conv_sparse<<<
          dim3(fp->layers[layerno]->output_size[0],
               fp->layers[layerno]->output_size[1],
               fp->layers[layerno]->output_size[2] / num_chunks),
          1>>>(uinf_coeff, usup_coeff, uinf_cst, usup_cst, aux_lb_array,
               aux_ub_array, fp->layers[k]->output_size[0],
               fp->layers[k]->output_size[1], fp->layers[k]->output_size[2],
               offset_x, offset_y, length_x, length_y, shift_x, shift_y);
    }

    size_t missing_length = ((length_x - 1) * fp->layers[k]->strides[0] +
                             fp->layers[k]->filter_size[0]) *
                            ((length_y - 1) * fp->layers[k]->strides[1] +
                             fp->layers[k]->filter_size[1]) *
                            fp->layers[k]->input_size[2];

    cudaMalloc((void **)&linf_coeff_tmp,
               x_y_size_last_layer * num_filters_last_layer * missing_length *
                   sizeof(float_type));
    cudaMalloc((void **)&lsup_coeff_tmp,
               x_y_size_last_layer * num_filters_last_layer * missing_length *
                   sizeof(float_type));
    cudaMalloc((void **)&uinf_coeff_tmp,
               x_y_size_last_layer * num_filters_last_layer * missing_length *
                   sizeof(float_type));
    cudaMalloc((void **)&usup_coeff_tmp,
               x_y_size_last_layer * num_filters_last_layer * missing_length *
                   sizeof(float_type));

    cudaMemset(linf_coeff_tmp, 0,
               x_y_size_last_layer * num_filters_last_layer * missing_length *
                   sizeof(float_type));
    cudaMemset(lsup_coeff_tmp, 0,
               x_y_size_last_layer * num_filters_last_layer * missing_length *
                   sizeof(float_type));
    cudaMemset(uinf_coeff_tmp, 0,
               x_y_size_last_layer * num_filters_last_layer * missing_length *
                   sizeof(float_type));
    cudaMemset(usup_coeff_tmp, 0,
               x_y_size_last_layer * num_filters_last_layer * missing_length *
                   sizeof(float_type));

    if (fp->layers[k]->input_size[2] >= 128) {
      coeffs_from_previous_layer_conv_sparse_filter_serial<<<
          dim3(fp->layers[layerno]->output_size[0],
               fp->layers[layerno]->output_size[1],
               fp->layers[layerno]->output_size[2] / num_chunks),
          fp->layers[k]->input_size[2]>>>(
          linf_coeff, lsup_coeff, linf_coeff_tmp, lsup_coeff_tmp, aux_coeffs,
          fp->layers[k]->output_size[0], fp->layers[k]->output_size[1],
          fp->layers[k]->output_size[2], fp->layers[k]->input_size[0],
          fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], offset_x,
          offset_y, length_x, length_y, shift_x, shift_y,
          fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1],
          fp->layers[k]->strides[0], fp->layers[k]->strides[1],
          fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
      coeffs_from_previous_layer_conv_sparse_filter_serial<<<
          dim3(fp->layers[layerno]->output_size[0],
               fp->layers[layerno]->output_size[1],
               fp->layers[layerno]->output_size[2] / num_chunks),
          fp->layers[k]->input_size[2]>>>(
          uinf_coeff, usup_coeff, uinf_coeff_tmp, usup_coeff_tmp, aux_coeffs,
          fp->layers[k]->output_size[0], fp->layers[k]->output_size[1],
          fp->layers[k]->output_size[2], fp->layers[k]->input_size[0],
          fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], offset_x,
          offset_y, length_x, length_y, shift_x, shift_y,
          fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1],
          fp->layers[k]->strides[0], fp->layers[k]->strides[1],
          fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
    } else {
      coeffs_from_previous_layer_conv_sparse<<<
          dim3(fp->layers[layerno]->output_size[0],
               fp->layers[layerno]->output_size[1],
               fp->layers[layerno]->output_size[2] / num_chunks),
          dim3(fp->layers[k]->input_size[2], fp->layers[k]->filter_size[1],
               fp->layers[k]->filter_size[0])>>>(
          linf_coeff, lsup_coeff, linf_coeff_tmp, lsup_coeff_tmp, aux_coeffs,
          fp->layers[k]->output_size[0], fp->layers[k]->output_size[1],
          fp->layers[k]->output_size[2], fp->layers[k]->input_size[0],
          fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], offset_x,
          offset_y, length_x, length_y, shift_x, shift_y,
          fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1],
          fp->layers[k]->strides[0], fp->layers[k]->strides[1],
          fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
      coeffs_from_previous_layer_conv_sparse<<<
          dim3(fp->layers[layerno]->output_size[0],
               fp->layers[layerno]->output_size[1],
               fp->layers[layerno]->output_size[2] / num_chunks),
          dim3(fp->layers[k]->input_size[2], fp->layers[k]->filter_size[1],
               fp->layers[k]->filter_size[0])>>>(
          uinf_coeff, usup_coeff, uinf_coeff_tmp, usup_coeff_tmp, aux_coeffs,
          fp->layers[k]->output_size[0], fp->layers[k]->output_size[1],
          fp->layers[k]->output_size[2], fp->layers[k]->input_size[0],
          fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], offset_x,
          offset_y, length_x, length_y, shift_x, shift_y,
          fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1],
          fp->layers[k]->strides[0], fp->layers[k]->strides[1],
          fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
    }

    csts_from_previous_layer_conv_sparse<<<
        dim3(fp->layers[layerno]->output_size[0],
             fp->layers[layerno]->output_size[1],
             fp->layers[layerno]->output_size[2] / num_chunks),
        1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, linf_cst_tmp,
             lsup_cst_tmp, aux_csts, fp->layers[k]->output_size[0],
             fp->layers[k]->output_size[1], fp->layers[k]->output_size[2],
             offset_x, offset_y, length_x, length_y, shift_x, shift_y);
    csts_from_previous_layer_conv_sparse<<<
        dim3(fp->layers[layerno]->output_size[0],
             fp->layers[layerno]->output_size[1],
             fp->layers[layerno]->output_size[2] / num_chunks),
        1>>>(uinf_coeff, usup_coeff, uinf_cst, usup_cst, uinf_cst_tmp,
             usup_cst_tmp, aux_csts, fp->layers[k]->output_size[0],
             fp->layers[k]->output_size[1], fp->layers[k]->output_size[2],
             offset_x, offset_y, length_x, length_y, shift_x, shift_y);

    offset_x = fp->layers[k]->strides[0] * offset_x - fp->layers[k]->pad[0];
    offset_y = fp->layers[k]->strides[1] * offset_y - fp->layers[k]->pad[1];

    length_x = (length_x - 1) * fp->layers[k]->strides[0] +
               fp->layers[k]->filter_size[0];
    length_y = (length_y - 1) * fp->layers[k]->strides[1] +
               fp->layers[k]->filter_size[1];

    shift_x = fp->layers[k]->strides[0] * shift_x;
    shift_y = fp->layers[k]->strides[1] * shift_y;

    std::swap(linf_coeff, linf_coeff_tmp);
    std::swap(lsup_coeff, lsup_coeff_tmp);
    std::swap(linf_cst, linf_cst_tmp);
    std::swap(lsup_cst, lsup_cst_tmp);

    std::swap(uinf_coeff, uinf_coeff_tmp);
    std::swap(usup_coeff, usup_coeff_tmp);
    std::swap(uinf_cst, uinf_cst_tmp);
    std::swap(usup_cst, usup_cst_tmp);

    cudaFree(linf_coeff_tmp);
    cudaFree(lsup_coeff_tmp);
    cudaFree(uinf_coeff_tmp);
    cudaFree(usup_coeff_tmp);
  }

  compute_lb_from_expr_conv_sparse<<<dim3(fp->layers[layerno]->output_size[0],
                                          fp->layers[layerno]->output_size[1],
                                          fp->layers[layerno]->output_size[2] /
                                              num_chunks),
                                     1>>>(
      lb_array, linf_coeff, lsup_coeff, linf_cst, fp->input_inf, fp->input_sup,
      num_chunks, chunk_counter, fp->layers[0]->input_size[0],
      fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], offset_x,
      offset_y, length_x, length_y, shift_x, shift_y);
  compute_ub_from_expr_conv_sparse<<<dim3(fp->layers[layerno]->output_size[0],
                                          fp->layers[layerno]->output_size[1],
                                          fp->layers[layerno]->output_size[2] /
                                              num_chunks),
                                     1>>>(
      ub_array, uinf_coeff, usup_coeff, usup_cst, fp->input_inf, fp->input_sup,
      num_chunks, chunk_counter, fp->layers[0]->input_size[0],
      fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], offset_x,
      offset_y, length_x, length_y, shift_x, shift_y);

  cudaFree(linf_coeff);
  cudaFree(lsup_coeff);
  cudaFree(linf_cst);
  cudaFree(lsup_cst);

  cudaFree(uinf_coeff);
  cudaFree(usup_coeff);
  cudaFree(uinf_cst);
  cudaFree(usup_cst);

  cudaFree(linf_cst_tmp);
  cudaFree(lsup_cst_tmp);

  cudaFree(uinf_cst_tmp);
  cudaFree(usup_cst_tmp);

  cudaDeviceSynchronize();

  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl
            << std::endl;
}

void update_state_using_previous_layers_conv(elina_manager_t *man, fppoly_t *fp,
                                             const size_t layerno) {
  const size_t num_chunks = predict_size(fp, layerno);

  for (size_t chunk_counter = 0; chunk_counter < num_chunks; chunk_counter++) {
    update_state_using_previous_layers_conv_chunk(man, fp, fp->numlayers - 1,
                                                  num_chunks, chunk_counter);
  }
}

void ffn_handle_intermediate_layer(elina_manager_t *man,
                                   elina_abstract0_t *element,
                                   const double **weights, const double *bias,
                                   const size_t num_out_neurons,
                                   const size_t num_in_neurons,
                                   const activation_type_t activation) {
  fppoly_t *fp = fppoly_of_abstract0(element);
  ffn_add_layer(fp, num_out_neurons, num_in_neurons, FFN, activation);

  float_type *coeffs = fp->layers[fp->numlayers - 1]->coeffs;
  float_type *csts = fp->layers[fp->numlayers - 1]->csts;

  layer_create_dense_exprs(coeffs, csts, weights, bias, num_out_neurons,
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

__global__ void print_bounds(const float_type *__restrict__ bounds_array,
                             const size_t num_out_neurons) {}

void ffn_handle_last_layer(elina_manager_t *man, elina_abstract0_t *element,
                           const double **weights, const double *bias,
                           const size_t num_out_neurons,
                           const size_t num_in_neurons,
                           const bool has_activation,
                           const activation_type_t activation) {
  fppoly_t *fp = fppoly_of_abstract0(element);
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  if (has_activation) {
    ffn_add_layer(fp, num_out_neurons, num_in_neurons, FFN, activation);
  } else {
    ffn_add_layer(fp, num_out_neurons, num_in_neurons, FFN, NONE);
  }

  float_type *coeffs = fp->layers[fp->numlayers - 1]->coeffs;
  float_type *csts = fp->layers[fp->numlayers - 1]->csts;

  layer_create_dense_exprs(coeffs, csts, weights, bias, num_out_neurons,
                           num_in_neurons);

  update_state_using_previous_layers(man, fp, fp->numlayers - 1);

  float_type *lb_array = fp->layers[fp->numlayers - 1]->lb_array;
  float_type *ub_array = fp->layers[fp->numlayers - 1]->ub_array;

  float_type *lb_array_host = (float_type *)malloc(
      fp->layers[fp->numlayers - 1]->num_out_neurons * sizeof(float_type));
  float_type *ub_array_host = (float_type *)malloc(
      fp->layers[fp->numlayers - 1]->num_out_neurons * sizeof(float_type));

  cudaMemcpy(lb_array_host, lb_array,
             fp->layers[fp->numlayers - 1]->num_out_neurons *
                 sizeof(float_type),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(ub_array_host, ub_array,
             fp->layers[fp->numlayers - 1]->num_out_neurons *
                 sizeof(float_type),
             cudaMemcpyDeviceToHost);

  for (size_t i = 0; i < num_out_neurons; i++) {
    printf("out inf number %lu is: %.*e\n", i, DIGS, lb_array_host[i]);
  }

  for (size_t i = 0; i < num_out_neurons; i++) {
    printf("out sup number %lu is: %.*e\n", i, DIGS, ub_array_host[i]);
  }

  std::cout << std::endl;

  free(lb_array_host);
  free(ub_array_host);
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

__global__ void create_sub_expr(float_type *__restrict__ inf_coeff,
                                float_type *__restrict__ sup_coeff,
                                float_type *__restrict__ inf_cst,
                                float_type *__restrict__ sup_cst,
                                const size_t index, const elina_dim_t y,
                                const elina_dim_t x) {
  inf_cst[index] = 0;
  sup_cst[index] = 0;

  for (size_t i = 0; i < 10; i++) {
    inf_coeff[index * 10 + i] = 0.;
    sup_coeff[index * 10 + i] = 0.;
  }

  inf_coeff[index * 10 + y] = 1.;
  sup_coeff[index * 10 + y] = 1.;

  inf_coeff[index * 10 + x] = -1.;
  sup_coeff[index * 10 + x] = -1.;
}

void get_lb_using_previous_layers(elina_manager_t *man,
                                  const fppoly_t *const fp) {
  const size_t numlayers = fp->numlayers;
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

  const size_t num_out_neurons_last_layer = 90;

  const size_t num_in_neurons_first_layer = fp->layers[0]->num_in_neurons;

  float_type *lb_dev;
  cudaMalloc((void **)&lb_dev, num_out_neurons_last_layer * sizeof(float_type));

  float_type *linf_coeff;
  float_type *lsup_coeff;
  float_type *linf_cst;
  float_type *lsup_cst;

  cudaMalloc((void **)&linf_coeff,
             num_out_neurons_last_layer * 10 * sizeof(float_type *));
  cudaMalloc((void **)&lsup_coeff,
             num_out_neurons_last_layer * 10 * sizeof(float_type *));
  cudaMalloc((void **)&linf_cst,
             num_out_neurons_last_layer * 10 * sizeof(float_type));
  cudaMalloc((void **)&lsup_cst,
             num_out_neurons_last_layer * 10 * sizeof(float_type));

  size_t index = 0;

  for (elina_dim_t y = 0; y < 10; y++) {
    for (elina_dim_t x = 0; x < 10; x++) {
      if (y != x) {
        create_sub_expr<<<1, 1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst,
                                  index, y, x);
        index++;
      }
    }
  }

  float_type *linf_coeff_tmp;
  float_type *lsup_coeff_tmp;
  float_type *linf_cst_tmp;
  float_type *lsup_cst_tmp;

  cudaMalloc((void **)&linf_coeff_tmp,
             num_out_neurons_last_layer * sizeof(float_type *));
  cudaMalloc((void **)&lsup_coeff_tmp,
             num_out_neurons_last_layer * sizeof(float_type *));
  cudaMalloc((void **)&linf_cst_tmp,
             num_out_neurons_last_layer * sizeof(float_type));
  cudaMalloc((void **)&lsup_cst_tmp,
             num_out_neurons_last_layer * sizeof(float_type));

  for (int k = numlayers - 1; k >= 0; k--) {
    const size_t num_out_neurons_current_layer = fp->layers[k]->num_out_neurons;
    const size_t num_in_neurons_current_layer = fp->layers[k]->num_in_neurons;

    const dim3 num_blocks_relu(num_out_neurons_last_layer,
                               num_out_neurons_current_layer / num_threads + 1,
                               1);
    const dim3 num_blocks_linear(num_out_neurons_last_layer,
                                 num_in_neurons_current_layer / num_threads + 1,
                                 1);

    float_type *aux_coeffs;
    float_type *aux_csts;

    if (fp->layers[k]->type == CONV) {
      aux_coeffs = fp->layers[k]->filter_weights;
      aux_csts = fp->layers[k]->filter_bias;
    } else {
      aux_coeffs = fp->layers[k]->coeffs;
      aux_csts = fp->layers[k]->csts;
    }

    float_type *aux_lb_array = fp->layers[k]->lb_array;
    float_type *aux_ub_array = fp->layers[k]->ub_array;

    if (fp->layers[k]->activation == RELU) {
      lexpr_replace_relu_bounds<<<num_blocks_relu, num_threads>>>(
          linf_coeff, lsup_coeff, linf_cst, lsup_cst, aux_lb_array,
          aux_ub_array, num_out_neurons_current_layer);
    }

    cudaMalloc((void **)&linf_coeff_tmp, num_out_neurons_last_layer *
                                             num_in_neurons_current_layer *
                                             sizeof(float_type));
    cudaMalloc((void **)&lsup_coeff_tmp, num_out_neurons_last_layer *
                                             num_in_neurons_current_layer *
                                             sizeof(float_type));

    cudaMemset(linf_coeff_tmp, 0,
               num_out_neurons_last_layer * num_in_neurons_current_layer *
                   sizeof(float_type));
    cudaMemset(lsup_coeff_tmp, 0,
               num_out_neurons_last_layer * num_in_neurons_current_layer *
                   sizeof(float_type));

    if (fp->layers[k]->type == CONV) {
      coeffs_from_previous_layer_conv<<<num_out_neurons_last_layer,
                                        dim3(fp->layers[k]->input_size[2],
                                             fp->layers[k]->filter_size[1],
                                             fp->layers[k]->filter_size[0])>>>(
          linf_coeff, lsup_coeff, linf_coeff_tmp, lsup_coeff_tmp, aux_coeffs,
          fp->layers[k]->output_size[0], fp->layers[k]->output_size[1],
          fp->layers[k]->output_size[2], fp->layers[k]->input_size[0],
          fp->layers[k]->input_size[1], fp->layers[k]->input_size[2],
          fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1],
          fp->layers[k]->strides[0], fp->layers[k]->strides[1],
          fp->layers[k]->pad[0], fp->layers[k]->pad[1]);

      csts_from_previous_layer_conv<<<num_out_neurons_last_layer, 1>>>(
          linf_coeff, lsup_coeff, linf_cst, lsup_cst, linf_cst_tmp,
          lsup_cst_tmp, aux_csts, fp->layers[k]->output_size[0],
          fp->layers[k]->output_size[1], fp->layers[k]->output_size[2]);
    } else {
      coeffs_from_previous_layer<<<num_blocks_linear, num_threads>>>(
          linf_coeff, lsup_coeff, linf_coeff_tmp, lsup_coeff_tmp, aux_coeffs,
          num_out_neurons_current_layer, num_in_neurons_current_layer);

      csts_from_previous_layer<<<num_out_neurons_last_layer, 1>>>(
          linf_coeff, lsup_coeff, linf_cst, lsup_cst, linf_cst_tmp,
          lsup_cst_tmp, aux_csts, num_out_neurons_current_layer);
    }

    std::swap(linf_coeff, linf_coeff_tmp);
    std::swap(lsup_coeff, lsup_coeff_tmp);
    std::swap(linf_cst, linf_cst_tmp);
    std::swap(lsup_cst, lsup_cst_tmp);

    cudaFree(linf_coeff_tmp);
    cudaFree(lsup_coeff_tmp);

    if (fp->layers[k]->type == CONV) {
      cudaFree(aux_coeffs);
      cudaFree(aux_csts);
    }
  }

  compute_lb_from_expr<<<num_out_neurons_last_layer, 1>>>(
      lb_dev, linf_coeff, lsup_coeff, linf_cst, fp->input_inf, fp->input_sup,
      num_in_neurons_first_layer);

  cudaFree(linf_coeff);
  cudaFree(lsup_coeff);
  cudaFree(linf_cst);
  cudaFree(lsup_cst);

  cudaFree(linf_cst_tmp);
  cudaFree(lsup_cst_tmp);

  float_type lb[num_out_neurons_last_layer];
  cudaMemcpy(&lb, lb_dev, num_out_neurons_last_layer * sizeof(float_type),
             cudaMemcpyDeviceToHost);

  cudaFree(lb_dev);

  for (size_t i = 0; i < num_out_neurons_last_layer; i++) {
    if (lb[i] > 0) {
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

void conv_add_layer(fppoly_t *const fp, const size_t num_out_neurons,
                    const size_t num_in_neurons,
                    const size_t num_nonzero_weights, const size_t num_biases,
                    const layertype_t type,
                    const activation_type_t activation) {
  layer_t *layer = (layer_t *)malloc(sizeof(layer_t));

  layer->num_out_neurons = num_out_neurons;
  layer->num_in_neurons = num_in_neurons;

  layer->type = type;
  layer->activation = activation;

  cudaMalloc((void **)&layer->lb_array, num_out_neurons * sizeof(float_type));
  cudaMalloc((void **)&layer->ub_array, num_out_neurons * sizeof(float_type));

  layer->coeffs = nullptr;
  layer->csts = nullptr;

  cudaMalloc((void **)&layer->filter_weights,
             num_nonzero_weights * sizeof(float_type));
  cudaMalloc((void **)&layer->filter_bias, num_biases * sizeof(float_type));
  cudaMemset(layer->filter_bias, 0, num_biases * sizeof(float_type));

  layer->input_size = (size_t *)malloc(3 * sizeof(size_t));
  layer->output_size = (size_t *)malloc(3 * sizeof(size_t));
  layer->filter_size = (size_t *)malloc(2 * sizeof(size_t));
  layer->strides = (size_t *)malloc(2 * sizeof(size_t));
  layer->pad = (long int *)malloc(2 * sizeof(long int));

  fp->layers[fp->numlayers] = layer;

  fp->numlayers++;
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
  const size_t size =
      filter_size[0] * filter_size[1] * input_size[2] * output_size[2];

  conv_add_layer(fp, num_out_neurons, num_pixels, size, output_size[2], CONV,
                 RELU);

  layer_t *current_layer = fp->layers[fp->numlayers - 1];

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

  const long int pad[2] = {pad_top, pad_left};

  float_type *filter_weights_tmp =
      (float_type *)malloc(size * sizeof(float_type));
  float_type *filter_bias_tmp = (float_type *)malloc(size * sizeof(float_type));

  for (size_t out_z = 0; out_z < output_size[2]; out_z++) {
    for (size_t inp_z = 0; inp_z < input_size[2]; inp_z++) {
      for (size_t x_shift = 0; x_shift < filter_size[0]; x_shift++) {
        for (size_t y_shift = 0; y_shift < filter_size[1]; y_shift++) {
          const size_t read_index =
              x_shift * filter_size[1] * input_size[2] * output_size[2] +
              y_shift * input_size[2] * output_size[2] +
              inp_z * output_size[2] + out_z;
          const size_t write_index =
              out_z * filter_size[0] * filter_size[1] * input_size[2] +
              x_shift * filter_size[1] * input_size[2] +
              y_shift * input_size[2] + inp_z;
          filter_weights_tmp[write_index] = filter_weights[read_index];
        }
      }
    }

    filter_bias_tmp[out_z] = filter_bias[out_z];
  }

  cudaMemcpy(current_layer->filter_weights, filter_weights_tmp,
             size * sizeof(float_type), cudaMemcpyHostToDevice);

  if (has_bias) {
    cudaMemcpy(current_layer->filter_bias, filter_bias_tmp,
               output_size[2] * sizeof(float_type), cudaMemcpyHostToDevice);
  }

  cudaMemcpy(current_layer->input_size, input_size, 3 * sizeof(size_t),
             cudaMemcpyHostToHost);
  cudaMemcpy(current_layer->output_size, output_size, 3 * sizeof(size_t),
             cudaMemcpyHostToHost);
  cudaMemcpy(current_layer->filter_size, filter_size, 2 * sizeof(size_t),
             cudaMemcpyHostToHost);
  cudaMemcpy(current_layer->strides, strides, 2 * sizeof(size_t),
             cudaMemcpyHostToHost);
  cudaMemcpy(current_layer->pad, pad, 2 * sizeof(long int),
             cudaMemcpyHostToHost);

  free(filter_weights_tmp);
  free(filter_bias_tmp);
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

  layer_compute_bounds_from_exprs_conv<<<dim3(fp->layers[0]->output_size[0],
                                              fp->layers[0]->output_size[1],
                                              fp->layers[0]->output_size[2]),
                                         1>>>(
      fp->layers[0]->filter_weights, fp->layers[0]->filter_bias,
      fp->layers[0]->lb_array, fp->layers[0]->ub_array, fp->input_inf,
      fp->input_sup, fp->layers[0]->output_size[0],
      fp->layers[0]->output_size[1], fp->layers[0]->output_size[2],
      fp->layers[0]->input_size[0], fp->layers[0]->input_size[1],
      fp->layers[0]->input_size[2], fp->layers[0]->filter_size[0],
      fp->layers[0]->filter_size[1], fp->layers[0]->strides[0],
      fp->layers[0]->strides[1], fp->layers[0]->pad[0], fp->layers[0]->pad[1]);
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

  update_state_using_previous_layers_conv(man, fp, fp->numlayers - 1);
}

void free_layer(layer_t *layer) {
  cudaFree(layer->coeffs);
  cudaFree(layer->csts);

  layer->coeffs = nullptr;
  layer->csts = nullptr;

  cudaFree(layer->lb_array);
  cudaFree(layer->ub_array);

  layer->lb_array = nullptr;
  layer->ub_array = nullptr;

  cudaFree(layer->filter_weights);
  cudaFree(layer->filter_bias);

  layer->filter_weights = nullptr;
  layer->filter_bias = nullptr;

  free(layer->input_size);
  free(layer->output_size);
  free(layer->filter_size);
  free(layer->strides);
  free(layer->pad);

  layer->input_size = nullptr;
  layer->output_size = nullptr;
  layer->filter_size = nullptr;
  layer->strides = nullptr;
  layer->pad = nullptr;

  free(layer);
  layer = nullptr;
}

void fppoly_free(elina_manager_t *man, fppoly_t *fp) {
  for (size_t i = 0; i < fp->numlayers; i++) {
    free_layer(fp->layers[i]);
  }

  free(fp->layers);
  fp->layers = nullptr;

  cudaFree(fp->input_inf);
  fp->input_inf = nullptr;
  cudaFree(fp->input_sup);
  fp->input_sup = nullptr;

  free(fp);
  fp = nullptr;
}

void layer_print(const layer_t *layer) {
  // neurons_print<<<1, 1>>>(layer->neurons, layer->num_out_neurons);
}

void fppoly_fprint(FILE *const stream, elina_manager_t *man,
                   const fppoly_t *const fp, const char **name_of_dim) {
  for (size_t i = 0; i < fp->numlayers; i++) {
    printf("layer: %zu\n", i);
    layer_print(fp->layers[i]);
  }
}

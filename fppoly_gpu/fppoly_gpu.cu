/**
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2020 Department of Computer Science, ETH Zurich
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
 *  @file fppoly_gpu.cu
 *  @author Christoph Müller
 *  @brief Main file of GPUPoly, all control logic is here.
 */

#include <iostream>
#include <cfenv>
#include <cuda.h>

#include "fppoly_gpu.h"

#ifdef single
__constant__ const float_type min_denormal = 1.40129846e-45;
__constant__ const float_type ulp =          1.1920929e-07;
#else
__constant__ const float_type min_denormal = 4.940656458412465441766e-324;
__constant__ const float_type ulp =          2.220446049250313080848e-16;
#endif

#include "cuda_util.h"
#include "rounding.h"
#include "relu_approx.h"
#include "linear.h"
#include "maxpool_convex_hull.h"

size_t maximum_backstep = 2;

const size_t num_threads = 256;

#define FULL_MASK 0xffffffff

bool results[10][9];
bool results_calculated;

size_t output_counter_x;
size_t output_counter_y;

int num_neurons_training_layer = 0;

float_type* adv_host = nullptr;

void clean_training_data()
{
    free(adv_host);

    adv_host = nullptr;

    num_neurons_training_layer = 0;
}


fppoly_t* fppoly_of_abstract0(elina_abstract0_t* a)
{
    return (fppoly_t*) a->value;
}


elina_abstract0_t* abstract0_of_fppoly(elina_manager_t* man, fppoly_t* fp)
{
    elina_abstract0_t* r = (elina_abstract0_t*) malloc(sizeof(elina_abstract0_t));
    assert(r);
    r->value = fp;
    r->man = elina_manager_copy(man);

    return r;
}


static inline void fppoly_internal_free(fppoly_internal_t* pr)
{
    if (pr)
    {
        pr->funid = ELINA_FUNID_UNKNOWN;
        free(pr);
        pr = nullptr;
    }
}


static inline fppoly_internal_t* fppoly_internal_alloc()
{
    fppoly_internal_t* pr = (fppoly_internal_t*) malloc(sizeof(fppoly_internal_t));
    pr->funid = ELINA_FUNID_UNKNOWN;
    pr->man = nullptr;
    pr->funopt = nullptr;

    return pr;
}


/* back pointer to our internal structure from the manager */
fppoly_internal_t* fppoly_init_from_manager(elina_manager_t* man, elina_funid_t funid)
{
    fppoly_internal_t* pr = (fppoly_internal_t*) man->internal;
    pr->funid = funid;

    if (!(pr->man))
    {
        pr->man = man;
    }

    return pr;
}


elina_manager_t* fppoly_manager_alloc(const size_t gpu_number, const size_t backstep_depth)
{
    results_calculated = false;

    auto status = cudaSetDevice(gpu_number);

    if(status == cudaErrorInvalidDevice)
    {
        std::cerr << "error: CUDA API call : " << cudaGetErrorString(status) << std::endl;
        exit(1);
    }

    if(status == cudaErrorSetOnActiveProcess)
    {
        std::cerr << "error: CUDA API call : " << cudaGetErrorString(status) << std::endl;
        exit(1);
    }

    maximum_backstep = backstep_depth;

    output_counter_x = 0;
    output_counter_y = 0;

    void** funptr;
    fppoly_internal_t* pr = fppoly_internal_alloc();

    elina_manager_t* man = elina_manager_alloc("fppoly",/* Library name */
            "1.0", /* version */
            pr, /* internal structure */
            (void (*)(void*)) fppoly_internal_free /* free function for internal */
            );

    funptr = man->funptr;
    funptr[ELINA_FUNID_FREE] = (void*) &fppoly_free;
    /* 3.Printing */
    funptr[ELINA_FUNID_FPRINT] = (void*) &fppoly_fprint;

    return man;
}


template <typename T>
__inline__ __device__
void warp_reduce_sum(T& val)
{
    for(int offset = warpSize/2; offset > 0; offset /= 2)
    {
        val = val + __shfl_down_sync(FULL_MASK, val, offset);
    }
}


template <typename T>
__inline__ __device__
void block_reduce_sum(T& val, const int number_threads)
{
    if(number_threads <= warpSize)
    {
        warp_reduce_sum(val);
    }
    else
    {
        static __shared__ T shared[32];
        const int lane = threadIdx.x%warpSize;
        const int wid = threadIdx.x/warpSize;

        warp_reduce_sum(val);

        if(lane == 0) shared[wid] = val;

        __syncthreads();

        val = (threadIdx.x < (blockDim.x/warpSize)) ? shared[lane] : 0;

        if(wid == 0)
        {
            warp_reduce_sum(val);
        }
    }
}


template <typename T>
__inline__ __device__
void warp_reduce_sum_rd(T& val)
{
    for(int offset = warpSize/2; offset > 0; offset /= 2)
    {
        val = add_rd(val, __shfl_down_sync(FULL_MASK, val, offset));
    }
}


template <typename T>
__inline__ __device__
void block_reduce_sum_rd(T& val, const int number_threads)
{
    if(number_threads <= warpSize)
    {
        warp_reduce_sum_rd(val);
    }
    else
    {
        static __shared__ T shared[32];
        const int lane = threadIdx.x%warpSize;
        const int wid = threadIdx.x/warpSize;

        warp_reduce_sum_rd(val);

        if(lane == 0) shared[wid] = val;

        __syncthreads();

        val = (threadIdx.x < (blockDim.x/warpSize)) ? shared[lane] : 0;

        if(wid == 0)
        {
            warp_reduce_sum_rd(val);
        }
    }
}


template <typename T>
__inline__ __device__
void warp_reduce_sum_ru(T& val)
{
    for(int offset = warpSize/2; offset > 0; offset /= 2)
    {
        val = add_ru(val, __shfl_down_sync(FULL_MASK, val, offset));
    }
}


template <typename T>
__inline__ __device__
void block_reduce_sum_ru(T& val, const int number_threads)
{
    if(number_threads <= warpSize)
    {
        warp_reduce_sum_ru(val);
    }
    else
    {
        static __shared__ T shared[32];
        const int lane = threadIdx.x%warpSize;
        const int wid = threadIdx.x/warpSize;

        warp_reduce_sum_ru(val);

        if(lane == 0) shared[wid] = val;

        __syncthreads();

        val = (threadIdx.x < (blockDim.x/warpSize)) ? shared[lane] : 0;

        if(wid == 0)
        {
            warp_reduce_sum_ru(val);
        }
    }
}


template <typename T>
__inline__ __device__
void warp_reduce_sum_csts(T& lval, T& uval)
{
    for(int offset = warpSize/2; offset > 0; offset /= 2)
    {
        const T lval_shfl = __shfl_down_sync(FULL_MASK, lval, offset);
        const T uval_shfl = __shfl_down_sync(FULL_MASK, uval, offset);

        T maxVal;
        T maxShfl;

        add_cst(lval, uval, maxVal, maxShfl, lval_shfl, uval_shfl);
    }
}


template <typename T>
__inline__ __device__
void block_reduce_sum_csts(T& lval, T& uval, const int number_threads)
{
    if(number_threads <= warpSize)
    {
        warp_reduce_sum_csts(lval, uval);
    }
    else
    {
        static __shared__ T lshared[32];
        static __shared__ T ushared[32];
        const int lane = threadIdx.x%warpSize;
        const int wid = threadIdx.x/warpSize;

        warp_reduce_sum_csts(lval, uval);

        if(lane == 0)
        {
            lshared[wid] = lval;
            ushared[wid] = uval;
        }

        __syncthreads();

        lval = (threadIdx.x < (blockDim.x/warpSize)) ? lshared[lane] : 0;
        uval = (threadIdx.x < (blockDim.x/warpSize)) ? ushared[lane] : 0;

        if(wid == 0)
        {
            warp_reduce_sum_csts(lval, uval);
        }
    }
}


void fppoly_from_network_input_box(fppoly_t* const res, const size_t intdim, const size_t realdim, const float_type* inf_array, const float_type* sup_array)
{
    res->layers = nullptr;
    res->numlayers = 0;

    const size_t num_pixels = intdim + realdim;

    float_type* tmp_input_inf = (float_type*) malloc(num_pixels*sizeof(float_type));
    float_type* tmp_input_sup = (float_type*) malloc(num_pixels*sizeof(float_type));

    for(size_t i = 0; i < num_pixels; i++)
    {
        tmp_input_inf[i] = inf_array[i];
        tmp_input_sup[i] = sup_array[i];
    }

    res->input_inf = malloc_device<float_type>(num_pixels);
    res->input_sup = malloc_device<float_type>(num_pixels);

    copy_to_device(res->input_inf, tmp_input_inf, num_pixels);
    copy_to_device(res->input_sup, tmp_input_sup, num_pixels);

    free(tmp_input_inf);
    free(tmp_input_sup);

    res->num_pixels = num_pixels;
}


elina_abstract0_t* fppoly_from_network_input(elina_manager_t* man, const size_t intdim, const size_t realdim, const float_type* inf_array, const float_type* sup_array)
{
    fppoly_t* res = (fppoly_t*) malloc(sizeof(fppoly_t));

    fppoly_from_network_input_box(res, intdim, realdim, inf_array, sup_array);

    res->input_lweights = nullptr;
    res->input_uweights = nullptr;

    res->input_lcst = nullptr;
    res->input_ucst = nullptr;

    return abstract0_of_fppoly(man, res);
}


elina_abstract0_t* fppoly_from_network_input_poly(elina_manager_t* man, const size_t intdim, const size_t realdim, const float_type* inf_array, const float_type* sup_array, const float_type* lexpr_weights, const float_type* lexpr_cst, const size_t* lexpr_dim, const float_type* uexpr_weights, const float_type* uexpr_cst, const size_t* uexpr_dim, const size_t expr_size)
{
    fppoly_t* res = (fppoly_t*) malloc(sizeof(fppoly_t));

    fppoly_from_network_input_box(res, intdim, realdim, inf_array, sup_array);

    const size_t num_pixels = intdim + realdim - expr_size;
    res->mu = expr_size;

    float_type* input_lweights_tmp = (float_type*) malloc(num_pixels*expr_size*sizeof(float_type));
    float_type* input_uweights_tmp = (float_type*) malloc(num_pixels*expr_size*sizeof(float_type));

    float_type* input_lcst_tmp = (float_type*) malloc(num_pixels*sizeof(float_type));
    float_type* input_ucst_tmp = (float_type*) malloc(num_pixels*sizeof(float_type));

    for(size_t i = 0; i < num_pixels; i++)
    {
        input_lcst_tmp[i] = lexpr_cst[i];

        for(size_t j = 0; j < expr_size; j++)
        {
            const size_t address = i*expr_size + lexpr_dim[i*expr_size + j] - num_pixels;
            input_lweights_tmp[address] = lexpr_weights[i*expr_size + j];
        }
    }

    for(size_t i = 0; i < num_pixels; i++)
    {
        input_ucst_tmp[i] = uexpr_cst[i];

        for(size_t j = 0; j < expr_size; j++)
        {
            const size_t address = i*expr_size + uexpr_dim[i*expr_size + j] - num_pixels;
            input_uweights_tmp[address] = uexpr_weights[i*expr_size + j];
        }
    }

    res->input_lweights = malloc_device<float_type>(num_pixels*expr_size);
    res->input_uweights = malloc_device<float_type>(num_pixels*expr_size);

    res->input_lcst = malloc_device<float_type>(num_pixels);
    res->input_ucst = malloc_device<float_type>(num_pixels);

    copy_to_device(res->input_lweights, input_lweights_tmp, num_pixels*expr_size);
    copy_to_device(res->input_uweights, input_uweights_tmp, num_pixels*expr_size);

    copy_to_device(res->input_lcst, input_lcst_tmp, num_pixels);
    copy_to_device(res->input_ucst, input_ucst_tmp, num_pixels);

    free(input_lweights_tmp);
    free(input_uweights_tmp);

    free(input_lcst_tmp);
    free(input_ucst_tmp);

    return abstract0_of_fppoly(man, res);
}


void ffn_add_layer(fppoly_t* const fp, const size_t num_out_neurons, const size_t num_in_neurons, const layertype_t type, const activation_type_t activation)
{
    layer_t* layer = (layer_t*) malloc(sizeof(layer_t));

    layer->num_out_neurons = num_out_neurons;
    layer->num_in_neurons = num_in_neurons;

    layer->type = type;
    layer->activation = activation;

    layer->lb_array = malloc_device<float_type>(num_out_neurons);
    layer->ub_array = malloc_device<float_type>(num_out_neurons);

    layer->coeffs = malloc_device<float_type>(num_out_neurons*num_in_neurons);
    layer->csts = malloc_device<float_type>(num_out_neurons);

    layer->filter_weights = nullptr;
    layer->filter_bias = nullptr;

    layer->input_size = nullptr;
    layer->output_size = nullptr;
    layer->filter_size = nullptr;
    layer->strides = nullptr;
    layer->pad = nullptr;

    layer->pool_lcoeffs = nullptr;
    layer->pool_ucoeffs = nullptr;

    layer->pool_lcsts = nullptr;
    layer->pool_ucsts = nullptr;

    fp->layers[fp->numlayers] = layer;

    fp->numlayers++;
}


__global__
void compute_lb_from_expr(float_type* __restrict__ lb_array, const float_type* __restrict__ inf_coeffs, const float_type* __restrict__ sup_coeffs, const float_type* __restrict__ inf_csts, const float_type* __restrict__ input_inf, const float_type* __restrict__ input_sup, const int expr_size)
{
    const int n = blockIdx.x;
    int i = threadIdx.x;

    float_type res_inf = 0;

    float_type tmp1, tmp2;

    while(i < expr_size)
    {
        lb_component(res_inf, tmp1, tmp2, inf_coeffs[n*expr_size + i], sup_coeffs[n*expr_size + i], input_inf[i], input_sup[i]);

        i += blockDim.x;
    }

    block_reduce_sum_rd(res_inf, blockDim.x);

    if(threadIdx.x == 0)
    {
        lb_array[n] = add_rd(res_inf, inf_csts[n]);
    }
}


__global__
void compute_ub_from_expr(float_type* __restrict__ ub_array, const float_type* __restrict__ inf_coeffs, const float_type* __restrict__ sup_coeffs, const float_type* __restrict__ sup_csts, const float_type* __restrict__ input_inf, const float_type* __restrict__ input_sup, const int expr_size)
{
    const int n = blockIdx.x;
    int i = threadIdx.x;

    float_type res_sup = 0;

    float_type tmp1, tmp2;

    while(i < expr_size)
    {
        ub_component(res_sup, tmp1, tmp2, inf_coeffs[n*expr_size + i], sup_coeffs[n*expr_size + i], input_inf[i], input_sup[i]);

        i += blockDim.x;
    }

    block_reduce_sum_ru(res_sup, blockDim.x);

    if(threadIdx.x == 0)
    {
        ub_array[n] = add_ru(res_sup, sup_csts[n]);
    }
}


__global__
void compute_lb_from_expr_conv_sparse(float_type* __restrict__ lb_array, const float_type* __restrict__ inf_coeffs, const float_type* __restrict__ sup_coeffs, const float_type* __restrict__ inf_csts, const float_type* __restrict__ input_inf, const float_type* __restrict__ input_sup, const int num_chunks, const int chunk_counter, const int output_size_x, const int output_size_y, const int output_size_z, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    const int local_n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;
    const int global_n = last_x*gridDim.y*num_chunks*gridDim.z + last_y*num_chunks*gridDim.z + chunk_counter*gridDim.z + last_z;

    int out_z = threadIdx.x;

    const int min_out_x = offset_x + last_x*shift_x;
    const int min_out_y = offset_y + last_y*shift_y;

    const int min_x = (min_out_x < 0) ? -min_out_x : 0;
    const int min_y = (min_out_y < 0) ? -min_out_y : 0;

    const int max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
    const int max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;

    float_type res_inf = 0;

    float_type tmp1, tmp2;

    while(out_z < output_size_z)
    {
        for(int out_x = min_x; out_x < max_x; out_x++)
        {
            for(int out_y = min_y; out_y < max_y; out_y++)
            {
                const int i = (out_x + min_out_x)*output_size_y*output_size_z + (out_y + min_out_y)*output_size_z + out_z;
                const int j = out_x*length_y*output_size_z + out_y*output_size_z + out_z;

                const int mat_out = local_n*length_x*length_y*output_size_z + j;

                float_type inf_coeff = inf_coeffs[mat_out];
                float_type sup_coeff = sup_coeffs[mat_out];

                if((inf_coeff != 0) || (sup_coeff != 0))
                {
                    lb_component(res_inf, tmp1, tmp2, inf_coeff, sup_coeff, input_inf[i], input_sup[i]);
                }
           }
        }

        out_z += blockDim.x;
    }

    block_reduce_sum_rd(res_inf, blockDim.x);

    if(threadIdx.x == 0)
    {
        lb_array[global_n] = add_rd(res_inf, inf_csts[local_n]);
    }
}


__global__
void compute_ub_from_expr_conv_sparse(float_type* __restrict__ ub_array, const float_type* __restrict__ inf_coeffs, const float_type* __restrict__ sup_coeffs, const float_type* __restrict__ sup_csts, const float_type* __restrict__ input_inf, const float_type* __restrict__ input_sup, const int num_chunks, const int chunk_counter, const int output_size_x, const int output_size_y, const int output_size_z, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    const int local_n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;
    const int global_n = last_x*gridDim.y*num_chunks*gridDim.z + last_y*num_chunks*gridDim.z + chunk_counter*gridDim.z + last_z;

    int out_z = threadIdx.x;

    const int min_out_x = offset_x + last_x*shift_x;
    const int min_out_y = offset_y + last_y*shift_y;

    const int min_x = (min_out_x < 0) ? -min_out_x : 0;
    const int min_y = (min_out_y < 0) ? -min_out_y : 0;

    const int max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
    const int max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;

    float_type res_sup = 0;

    float_type tmp1, tmp2;

    while(out_z < output_size_z)
    {
        for(int out_x = min_x; out_x < max_x; out_x++)
        {
            for(int out_y = min_y; out_y < max_y; out_y++)
            {
                const int i = (out_x + min_out_x)*output_size_y*output_size_z + (out_y + min_out_y)*output_size_z + out_z;
                const int j = out_x*length_y*output_size_z + out_y*output_size_z + out_z;

                const int mat_out = local_n*length_x*length_y*output_size_z + j;

                float_type inf_coeff = inf_coeffs[mat_out];
                float_type sup_coeff = sup_coeffs[mat_out];

                if((inf_coeff != 0) || (sup_coeff != 0))
                {
                    ub_component(res_sup, tmp1, tmp2, inf_coeff, sup_coeff, input_inf[i], input_sup[i]);
                }
            }
        }

        out_z += blockDim.x;
    }

    block_reduce_sum_ru(res_sup, blockDim.x);

    if(threadIdx.x == 0)
    {
        ub_array[global_n] = add_ru(res_sup, sup_csts[local_n]);
    }
}


__global__
void compute_lb_from_expr_input_poly_sparse(float_type* __restrict__ lb_array, const float_type* __restrict__ inf_coeffs, const float_type* __restrict__ sup_coeffs, const float_type* __restrict__ inf_csts, const float_type* __restrict__ input_inf, const float_type* __restrict__ input_sup, const int num_chunks, const int chunk_counter, const int mu)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    const int local_n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;
    const int global_n = last_x*gridDim.y*num_chunks*gridDim.z + last_y*num_chunks*gridDim.z + chunk_counter*gridDim.z + last_z;

    float_type res_inf = inf_csts[local_n];

    float_type tmp1, tmp2;

    for(int i = 0; i < mu; i++)
    {
        const int mat_out = local_n*mu + i;

        float_type inf_coeff = inf_coeffs[mat_out];
        float_type sup_coeff = sup_coeffs[mat_out];

        if((inf_coeff != 0) || (sup_coeff != 0))
        {
            lb_component(res_inf, tmp1, tmp2, inf_coeff, sup_coeff, input_inf[i], input_sup[i]);
        }
    }

    lb_array[global_n] = res_inf;
}


__global__
void compute_ub_from_expr_input_poly_sparse(float_type* __restrict__ ub_array, const float_type* __restrict__ inf_coeffs, const float_type* __restrict__ sup_coeffs, const float_type* __restrict__ sup_csts, const float_type* __restrict__ input_inf, const float_type* __restrict__ input_sup, const int num_chunks, const int chunk_counter, const int mu)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    const int local_n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;
    const int global_n = last_x*gridDim.y*num_chunks*gridDim.z + last_y*num_chunks*gridDim.z + chunk_counter*gridDim.z + last_z;

    float_type res_sup = sup_csts[local_n];

    float_type tmp1, tmp2;

    for(int i = 0; i < mu; i++)
    {
        const int mat_out = local_n*mu + i;

        float_type inf_coeff = inf_coeffs[mat_out];
        float_type sup_coeff = sup_coeffs[mat_out];

        if((inf_coeff != 0) || (sup_coeff != 0))
        {
            ub_component(res_sup, tmp1, tmp2, inf_coeff, sup_coeff, input_inf[i], input_sup[i]);
        }
    }

    ub_array[global_n] = res_sup;
}


__global__
void device_layer_create_dense_expr(float_type* __restrict__ coeffs, float_type* __restrict__ csts, const float_type* __restrict__ weights, const float_type* __restrict__ bias, const int num_out_neurons, const int num_in_neurons)
{
    const int i = blockIdx.x;

    const float_type* weight_i = weights + i*num_in_neurons;
    const float_type bias_i = bias[i];

    csts[i] = bias_i;

    for(int j = 0; j < num_in_neurons; j++)
    {
        coeffs[i*num_in_neurons + j] = weight_i[j];
    }
}


void layer_create_dense_exprs(float_type* coeffs, float_type* csts, const float_type* const * weights, const float_type* bias, const size_t num_out_neurons, const size_t num_in_neurons)
{
    float_type* tmp_weights = malloc_device<float_type>(num_out_neurons*num_in_neurons);
    float_type* tmp_bias = malloc_device<float_type>(num_out_neurons);

    for(size_t i = 0; i < num_out_neurons; i++)
    {
        copy_to_device(tmp_weights + i*num_in_neurons, weights[i], num_in_neurons);
    }

    copy_to_device(tmp_bias, bias, num_out_neurons);

    device_layer_create_dense_expr<<<num_out_neurons, 1>>>(coeffs, csts, tmp_weights, tmp_bias, num_out_neurons, num_in_neurons);

    cudaFree(tmp_weights);
    cudaFree(tmp_bias);
}


__global__
void copy_coeffs_and_csts(float_type* __restrict__ target_coeff, float_type* __restrict__ target_cst, const float_type* __restrict__ source_coeffs, const float_type* __restrict__ source_csts, const int coeff_size)
{
    const int i = blockIdx.x;

    for(int j = 0; j < coeff_size; j++)
    {
        target_coeff[i*coeff_size + j] = source_coeffs[i*coeff_size + j];
    }

    target_cst[i] = source_csts[i];
}


__global__
void add_coeffs_and_csts(float_type* __restrict__ target_inf_coeff, float_type* __restrict__ target_sup_coeff, float_type* __restrict__ target_inf_cst, float_type* __restrict__ target_sup_cst, const float_type* __restrict__ source_inf_coeff, const float_type* __restrict__ source_sup_coeff, const float_type* __restrict__ source_inf_cst, const float_type* __restrict__ source_sup_cst, const int expr_size)
{
    const int i = blockIdx.x;

    float_type maxRes;
    float_type maxMul;

    for(int j = 0; j < expr_size; j++)
    {
        add_coeff(target_inf_coeff[i*expr_size + j], target_sup_coeff[i*expr_size + j], maxRes, maxMul, source_inf_coeff[i*expr_size + j], source_sup_coeff[i*expr_size + j]);
    }

    add_cst(target_inf_cst[i], target_sup_cst[i], maxRes, maxMul, source_inf_cst[i], source_sup_cst[i]);
}


__global__
void add_coeffs_and_csts_sparse(float_type* __restrict__ target_inf_coeff, float_type* __restrict__ target_sup_coeff, float_type* __restrict__ target_inf_cst, float_type* __restrict__ target_sup_cst, const float_type* __restrict__ source_inf_coeff, const float_type* __restrict__ source_sup_coeff, const float_type* __restrict__ source_inf_cst, const float_type* __restrict__ source_sup_cst, const int output_size_target_x, const int output_size_target_y, const int output_size_source_x, const int output_size_source_y, const int output_size_z, const int relative_offset_x, const int relative_offset_y)
{
    int last_x = blockIdx.x;
    int last_y = blockIdx.y;
    int last_z = blockIdx.z;

    const int n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;

    const int num_coeffs_target = output_size_target_x*output_size_target_y*output_size_z;
    const int num_coeffs_source = output_size_source_x*output_size_source_y*output_size_z;

    float_type maxRes;
    float_type maxMul;

    for(int out_x = 0; out_x < output_size_source_x; out_x++)
    {
        for(int out_y = 0; out_y < output_size_source_y; out_y++)
        {
            for(int out_z = 0; out_z < output_size_z; out_z++)
            {
                const int target_idx = (out_x + relative_offset_x)*output_size_target_y*output_size_z + (out_y + relative_offset_y)*output_size_z + out_z;
                const int source_idx = out_x*output_size_source_y*output_size_z + out_y*output_size_z + out_z;

                add_coeff(target_inf_coeff[n*num_coeffs_target + target_idx], target_sup_coeff[n*num_coeffs_target + target_idx], maxRes, maxMul, source_inf_coeff[n*num_coeffs_source + source_idx], source_sup_coeff[n*num_coeffs_source + source_idx]);
            }
        }
    }

    add_cst(target_inf_cst[n], target_sup_cst[n], maxRes, maxMul, source_inf_cst[n], source_sup_cst[n]);
}


__global__
void lcoeffs_from_previous_layer_pool(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, float_type* __restrict__ res_inf_coeff, float_type* __restrict__ res_sup_coeff, const float_type* __restrict__ aux_lcoeffs, const float_type* __restrict__ aux_ucoeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int input_size_x, const int input_size_y, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int n = blockIdx.x;

    int out_z = threadIdx.x;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(out_z < output_size_z)
    {
        for(int out_x = 0; out_x < output_size_x; out_x++)
        {
            for(int out_y = 0; out_y < output_size_y; out_y++)
            {
                for(size_t x_shift = 0; x_shift < filter_size_x; x_shift++)
                {
                    for(size_t y_shift = 0; y_shift < filter_size_y; y_shift++)
                    {
                        const int x_val = out_x*stride_x + x_shift - pad_x;
                        const int y_val = out_y*stride_y + y_shift - pad_y;

                        if(!((y_val < 0) || (y_val >= input_size_y)))
                        {
                            if(!((x_val < 0) || (x_val >= input_size_x)))
                            {
                                const int mat_in = x_val*input_size_y*output_size_z + y_val*output_size_z + out_z;
                                const int b = n*input_size_x*input_size_y*output_size_z + mat_in;

                                const int mat_out = out_x*output_size_y*output_size_z + out_y*output_size_z + out_z;
                                const int a = n*output_size_x*output_size_y*output_size_z + mat_out;

                                const int c = out_x*output_size_y*filter_size_x*filter_size_y*output_size_z + out_y*filter_size_x*filter_size_y*output_size_z + (x_shift*filter_size_y + y_shift)*output_size_z + out_z;
                                const float_type prev_inf_coeff = expr_inf_coeff[a];
                                const float_type prev_sup_coeff = expr_sup_coeff[a];

                                if(prev_inf_coeff > 0)
                                {
                                    affine_trans_coeff(res_inf_coeff[b], res_sup_coeff[b], tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, aux_lcoeffs[c]);
                                }

                                if(prev_sup_coeff < 0)
                                {
                                    affine_trans_coeff(res_inf_coeff[b], res_sup_coeff[b], tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, aux_ucoeffs[c]);
                                }
                            }
                        }
                    }
                }
            }
        }

        out_z += blockDim.x;
    }
}


__global__
void ucoeffs_from_previous_layer_pool(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, float_type* __restrict__ res_inf_coeff, float_type* __restrict__ res_sup_coeff, const float_type* __restrict__ aux_lcoeffs, const float_type* __restrict__ aux_ucoeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int input_size_x, const int input_size_y, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int n = blockIdx.x;

    int out_z = threadIdx.x;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(out_z < output_size_z)
    {
        for(int out_x = 0; out_x < output_size_x; out_x++)
        {
            for(int out_y = 0; out_y < output_size_y; out_y++)
            {
                for(size_t x_shift = 0; x_shift < filter_size_x; x_shift++)
                {
                    for(size_t y_shift = 0; y_shift < filter_size_y; y_shift++)
                    {
                        const int x_val = out_x*stride_x + x_shift - pad_x;
                        const int y_val = out_y*stride_y + y_shift - pad_y;

                        if(!((y_val < 0) || (y_val >= input_size_y)))
                        {
                            if(!((x_val < 0) || (x_val >= input_size_x)))
                            {
                                const int mat_in = x_val*input_size_y*output_size_z + y_val*output_size_z + out_z;
                                const int b = n*input_size_x*input_size_y*output_size_z + mat_in;

                                const int mat_out = out_x*output_size_y*output_size_z + out_y*output_size_z + out_z;
                                const int a = n*output_size_x*output_size_y*output_size_z + mat_out;

                                const int c = out_x*output_size_y*filter_size_x*filter_size_y*output_size_z + out_y*filter_size_x*filter_size_y*output_size_z + (x_shift*filter_size_y + y_shift)*output_size_z + out_z;
                                const float_type prev_inf_coeff = expr_inf_coeff[a];
                                const float_type prev_sup_coeff = expr_sup_coeff[a];

                                if(prev_inf_coeff > 0)
                                {
                                    affine_trans_coeff(res_inf_coeff[b], res_sup_coeff[b], tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, aux_ucoeffs[c]);
                                }

                                if(prev_sup_coeff < 0)
                                {
                                    affine_trans_coeff(res_inf_coeff[b], res_sup_coeff[b], tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, aux_lcoeffs[c]);
                                }
                            }
                        }
                    }
                }
            }
        }

        out_z += blockDim.x;
    }
}


__global__
void lcoeffs_from_input_poly(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, float_type* __restrict__ res_inf_coeff, float_type* __restrict__ res_sup_coeff, const float_type* __restrict__ input_lcoeffs, const float_type* __restrict__ input_ucoeffs, const int num_out_neurons_current_layer, const int num_in_neurons_current_layer)
{
    const int n = blockIdx.x;

    int j = threadIdx.x;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(j < num_in_neurons_current_layer)
    {
        const int b = n*num_in_neurons_current_layer + j;

        float_type inf_coeff = 0;
        float_type sup_coeff = 0;

        for(int i = 0; i < num_out_neurons_current_layer; i++)
        {
            const int a = n*num_out_neurons_current_layer + i;
            const int c = i*num_in_neurons_current_layer + j;

            const float_type prev_inf_coeff = expr_inf_coeff[a];
            const float_type prev_sup_coeff = expr_sup_coeff[a];

            if(prev_inf_coeff > 0)
            {
                affine_trans_coeff(inf_coeff, sup_coeff, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_lcoeffs[c]);
            }

            if(prev_sup_coeff < 0)
            {
                affine_trans_coeff(inf_coeff, sup_coeff, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_ucoeffs[c]);
            }
        }

        res_inf_coeff[b] = inf_coeff;
        res_sup_coeff[b] = sup_coeff;

        j += blockDim.x;
    }
}


__global__
void ucoeffs_from_input_poly(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, float_type* __restrict__ res_inf_coeff, float_type* __restrict__ res_sup_coeff, const float_type* __restrict__ input_lcoeffs, const float_type* __restrict__ input_ucoeffs, const int num_out_neurons_current_layer, const int num_in_neurons_current_layer)
{
    const int n = blockIdx.x;

    int j = threadIdx.x;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(j < num_in_neurons_current_layer)
    {
        const int b = n*num_in_neurons_current_layer + j;

        float_type inf_coeff = 0;
        float_type sup_coeff = 0;

        for(int i = 0; i < num_out_neurons_current_layer; i++)
        {
            const int a = n*num_out_neurons_current_layer + i;
            const int c = i*num_in_neurons_current_layer + j;

            const float_type prev_inf_coeff = expr_inf_coeff[a];
            const float_type prev_sup_coeff = expr_sup_coeff[a];

            if(prev_inf_coeff > 0)
            {
                affine_trans_coeff(inf_coeff, sup_coeff, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_ucoeffs[c]);
            }

            if(prev_sup_coeff < 0)
            {
                affine_trans_coeff(inf_coeff, sup_coeff, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_lcoeffs[c]);
            }
        }

        res_inf_coeff[b] = inf_coeff;
        res_sup_coeff[b] = sup_coeff;

        j += blockDim.x;
    }
}


__global__
void lcsts_from_input_poly(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, const float_type* __restrict__ expr_inf_cst, const float_type* __restrict__ expr_sup_cst, float_type* __restrict__ res_inf_cst, float_type* __restrict__ res_sup_cst, const float_type* __restrict__ input_lcsts, const float_type* __restrict__ input_ucsts, const float_type* __restrict__ input_inf, const float_type* __restrict__ input_sup, const int num_out_neurons_current_layer)
{
    const int n = blockIdx.x;

    float_type inf_cst = expr_inf_cst[n];
    float_type sup_cst = expr_sup_cst[n];

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    for(int i = 0; i < num_out_neurons_current_layer; i++)
    {
        const int a = n*num_out_neurons_current_layer + i;

        const float_type prev_inf_coeff = expr_inf_coeff[a];
        const float_type prev_sup_coeff = expr_sup_coeff[a];

        if(prev_inf_coeff > 0)
        {
            affine_trans_cst(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_lcsts[i]);
        }
        else if(prev_sup_coeff < 0)
        {
            affine_trans_cst(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_ucsts[i]);
        }
        else
        {
            lcst_input_poly_neutral(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_inf[i], input_sup[i]);
        }
    }

    res_inf_cst[n] = inf_cst;
    res_sup_cst[n] = sup_cst;
}


__global__
void ucsts_from_input_poly(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, const float_type* __restrict__ expr_inf_cst, const float_type* __restrict__ expr_sup_cst, float_type* __restrict__ res_inf_cst, float_type* __restrict__ res_sup_cst, const float_type* __restrict__ input_lcsts, const float_type* __restrict__ input_ucsts, const float_type* __restrict__ input_inf, const float_type* __restrict__ input_sup, const int num_out_neurons_current_layer)
{
    const int n = blockIdx.x;

    float_type inf_cst = expr_inf_cst[n];
    float_type sup_cst = expr_sup_cst[n];

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    for(int i = 0; i < num_out_neurons_current_layer; i++)
    {
        const int a = n*num_out_neurons_current_layer + i;

        const float_type prev_inf_coeff = expr_inf_coeff[a];
        const float_type prev_sup_coeff = expr_sup_coeff[a];

        if(prev_inf_coeff > 0)
        {
            affine_trans_cst(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_ucsts[i]);
        }
        else if(prev_sup_coeff < 0)
        {
            affine_trans_cst(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_lcsts[i]);
        }
        else
        {
            ucst_input_poly_neutral(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_inf[i], input_sup[i]);
        }
    }

    res_inf_cst[n] = inf_cst;
    res_sup_cst[n] = sup_cst;
}


void ffn_handle_first_layer(elina_manager_t* man, elina_abstract0_t* abs, const float_type* const * weights, const float_type* cst, const size_t size, const size_t num_pixels, size_t* predecessors, const activation_type_t activation, const bool alloc, fnn_op OP)
{
    fppoly_t* res = fppoly_of_abstract0(abs);
    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    res->layers = (layer_t**) malloc(2000*sizeof(layer_t*));
    ffn_add_layer(res, size, num_pixels, FFN, activation);
	res->layers[0]->predecessors = predecessors;

    float_type* coeffs = res->layers[0]->coeffs;
    float_type* csts = res->layers[0]->csts;

    if(OP == MUL || OP == SUB1 || OP == SUB2)
    {
        float_type** dense_weights = (float_type**) malloc(num_pixels*sizeof(float_type*));

        for(size_t j = 0; j < num_pixels; j++)
        {
            dense_weights[j] = (float_type*) calloc(num_pixels, sizeof(float_type));
        }

        if(OP == MUL)
        {
            float_type* dense_cst = (float_type*) calloc(num_pixels, sizeof(float_type));

            for(size_t j = 0; j < num_pixels; j++)
            {
                dense_weights[j][j] = cst[j];
            }

            layer_create_dense_exprs(coeffs, csts, dense_weights, dense_cst, num_pixels, num_pixels);

            free(dense_cst);
        }
        else if(OP == SUB1)
        {
            for(size_t j = 0; j < num_pixels; j++)
            {
                dense_weights[j][j] = -1;
            }

            layer_create_dense_exprs(coeffs, csts, dense_weights, cst, num_pixels, num_pixels);
        }
        else if(OP == SUB2)
        {
            float_type* dense_cst = (float_type*) calloc(num_pixels, sizeof(float_type));

            for(size_t j = 0; j < num_pixels; j++)
            {
                dense_weights[j][j] = 1;
            }

            for(size_t j = 0; j < num_pixels; j++)
            {
                dense_cst[j] = -cst[j];
            }

            layer_create_dense_exprs(coeffs, csts, dense_weights, dense_cst, num_pixels, num_pixels);

            free(dense_cst);
        }

        for(size_t j = 0; j < num_pixels; j++)
        {
            free(dense_weights[j]);
        }

        free(dense_weights);
    }
    else
    {
        layer_create_dense_exprs(coeffs, csts, weights, cst, size, num_pixels);
    }

    if((res->input_lweights != nullptr) && (res->input_uweights != nullptr) && (res->input_lcst != nullptr) && (res->input_ucst != nullptr))
    {
        const size_t num_out_neurons_0_layer = res->layers[0]->num_out_neurons;
        const size_t num_in_neurons_0_layer = res->layers[0]->num_in_neurons;
        const size_t mu = res->mu;

        float_type* linf_coeff = malloc_device<float_type>(num_out_neurons_0_layer*mu);
        float_type* lsup_coeff = malloc_device<float_type>(num_out_neurons_0_layer*mu);
        float_type* uinf_coeff = malloc_device<float_type>(num_out_neurons_0_layer*mu);
        float_type* usup_coeff = malloc_device<float_type>(num_out_neurons_0_layer*mu);

        float_type* linf_cst = malloc_device<float_type>(num_out_neurons_0_layer);
        float_type* lsup_cst = malloc_device<float_type>(num_out_neurons_0_layer);
        float_type* uinf_cst = malloc_device<float_type>(num_out_neurons_0_layer);
        float_type* usup_cst = malloc_device<float_type>(num_out_neurons_0_layer);

        cudaMemset(linf_coeff, 0, num_out_neurons_0_layer*mu*sizeof(float_type));
        cudaMemset(lsup_coeff, 0, num_out_neurons_0_layer*mu*sizeof(float_type));
        cudaMemset(uinf_coeff, 0, num_out_neurons_0_layer*mu*sizeof(float_type));
        cudaMemset(usup_coeff, 0, num_out_neurons_0_layer*mu*sizeof(float_type));

        lcoeffs_from_input_poly<<<res->layers[0]->num_out_neurons, num_threads>>>(coeffs, coeffs, linf_coeff, lsup_coeff, res->input_lweights, res->input_uweights, num_in_neurons_0_layer, mu);
        ucoeffs_from_input_poly<<<res->layers[0]->num_out_neurons, num_threads>>>(coeffs, coeffs, uinf_coeff, usup_coeff, res->input_lweights, res->input_uweights, num_in_neurons_0_layer, mu);

        lcsts_from_input_poly<<<res->layers[0]->num_out_neurons, 1>>>(coeffs, coeffs, csts, csts, linf_cst, lsup_cst, res->input_lcst, res->input_ucst, res->input_inf, res->input_sup, num_in_neurons_0_layer);
        ucsts_from_input_poly<<<res->layers[0]->num_out_neurons, 1>>>(coeffs, coeffs, csts, csts, uinf_cst, usup_cst, res->input_lcst, res->input_ucst, res->input_inf, res->input_sup, num_in_neurons_0_layer);

        compute_lb_from_expr<<<res->layers[0]->num_out_neurons, num_threads>>>(res->layers[0]->lb_array, linf_coeff, lsup_coeff, linf_cst, res->input_inf + num_in_neurons_0_layer, res->input_sup + num_in_neurons_0_layer, mu);
        compute_ub_from_expr<<<res->layers[0]->num_out_neurons, num_threads>>>(res->layers[0]->ub_array, uinf_coeff, usup_coeff, usup_cst, res->input_inf + num_in_neurons_0_layer, res->input_sup + num_in_neurons_0_layer, mu);

        cudaFree(linf_coeff);
        cudaFree(lsup_coeff);
        cudaFree(uinf_coeff);
        cudaFree(usup_coeff);

        linf_coeff = nullptr;
        lsup_coeff = nullptr;
        uinf_coeff = nullptr;
        usup_coeff = nullptr;
    }
    else
    {
        compute_lb_from_expr<<<res->layers[0]->num_out_neurons, num_threads>>>(res->layers[0]->lb_array, coeffs, coeffs, csts, res->input_inf, res->input_sup, res->layers[0]->num_in_neurons);
        compute_ub_from_expr<<<res->layers[0]->num_out_neurons, num_threads>>>(res->layers[0]->ub_array, coeffs, coeffs, csts, res->input_inf, res->input_sup, res->layers[0]->num_in_neurons);
    }
}


void ffn_handle_first_relu_layer(elina_manager_t* man, elina_abstract0_t* abs, const float_type* const * weights, const float_type* bias, const size_t size, const size_t num_pixels, size_t* predecessors)
{
    ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels, predecessors, RELU, true, MATMULT);
}


void ffn_handle_first_sub_layer(elina_manager_t* man, elina_abstract0_t* abs, const float_type* cst, const bool is_minuend, const size_t size, size_t* predecessors)
{
	if(is_minuend==true)
    {
        ffn_handle_first_layer(man, abs, NULL, cst, size, size, predecessors, NONE, true, SUB1);
	}
	else
    {
        ffn_handle_first_layer(man, abs, NULL, cst, size, size, predecessors, NONE, true, SUB2);
	}
}


void ffn_handle_first_mul_layer(elina_manager_t* man, elina_abstract0_t* abs, const float_type* bias, const size_t size, size_t* predecessors)
{
    ffn_handle_first_layer(man, abs, NULL, bias, size, size, predecessors, NONE, true, MUL);
}


__global__
void expr_replace_relu_bounds(float_type* __restrict__ linf_coeff, float_type* __restrict__ lsup_coeff, float_type* __restrict__ uinf_coeff, float_type* __restrict__ usup_coeff, float_type* __restrict__ linf_cst, float_type* __restrict__ lsup_cst, float_type* __restrict__ uinf_cst, float_type* __restrict__ usup_cst, const float_type* __restrict__ lb_array, const float_type* __restrict__ ub_array, const int num_out_neurons_current_layer, const bool use_area_heuristic)
{
    const int n = blockIdx.x;

    int i = threadIdx.x;

    float_type res_linf_cst = 0;
    float_type res_lsup_cst = 0;
    float_type res_uinf_cst = 0;
    float_type res_usup_cst = 0;

    while(i < num_out_neurons_current_layer)
    {
        const int a = n*num_out_neurons_current_layer + i;

        lcoeff_replace_relu_bounds(linf_coeff[a], lsup_coeff[a], res_linf_cst, res_lsup_cst, lb_array[i], ub_array[i], use_area_heuristic);
        ucoeff_replace_relu_bounds(uinf_coeff[a], usup_coeff[a], res_uinf_cst, res_usup_cst, lb_array[i], ub_array[i], use_area_heuristic);

        i += blockDim.x;
    }

    block_reduce_sum_rd(res_linf_cst, blockDim.x);
    block_reduce_sum_ru(res_lsup_cst, blockDim.x);
    block_reduce_sum_rd(res_uinf_cst, blockDim.x);
    block_reduce_sum_ru(res_usup_cst, blockDim.x);

    if(threadIdx.x == 0)
    {
        linf_cst[n] = add_rd(res_linf_cst, linf_cst[n]);
        lsup_cst[n] = add_ru(res_lsup_cst, lsup_cst[n]);
        uinf_cst[n] = add_rd(res_uinf_cst, uinf_cst[n]);
        usup_cst[n] = add_ru(res_usup_cst, usup_cst[n]);
    }
}


__global__
void expr_replace_relu_bounds_conv_sparse(float_type* __restrict__ linf_coeff, float_type* __restrict__ lsup_coeff, float_type* __restrict__ uinf_coeff, float_type* __restrict__ usup_coeff, float_type* __restrict__ linf_cst, float_type* __restrict__ lsup_cst, float_type* __restrict__ uinf_cst, float_type* __restrict__ usup_cst, const float_type* __restrict__ lb_array, const float_type* __restrict__ ub_array, const int output_size_x, const int output_size_y, const int output_size_z, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y, const bool use_area_heuristic)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    int out_z = threadIdx.x;

    __shared__ int n;

    __shared__ int min_out_x;
    __shared__ int min_out_y;

    __shared__ int min_x;
    __shared__ int min_y;

    __shared__ int max_x;
    __shared__ int max_y;

    if(out_z == 0)
    {
        n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;

        min_out_x = offset_x + last_x*shift_x;
        min_out_y = offset_y + last_y*shift_y;

        min_x = (min_out_x < 0) ? -min_out_x : 0;
        min_y = (min_out_y < 0) ? -min_out_y : 0;

        max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
        max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;
    }

    __syncthreads();

    float_type res_linf_cst = 0;
    float_type res_lsup_cst = 0;
    float_type res_uinf_cst = 0;
    float_type res_usup_cst = 0;

    while(out_z < output_size_z)
    {
        for(int out_x = min_x; out_x < max_x; out_x++)
        {
            for(int out_y = min_y; out_y < max_y; out_y++)
            {
                const int i = (out_x + min_out_x)*output_size_y*output_size_z + (out_y + min_out_y)*output_size_z + out_z;
                const int j = out_x*length_y*output_size_z + out_y*output_size_z + out_z;

                const int a = n*length_x*length_y*output_size_z + j;

                lcoeff_replace_relu_bounds(linf_coeff[a], lsup_coeff[a], res_linf_cst, res_lsup_cst, lb_array[i], ub_array[i], use_area_heuristic);
                ucoeff_replace_relu_bounds(uinf_coeff[a], usup_coeff[a], res_uinf_cst, res_usup_cst, lb_array[i], ub_array[i], use_area_heuristic);
            }
        }

        out_z += blockDim.x;
    }

    block_reduce_sum_rd(res_linf_cst, blockDim.x);
    block_reduce_sum_ru(res_lsup_cst, blockDim.x);
    block_reduce_sum_rd(res_uinf_cst, blockDim.x);
    block_reduce_sum_ru(res_usup_cst, blockDim.x);

    if(threadIdx.x == 0)
    {
        linf_cst[n] = add_rd(res_linf_cst, linf_cst[n]);
        lsup_cst[n] = add_ru(res_lsup_cst, lsup_cst[n]);
        uinf_cst[n] = add_rd(res_uinf_cst, uinf_cst[n]);
        usup_cst[n] = add_ru(res_usup_cst, usup_cst[n]);
    }
}


__global__
void coeffs_from_previous_layer(const float_type* __restrict__ expr_linf_coeff, const float_type* __restrict__ expr_lsup_coeff, const float_type* __restrict__ expr_uinf_coeff, const float_type* __restrict__ expr_usup_coeff, float_type* __restrict__ res_linf_coeff, float_type* __restrict__ res_lsup_coeff, float_type* __restrict__ res_uinf_coeff, float_type* __restrict__ res_usup_coeff, const float_type* __restrict__ aux_coeffs, const int num_out_neurons_current_layer, const int num_in_neurons_current_layer)
{
    const int n = blockIdx.x;

    int j = threadIdx.x;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(j < num_in_neurons_current_layer)
    {
        const int b = n*num_in_neurons_current_layer + j;

        float_type linf_coeff = 0;
        float_type lsup_coeff = 0;
        float_type uinf_coeff = 0;
        float_type usup_coeff = 0;

        for(int i = 0; i < num_out_neurons_current_layer; i++)
        {
            const int a = n*num_out_neurons_current_layer + i;
            const int c = i*num_in_neurons_current_layer + j;

            const float_type prev_linf_coeff = expr_linf_coeff[a];
            const float_type prev_lsup_coeff = expr_lsup_coeff[a];
            const float_type prev_uinf_coeff = expr_uinf_coeff[a];
            const float_type prev_usup_coeff = expr_usup_coeff[a];

            const bool lnonnull = (prev_linf_coeff != 0) || (prev_lsup_coeff != 0);
            const bool unonnull = (prev_uinf_coeff != 0) || (prev_usup_coeff != 0);

            if(lnonnull || unonnull)
            {
                const float_type aux_coeff = aux_coeffs[c];

                if(lnonnull)
                {
                    affine_trans_coeff(linf_coeff, lsup_coeff, tmp1, tmp2, maxRes, maxMul, prev_linf_coeff, prev_lsup_coeff, aux_coeff);
                }

                if(unonnull)
                {
                    affine_trans_coeff(uinf_coeff, usup_coeff, tmp1, tmp2, maxRes, maxMul, prev_uinf_coeff, prev_usup_coeff, aux_coeff);
                }
            }
        }

        res_linf_coeff[b] = linf_coeff;
        res_lsup_coeff[b] = lsup_coeff;
        res_uinf_coeff[b] = uinf_coeff;
        res_usup_coeff[b] = usup_coeff;

        j += blockDim.x;
    }
}


__global__
void coeffs_from_previous_layer_conv(const float_type* __restrict__ expr_linf_coeff, const float_type* __restrict__ expr_lsup_coeff, const float_type* __restrict__ expr_uinf_coeff, const float_type* __restrict__ expr_usup_coeff, float_type* __restrict__ res_linf_coeff, float_type* __restrict__ res_lsup_coeff, float_type* __restrict__ res_uinf_coeff, float_type* __restrict__ res_usup_coeff, const float_type* __restrict__ aux_coeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int input_size_x, const int input_size_y, const int input_size_z, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int n = blockIdx.x;

    const int inp_z = threadIdx.x;
    const int y_shift = threadIdx.y;
    const int x_shift = threadIdx.z;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    if(x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z < filter_size_x*filter_size_y*input_size_z)
    {
        for(int out_x = 0; out_x < output_size_x; out_x++)
        {
            for(int out_y = 0; out_y < output_size_y; out_y++)
            {
                const int x_val = out_x*stride_x + x_shift - pad_x;
                const int y_val = out_y*stride_y + y_shift - pad_y;

                if(!((y_val < 0) || (y_val >= input_size_y)))
                {
                    if(!((x_val < 0) || (x_val >= input_size_x)))
                    {
                        const int mat_in = x_val*input_size_y*input_size_z + y_val*input_size_z + inp_z;
                        const int b = n*input_size_x*input_size_y*input_size_z + mat_in;

                        float_type linf_coeff = res_linf_coeff[b];
                        float_type lsup_coeff = res_lsup_coeff[b];
                        float_type uinf_coeff = res_uinf_coeff[b];
                        float_type usup_coeff = res_usup_coeff[b];

                        for(int out_z = 0; out_z < output_size_z; out_z++)
                        {
                            const int mat_out = out_x*output_size_y*output_size_z + out_y*output_size_z + out_z;

                            const int a = n*output_size_x*output_size_y*output_size_z + mat_out;

                            const float_type prev_linf_coeff = expr_linf_coeff[a];
                            const float_type prev_lsup_coeff = expr_lsup_coeff[a];
                            const float_type prev_uinf_coeff = expr_uinf_coeff[a];
                            const float_type prev_usup_coeff = expr_usup_coeff[a];

                            const bool lnonnull = (prev_linf_coeff != 0) || (prev_lsup_coeff != 0);
                            const bool unonnull = (prev_uinf_coeff != 0) || (prev_usup_coeff != 0);

                            if(lnonnull || unonnull)
                            {
                                const int filter_index = out_z*filter_size_x*filter_size_y*input_size_z + x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z;
                                const float_type aux_coeff = aux_coeffs[filter_index];

                                if(lnonnull)
                                {
                                    affine_trans_coeff(linf_coeff, lsup_coeff, tmp1, tmp2, maxRes, maxMul, prev_linf_coeff, prev_lsup_coeff, aux_coeff);
                                }

                                if(unonnull)
                                {
                                    affine_trans_coeff(uinf_coeff, usup_coeff, tmp1, tmp2, maxRes, maxMul, prev_uinf_coeff, prev_usup_coeff, aux_coeff);
                                }
                            }
                        }

                        res_linf_coeff[b] = linf_coeff;
                        res_lsup_coeff[b] = lsup_coeff;
                        res_uinf_coeff[b] = uinf_coeff;
                        res_usup_coeff[b] = usup_coeff;
                    }
                }

                __syncthreads();
            }
        }
    }
}


__global__
void coeffs_from_previous_layer_conv_x_filter_serial(const float_type* __restrict__ expr_linf_coeff, const float_type* __restrict__ expr_lsup_coeff, const float_type* __restrict__ expr_uinf_coeff, const float_type* __restrict__ expr_usup_coeff, float_type* __restrict__ res_linf_coeff, float_type* __restrict__ res_lsup_coeff, float_type* __restrict__ res_uinf_coeff, float_type* __restrict__ res_usup_coeff, const float_type* __restrict__ aux_coeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int input_size_x, const int input_size_y, const int input_size_z, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int n = blockIdx.x;

    const int inp_z = threadIdx.x;
    const int y_shift = threadIdx.y;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    if(y_shift*input_size_z + inp_z < filter_size_y*input_size_z)
    {
        for(int out_x = 0; out_x < output_size_x; out_x++)
        {
            for(int out_y = 0; out_y < output_size_y; out_y++)
            {
                for(int x_shift = 0; x_shift < filter_size_x; x_shift++)
                {
                    const int x_val = out_x*stride_x + x_shift - pad_x;
                    const int y_val = out_y*stride_y + y_shift - pad_y;

                    if(!((y_val < 0) || (y_val >= input_size_y)))
                    {
                        if(!((x_val < 0) || (x_val >= input_size_x)))
                        {
                            const int mat_in = x_val*input_size_y*input_size_z + y_val*input_size_z + inp_z;
                            const int b = n*input_size_x*input_size_y*input_size_z + mat_in;

                            float_type linf_coeff = res_linf_coeff[b];
                            float_type lsup_coeff = res_lsup_coeff[b];
                            float_type uinf_coeff = res_uinf_coeff[b];
                            float_type usup_coeff = res_usup_coeff[b];

                            for(int out_z = 0; out_z < output_size_z; out_z++)
                            {
                                const int mat_out = out_x*output_size_y*output_size_z + out_y*output_size_z + out_z;

                                const int a = n*output_size_x*output_size_y*output_size_z + mat_out;

                                const float_type prev_linf_coeff = expr_linf_coeff[a];
                                const float_type prev_lsup_coeff = expr_lsup_coeff[a];
                                const float_type prev_uinf_coeff = expr_uinf_coeff[a];
                                const float_type prev_usup_coeff = expr_usup_coeff[a];

                                const bool lnonnull = (prev_linf_coeff != 0) || (prev_lsup_coeff != 0);
                                const bool unonnull = (prev_uinf_coeff != 0) || (prev_usup_coeff != 0);

                                if(lnonnull || unonnull)
                                {
                                    const int filter_index = out_z*filter_size_x*filter_size_y*input_size_z + x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z;
                                    const float_type aux_coeff = aux_coeffs[filter_index];

                                    if(lnonnull)
                                    {
                                        affine_trans_coeff(linf_coeff, lsup_coeff, tmp1, tmp2, maxRes, maxMul, prev_linf_coeff, prev_lsup_coeff, aux_coeff);
                                    }

                                    if(unonnull)
                                    {
                                        affine_trans_coeff(uinf_coeff, usup_coeff, tmp1, tmp2, maxRes, maxMul, prev_uinf_coeff, prev_usup_coeff, aux_coeff);
                                    }
                                }
                            }

                            res_linf_coeff[b] = linf_coeff;
                            res_lsup_coeff[b] = lsup_coeff;
                            res_uinf_coeff[b] = uinf_coeff;
                            res_usup_coeff[b] = usup_coeff;
                        }
                    }

                    __syncthreads();
                }
            }
        }
    }
}


__global__
void coeffs_from_previous_layer_conv_filter_serial(const float_type* __restrict__ expr_linf_coeff, const float_type* __restrict__ expr_lsup_coeff, const float_type* __restrict__ expr_uinf_coeff, const float_type* __restrict__ expr_usup_coeff, float_type* __restrict__ res_linf_coeff, float_type* __restrict__ res_lsup_coeff, float_type* __restrict__ res_uinf_coeff, float_type* __restrict__ res_usup_coeff, const float_type* __restrict__ aux_coeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int input_size_x, const int input_size_y, const int input_size_z, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int n = blockIdx.x;

    int inp_z = threadIdx.x;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(inp_z < input_size_z)
    {
        for(int out_x = 0; out_x < output_size_x; out_x++)
        {
            for(int out_y = 0; out_y < output_size_y; out_y++)
            {
                for(int x_shift = 0; x_shift < filter_size_x; x_shift++)
                {
                    for(int y_shift = 0; y_shift < filter_size_y; y_shift++)
                    {
                        const int x_val = out_x*stride_x + x_shift - pad_x;
                        const int y_val = out_y*stride_y + y_shift - pad_y;

                        if(!((y_val < 0) || (y_val >= input_size_y)))
                        {
                            if(!((x_val < 0) || (x_val >= input_size_x)))
                            {
                                const int mat_in = x_val*input_size_y*input_size_z + y_val*input_size_z + inp_z;
                                const int b = n*input_size_x*input_size_y*input_size_z + mat_in;

                                float_type linf_coeff = res_linf_coeff[b];
                                float_type lsup_coeff = res_lsup_coeff[b];
                                float_type uinf_coeff = res_uinf_coeff[b];
                                float_type usup_coeff = res_usup_coeff[b];

                                for(int out_z = 0; out_z < output_size_z; out_z++)
                                {
                                    const int mat_out = out_x*output_size_y*output_size_z + out_y*output_size_z + out_z;

                                    const int a = n*output_size_x*output_size_y*output_size_z + mat_out;

                                    const float_type prev_linf_coeff = expr_linf_coeff[a];
                                    const float_type prev_lsup_coeff = expr_lsup_coeff[a];
                                    const float_type prev_uinf_coeff = expr_uinf_coeff[a];
                                    const float_type prev_usup_coeff = expr_usup_coeff[a];

                                    const bool lnonnull = (prev_linf_coeff != 0) || (prev_lsup_coeff != 0);
                                    const bool unonnull = (prev_uinf_coeff != 0) || (prev_usup_coeff != 0);

                                    if(lnonnull || unonnull)
                                    {
                                        const int filter_index = out_z*filter_size_x*filter_size_y*input_size_z + x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z;
                                        const float_type aux_coeff = aux_coeffs[filter_index];

                                        if(lnonnull)
                                        {
                                            affine_trans_coeff(linf_coeff, lsup_coeff, tmp1, tmp2, maxRes, maxMul, prev_linf_coeff, prev_lsup_coeff, aux_coeff);
                                        }

                                        if(unonnull)
                                        {
                                            affine_trans_coeff(uinf_coeff, usup_coeff, tmp1, tmp2, maxRes, maxMul, prev_uinf_coeff, prev_usup_coeff, aux_coeff);
                                        }
                                    }
                                }

                                res_linf_coeff[b] = linf_coeff;
                                res_lsup_coeff[b] = lsup_coeff;
                                res_uinf_coeff[b] = uinf_coeff;
                                res_usup_coeff[b] = usup_coeff;
                            }
                        }
                    }
                }
            }
        }

        inp_z += blockDim.x;
    }
}


__global__
void coeffs_from_previous_layer_conv_sparse(const float_type* __restrict__ expr_linf_coeff, const float_type* __restrict__ expr_lsup_coeff, const float_type* __restrict__ expr_uinf_coeff, const float_type* __restrict__ expr_usup_coeff, float_type* __restrict__ res_linf_coeff, float_type* __restrict__ res_lsup_coeff, float_type* __restrict__ res_uinf_coeff, float_type* __restrict__ res_usup_coeff, const float_type* __restrict__ aux_coeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int input_size_x, const int input_size_y, const int input_size_z, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y, const int new_length_x, const int new_length_y, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    const int inp_z = threadIdx.x;
    const int y_shift = threadIdx.y;
    const int x_shift = threadIdx.z;

    __shared__ int n;

    __shared__ int min_out_x;
    __shared__ int min_out_y;

    __shared__ int min_x;
    __shared__ int min_y;

    __shared__ int max_x;
    __shared__ int max_y;

    if(x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z == 0)
    {
        n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;

        min_out_x = offset_x + last_x*shift_x;
        min_out_y = offset_y + last_y*shift_y;

        min_x = (min_out_x < 0) ? -min_out_x : 0;
        min_y = (min_out_y < 0) ? -min_out_y : 0;

        max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
        max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;
    }

    __syncthreads();

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    if(x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z < filter_size_x*filter_size_y*input_size_z)
    {
        for(int out_x = min_x; out_x < max_x; out_x++)
        {
            for(int out_y = min_y; out_y < max_y; out_y++)
            {
                const int x_val = out_x*stride_x + min_out_x*stride_x + x_shift - pad_x;
                const int y_val = out_y*stride_y + min_out_y*stride_y + y_shift - pad_y;

                if(!((x_val < 0) || (x_val >= input_size_x)))
                {
                    if(!((y_val < 0) || (y_val >= input_size_y)))
                    {
                        const int mat_in = (out_x*stride_x + x_shift)*new_length_y*input_size_z + (out_y*stride_y + y_shift)*input_size_z + inp_z;
                        const int b = n*new_length_x*new_length_y*input_size_z + mat_in;

                        float_type linf_coeff = res_linf_coeff[b];
                        float_type lsup_coeff = res_lsup_coeff[b];
                        float_type uinf_coeff = res_uinf_coeff[b];
                        float_type usup_coeff = res_usup_coeff[b];

                        for(int out_z = 0; out_z < output_size_z; out_z++)
                        {
                            const int mat_out = out_x*length_y*output_size_z + out_y*output_size_z + out_z;
                            const int a = n*length_x*length_y*output_size_z + mat_out;

                            const float_type prev_linf_coeff = expr_linf_coeff[a];
                            const float_type prev_lsup_coeff = expr_lsup_coeff[a];
                            const float_type prev_uinf_coeff = expr_uinf_coeff[a];
                            const float_type prev_usup_coeff = expr_usup_coeff[a];

                            const bool lnonnull = (prev_linf_coeff != 0) || (prev_lsup_coeff != 0);
                            const bool unonnull = (prev_uinf_coeff != 0) || (prev_usup_coeff != 0);

                            if(lnonnull || unonnull)
                            {
                                const int filter_index = out_z*filter_size_x*filter_size_y*input_size_z + x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z;
                                const float_type aux_coeff = aux_coeffs[filter_index];

                                if(lnonnull)
                                {
                                    affine_trans_coeff(linf_coeff, lsup_coeff, tmp1, tmp2, maxRes, maxMul, prev_linf_coeff, prev_lsup_coeff, aux_coeff);
                                }

                                if(unonnull)
                                {
                                    affine_trans_coeff(uinf_coeff, usup_coeff, tmp1, tmp2, maxRes, maxMul, prev_uinf_coeff, prev_usup_coeff, aux_coeff);
                                }
                            }
                        }

                        res_linf_coeff[b] = linf_coeff;
                        res_lsup_coeff[b] = lsup_coeff;
                        res_uinf_coeff[b] = uinf_coeff;
                        res_usup_coeff[b] = usup_coeff;
                    }
                }

                __syncthreads();
            }
        }
    }
}


__global__
void coeffs_from_previous_layer_conv_sparse_x_filter_serial(const float_type* __restrict__ expr_linf_coeff, const float_type* __restrict__ expr_lsup_coeff, const float_type* __restrict__ expr_uinf_coeff, const float_type* __restrict__ expr_usup_coeff, float_type* __restrict__ res_linf_coeff, float_type* __restrict__ res_lsup_coeff, float_type* __restrict__ res_uinf_coeff, float_type* __restrict__ res_usup_coeff, const float_type* __restrict__ aux_coeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int input_size_x, const int input_size_y, const int input_size_z, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y, const int new_length_x, const int new_length_y, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    const int inp_z = threadIdx.x;
    const int y_shift = threadIdx.y;

    __shared__ int n;

    __shared__ int min_out_x;
    __shared__ int min_out_y;

    __shared__ int min_x;
    __shared__ int min_y;

    __shared__ int max_x;
    __shared__ int max_y;

    if(y_shift*input_size_z + inp_z == 0)
    {
        n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;

        min_out_x = offset_x + last_x*shift_x;
        min_out_y = offset_y + last_y*shift_y;

        min_x = (min_out_x < 0) ? -min_out_x : 0;
        min_y = (min_out_y < 0) ? -min_out_y : 0;

        max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
        max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;
    }

    __syncthreads();

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    if(y_shift*input_size_z + inp_z < filter_size_y*input_size_z)
    {
        for(int out_x = min_x; out_x < max_x; out_x++)
        {
            for(int out_y = min_y; out_y < max_y; out_y++)
            {
                for(int x_shift = 0; x_shift < filter_size_x; x_shift++)
                {
                    const int x_val = out_x*stride_x + min_out_x*stride_x + x_shift - pad_x;
                    const int y_val = out_y*stride_y + min_out_y*stride_y + y_shift - pad_y;

                    if(!((x_val < 0) || (x_val >= input_size_x)))
                    {
                        if(!((y_val < 0) || (y_val >= input_size_y)))
                        {
                            const int mat_in = (out_x*stride_x + x_shift)*new_length_y*input_size_z + (out_y*stride_y + y_shift)*input_size_z + inp_z;
                            const int b = n*new_length_x*new_length_y*input_size_z + mat_in;

                            float_type linf_coeff = res_linf_coeff[b];
                            float_type lsup_coeff = res_lsup_coeff[b];
                            float_type uinf_coeff = res_uinf_coeff[b];
                            float_type usup_coeff = res_usup_coeff[b];

                            for(int out_z = 0; out_z < output_size_z; out_z++)
                            {
                                const int mat_out = out_x*length_y*output_size_z + out_y*output_size_z + out_z;
                                const int a = n*length_x*length_y*output_size_z + mat_out;

                                const float_type prev_linf_coeff = expr_linf_coeff[a];
                                const float_type prev_lsup_coeff = expr_lsup_coeff[a];
                                const float_type prev_uinf_coeff = expr_uinf_coeff[a];
                                const float_type prev_usup_coeff = expr_usup_coeff[a];

                                const bool lnonnull = (prev_linf_coeff != 0) || (prev_lsup_coeff != 0);
                                const bool unonnull = (prev_uinf_coeff != 0) || (prev_usup_coeff != 0);

                                if(lnonnull || unonnull)
                                {
                                    const int filter_index = out_z*filter_size_x*filter_size_y*input_size_z + x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z;
                                    const float_type aux_coeff = aux_coeffs[filter_index];

                                    if(lnonnull)
                                    {
                                        affine_trans_coeff(linf_coeff, lsup_coeff, tmp1, tmp2, maxRes, maxMul, prev_linf_coeff, prev_lsup_coeff, aux_coeff);
                                    }

                                    if(unonnull)
                                    {
                                        affine_trans_coeff(uinf_coeff, usup_coeff, tmp1, tmp2, maxRes, maxMul, prev_uinf_coeff, prev_usup_coeff, aux_coeff);
                                    }
                                }
                            }

                            res_linf_coeff[b] = linf_coeff;
                            res_lsup_coeff[b] = lsup_coeff;
                            res_uinf_coeff[b] = uinf_coeff;
                            res_usup_coeff[b] = usup_coeff;
                        }
                    }

                    __syncthreads();
                }
            }
        }
    }
}


__global__
void coeffs_from_previous_layer_conv_sparse_filter_serial(const float_type* __restrict__ expr_linf_coeff, const float_type* __restrict__ expr_lsup_coeff, const float_type* __restrict__ expr_uinf_coeff, const float_type* __restrict__ expr_usup_coeff, float_type* __restrict__ res_linf_coeff, float_type* __restrict__ res_lsup_coeff, float_type* __restrict__ res_uinf_coeff, float_type* __restrict__ res_usup_coeff, const float_type* __restrict__ aux_coeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int input_size_x, const int input_size_y, const int input_size_z, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y, const int new_length_x, const int new_length_y, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    int inp_z = threadIdx.x;

    __shared__ int n;

    __shared__ int min_out_x;
    __shared__ int min_out_y;

    __shared__ int min_x;
    __shared__ int min_y;

    __shared__ int max_x;
    __shared__ int max_y;

    if(inp_z == 0)
    {
        n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;

        min_out_x = offset_x + last_x*shift_x;
        min_out_y = offset_y + last_y*shift_y;

        min_x = (min_out_x < 0) ? -min_out_x : 0;
        min_y = (min_out_y < 0) ? -min_out_y : 0;

        max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
        max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;
    }

    __syncthreads();

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(inp_z < input_size_z)
    {
        for(int out_x = min_x; out_x < max_x; out_x++)
        {
            for(int out_y = min_y; out_y < max_y; out_y++)
            {
                for(int x_shift = 0; x_shift < filter_size_x; x_shift++)
                {
                    for(int y_shift = 0; y_shift < filter_size_y; y_shift++)
                    {
                        const int x_val = out_x*stride_x + min_out_x*stride_x + x_shift - pad_x;
                        const int y_val = out_y*stride_y + min_out_y*stride_y + y_shift - pad_y;

                        if(!((x_val < 0) || (x_val >= input_size_x)))
                        {
                            if(!((y_val < 0) || (y_val >= input_size_y)))
                            {
                                const int mat_in = (out_x*stride_x + x_shift)*new_length_y*input_size_z + (out_y*stride_y + y_shift)*input_size_z + inp_z;
                                const int b = n*new_length_x*new_length_y*input_size_z + mat_in;

                                float_type linf_coeff = res_linf_coeff[b];
                                float_type lsup_coeff = res_lsup_coeff[b];
                                float_type uinf_coeff = res_uinf_coeff[b];
                                float_type usup_coeff = res_usup_coeff[b];

                                for(int out_z = 0; out_z < output_size_z; out_z++)
                                {
                                    const int mat_out = out_x*length_y*output_size_z + out_y*output_size_z + out_z;
                                    const int a = n*length_x*length_y*output_size_z + mat_out;

                                    const float_type prev_linf_coeff = expr_linf_coeff[a];
                                    const float_type prev_lsup_coeff = expr_lsup_coeff[a];
                                    const float_type prev_uinf_coeff = expr_uinf_coeff[a];
                                    const float_type prev_usup_coeff = expr_usup_coeff[a];

                                    const bool lnonnull = (prev_linf_coeff != 0) || (prev_lsup_coeff != 0);
                                    const bool unonnull = (prev_uinf_coeff != 0) || (prev_usup_coeff != 0);

                                    if(lnonnull || unonnull)
                                    {
                                        const int filter_index = out_z*filter_size_x*filter_size_y*input_size_z + x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z;
                                        const float_type aux_coeff = aux_coeffs[filter_index];

                                        if(lnonnull)
                                        {
                                            affine_trans_coeff(linf_coeff, lsup_coeff, tmp1, tmp2, maxRes, maxMul, prev_linf_coeff, prev_lsup_coeff, aux_coeff);
                                        }

                                        if(unonnull)
                                        {
                                            affine_trans_coeff(uinf_coeff, usup_coeff, tmp1, tmp2, maxRes, maxMul, prev_uinf_coeff, prev_usup_coeff, aux_coeff);
                                        }
                                    }
                                }

                                res_linf_coeff[b] = linf_coeff;
                                res_lsup_coeff[b] = lsup_coeff;
                                res_uinf_coeff[b] = uinf_coeff;
                                res_usup_coeff[b] = usup_coeff;
                            }
                        }
                    }
                }
            }
        }

        inp_z += blockDim.x;
    }
}


__global__
void lcoeffs_from_previous_layer_pool_sparse(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, float_type* __restrict__ res_inf_coeff, float_type* __restrict__ res_sup_coeff, const float_type* __restrict__ aux_lcoeffs, const float_type* __restrict__ aux_ucoeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int input_size_x, const int input_size_y, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y, const int new_length_x, const int new_length_y, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    int out_z = threadIdx.x;

    __shared__ int n;

    __shared__ int min_out_x;
    __shared__ int min_out_y;

    __shared__ int min_x;
    __shared__ int min_y;

    __shared__ int max_x;
    __shared__ int max_y;

    if(out_z == 0)
    {
        n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;

        min_out_x = offset_x + last_x*shift_x;
        min_out_y = offset_y + last_y*shift_y;

        min_x = (min_out_x < 0) ? -min_out_x : 0;
        min_y = (min_out_y < 0) ? -min_out_y : 0;

        max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
        max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;
    }

    __syncthreads();

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(out_z < output_size_z)
    {
        for(int out_x = min_x; out_x < max_x; out_x++)
        {
            for(int out_y = min_y; out_y < max_y; out_y++)
            {
                for(size_t x_shift = 0; x_shift < filter_size_x; x_shift++)
                {
                    for(size_t y_shift = 0; y_shift < filter_size_y; y_shift++)
                    {
                        const int x_val = out_x*stride_x + min_out_x*stride_x + x_shift - pad_x;
                        const int y_val = out_y*stride_y + min_out_y*stride_y + y_shift - pad_y;

                        if(!((x_val < 0) || (x_val >= input_size_x)))
                        {
                            if(!((y_val < 0) || (y_val >= input_size_y)))
                            {
                                const int mat_in = (out_x*stride_x + x_shift)*new_length_y*output_size_z + (out_y*stride_y + y_shift)*output_size_z + out_z;
                                const int b = n*new_length_x*new_length_y*output_size_z + mat_in;

                                const int mat_out = out_x*length_y*output_size_z + out_y*output_size_z + out_z;
                                const int a = n*length_x*length_y*output_size_z + mat_out;

                                const int c = (out_x + min_out_x)*output_size_y*filter_size_x*filter_size_y*output_size_z + (out_y + min_out_y)*filter_size_x*filter_size_y*output_size_z + (x_shift*filter_size_y + y_shift)*output_size_z + out_z;

                                const float_type prev_inf_coeff = expr_inf_coeff[a];
                                const float_type prev_sup_coeff = expr_sup_coeff[a];

                                if(prev_inf_coeff > 0)
                                {
                                    affine_trans_coeff(res_inf_coeff[b], res_sup_coeff[b], tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, aux_lcoeffs[c]);
                                }

                                if(prev_sup_coeff < 0)
                                {
                                    affine_trans_coeff(res_inf_coeff[b], res_sup_coeff[b], tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, aux_ucoeffs[c]);
                                }
                            }
                        }
                    }
                }
            }
        }

        out_z += blockDim.x;
    }
}


__global__
void ucoeffs_from_previous_layer_pool_sparse(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, float_type* __restrict__ res_inf_coeff, float_type* __restrict__ res_sup_coeff, const float_type* __restrict__ aux_lcoeffs, const float_type* __restrict__ aux_ucoeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int input_size_x, const int input_size_y, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y, const int new_length_x, const int new_length_y, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    int out_z = threadIdx.x;

    __shared__ int n;

    __shared__ int min_out_x;
    __shared__ int min_out_y;

    __shared__ int min_x;
    __shared__ int min_y;

    __shared__ int max_x;
    __shared__ int max_y;

    if(out_z == 0)
    {
        n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;

        min_out_x = offset_x + last_x*shift_x;
        min_out_y = offset_y + last_y*shift_y;

        min_x = (min_out_x < 0) ? -min_out_x : 0;
        min_y = (min_out_y < 0) ? -min_out_y : 0;

        max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
        max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;
    }

    __syncthreads();

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(out_z < output_size_z)
    {
        for(int out_x = min_x; out_x < max_x; out_x++)
        {
            for(int out_y = min_y; out_y < max_y; out_y++)
            {
                for(size_t x_shift = 0; x_shift < filter_size_x; x_shift++)
                {
                    for(size_t y_shift = 0; y_shift < filter_size_y; y_shift++)
                    {
                        const int x_val = out_x*stride_x + min_out_x*stride_x + x_shift - pad_x;
                        const int y_val = out_y*stride_y + min_out_y*stride_y + y_shift - pad_y;

                        if(!((x_val < 0) || (x_val >= input_size_x)))
                        {
                            if(!((y_val < 0) || (y_val >= input_size_y)))
                            {
                                const int mat_in = (out_x*stride_x + x_shift)*new_length_y*output_size_z + (out_y*stride_y + y_shift)*output_size_z + out_z;
                                const int b = n*new_length_x*new_length_y*output_size_z + mat_in;

                                const int mat_out = out_x*length_y*output_size_z + out_y*output_size_z + out_z;
                                const int a = n*length_x*length_y*output_size_z + mat_out;

                                const int c = (out_x + min_out_x)*output_size_y*filter_size_x*filter_size_y*output_size_z + (out_y + min_out_y)*filter_size_x*filter_size_y*output_size_z + (x_shift*filter_size_y + y_shift)*output_size_z + out_z;

                                const float_type prev_inf_coeff = expr_inf_coeff[a];
                                const float_type prev_sup_coeff = expr_sup_coeff[a];

                                if(prev_inf_coeff > 0)
                                {
                                    affine_trans_coeff(res_inf_coeff[b], res_sup_coeff[b], tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, aux_ucoeffs[c]);
                                }

                                if(prev_sup_coeff < 0)
                                {
                                    affine_trans_coeff(res_inf_coeff[b], res_sup_coeff[b], tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, aux_lcoeffs[c]);
                                }
                            }
                        }
                    }
                }
            }
        }

        out_z += blockDim.x;
    }
}


__global__
void lcoeffs_from_input_poly_sparse(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, float_type* __restrict__ res_inf_coeff, float_type* __restrict__ res_sup_coeff, const float_type* __restrict__ input_lcoeffs, const float_type* __restrict__ input_ucoeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int num_in_neurons, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    int i = threadIdx.x;

    __shared__ int n;

    __shared__ int min_out_x;
    __shared__ int min_out_y;

    __shared__ int min_x;
    __shared__ int min_y;

    __shared__ int max_x;
    __shared__ int max_y;

    if(i == 0)
    {
        n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;

        min_out_x = offset_x + last_x*shift_x;
        min_out_y = offset_y + last_y*shift_y;

        min_x = (min_out_x < 0) ? -min_out_x : 0;
        min_y = (min_out_y < 0) ? -min_out_y : 0;

        max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
        max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;
    }

    __syncthreads();

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(i < num_in_neurons)
    {
        const int b = n*num_in_neurons + i;

        float_type inf_coeff = res_inf_coeff[b];
        float_type sup_coeff = res_sup_coeff[b];

        for(int out_x = min_x; out_x < max_x; out_x++)
        {
            for(int out_y = min_y; out_y < max_y; out_y++)
            {
                for(int out_z = 0; out_z < output_size_z; out_z++)
                {
                    const int mat_out = out_x*length_y*output_size_z + out_y*output_size_z + out_z;
                    const int a = n*length_x*length_y*output_size_z + mat_out;

                    const float_type prev_inf_coeff = expr_inf_coeff[a];
                    const float_type prev_sup_coeff = expr_sup_coeff[a];

                    const int c = ((out_x + min_out_x)*output_size_y*output_size_z + (out_y + min_out_y)*output_size_z + out_z)*num_in_neurons + i;

                    if(prev_inf_coeff > 0)
                    {
                        affine_trans_coeff(inf_coeff, sup_coeff, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_lcoeffs[c]);
                    }

                    if(prev_sup_coeff < 0)
                    {
                        affine_trans_coeff(inf_coeff, sup_coeff, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_ucoeffs[c]);
                    }
                }
            }
        }

        res_inf_coeff[b] = inf_coeff;
        res_sup_coeff[b] = sup_coeff;

        i += blockDim.x;
    }
}


__global__
void ucoeffs_from_input_poly_sparse(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, float_type* __restrict__ res_inf_coeff, float_type* __restrict__ res_sup_coeff, const float_type* __restrict__ input_lcoeffs, const float_type* __restrict__ input_ucoeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int num_in_neurons, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    int i = threadIdx.x;

    __shared__ int n;

    __shared__ int min_out_x;
    __shared__ int min_out_y;

    __shared__ int min_x;
    __shared__ int min_y;

    __shared__ int max_x;
    __shared__ int max_y;

    if(i == 0)
    {
        n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;

        min_out_x = offset_x + last_x*shift_x;
        min_out_y = offset_y + last_y*shift_y;

        min_x = (min_out_x < 0) ? -min_out_x : 0;
        min_y = (min_out_y < 0) ? -min_out_y : 0;

        max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
        max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;
    }

    __syncthreads();

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(i < num_in_neurons)
    {
        const int b = n*num_in_neurons + i;

        float_type inf_coeff = res_inf_coeff[b];
        float_type sup_coeff = res_sup_coeff[b];

        for(int out_x = min_x; out_x < max_x; out_x++)
        {
            for(int out_y = min_y; out_y < max_y; out_y++)
            {
                for(int out_z = 0; out_z < output_size_z; out_z++)
                {
                    const int mat_out = out_x*length_y*output_size_z + out_y*output_size_z + out_z;
                    const int a = n*length_x*length_y*output_size_z + mat_out;

                    const float_type prev_inf_coeff = expr_inf_coeff[a];
                    const float_type prev_sup_coeff = expr_sup_coeff[a];

                    const int c = ((out_x + min_out_x)*output_size_y*output_size_z + (out_y + min_out_y)*output_size_z + out_z)*num_in_neurons + i;

                    if(prev_inf_coeff > 0)
                    {
                        affine_trans_coeff(inf_coeff, sup_coeff, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_ucoeffs[c]);
                    }

                    if(prev_sup_coeff < 0)
                    {
                        affine_trans_coeff(inf_coeff, sup_coeff, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_lcoeffs[c]);
                    }
                }
            }
        }

        res_inf_coeff[b] = inf_coeff;
        res_sup_coeff[b] = sup_coeff;

        i += blockDim.x;
    }
}


__global__
void csts_from_previous_layer(const float_type* __restrict__ expr_linf_coeff, const float_type* __restrict__ expr_lsup_coeff, const float_type* __restrict__ expr_uinf_coeff, const float_type* __restrict__ expr_usup_coeff, const float_type* __restrict__ expr_linf_cst, const float_type* __restrict__ expr_lsup_cst, const float_type* __restrict__ expr_uinf_cst, const float_type* __restrict__ expr_usup_cst, float_type* __restrict__ res_linf_cst, float_type* __restrict__ res_lsup_cst, float_type* __restrict__ res_uinf_cst, float_type* __restrict__ res_usup_cst, const float_type* __restrict__ aux_csts, const int num_out_neurons_current_layer)
{
    const int n = blockIdx.x;
    int i = threadIdx.x;

    float_type linf_cst = 0;
    float_type lsup_cst = 0;
    float_type uinf_cst = 0;
    float_type usup_cst = 0;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(i < num_out_neurons_current_layer)
    {
        const int a = n*num_out_neurons_current_layer + i;
        const float_type aux_cst = aux_csts[i];

        affine_trans_cst(linf_cst, lsup_cst, tmp1, tmp2, maxRes, maxMul, expr_linf_coeff[a], expr_lsup_coeff[a], aux_cst);
        affine_trans_cst(uinf_cst, usup_cst, tmp1, tmp2, maxRes, maxMul, expr_uinf_coeff[a], expr_usup_coeff[a], aux_cst);

        i += blockDim.x;
    }

    block_reduce_sum_csts(linf_cst, lsup_cst, blockDim.x);
    block_reduce_sum_csts(uinf_cst, usup_cst, blockDim.x);

    if(threadIdx.x == 0)
    {
        res_linf_cst[n] = add_rd(linf_cst, expr_linf_cst[n]);
        res_lsup_cst[n] = add_ru(lsup_cst, expr_lsup_cst[n]);
        res_uinf_cst[n] = add_rd(uinf_cst, expr_uinf_cst[n]);
        res_usup_cst[n] = add_ru(usup_cst, expr_usup_cst[n]);
    }
}


__global__
void csts_from_previous_layer_conv(const float_type* __restrict__ expr_linf_coeff, const float_type* __restrict__ expr_lsup_coeff, const float_type* __restrict__ expr_uinf_coeff, const float_type* __restrict__ expr_usup_coeff, const float_type* __restrict__ expr_linf_cst, const float_type* __restrict__ expr_lsup_cst, const float_type* __restrict__ expr_uinf_cst, const float_type* __restrict__ expr_usup_cst, float_type* __restrict__ res_linf_cst, float_type* __restrict__ res_lsup_cst, float_type* __restrict__ res_uinf_cst, float_type* __restrict__ res_usup_cst, const float_type* __restrict__ aux_csts, const int current_layer_out_size_x, const int current_layer_out_size_y, const int current_layer_out_size_z)
{
    const int n = blockIdx.x;
    int j = threadIdx.x;

    float_type linf_cst = 0;
    float_type lsup_cst = 0;
    float_type uinf_cst = 0;
    float_type usup_cst = 0;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(j < current_layer_out_size_z)
    {
        for(int i = 0; i < current_layer_out_size_x*current_layer_out_size_y; i++)
        {
            const int a = n*current_layer_out_size_x*current_layer_out_size_y*current_layer_out_size_z + i*current_layer_out_size_z + j;
            const float_type aux_cst = aux_csts[j];

            affine_trans_cst(linf_cst, lsup_cst, tmp1, tmp2, maxRes, maxMul, expr_linf_coeff[a], expr_lsup_coeff[a], aux_cst);
            affine_trans_cst(uinf_cst, usup_cst, tmp1, tmp2, maxRes, maxMul, expr_uinf_coeff[a], expr_usup_coeff[a], aux_cst);
        }

        j += blockDim.x;
    }

    block_reduce_sum_csts(linf_cst, lsup_cst, blockDim.x);
    block_reduce_sum_csts(uinf_cst, usup_cst, blockDim.x);

    if(threadIdx.x == 0)
    {
        res_linf_cst[n] = add_rd(linf_cst, expr_linf_cst[n]);
        res_lsup_cst[n] = add_ru(lsup_cst, expr_lsup_cst[n]);
        res_uinf_cst[n] = add_rd(uinf_cst, expr_uinf_cst[n]);
        res_usup_cst[n] = add_ru(usup_cst, expr_usup_cst[n]);
    }
}


__global__
void csts_from_previous_layer_conv_sparse(const float_type* __restrict__ expr_linf_coeff, const float_type* __restrict__ expr_lsup_coeff, const float_type* __restrict__ expr_uinf_coeff, const float_type* __restrict__ expr_usup_coeff, const float_type* __restrict__ expr_linf_cst, const float_type* __restrict__ expr_lsup_cst, const float_type* __restrict__ expr_uinf_cst, const float_type* __restrict__ expr_usup_cst, float_type* __restrict__ res_linf_cst, float_type* __restrict__ res_lsup_cst, float_type* __restrict__ res_uinf_cst, float_type* __restrict__ res_usup_cst, const float_type* __restrict__ aux_csts, const int output_size_x, const int output_size_y, const int output_size_z, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    int out_z = threadIdx.x;

    __shared__ int n;

    __shared__ int min_out_x;
    __shared__ int min_out_y;

    __shared__ int min_x;
    __shared__ int min_y;

    __shared__ int max_x;
    __shared__ int max_y;

    if(out_z == 0)
    {
        n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;

        min_out_x = offset_x + last_x*shift_x;
        min_out_y = offset_y + last_y*shift_y;

        min_x = (min_out_x < 0) ? -min_out_x : 0;
        min_y = (min_out_y < 0) ? -min_out_y : 0;

        max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
        max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;
    }

    __syncthreads();

    float_type linf_cst = 0;
    float_type lsup_cst = 0;
    float_type uinf_cst = 0;
    float_type usup_cst = 0;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(out_z < output_size_z)
    {
        for(int out_x = min_x; out_x < max_x; out_x++)
        {
            for(int out_y = min_y; out_y < max_y; out_y++)
            {
                const int mat_out = n*length_x*length_y*output_size_z + out_x*length_y*output_size_z + out_y*output_size_z + out_z;
                const float_type aux_cst = aux_csts[out_z];

                affine_trans_cst(linf_cst, lsup_cst, tmp1, tmp2, maxRes, maxMul, expr_linf_coeff[mat_out], expr_lsup_coeff[mat_out], aux_cst);
                affine_trans_cst(uinf_cst, usup_cst, tmp1, tmp2, maxRes, maxMul, expr_uinf_coeff[mat_out], expr_usup_coeff[mat_out], aux_cst);
            }
        }

        out_z += blockDim.x;
    }

    block_reduce_sum_csts(linf_cst, lsup_cst, blockDim.x);
    block_reduce_sum_csts(uinf_cst, usup_cst, blockDim.x);

    if(threadIdx.x == 0)
    {
        res_linf_cst[n] = add_rd(linf_cst, expr_linf_cst[n]);
        res_lsup_cst[n] = add_ru(lsup_cst, expr_lsup_cst[n]);
        res_uinf_cst[n] = add_rd(uinf_cst, expr_uinf_cst[n]);
        res_usup_cst[n] = add_ru(usup_cst, expr_usup_cst[n]);
    }
}


__global__
void lcsts_from_input_poly_sparse(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, const float_type* __restrict__ expr_inf_cst, const float_type* __restrict__ expr_sup_cst, float_type* __restrict__ res_inf_cst, float_type* __restrict__ res_sup_cst, const float_type* __restrict__ input_lcsts, const float_type* __restrict__ input_ucsts, const float_type* __restrict__ input_inf, const float_type* __restrict__ input_sup, const int output_size_x, const int output_size_y, const int output_size_z, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    const int n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;

    const int min_out_x = offset_x + last_x*shift_x;
    const int min_out_y = offset_y + last_y*shift_y;

    const int min_x = (min_out_x < 0) ? -min_out_x : 0;
    const int min_y = (min_out_y < 0) ? -min_out_y : 0;

    const int max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
    const int max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;

    float_type inf_cst = expr_inf_cst[n];
    float_type sup_cst = expr_sup_cst[n];

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    for(int out_x = min_x; out_x < max_x; out_x++)
    {
        for(int out_y = min_y; out_y < max_y; out_y++)
        {
            for(int out_z = 0; out_z < output_size_z; out_z++)
            {
                const int j = out_x*length_y*output_size_z + out_y*output_size_z + out_z;
                const int i = (out_x + min_out_x)*output_size_y*output_size_z + (out_y + min_out_y)*output_size_z + out_z;

                const int mat_out = n*length_x*length_y*output_size_z + j;

                const float_type prev_inf_coeff = expr_inf_coeff[mat_out];
                const float_type prev_sup_coeff = expr_sup_coeff[mat_out];

                if(prev_inf_coeff > 0)
                {
                    affine_trans_cst(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_lcsts[i]);
                }
                else if(prev_sup_coeff < 0)
                {
                    affine_trans_cst(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_ucsts[i]);
                }
                else
                {
                    lcst_input_poly_neutral(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_inf[i], input_sup[i]);
                }
            }
        }
    }

    res_inf_cst[n] = inf_cst;
    res_sup_cst[n] = sup_cst;
}


__global__
void ucsts_from_input_poly_sparse(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, const float_type* __restrict__ expr_inf_cst, const float_type* __restrict__ expr_sup_cst, float_type* __restrict__ res_inf_cst, float_type* __restrict__ res_sup_cst, const float_type* __restrict__ input_lcsts, const float_type* __restrict__ input_ucsts, const float_type* __restrict__ input_inf, const float_type* __restrict__ input_sup, const int output_size_x, const int output_size_y, const int output_size_z, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    const int n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;

    const int min_out_x = offset_x + last_x*shift_x;
    const int min_out_y = offset_y + last_y*shift_y;

    const int min_x = (min_out_x < 0) ? -min_out_x : 0;
    const int min_y = (min_out_y < 0) ? -min_out_y : 0;

    const int max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
    const int max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;

    float_type inf_cst = expr_inf_cst[n];
    float_type sup_cst = expr_sup_cst[n];

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    for(int out_x = min_x; out_x < max_x; out_x++)
    {
        for(int out_y = min_y; out_y < max_y; out_y++)
        {
            for(int out_z = 0; out_z < output_size_z; out_z++)
            {
                const int j = out_x*length_y*output_size_z + out_y*output_size_z + out_z;
                const int i = (out_x + min_out_x)*output_size_y*output_size_z + (out_y + min_out_y)*output_size_z + out_z;

                const int mat_out = n*length_x*length_y*output_size_z + j;

                const float_type prev_inf_coeff = expr_inf_coeff[mat_out];
                const float_type prev_sup_coeff = expr_sup_coeff[mat_out];

                if(prev_inf_coeff > 0)
                {
                    affine_trans_cst(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_ucsts[i]);
                }
                else if(prev_sup_coeff < 0)
                {
                    affine_trans_cst(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_lcsts[i]);
                }
                else
                {
                    ucst_input_poly_neutral(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, input_inf[i], input_sup[i]);
                }
            }
        }
    }

    res_inf_cst[n] = inf_cst;
    res_sup_cst[n] = sup_cst;
}


__global__
void device_layer_create_sparse_exprs(float_type* dense_coeff, float_type* bias, const float_type* filter_weights, const float_type* filter_bias, const int chunk_counter, const int input_size_x, const int input_size_y, const int input_size_z, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int out_x = blockIdx.x;
    const int out_y = blockIdx.y;
    const int out_z = blockIdx.z;

    const int local_mat_x = out_x*gridDim.y*gridDim.z + out_y*gridDim.z + out_z;

    for(int x_shift = 0; x_shift < filter_size_x; x_shift++)
    {
        for(int y_shift = 0; y_shift < filter_size_y; y_shift++)
        {
            for(int inp_z = 0; inp_z < input_size_z; inp_z++)
            {
                const int x_val = out_x*stride_x + x_shift - pad_x;
                const int y_val = out_y*stride_y + y_shift - pad_y;

                if((y_val < 0) || (y_val >= input_size_y))
                {
                    continue;
                }

                if((x_val < 0) || (x_val >= input_size_x))
                {
                    continue;
                }

                const int mat_y = x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z;
                const int filter_index = (chunk_counter*gridDim.z + out_z)*filter_size_x*filter_size_y*input_size_z + x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z;

                dense_coeff[local_mat_x*filter_size_x*filter_size_y*input_size_z + mat_y] = filter_weights[filter_index];
            }
        }
    }

    bias[local_mat_x] = filter_bias[chunk_counter*gridDim.z + out_z];
}


void update_state_using_predecessor_layer(fppoly_t* fp, float_type** linf_coeff, float_type** lsup_coeff, float_type** linf_cst, float_type** lsup_cst, float_type** uinf_coeff, float_type** usup_coeff, float_type** uinf_cst, float_type** usup_cst, float_type** linf_coeff_tmp, float_type** lsup_coeff_tmp, float_type** linf_cst_tmp, float_type** lsup_cst_tmp, float_type** uinf_coeff_tmp, float_type** usup_coeff_tmp, float_type** uinf_cst_tmp, float_type** usup_cst_tmp, const size_t k, const bool use_area_heuristic, const size_t num_out_neurons_last_layer)
{
    const size_t num_out_neurons_current_layer = fp->layers[k]->num_out_neurons;
    const size_t num_in_neurons_current_layer  = fp->layers[k]->num_in_neurons;

    if(fp->layers[k]->activation == RELU)
    {
        expr_replace_relu_bounds<<<num_out_neurons_last_layer, num_threads>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_cst, *lsup_cst, *uinf_cst, *usup_cst, fp->layers[k]->lb_array, fp->layers[k]->ub_array, num_out_neurons_current_layer, use_area_heuristic);
    }

    *linf_coeff_tmp = malloc_device<float_type>(num_out_neurons_last_layer*num_in_neurons_current_layer);
    *lsup_coeff_tmp = malloc_device<float_type>(num_out_neurons_last_layer*num_in_neurons_current_layer);
    *uinf_coeff_tmp = malloc_device<float_type>(num_out_neurons_last_layer*num_in_neurons_current_layer);
    *usup_coeff_tmp = malloc_device<float_type>(num_out_neurons_last_layer*num_in_neurons_current_layer);

    cudaMemset(*linf_coeff_tmp, 0, num_out_neurons_last_layer*num_in_neurons_current_layer*sizeof(float_type));
    cudaMemset(*lsup_coeff_tmp, 0, num_out_neurons_last_layer*num_in_neurons_current_layer*sizeof(float_type));
    cudaMemset(*uinf_coeff_tmp, 0, num_out_neurons_last_layer*num_in_neurons_current_layer*sizeof(float_type));
    cudaMemset(*usup_coeff_tmp, 0, num_out_neurons_last_layer*num_in_neurons_current_layer*sizeof(float_type));

    if(fp->layers[k]->type == CONV)
    {
        if(fp->layers[k]->input_size[0]*fp->layers[k]->input_size[1]*fp->layers[k]->input_size[2] > num_threads)
        {
            if(fp->layers[k]->input_size[1]*fp->layers[k]->input_size[2] > num_threads)
            {
                if(fp->layers[k]->input_size[2] > num_threads)
                {
                    coeffs_from_previous_layer_conv_filter_serial<<<num_out_neurons_last_layer, num_threads>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, *uinf_coeff_tmp, *usup_coeff_tmp, fp->layers[k]->filter_weights, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
                }
                else
                {
                    coeffs_from_previous_layer_conv_filter_serial<<<num_out_neurons_last_layer, fp->layers[k]->input_size[2]>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, *uinf_coeff_tmp, *usup_coeff_tmp, fp->layers[k]->filter_weights, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
                }
            }
            else
            {
                coeffs_from_previous_layer_conv_x_filter_serial<<<num_out_neurons_last_layer, dim3(fp->layers[k]->input_size[2], fp->layers[k]->filter_size[1], 1)>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, *uinf_coeff_tmp, *usup_coeff_tmp, fp->layers[k]->filter_weights, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
            }
        }
        else
        {
            coeffs_from_previous_layer_conv<<<num_out_neurons_last_layer, dim3(fp->layers[k]->input_size[2], fp->layers[k]->filter_size[1], fp->layers[k]->filter_size[0])>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, *uinf_coeff_tmp, *usup_coeff_tmp, fp->layers[k]->filter_weights, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
        }

        csts_from_previous_layer_conv<<<num_out_neurons_last_layer, fp->layers[k]->output_size[2]>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_cst, *lsup_cst, *uinf_cst, *usup_cst, *linf_cst_tmp, *lsup_cst_tmp, *uinf_cst_tmp, *usup_cst_tmp, fp->layers[k]->filter_bias, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2]);
    }
    else
    {
        if(fp->layers[k]->type == FFN)
        {
            coeffs_from_previous_layer<<<num_out_neurons_last_layer, num_threads>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, *uinf_coeff_tmp, *usup_coeff_tmp, fp->layers[k]->coeffs, num_out_neurons_current_layer, num_in_neurons_current_layer);

            csts_from_previous_layer<<<num_out_neurons_last_layer, num_threads>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_cst, *lsup_cst, *uinf_cst, *usup_cst, *linf_cst_tmp, *lsup_cst_tmp, *uinf_cst_tmp, *usup_cst_tmp, fp->layers[k]->csts, num_out_neurons_current_layer);
        }
        else
        {
            lcoeffs_from_previous_layer_pool<<<num_out_neurons_last_layer, fp->layers[k]->output_size[2]>>>(*linf_coeff, *lsup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, fp->layers[k]->pool_lcoeffs, fp->layers[k]->pool_ucoeffs, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
            ucoeffs_from_previous_layer_pool<<<num_out_neurons_last_layer, fp->layers[k]->output_size[2]>>>(*uinf_coeff, *usup_coeff, *uinf_coeff_tmp, *usup_coeff_tmp, fp->layers[k]->pool_lcoeffs, fp->layers[k]->pool_ucoeffs, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);

            lcsts_from_input_poly<<<num_out_neurons_last_layer, 1>>>(*linf_coeff, *lsup_coeff, *linf_cst, *lsup_cst, *linf_cst_tmp, *lsup_cst_tmp, fp->layers[k]->pool_lcsts, fp->layers[k]->pool_ucsts, fp->layers[k]->lb_array, fp->layers[k]->ub_array, num_out_neurons_current_layer);
            ucsts_from_input_poly<<<num_out_neurons_last_layer, 1>>>(*uinf_coeff, *usup_coeff, *uinf_cst, *usup_cst, *uinf_cst_tmp, *usup_cst_tmp, fp->layers[k]->pool_lcsts, fp->layers[k]->pool_ucsts, fp->layers[k]->lb_array, fp->layers[k]->ub_array, num_out_neurons_current_layer);
        }
    }

    std::swap(*linf_coeff, *linf_coeff_tmp);
    std::swap(*lsup_coeff, *lsup_coeff_tmp);
    std::swap(*linf_cst, *linf_cst_tmp);
    std::swap(*lsup_cst, *lsup_cst_tmp);

    std::swap(*uinf_coeff, *uinf_coeff_tmp);
    std::swap(*usup_coeff, *usup_coeff_tmp);
    std::swap(*uinf_cst, *uinf_cst_tmp);
    std::swap(*usup_cst, *usup_cst_tmp);

    cudaFree(*linf_coeff_tmp);
    cudaFree(*lsup_coeff_tmp);
    cudaFree(*uinf_coeff_tmp);
    cudaFree(*usup_coeff_tmp);

    *linf_coeff_tmp = nullptr;
    *lsup_coeff_tmp = nullptr;
    *uinf_coeff_tmp = nullptr;
    *usup_coeff_tmp = nullptr;
}


void update_state_using_previous_layers(elina_manager_t* man, fppoly_t* fp, const size_t layerno, const bool use_area_heuristic)
{
    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    const size_t num_out_neurons_last_layer = fp->layers[layerno]->num_out_neurons;
    const size_t num_in_neurons_last_layer = fp->layers[layerno]->num_in_neurons;

    const size_t num_in_neurons_0_layer = fp->layers[0]->num_in_neurons;

    float_type* coeffs = fp->layers[layerno]->coeffs;
    float_type* csts = fp->layers[layerno]->csts;

    float_type* lb_array = fp->layers[layerno]->lb_array;
    float_type* ub_array = fp->layers[layerno]->ub_array;

    float_type* linf_coeff = malloc_device<float_type>(num_out_neurons_last_layer*num_in_neurons_last_layer);
    float_type* lsup_coeff = malloc_device<float_type>(num_out_neurons_last_layer*num_in_neurons_last_layer);
    float_type* uinf_coeff = malloc_device<float_type>(num_out_neurons_last_layer*num_in_neurons_last_layer);
    float_type* usup_coeff = malloc_device<float_type>(num_out_neurons_last_layer*num_in_neurons_last_layer);

    float_type* linf_cst = malloc_device<float_type>(num_out_neurons_last_layer);
    float_type* lsup_cst = malloc_device<float_type>(num_out_neurons_last_layer);
    float_type* uinf_cst = malloc_device<float_type>(num_out_neurons_last_layer);
    float_type* usup_cst = malloc_device<float_type>(num_out_neurons_last_layer);

    copy_coeffs_and_csts<<<num_out_neurons_last_layer, 1>>>(linf_coeff, linf_cst, coeffs, csts, num_in_neurons_last_layer);
    copy_coeffs_and_csts<<<num_out_neurons_last_layer, 1>>>(lsup_coeff, lsup_cst, coeffs, csts, num_in_neurons_last_layer);
    copy_coeffs_and_csts<<<num_out_neurons_last_layer, 1>>>(uinf_coeff, uinf_cst, coeffs, csts, num_in_neurons_last_layer);
    copy_coeffs_and_csts<<<num_out_neurons_last_layer, 1>>>(usup_coeff, usup_cst, coeffs, csts, num_in_neurons_last_layer);

    float_type* linf_coeff_tmp;
    float_type* lsup_coeff_tmp;
    float_type* uinf_coeff_tmp;
    float_type* usup_coeff_tmp;

    float_type* linf_cst_tmp = malloc_device<float_type>(num_out_neurons_last_layer);
    float_type* lsup_cst_tmp = malloc_device<float_type>(num_out_neurons_last_layer);
    float_type* uinf_cst_tmp = malloc_device<float_type>(num_out_neurons_last_layer);
    float_type* usup_cst_tmp = malloc_device<float_type>(num_out_neurons_last_layer);

    int k = fp->layers[layerno]->predecessors[0] - 1;

    while(k >= 0)
    {
        if(fp->layers[k]->type == RESIDUAL)
        {
            if(fp->layers[k]->activation==RELU)
            {
                expr_replace_relu_bounds<<<num_out_neurons_last_layer, num_threads>>>(linf_coeff, lsup_coeff, uinf_coeff, usup_coeff, linf_cst, lsup_cst, uinf_cst, usup_cst, fp->layers[k]->lb_array, fp->layers[k]->ub_array, fp->layers[k]->num_out_neurons, use_area_heuristic);
            }

            const size_t reslayer_size = fp->layers[k]->num_out_neurons;

            const size_t predecessor1 = fp->layers[k]->predecessors[0] - 1;
            const size_t predecessor2 = fp->layers[k]->predecessors[1] - 1;

            char* predecessor_map = (char*) calloc(k, sizeof(char));
            int iter = predecessor1;

            while(iter >= 0)
            {
                predecessor_map[iter] = 1;
                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            iter =  predecessor2;
            int common_predecessor = 0;

            while(iter >= 0)
            {
                if(predecessor_map[iter] == 1)
                {
                    common_predecessor = iter;

                    break;
                }

                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            free(predecessor_map);

            float_type* copy_linf_coeff = malloc_device<float_type>(num_out_neurons_last_layer*reslayer_size);
            float_type* copy_lsup_coeff = malloc_device<float_type>(num_out_neurons_last_layer*reslayer_size);
            float_type* copy_uinf_coeff = malloc_device<float_type>(num_out_neurons_last_layer*reslayer_size);
            float_type* copy_usup_coeff = malloc_device<float_type>(num_out_neurons_last_layer*reslayer_size);

            float_type* copy_linf_cst = malloc_device<float_type>(num_out_neurons_last_layer);
            float_type* copy_lsup_cst = malloc_device<float_type>(num_out_neurons_last_layer);
            float_type* copy_uinf_cst = malloc_device<float_type>(num_out_neurons_last_layer);
            float_type* copy_usup_cst = malloc_device<float_type>(num_out_neurons_last_layer);

            copy_coeffs_and_csts<<<num_out_neurons_last_layer, 1>>>(copy_linf_coeff, copy_linf_cst, linf_coeff, linf_cst, reslayer_size);
            copy_coeffs_and_csts<<<num_out_neurons_last_layer, 1>>>(copy_lsup_coeff, copy_lsup_cst, lsup_coeff, lsup_cst, reslayer_size);
            copy_coeffs_and_csts<<<num_out_neurons_last_layer, 1>>>(copy_uinf_coeff, copy_uinf_cst, uinf_coeff, uinf_cst, reslayer_size);
            copy_coeffs_and_csts<<<num_out_neurons_last_layer, 1>>>(copy_usup_coeff, copy_usup_cst, usup_coeff, usup_cst, reslayer_size);

            cudaMemset(copy_linf_cst, 0, num_out_neurons_last_layer*sizeof(float_type));
            cudaMemset(copy_lsup_cst, 0, num_out_neurons_last_layer*sizeof(float_type));
            cudaMemset(copy_uinf_cst, 0, num_out_neurons_last_layer*sizeof(float_type));
            cudaMemset(copy_usup_cst, 0, num_out_neurons_last_layer*sizeof(float_type));

            float_type* copy_linf_coeff_tmp;
            float_type* copy_lsup_coeff_tmp;
            float_type* copy_uinf_coeff_tmp;
            float_type* copy_usup_coeff_tmp;

            float_type* copy_linf_cst_tmp = malloc_device<float_type>(num_out_neurons_last_layer);
            float_type* copy_lsup_cst_tmp = malloc_device<float_type>(num_out_neurons_last_layer);
            float_type* copy_uinf_cst_tmp = malloc_device<float_type>(num_out_neurons_last_layer);
            float_type* copy_usup_cst_tmp = malloc_device<float_type>(num_out_neurons_last_layer);

            iter = predecessor1;

            while(iter != common_predecessor)
            {
                update_state_using_predecessor_layer(fp, &linf_coeff, &lsup_coeff, &linf_cst, &lsup_cst, &uinf_coeff, &usup_coeff, &uinf_cst, &usup_cst, &linf_coeff_tmp, &lsup_coeff_tmp, &linf_cst_tmp, &lsup_cst_tmp, &uinf_coeff_tmp, &usup_coeff_tmp, &uinf_cst_tmp, &usup_cst_tmp, iter, use_area_heuristic, num_out_neurons_last_layer);

                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            iter = predecessor2;

            while(iter != common_predecessor)
            {
                update_state_using_predecessor_layer(fp, &copy_linf_coeff, &copy_lsup_coeff, &copy_linf_cst, &copy_lsup_cst, &copy_uinf_coeff, &copy_usup_coeff, &copy_uinf_cst, &copy_usup_cst, &copy_linf_coeff_tmp, &copy_lsup_coeff_tmp, &copy_linf_cst_tmp, &copy_lsup_cst_tmp, &copy_uinf_coeff_tmp, &copy_usup_coeff_tmp, &copy_uinf_cst_tmp, &copy_usup_cst_tmp, iter, use_area_heuristic, num_out_neurons_last_layer);

                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            add_coeffs_and_csts<<<num_out_neurons_last_layer, 1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, copy_linf_coeff, copy_lsup_coeff, copy_linf_cst, copy_lsup_cst, fp->layers[common_predecessor]->num_out_neurons);
            add_coeffs_and_csts<<<num_out_neurons_last_layer, 1>>>(uinf_coeff, usup_coeff, uinf_cst, usup_cst, copy_uinf_coeff, copy_usup_coeff, copy_uinf_cst, copy_usup_cst, fp->layers[common_predecessor]->num_out_neurons);

            cudaFree(copy_linf_coeff);
            cudaFree(copy_lsup_coeff);
            cudaFree(copy_linf_cst);
            cudaFree(copy_lsup_cst);

            cudaFree(copy_uinf_coeff);
            cudaFree(copy_usup_coeff);
            cudaFree(copy_uinf_cst);
            cudaFree(copy_usup_cst);

            cudaFree(copy_linf_cst_tmp);
            cudaFree(copy_lsup_cst_tmp);

            cudaFree(copy_uinf_cst_tmp);
            cudaFree(copy_usup_cst_tmp);

            copy_linf_coeff = nullptr;
            copy_lsup_coeff = nullptr;
            copy_linf_cst = nullptr;
            copy_lsup_cst = nullptr;

            copy_uinf_coeff = nullptr;
            copy_usup_coeff = nullptr;
            copy_uinf_cst = nullptr;
            copy_usup_cst = nullptr;

            copy_linf_cst_tmp = nullptr;
            copy_lsup_cst_tmp = nullptr;

            copy_uinf_cst_tmp = nullptr;
            copy_usup_cst_tmp = nullptr;

            k = common_predecessor;
        }
        else
        {
            update_state_using_predecessor_layer(fp, &linf_coeff, &lsup_coeff, &linf_cst, &lsup_cst, &uinf_coeff, &usup_coeff, &uinf_cst, &usup_cst, &linf_coeff_tmp, &lsup_coeff_tmp, &linf_cst_tmp, &lsup_cst_tmp, &uinf_coeff_tmp, &usup_coeff_tmp, &uinf_cst_tmp, &usup_cst_tmp, k, use_area_heuristic, num_out_neurons_last_layer);

            k = fp->layers[k]->predecessors[0] - 1;
        }
    }

    if((fp->input_lweights != nullptr) && (fp->input_uweights != nullptr) && (fp->input_lcst != nullptr) && (fp->input_ucst != nullptr))
    {
        const size_t mu = fp->mu;

        linf_coeff_tmp = malloc_device<float_type>(num_out_neurons_last_layer*mu);
        lsup_coeff_tmp = malloc_device<float_type>(num_out_neurons_last_layer*mu);
        uinf_coeff_tmp = malloc_device<float_type>(num_out_neurons_last_layer*mu);
        usup_coeff_tmp = malloc_device<float_type>(num_out_neurons_last_layer*mu);

        cudaMemset(linf_coeff_tmp, 0, num_out_neurons_last_layer*mu*sizeof(float_type));
        cudaMemset(lsup_coeff_tmp, 0, num_out_neurons_last_layer*mu*sizeof(float_type));
        cudaMemset(uinf_coeff_tmp, 0, num_out_neurons_last_layer*mu*sizeof(float_type));
        cudaMemset(usup_coeff_tmp, 0, num_out_neurons_last_layer*mu*sizeof(float_type));

        lcoeffs_from_input_poly<<<num_out_neurons_last_layer, num_threads>>>(linf_coeff, lsup_coeff, linf_coeff_tmp, lsup_coeff_tmp, fp->input_lweights, fp->input_uweights, num_in_neurons_0_layer, mu);
        ucoeffs_from_input_poly<<<num_out_neurons_last_layer, num_threads>>>(uinf_coeff, usup_coeff, uinf_coeff_tmp, usup_coeff_tmp, fp->input_lweights, fp->input_uweights, num_in_neurons_0_layer, mu);

        lcsts_from_input_poly<<<num_out_neurons_last_layer, 1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, linf_cst_tmp, lsup_cst_tmp, fp->input_lcst, fp->input_ucst, fp->input_inf, fp->input_sup, num_in_neurons_0_layer);
        ucsts_from_input_poly<<<num_out_neurons_last_layer, 1>>>(uinf_coeff, usup_coeff, uinf_cst, usup_cst, uinf_cst_tmp, usup_cst_tmp, fp->input_lcst, fp->input_ucst, fp->input_inf, fp->input_sup, num_in_neurons_0_layer);

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

        linf_coeff_tmp = nullptr;
        lsup_coeff_tmp = nullptr;
        uinf_coeff_tmp = nullptr;
        usup_coeff_tmp = nullptr;

        compute_lb_from_expr<<<num_out_neurons_last_layer, num_threads>>>(lb_array, linf_coeff, lsup_coeff, linf_cst, fp->input_inf + num_in_neurons_0_layer, fp->input_sup + num_in_neurons_0_layer, mu);
        compute_ub_from_expr<<<num_out_neurons_last_layer, num_threads>>>(ub_array, uinf_coeff, usup_coeff, usup_cst, fp->input_inf + num_in_neurons_0_layer, fp->input_sup + num_in_neurons_0_layer, mu);
    }
    else
    {
        compute_lb_from_expr<<<num_out_neurons_last_layer, num_threads>>>(lb_array, linf_coeff, lsup_coeff, linf_cst, fp->input_inf, fp->input_sup, num_in_neurons_0_layer);
        compute_ub_from_expr<<<num_out_neurons_last_layer, num_threads>>>(ub_array, uinf_coeff, usup_coeff, usup_cst, fp->input_inf, fp->input_sup, num_in_neurons_0_layer);
    }

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

    linf_coeff = nullptr;
    lsup_coeff = nullptr;
    linf_cst = nullptr;
    lsup_cst = nullptr;

    uinf_coeff = nullptr;
    usup_coeff = nullptr;
    uinf_cst = nullptr;
    usup_cst = nullptr;

    linf_cst_tmp = nullptr;
    lsup_cst_tmp = nullptr;

    uinf_cst_tmp = nullptr;
    usup_cst_tmp = nullptr;
}

void predict_size_of_conv_layer(fppoly_t* fp, size_t& current_size, size_t& last_size, long int& offset_x, long int& offset_y, long int& length_x, long int& length_y, long int& shift_x, long int& shift_y, const size_t k, const size_t num_out_neurons_last_layer)
{
    offset_x = fp->layers[k]->strides[0]*offset_x - fp->layers[k]->pad[0];
    offset_y = fp->layers[k]->strides[1]*offset_y - fp->layers[k]->pad[1];

    length_x = (length_x - 1)*fp->layers[k]->strides[0] + fp->layers[k]->filter_size[0];
    length_y = (length_y - 1)*fp->layers[k]->strides[1] + fp->layers[k]->filter_size[1];

    shift_x = fp->layers[k]->strides[0]*shift_x;
    shift_y = fp->layers[k]->strides[1]*shift_y;

    last_size = current_size;
    current_size = 4*num_out_neurons_last_layer*length_x*length_y*fp->layers[k]->input_size[2]*sizeof(float_type);
}


size_t predict_size(fppoly_t* fp, const size_t layerno)
{
    size_t backstep_counter = 0;

    size_t free_space;
    size_t total_space;

    cudaMemGetInfo(&free_space, &total_space);

    const size_t num_out_neurons_last_layer = fp->layers[layerno]->num_out_neurons;

    long int offset_x = -fp->layers[layerno]->pad[0];
    long int offset_y = -fp->layers[layerno]->pad[1];

    long int length_x = fp->layers[layerno]->filter_size[0];
    long int length_y = fp->layers[layerno]->filter_size[1];

    long int shift_x = fp->layers[layerno]->strides[0];
    long int shift_y = fp->layers[layerno]->strides[1];

    size_t current_size = 4*num_out_neurons_last_layer*length_x*length_y*fp->layers[layerno]->input_size[2]*sizeof(float_type);
    size_t last_size = 0;

    size_t maximum_size = 0;

    int k;

    if(fp->layers[layerno]->type == RESIDUAL)
    {
        const size_t predecessor1 = fp->layers[layerno]->predecessors[0] - 1;
        const size_t predecessor2 = fp->layers[layerno]->predecessors[1] - 1;

        char* predecessor_map = (char*) calloc(layerno, sizeof(char));
        int iter = predecessor1;

        while(iter >= 0)
        {
            predecessor_map[iter] = 1;
            iter = fp->layers[iter]->predecessors[0] - 1;
        }

        iter =  predecessor2;
        int common_predecessor = 0;

        while(iter >= 0)
        {
            if(predecessor_map[iter] == 1)
            {
                common_predecessor = iter;

                break;
            }

            iter = fp->layers[iter]->predecessors[0] - 1;
        }

        free(predecessor_map);

        long int copy_offset_x = offset_x;
        long int copy_offset_y = offset_y;

        long int copy_length_x = length_x;
        long int copy_length_y = length_y;

        long int copy_shift_x = shift_x;
        long int copy_shift_y = shift_y;

        size_t copy_current_size = current_size;
        size_t copy_last_size = last_size;

        iter = predecessor1;

        while(iter != common_predecessor)
        {
            predict_size_of_conv_layer(fp, copy_current_size, copy_last_size, copy_offset_x, copy_offset_y, copy_length_x, copy_length_y, copy_shift_x, copy_shift_y, iter, num_out_neurons_last_layer);

            if(copy_last_size + copy_current_size + current_size > maximum_size)
            {
                maximum_size = copy_last_size + copy_current_size + current_size;
            }

            iter = fp->layers[iter]->predecessors[0] - 1;
        }

        iter = predecessor2;

        while(iter != common_predecessor)
        {
            predict_size_of_conv_layer(fp, current_size, last_size, offset_x, offset_y, length_x, length_y, shift_x, shift_y, iter, num_out_neurons_last_layer);

            if(last_size + current_size + copy_current_size > maximum_size)
            {
                maximum_size = last_size + current_size + copy_current_size;
            }

            iter = fp->layers[iter]->predecessors[0] - 1;
        }

        if(copy_length_x > length_x)
        {
            std::swap(copy_current_size, current_size);

            std::swap(offset_x, copy_offset_x);
            std::swap(offset_y, copy_offset_y);

            std::swap(length_x, copy_length_x);
            std::swap(length_y, copy_length_y);

            std::swap(shift_x, copy_shift_x);
            std::swap(shift_y, copy_shift_y);
        }

        k = common_predecessor;
    }
    else
    {
        k = fp->layers[layerno]->predecessors[0] - 1;
    }

    while(k >= 0)
    {
        if(fp->layers[k]->type == RESIDUAL)
        {
            const size_t predecessor1 = fp->layers[k]->predecessors[0] - 1;
            const size_t predecessor2 = fp->layers[k]->predecessors[1] - 1;

            char* predecessor_map = (char*) calloc(k, sizeof(char));
            int iter = predecessor1;

            while(iter >= 0)
            {
                predecessor_map[iter] = 1;
                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            iter =  predecessor2;
            int common_predecessor = 0;

            while(iter >= 0)
            {
                if(predecessor_map[iter] == 1)
                {
                    common_predecessor = iter;

                    break;
                }

                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            free(predecessor_map);

            long int copy_offset_x = offset_x;
            long int copy_offset_y = offset_y;

            long int copy_length_x = length_x;
            long int copy_length_y = length_y;

            long int copy_shift_x = shift_x;
            long int copy_shift_y = shift_y;

            size_t copy_current_size = current_size;
            size_t copy_last_size = last_size;

            iter = predecessor1;

            while(iter != common_predecessor)
            {
                predict_size_of_conv_layer(fp, copy_current_size, copy_last_size, copy_offset_x, copy_offset_y, copy_length_x, copy_length_y, copy_shift_x, copy_shift_y, iter, num_out_neurons_last_layer);

                if(copy_last_size + copy_current_size + current_size > maximum_size)
                {
                    maximum_size = copy_last_size + copy_current_size + current_size;
                }

                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            iter = predecessor2;

            while(iter != common_predecessor)
            {
                predict_size_of_conv_layer(fp, current_size, last_size, offset_x, offset_y, length_x, length_y, shift_x, shift_y, iter, num_out_neurons_last_layer);

                if(last_size + current_size + copy_current_size > maximum_size)
                {
                    maximum_size = last_size + current_size + copy_current_size;
                }

                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            if(copy_length_x > length_x)
            {
                std::swap(copy_current_size, current_size);

                std::swap(offset_x, copy_offset_x);
                std::swap(offset_y, copy_offset_y);

                std::swap(length_x, copy_length_x);
                std::swap(length_y, copy_length_y);

                std::swap(shift_x, copy_shift_x);
                std::swap(shift_y, copy_shift_y);
            }

            k = common_predecessor;

            backstep_counter++;

            if(backstep_counter >= maximum_backstep)
            {
                break;
            }
        }
        else
        {
            predict_size_of_conv_layer(fp, current_size, last_size, offset_x, offset_y, length_x, length_y, shift_x, shift_y, k, num_out_neurons_last_layer);

            if(last_size + current_size > maximum_size)
            {
                maximum_size = last_size + current_size;
            }

            k = fp->layers[k]->predecessors[0] - 1;
        }
    }

    if((fp->input_lweights != nullptr) && (fp->input_uweights != nullptr) && (fp->input_lcst != nullptr) && (fp->input_ucst != nullptr))
    {
        last_size = current_size;
        current_size = 4*num_out_neurons_last_layer*fp->mu*sizeof(float_type);

        if(last_size + current_size > maximum_size)
        {
            maximum_size = last_size + current_size;
        }
    }

    size_t num_chunks = 1;

    while(maximum_size > (free_space - (1 << 27)))
    {
        maximum_size /= 2;
        num_chunks *= 2;
    }

    return num_chunks;
}


void update_state_using_predecessor_layer_sparse(fppoly_t* fp, float_type** linf_coeff, float_type** lsup_coeff, float_type** linf_cst, float_type** lsup_cst, float_type** uinf_coeff, float_type** usup_coeff, float_type** uinf_cst, float_type** usup_cst, float_type** linf_coeff_tmp, float_type** lsup_coeff_tmp, float_type** linf_cst_tmp, float_type** lsup_cst_tmp, float_type** uinf_coeff_tmp, float_type** usup_coeff_tmp, float_type** uinf_cst_tmp, float_type** usup_cst_tmp, const size_t layerno, const size_t k, const bool use_area_heuristic, const size_t x_y_size_last_layer, const size_t num_filters_last_layer, const size_t num_chunks, long int& offset_x, long int& offset_y, long int& length_x, long int& length_y, long int& shift_x, long int& shift_y)
{
    if(fp->layers[k]->activation == RELU)
    {
        expr_replace_relu_bounds_conv_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), fp->layers[k]->output_size[2]>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_cst, *lsup_cst, *uinf_cst, *usup_cst, fp->layers[k]->lb_array, fp->layers[k]->ub_array, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y, use_area_heuristic);
    }

    const long int new_offset_x = fp->layers[k]->strides[0]*offset_x - fp->layers[k]->pad[0];
    const long int new_offset_y = fp->layers[k]->strides[1]*offset_y - fp->layers[k]->pad[1];

    const long int new_length_x = (length_x - 1)*fp->layers[k]->strides[0] + fp->layers[k]->filter_size[0];
    const long int new_length_y = (length_y - 1)*fp->layers[k]->strides[1] + fp->layers[k]->filter_size[1];

    const long int new_shift_x = fp->layers[k]->strides[0]*shift_x;
    const long int new_shift_y = fp->layers[k]->strides[1]*shift_y;

    const size_t missing_length = new_length_x*new_length_y*fp->layers[k]->input_size[2];

    *linf_coeff_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*missing_length);
    *lsup_coeff_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*missing_length);
    *uinf_coeff_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*missing_length);
    *usup_coeff_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*missing_length);

    cudaMemset(*linf_coeff_tmp, 0, x_y_size_last_layer*num_filters_last_layer*missing_length*sizeof(float_type));
    cudaMemset(*lsup_coeff_tmp, 0, x_y_size_last_layer*num_filters_last_layer*missing_length*sizeof(float_type));
    cudaMemset(*uinf_coeff_tmp, 0, x_y_size_last_layer*num_filters_last_layer*missing_length*sizeof(float_type));
    cudaMemset(*usup_coeff_tmp, 0, x_y_size_last_layer*num_filters_last_layer*missing_length*sizeof(float_type));

    if(fp->layers[k]->type == CONV)
    {
        if(fp->layers[k]->input_size[0]*fp->layers[k]->input_size[1]*fp->layers[k]->input_size[2] > num_threads)
        {
            if(fp->layers[k]->input_size[1]*fp->layers[k]->input_size[2] > num_threads)
            {
                if(fp->layers[k]->input_size[2] > num_threads)
                {
                    coeffs_from_previous_layer_conv_sparse_filter_serial<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), num_threads>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, *uinf_coeff_tmp, *usup_coeff_tmp, fp->layers[k]->filter_weights, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y, new_length_x, new_length_y, fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
                }
                else
                {
                    coeffs_from_previous_layer_conv_sparse_filter_serial<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), fp->layers[k]->input_size[2]>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, *uinf_coeff_tmp, *usup_coeff_tmp, fp->layers[k]->filter_weights, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y, new_length_x, new_length_y, fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
                }
            }
            else
            {
                coeffs_from_previous_layer_conv_sparse_x_filter_serial<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), dim3(fp->layers[k]->input_size[2], fp->layers[k]->filter_size[1], 1)>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, *uinf_coeff_tmp, *usup_coeff_tmp, fp->layers[k]->filter_weights, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y, new_length_x, new_length_y, fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
            }
        }
        else
        {
            coeffs_from_previous_layer_conv_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), dim3(fp->layers[k]->input_size[2], fp->layers[k]->filter_size[1], fp->layers[k]->filter_size[0])>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, *uinf_coeff_tmp, *usup_coeff_tmp, fp->layers[k]->filter_weights, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y, new_length_x, new_length_y, fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
        }

        csts_from_previous_layer_conv_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), fp->layers[k]->output_size[2]>>>(*linf_coeff, *lsup_coeff, *uinf_coeff, *usup_coeff, *linf_cst, *lsup_cst, *uinf_cst, *usup_cst, *linf_cst_tmp, *lsup_cst_tmp, *uinf_cst_tmp, *usup_cst_tmp, fp->layers[k]->filter_bias, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y);
    }
    else
    {
        lcoeffs_from_previous_layer_pool_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), fp->layers[k]->output_size[2]>>>(*linf_coeff, *lsup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, fp->layers[k]->pool_lcoeffs, fp->layers[k]->pool_ucoeffs, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], offset_x, offset_y, length_x, length_y, shift_x, shift_y, new_length_x, new_length_y, fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
        ucoeffs_from_previous_layer_pool_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), fp->layers[k]->output_size[2]>>>(*uinf_coeff, *usup_coeff, *uinf_coeff_tmp, *usup_coeff_tmp, fp->layers[k]->pool_lcoeffs, fp->layers[k]->pool_ucoeffs, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], offset_x, offset_y, length_x, length_y, shift_x, shift_y, new_length_x, new_length_y, fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);

        lcsts_from_input_poly_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), 1>>>(*linf_coeff, *lsup_coeff, *linf_cst, *lsup_cst, *linf_cst_tmp, *lsup_cst_tmp, fp->layers[k]->pool_lcsts, fp->layers[k]->pool_ucsts, fp->layers[k]->lb_array, fp->layers[k]->ub_array, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y);
        ucsts_from_input_poly_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), 1>>>(*uinf_coeff, *usup_coeff, *uinf_cst, *usup_cst, *uinf_cst_tmp, *usup_cst_tmp, fp->layers[k]->pool_lcsts, fp->layers[k]->pool_ucsts, fp->layers[k]->lb_array, fp->layers[k]->ub_array, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y);
    }

    offset_x = new_offset_x;
    offset_y = new_offset_y;

    length_x = new_length_x;
    length_y = new_length_y;

    shift_x = new_shift_x;
    shift_y = new_shift_y;

    std::swap(*linf_coeff, *linf_coeff_tmp);
    std::swap(*lsup_coeff, *lsup_coeff_tmp);
    std::swap(*linf_cst, *linf_cst_tmp);
    std::swap(*lsup_cst, *lsup_cst_tmp);

    std::swap(*uinf_coeff, *uinf_coeff_tmp);
    std::swap(*usup_coeff, *usup_coeff_tmp);
    std::swap(*uinf_cst, *uinf_cst_tmp);
    std::swap(*usup_cst, *usup_cst_tmp);

    cudaFree(*linf_coeff_tmp);
    cudaFree(*lsup_coeff_tmp);
    cudaFree(*uinf_coeff_tmp);
    cudaFree(*usup_coeff_tmp);

    *linf_coeff_tmp = nullptr;
    *lsup_coeff_tmp = nullptr;
    *uinf_coeff_tmp = nullptr;
    *usup_coeff_tmp = nullptr;
}


__global__
void create_res_coeffs_csts(float_type* coeffs, float_type* bias, const int num_chunks, const int chunk_counter, const int input_size_z)
{
    const int out_x = blockIdx.x;
    const int out_y = blockIdx.y;
    const int out_z = blockIdx.z;

    const int local_mat_x = out_x*gridDim.y*gridDim.z + out_y*gridDim.z + out_z;

    coeffs[local_mat_x*input_size_z + chunk_counter*gridDim.z + out_z] = 1;

    bias[local_mat_x] = 0;
}


__global__
void compute_c_array(float_type* __restrict__ c_arrays, const float_type* __restrict__ gradient, const float_type* __restrict__ linf_coeffs, const float_type* __restrict__ usup_coeffs, const int num_chunks, const int chunk_counter, const int output_size_x, const int output_size_y, const int output_size_z, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    int out_z = threadIdx.x;

    const int local_n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;
    const int global_n = last_x*gridDim.y*num_chunks*gridDim.z + last_y*num_chunks*gridDim.z + chunk_counter*gridDim.z + last_z;

    const int min_out_x = offset_x + last_x*shift_x;
    const int min_out_y = offset_y + last_y*shift_y;

    const int min_x = (min_out_x < 0) ? -min_out_x : 0;
    const int min_y = (min_out_y < 0) ? -min_out_y : 0;

    const int max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
    const int max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;

    const float_type grad_entry = gradient[global_n];

    while(out_z < output_size_z)
    {
        for(int out_x = min_x; out_x < max_x; out_x++)
        {
            for(int out_y = min_y; out_y < max_y; out_y++)
            {
                const int j = out_x*length_y*output_size_z + out_y*output_size_z + out_z;
                const int mat_out = local_n*length_x*length_y*output_size_z + j;

                const int i = (out_x + min_out_x)*output_size_y*output_size_z + (out_y + min_out_y)*output_size_z + out_z;
                const int global_out = global_n*output_size_x*output_size_y*output_size_z + i;

                float_type linf_coeff = linf_coeffs[mat_out];
                float_type usup_coeff = usup_coeffs[mat_out];

                if((linf_coeff != 0) || (usup_coeff != 0))
                {
                    if(grad_entry < 0)
                    {
                        c_arrays[global_out] = grad_entry*linf_coeff;
                    }
                    else
                    {
                        c_arrays[global_out] = grad_entry*usup_coeff;
                    }
                }
            }
        }

        out_z += num_threads;
    }
}


__global__
void reduce_c_compute_xval(float_type* __restrict__ xval, const float_type* __restrict__ c_arrays, const float_type* __restrict__ input_inf, const float_type* __restrict__ input_sup, const int num_chunks, const int chunk_counter, const int num_neurons_last_layer, const int num_pixels)
{
    int j = threadIdx.x;

    while(j < num_pixels)
    {
        for(int i = 0; i < num_neurons_last_layer; i++)
        {
            xval[j] += c_arrays[i*num_pixels + j];
        }

        j += num_threads;
    }

    j = threadIdx.x;

    while(j < num_pixels)
    {
        if(xval[j] < 0)
        {
            xval[j] = input_inf[j];
        }
        else
        {
            xval[j] = input_sup[j];
        }

        j += num_threads;
    }
}


__global__
void compute_adv(float_type* __restrict__ adv, const float_type* __restrict__ gradient, const float_type* __restrict__ xval, const float_type* __restrict__ linf_coeffs, const float_type* __restrict__ usup_coeffs, const float_type* __restrict__ linf_csts, const float_type* __restrict__ usup_csts, const int num_chunks, const int chunk_counter, const int output_size_x, const int output_size_y, const int output_size_z, const int offset_x, const int offset_y, const int length_x, const int length_y, const int shift_x, const int shift_y)
{
    const int last_x = blockIdx.x;
    const int last_y = blockIdx.y;
    const int last_z = blockIdx.z;

    int out_z = threadIdx.x;

    const int local_n = last_x*gridDim.y*gridDim.z + last_y*gridDim.z + last_z;
    const int global_n = last_x*gridDim.y*num_chunks*gridDim.z + last_y*num_chunks*gridDim.z + chunk_counter*gridDim.z + last_z;

    const int min_out_x = offset_x + last_x*shift_x;
    const int min_out_y = offset_y + last_y*shift_y;

    const int min_x = (min_out_x < 0) ? -min_out_x : 0;
    const int min_y = (min_out_y < 0) ? -min_out_y : 0;

    const int max_x = (length_x + min_out_x > output_size_x) ? output_size_x - min_out_x : length_x;
    const int max_y = (length_y + min_out_y > output_size_y) ? output_size_y - min_out_y : length_y;

    if(threadIdx.x == 0)
    {
        if(gradient[global_n] < 0)
        {
            adv[global_n] = linf_csts[global_n];
        }
        else
        {
            adv[global_n] = usup_csts[global_n];
        }
    }

    __syncthreads();

    float_type advec = 0;

    while(out_z < output_size_z)
    {
        for(int out_x = min_x; out_x < max_x; out_x++)
        {
            for(int out_y = min_y; out_y < max_y; out_y++)
            {
                const int j = out_x*length_y*output_size_z + out_y*output_size_z + out_z;
                const int mat_out = local_n*length_x*length_y*output_size_z + j;

                const int i = (out_x + min_out_x)*output_size_y*output_size_z + (out_y + min_out_y)*output_size_z + out_z;

                float_type linf_coeff = linf_coeffs[mat_out];
                float_type usup_coeff = usup_coeffs[mat_out];

                if((linf_coeff != 0) || (usup_coeff != 0))
                {
                    if(gradient[global_n] < 0)
                    {
                        advec += linf_coeff*xval[i];
                    }
                    else
                    {
                        advec += usup_coeff*xval[i];
                    }
                }
            }
        }

        out_z += num_threads;
    }

    block_reduce_sum(advec, blockDim.x);

    if(threadIdx.x == 0)
    {
        adv[global_n] += advec;
    }
}


void update_state_using_previous_layers_sparse(elina_manager_t* man, fppoly_t* fp, const size_t layerno, const size_t num_chunks, const size_t chunk_counter, const bool use_area_heuristic, const bool retain_training_data, const float_type* gradient)
{
    size_t backstep_counter = 0;

    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    const size_t x_y_size_last_layer = fp->layers[layerno]->output_size[0]*fp->layers[layerno]->output_size[1];
    const size_t num_filters_last_layer = fp->layers[layerno]->output_size[2]/num_chunks;

    const size_t num_in_neurons_0_layer = fp->layers[0]->num_in_neurons;

    long int offset_x = -fp->layers[layerno]->pad[0];
    long int offset_y = -fp->layers[layerno]->pad[1];

    long int length_x = fp->layers[layerno]->filter_size[0];
    long int length_y = fp->layers[layerno]->filter_size[1];

    long int shift_x = fp->layers[layerno]->strides[0];
    long int shift_y = fp->layers[layerno]->strides[1];

    float_type* coeffs = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*length_x*length_y*fp->layers[layerno]->input_size[2]);
    float_type* csts = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);

    cudaMemset(coeffs, 0, x_y_size_last_layer*num_filters_last_layer*length_x*length_y*fp->layers[layerno]->input_size[2]*sizeof(float_type));
    cudaMemset(csts, 0, x_y_size_last_layer*num_filters_last_layer*sizeof(float_type));

    if(fp->layers[layerno]->type == RESIDUAL)
    {
        create_res_coeffs_csts<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), 1>>>(coeffs, csts, num_chunks, chunk_counter, fp->layers[layerno]->input_size[2]);
    }
    else
    {
        device_layer_create_sparse_exprs<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), 1>>>(coeffs, csts, fp->layers[layerno]->filter_weights, fp->layers[layerno]->filter_bias, chunk_counter, fp->layers[layerno]->input_size[0], fp->layers[layerno]->input_size[1], fp->layers[layerno]->input_size[2], fp->layers[layerno]->filter_size[0], fp->layers[layerno]->filter_size[1], fp->layers[layerno]->strides[0], fp->layers[layerno]->strides[1], fp->layers[layerno]->pad[0], fp->layers[layerno]->pad[1]);
    }

    float_type* lb_array = fp->layers[layerno]->lb_array;
    float_type* ub_array = fp->layers[layerno]->ub_array;

    float_type* linf_coeff = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*length_x*length_y*fp->layers[layerno]->input_size[2]);
    float_type* lsup_coeff = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*length_x*length_y*fp->layers[layerno]->input_size[2]);
    float_type* uinf_coeff = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*length_x*length_y*fp->layers[layerno]->input_size[2]);
    float_type* usup_coeff = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*length_x*length_y*fp->layers[layerno]->input_size[2]);

    float_type* linf_cst = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
    float_type* lsup_cst = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
    float_type* uinf_cst = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
    float_type* usup_cst = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);

    copy_coeffs_and_csts<<<x_y_size_last_layer*num_filters_last_layer, 1>>>(linf_coeff, linf_cst, coeffs, csts, length_x*length_y*fp->layers[layerno]->input_size[2]);
    copy_coeffs_and_csts<<<x_y_size_last_layer*num_filters_last_layer, 1>>>(lsup_coeff, lsup_cst, coeffs, csts, length_x*length_y*fp->layers[layerno]->input_size[2]);
    copy_coeffs_and_csts<<<x_y_size_last_layer*num_filters_last_layer, 1>>>(uinf_coeff, uinf_cst, coeffs, csts, length_x*length_y*fp->layers[layerno]->input_size[2]);
    copy_coeffs_and_csts<<<x_y_size_last_layer*num_filters_last_layer, 1>>>(usup_coeff, usup_cst, coeffs, csts, length_x*length_y*fp->layers[layerno]->input_size[2]);

    cudaFree(coeffs);
    cudaFree(csts);

    coeffs = nullptr;
    csts = nullptr;

    float_type* linf_coeff_tmp;
    float_type* lsup_coeff_tmp;
    float_type* uinf_coeff_tmp;
    float_type* usup_coeff_tmp;

    float_type* linf_cst_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
    float_type* lsup_cst_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
    float_type* uinf_cst_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
    float_type* usup_cst_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);

    int k;

    if(fp->layers[layerno]->type == RESIDUAL)
    {
        const size_t predecessor1 = fp->layers[layerno]->predecessors[0] - 1;
        const size_t predecessor2 = fp->layers[layerno]->predecessors[1] - 1;

        char* predecessor_map = (char*) calloc(layerno, sizeof(char));
        int iter = predecessor1;

        while(iter >= 0)
        {
            predecessor_map[iter] = 1;
            iter = fp->layers[iter]->predecessors[0] - 1;
        }

        iter =  predecessor2;
        int common_predecessor = 0;

        while(iter >= 0)
        {
            if(predecessor_map[iter] == 1)
            {
                common_predecessor = iter;

                break;
            }

            iter = fp->layers[iter]->predecessors[0] - 1;
        }

        free(predecessor_map);

        const size_t num_filters_current_layer = fp->layers[layerno]->output_size[2];

        float_type* copy_linf_coeff = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*length_x*length_y*num_filters_current_layer);
        float_type* copy_lsup_coeff = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*length_x*length_y*num_filters_current_layer);
        float_type* copy_uinf_coeff = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*length_x*length_y*num_filters_current_layer);
        float_type* copy_usup_coeff = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*length_x*length_y*num_filters_current_layer);

        float_type* copy_linf_cst = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
        float_type* copy_lsup_cst = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
        float_type* copy_uinf_cst = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
        float_type* copy_usup_cst = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);

        copy_coeffs_and_csts<<<x_y_size_last_layer*num_filters_last_layer, 1>>>(copy_linf_coeff, copy_linf_cst, linf_coeff, linf_cst, length_x*length_y*num_filters_current_layer);
        copy_coeffs_and_csts<<<x_y_size_last_layer*num_filters_last_layer, 1>>>(copy_lsup_coeff, copy_lsup_cst, lsup_coeff, lsup_cst, length_x*length_y*num_filters_current_layer);
        copy_coeffs_and_csts<<<x_y_size_last_layer*num_filters_last_layer, 1>>>(copy_uinf_coeff, copy_uinf_cst, uinf_coeff, uinf_cst, length_x*length_y*num_filters_current_layer);
        copy_coeffs_and_csts<<<x_y_size_last_layer*num_filters_last_layer, 1>>>(copy_usup_coeff, copy_usup_cst, usup_coeff, usup_cst, length_x*length_y*num_filters_current_layer);

        cudaMemset(copy_linf_cst, 0, x_y_size_last_layer*num_filters_last_layer*sizeof(float_type));
        cudaMemset(copy_lsup_cst, 0, x_y_size_last_layer*num_filters_last_layer*sizeof(float_type));
        cudaMemset(copy_uinf_cst, 0, x_y_size_last_layer*num_filters_last_layer*sizeof(float_type));
        cudaMemset(copy_usup_cst, 0, x_y_size_last_layer*num_filters_last_layer*sizeof(float_type));

        float_type* copy_linf_coeff_tmp;
        float_type* copy_lsup_coeff_tmp;
        float_type* copy_uinf_coeff_tmp;
        float_type* copy_usup_coeff_tmp;

        float_type* copy_linf_cst_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
        float_type* copy_lsup_cst_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
        float_type* copy_uinf_cst_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
        float_type* copy_usup_cst_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);

        long int copy_offset_x = offset_x;
        long int copy_offset_y = offset_y;

        long int copy_length_x = length_x;
        long int copy_length_y = length_y;

        long int copy_shift_x = shift_x;
        long int copy_shift_y = shift_y;

        iter = predecessor1;

        while(iter != common_predecessor)
        {
            update_state_using_predecessor_layer_sparse(fp, &copy_linf_coeff, &copy_lsup_coeff, &copy_linf_cst, &copy_lsup_cst, &copy_uinf_coeff, &copy_usup_coeff, &copy_uinf_cst, &copy_usup_cst, &copy_linf_coeff_tmp, &copy_lsup_coeff_tmp, &copy_linf_cst_tmp, &copy_lsup_cst_tmp, &copy_uinf_coeff_tmp, &copy_usup_coeff_tmp, &copy_uinf_cst_tmp, &copy_usup_cst_tmp, layerno, iter, use_area_heuristic, x_y_size_last_layer, num_filters_last_layer, num_chunks, copy_offset_x, copy_offset_y, copy_length_x, copy_length_y, copy_shift_x, copy_shift_y);

            iter = fp->layers[iter]->predecessors[0] - 1;
        }

        iter = predecessor2;

        while(iter != common_predecessor)
        {
            update_state_using_predecessor_layer_sparse(fp, &linf_coeff, &lsup_coeff, &linf_cst, &lsup_cst, &uinf_coeff, &usup_coeff, &uinf_cst, &usup_cst, &linf_coeff_tmp, &lsup_coeff_tmp, &linf_cst_tmp, &lsup_cst_tmp, &uinf_coeff_tmp, &usup_coeff_tmp, &uinf_cst_tmp, &usup_cst_tmp, layerno, iter, use_area_heuristic, x_y_size_last_layer, num_filters_last_layer, num_chunks, offset_x, offset_y, length_x, length_y, shift_x, shift_y);

            iter = fp->layers[iter]->predecessors[0] - 1;
        }

        if(copy_length_x > length_x)
        {
            std::swap(linf_coeff, copy_linf_coeff);
            std::swap(lsup_coeff, copy_lsup_coeff);
            std::swap(linf_cst, copy_linf_cst);
            std::swap(lsup_cst, copy_lsup_cst);

            std::swap(uinf_coeff, copy_uinf_coeff);
            std::swap(usup_coeff, copy_usup_coeff);
            std::swap(uinf_cst, copy_uinf_cst);
            std::swap(usup_cst, copy_usup_cst);

            std::swap(offset_x, copy_offset_x);
            std::swap(offset_y, copy_offset_y);

            std::swap(length_x, copy_length_x);
            std::swap(length_y, copy_length_y);

            std::swap(shift_x, copy_shift_x);
            std::swap(shift_y, copy_shift_y);
        }

        add_coeffs_and_csts_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), 1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, copy_linf_coeff, copy_lsup_coeff, copy_linf_cst, copy_lsup_cst, length_x, length_y, copy_length_x, copy_length_y,  fp->layers[common_predecessor]->output_size[2], copy_offset_x - offset_x, copy_offset_y - offset_y);
        add_coeffs_and_csts_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), 1>>>(uinf_coeff, usup_coeff, uinf_cst, usup_cst, copy_uinf_coeff, copy_usup_coeff, copy_uinf_cst, copy_usup_cst, length_x, length_y, copy_length_x, copy_length_y,  fp->layers[common_predecessor]->output_size[2], copy_offset_x - offset_x, copy_offset_y - offset_y);

        cudaFree(copy_linf_coeff);
        cudaFree(copy_lsup_coeff);
        cudaFree(copy_linf_cst);
        cudaFree(copy_lsup_cst);

        cudaFree(copy_uinf_coeff);
        cudaFree(copy_usup_coeff);
        cudaFree(copy_uinf_cst);
        cudaFree(copy_usup_cst);

        cudaFree(copy_linf_cst_tmp);
        cudaFree(copy_lsup_cst_tmp);

        cudaFree(copy_uinf_cst_tmp);
        cudaFree(copy_usup_cst_tmp);

        copy_linf_coeff = nullptr;
        copy_lsup_coeff = nullptr;
        copy_linf_cst = nullptr;
        copy_lsup_cst = nullptr;

        copy_uinf_coeff = nullptr;
        copy_usup_coeff = nullptr;
        copy_uinf_cst = nullptr;
        copy_usup_cst = nullptr;

        copy_linf_cst_tmp = nullptr;
        copy_lsup_cst_tmp = nullptr;

        copy_uinf_cst_tmp = nullptr;
        copy_usup_cst_tmp = nullptr;

        k = common_predecessor;
    }
    else
    {
        k = fp->layers[layerno]->predecessors[0] - 1;
    }

    while(k >= 0)
    {
        if(fp->layers[k]->type == RESIDUAL)
        {
            if(fp->layers[k]->activation == RELU)
            {
                expr_replace_relu_bounds_conv_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), fp->layers[k]->output_size[2]>>>(linf_coeff, lsup_coeff, uinf_coeff, usup_coeff, linf_cst, lsup_cst, uinf_cst, usup_cst, fp->layers[k]->lb_array, fp->layers[k]->ub_array, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y, use_area_heuristic);
            }

            const size_t predecessor1 = fp->layers[k]->predecessors[0] - 1;
            const size_t predecessor2 = fp->layers[k]->predecessors[1] - 1;

            char* predecessor_map = (char*) calloc(k, sizeof(char));
            int iter = predecessor1;

            while(iter >= 0)
            {
                predecessor_map[iter] = 1;
                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            iter =  predecessor2;
            int common_predecessor = 0;

            while(iter >= 0)
            {
                if(predecessor_map[iter] == 1)
                {
                    common_predecessor = iter;

                    break;
                }

                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            free(predecessor_map);

            const size_t num_filters_current_layer = fp->layers[k]->output_size[2];

            float_type* copy_linf_coeff = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*length_x*length_y*num_filters_current_layer);
            float_type* copy_lsup_coeff = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*length_x*length_y*num_filters_current_layer);
            float_type* copy_uinf_coeff = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*length_x*length_y*num_filters_current_layer);
            float_type* copy_usup_coeff = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*length_x*length_y*num_filters_current_layer);

            float_type* copy_linf_cst = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
            float_type* copy_lsup_cst = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
            float_type* copy_uinf_cst = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
            float_type* copy_usup_cst = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);

            copy_coeffs_and_csts<<<x_y_size_last_layer*num_filters_last_layer, 1>>>(copy_linf_coeff, copy_linf_cst, linf_coeff, linf_cst, length_x*length_y*num_filters_current_layer);
            copy_coeffs_and_csts<<<x_y_size_last_layer*num_filters_last_layer, 1>>>(copy_lsup_coeff, copy_lsup_cst, lsup_coeff, lsup_cst, length_x*length_y*num_filters_current_layer);
            copy_coeffs_and_csts<<<x_y_size_last_layer*num_filters_last_layer, 1>>>(copy_uinf_coeff, copy_uinf_cst, uinf_coeff, uinf_cst, length_x*length_y*num_filters_current_layer);
            copy_coeffs_and_csts<<<x_y_size_last_layer*num_filters_last_layer, 1>>>(copy_usup_coeff, copy_usup_cst, usup_coeff, usup_cst, length_x*length_y*num_filters_current_layer);

            cudaMemset(copy_linf_cst, 0, x_y_size_last_layer*num_filters_last_layer*sizeof(float_type));
            cudaMemset(copy_lsup_cst, 0, x_y_size_last_layer*num_filters_last_layer*sizeof(float_type));
            cudaMemset(copy_uinf_cst, 0, x_y_size_last_layer*num_filters_last_layer*sizeof(float_type));
            cudaMemset(copy_usup_cst, 0, x_y_size_last_layer*num_filters_last_layer*sizeof(float_type));

            float_type* copy_linf_coeff_tmp;
            float_type* copy_lsup_coeff_tmp;
            float_type* copy_uinf_coeff_tmp;
            float_type* copy_usup_coeff_tmp;

            float_type* copy_linf_cst_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
            float_type* copy_lsup_cst_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
            float_type* copy_uinf_cst_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
            float_type* copy_usup_cst_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);

            long int copy_offset_x = offset_x;
            long int copy_offset_y = offset_y;

            long int copy_length_x = length_x;
            long int copy_length_y = length_y;

            long int copy_shift_x = shift_x;
            long int copy_shift_y = shift_y;

            iter = predecessor1;

            while(iter != common_predecessor)
            {
                update_state_using_predecessor_layer_sparse(fp, &copy_linf_coeff, &copy_lsup_coeff, &copy_linf_cst, &copy_lsup_cst, &copy_uinf_coeff, &copy_usup_coeff, &copy_uinf_cst, &copy_usup_cst, &copy_linf_coeff_tmp, &copy_lsup_coeff_tmp, &copy_linf_cst_tmp, &copy_lsup_cst_tmp, &copy_uinf_coeff_tmp, &copy_usup_coeff_tmp, &copy_uinf_cst_tmp, &copy_usup_cst_tmp, layerno, iter, use_area_heuristic, x_y_size_last_layer, num_filters_last_layer, num_chunks, copy_offset_x, copy_offset_y, copy_length_x, copy_length_y, copy_shift_x, copy_shift_y);

                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            iter = predecessor2;

            while(iter != common_predecessor)
            {
                update_state_using_predecessor_layer_sparse(fp, &linf_coeff, &lsup_coeff, &linf_cst, &lsup_cst, &uinf_coeff, &usup_coeff, &uinf_cst, &usup_cst, &linf_coeff_tmp, &lsup_coeff_tmp, &linf_cst_tmp, &lsup_cst_tmp, &uinf_coeff_tmp, &usup_coeff_tmp, &uinf_cst_tmp, &usup_cst_tmp, layerno, iter, use_area_heuristic, x_y_size_last_layer, num_filters_last_layer, num_chunks, offset_x, offset_y, length_x, length_y, shift_x, shift_y);

                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            if(copy_length_x > length_x)
            {
                std::swap(linf_coeff, copy_linf_coeff);
                std::swap(lsup_coeff, copy_lsup_coeff);
                std::swap(linf_cst, copy_linf_cst);
                std::swap(lsup_cst, copy_lsup_cst);

                std::swap(uinf_coeff, copy_uinf_coeff);
                std::swap(usup_coeff, copy_usup_coeff);
                std::swap(uinf_cst, copy_uinf_cst);
                std::swap(usup_cst, copy_usup_cst);

                std::swap(offset_x, copy_offset_x);
                std::swap(offset_y, copy_offset_y);

                std::swap(length_x, copy_length_x);
                std::swap(length_y, copy_length_y);

                std::swap(shift_x, copy_shift_x);
                std::swap(shift_y, copy_shift_y);
            }

            add_coeffs_and_csts_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), 1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, copy_linf_coeff, copy_lsup_coeff, copy_linf_cst, copy_lsup_cst, length_x, length_y, copy_length_x, copy_length_y,  fp->layers[common_predecessor]->output_size[2], copy_offset_x - offset_x, copy_offset_y - offset_y);
            add_coeffs_and_csts_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), 1>>>(uinf_coeff, usup_coeff, uinf_cst, usup_cst, copy_uinf_coeff, copy_usup_coeff, copy_uinf_cst, copy_usup_cst, length_x, length_y, copy_length_x, copy_length_y,  fp->layers[common_predecessor]->output_size[2], copy_offset_x - offset_x, copy_offset_y - offset_y);

            cudaFree(copy_linf_coeff);
            cudaFree(copy_lsup_coeff);
            cudaFree(copy_linf_cst);
            cudaFree(copy_lsup_cst);

            cudaFree(copy_uinf_coeff);
            cudaFree(copy_usup_coeff);
            cudaFree(copy_uinf_cst);
            cudaFree(copy_usup_cst);

            cudaFree(copy_linf_cst_tmp);
            cudaFree(copy_lsup_cst_tmp);

            cudaFree(copy_uinf_cst_tmp);
            cudaFree(copy_usup_cst_tmp);

            copy_linf_coeff = nullptr;
            copy_lsup_coeff = nullptr;
            copy_linf_cst = nullptr;
            copy_lsup_cst = nullptr;

            copy_uinf_coeff = nullptr;
            copy_usup_coeff = nullptr;
            copy_uinf_cst = nullptr;
            copy_usup_cst = nullptr;

            copy_linf_cst_tmp = nullptr;
            copy_lsup_cst_tmp = nullptr;

            copy_uinf_cst_tmp = nullptr;
            copy_usup_cst_tmp = nullptr;

            k = common_predecessor;

            backstep_counter++;

            if(backstep_counter >= maximum_backstep)
            {
                if(fp->layers[k]->activation == RELU)
                {
                    expr_replace_relu_bounds_conv_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), fp->layers[k]->output_size[2]>>>(linf_coeff, lsup_coeff, uinf_coeff, usup_coeff, linf_cst, lsup_cst, uinf_cst, usup_cst, fp->layers[k]->lb_array, fp->layers[k]->ub_array, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y, use_area_heuristic);
                }

                compute_lb_from_expr_conv_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), fp->layers[k]->output_size[2]>>>(lb_array, linf_coeff, lsup_coeff, linf_cst, fp->layers[k]->lb_array, fp->layers[k]->ub_array, num_chunks, chunk_counter, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y);
                compute_ub_from_expr_conv_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), fp->layers[k]->output_size[2]>>>(ub_array, uinf_coeff, usup_coeff, usup_cst, fp->layers[k]->lb_array, fp->layers[k]->ub_array, num_chunks, chunk_counter, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y);

                break;
            }
        }
        else
        {
            update_state_using_predecessor_layer_sparse(fp, &linf_coeff, &lsup_coeff, &linf_cst, &lsup_cst, &uinf_coeff, &usup_coeff, &uinf_cst, &usup_cst, &linf_coeff_tmp, &lsup_coeff_tmp, &linf_cst_tmp, &lsup_cst_tmp, &uinf_coeff_tmp, &usup_coeff_tmp, &uinf_cst_tmp, &usup_cst_tmp, layerno, k, use_area_heuristic, x_y_size_last_layer, num_filters_last_layer, num_chunks, offset_x, offset_y, length_x, length_y, shift_x, shift_y);

            k = fp->layers[k]->predecessors[0] - 1;
        }
    }

    if(backstep_counter < maximum_backstep)
    {
        if((fp->input_lweights != nullptr) && (fp->input_uweights != nullptr) && (fp->input_lcst != nullptr) && (fp->input_ucst != nullptr))
        {
            const size_t mu = fp->mu;

            linf_coeff_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*mu);
            lsup_coeff_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*mu);
            uinf_coeff_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*mu);
            usup_coeff_tmp = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*mu);

            cudaMemset(linf_coeff_tmp, 0, x_y_size_last_layer*num_filters_last_layer*mu*sizeof(float_type));
            cudaMemset(lsup_coeff_tmp, 0, x_y_size_last_layer*num_filters_last_layer*mu*sizeof(float_type));
            cudaMemset(uinf_coeff_tmp, 0, x_y_size_last_layer*num_filters_last_layer*mu*sizeof(float_type));
            cudaMemset(usup_coeff_tmp, 0, x_y_size_last_layer*num_filters_last_layer*mu*sizeof(float_type));

            lcoeffs_from_input_poly_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), num_threads>>>(linf_coeff, lsup_coeff, linf_coeff_tmp, lsup_coeff_tmp, fp->input_lweights, fp->input_uweights, fp->layers[0]->input_size[0], fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], mu, offset_x, offset_y, length_x, length_y, shift_x, shift_y);
            ucoeffs_from_input_poly_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), num_threads>>>(uinf_coeff, usup_coeff, uinf_coeff_tmp, usup_coeff_tmp, fp->input_lweights, fp->input_uweights, fp->layers[0]->input_size[0], fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], mu, offset_x, offset_y, length_x, length_y, shift_x, shift_y);

            lcsts_from_input_poly_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), 1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, linf_cst_tmp, lsup_cst_tmp, fp->input_lcst, fp->input_ucst, fp->input_inf, fp->input_sup, fp->layers[0]->input_size[0], fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y);
            ucsts_from_input_poly_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), 1>>>(uinf_coeff, usup_coeff, uinf_cst, usup_cst, uinf_cst_tmp, usup_cst_tmp, fp->input_lcst, fp->input_ucst, fp->input_inf, fp->input_sup, fp->layers[0]->input_size[0], fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y);

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

            linf_coeff_tmp = nullptr;
            lsup_coeff_tmp = nullptr;
            uinf_coeff_tmp = nullptr;
            usup_coeff_tmp = nullptr;

            compute_lb_from_expr_input_poly_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), 1>>>(lb_array, linf_coeff, lsup_coeff, linf_cst, fp->input_inf + num_in_neurons_0_layer, fp->input_sup + num_in_neurons_0_layer, num_chunks, chunk_counter, mu);
            compute_ub_from_expr_input_poly_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), 1>>>(ub_array, uinf_coeff, usup_coeff, usup_cst, fp->input_inf + num_in_neurons_0_layer, fp->input_sup + num_in_neurons_0_layer, num_chunks, chunk_counter, mu);
        }
        else
        {
            compute_lb_from_expr_conv_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), fp->layers[0]->input_size[2]>>>(lb_array, linf_coeff, lsup_coeff, linf_cst, fp->input_inf, fp->input_sup, num_chunks, chunk_counter, fp->layers[0]->input_size[0], fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y);
            compute_ub_from_expr_conv_sparse<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), fp->layers[0]->input_size[2]>>>(ub_array, uinf_coeff, usup_coeff, usup_cst, fp->input_inf, fp->input_sup, num_chunks, chunk_counter, fp->layers[0]->input_size[0], fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y);
        }
    }

    if(retain_training_data)
    {
        float_type* gradient_dev = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);

        copy_to_device(gradient_dev, gradient, x_y_size_last_layer*num_filters_last_layer);

        num_neurons_training_layer = x_y_size_last_layer*num_filters_last_layer;

        float_type* c_array = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer*num_in_neurons_0_layer);
        cudaMemset(c_array, 0, x_y_size_last_layer*num_filters_last_layer*num_in_neurons_0_layer*sizeof(float_type));

        const int input_size_x = fp->layers[0]->input_size[0];
        const int input_size_y = fp->layers[0]->input_size[1];
        const int input_size_z = fp->layers[0]->input_size[2];

        compute_c_array<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), num_threads>>>(c_array, gradient_dev, linf_coeff, usup_coeff, num_chunks, chunk_counter, input_size_x, input_size_y, input_size_z, offset_x, offset_y, length_x, length_y, shift_x, shift_y);

        float_type* xval = malloc_device<float_type>(num_in_neurons_0_layer);
        cudaMemset(xval, 0, num_in_neurons_0_layer*sizeof(float_type));

        reduce_c_compute_xval<<<1, num_threads>>>(xval, c_array, fp->input_inf, fp->input_sup, num_chunks, chunk_counter, x_y_size_last_layer*num_filters_last_layer, num_in_neurons_0_layer);

        float_type* adv = malloc_device<float_type>(x_y_size_last_layer*num_filters_last_layer);
        cudaMemset(adv, 0, x_y_size_last_layer*num_filters_last_layer*sizeof(float_type));

        compute_adv<<<dim3(fp->layers[layerno]->output_size[0], fp->layers[layerno]->output_size[1], fp->layers[layerno]->output_size[2]/num_chunks), num_threads>>>(adv, gradient_dev, xval, linf_coeff, usup_coeff, linf_cst, usup_cst, num_chunks, chunk_counter, input_size_x, input_size_y, input_size_z, offset_x, offset_y, length_x, length_y, shift_x, shift_y);

        adv_host = (float_type*) malloc(x_y_size_last_layer*num_filters_last_layer*sizeof(float_type));
        copy_to_host(adv_host, adv, x_y_size_last_layer*num_filters_last_layer);

        cudaFree(gradient_dev);
        cudaFree(c_array);
        cudaFree(xval);
        cudaFree(adv);

        gradient_dev = nullptr;
        c_array = nullptr;
        xval = nullptr;
        adv = nullptr;
    }

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

    linf_coeff = nullptr;
    lsup_coeff = nullptr;
    linf_cst = nullptr;
    lsup_cst = nullptr;

    uinf_coeff = nullptr;
    usup_coeff = nullptr;
    uinf_cst = nullptr;
    usup_cst = nullptr;

    linf_cst_tmp = nullptr;
    lsup_cst_tmp = nullptr;

    uinf_cst_tmp = nullptr;
    usup_cst_tmp = nullptr;
}


void update_state_using_previous_layers_sparse_full(elina_manager_t* man, fppoly_t* fp, const size_t layerno, const bool use_area_heuristic, const bool retain_training_data, const float_type* gradient)
{
    const size_t num_chunks = predict_size(fp, layerno);

    for(size_t chunk_counter = 0; chunk_counter < num_chunks; chunk_counter++)
    {
        update_state_using_previous_layers_sparse(man, fp, fp->numlayers - 1, num_chunks, chunk_counter, use_area_heuristic, retain_training_data, gradient);
    }
}


void ffn_handle_intermediate_layer(elina_manager_t* man, elina_abstract0_t* element, const float_type* const * weights, const float_type* cst, const size_t num_out_neurons, const size_t num_in_neurons, size_t* predecessors, const activation_type_t activation, const bool alloc, const bool use_area_heuristic, fnn_op OP)
{
    fppoly_t* fp = fppoly_of_abstract0(element);
    ffn_add_layer(fp, num_out_neurons, num_in_neurons, FFN, activation);

    fp->layers[fp->numlayers - 1]->predecessors = predecessors;

    float_type* coeffs = fp->layers[fp->numlayers - 1]->coeffs;
    float_type* csts = fp->layers[fp->numlayers - 1]->csts;

    if(OP == MUL || OP == SUB1 || OP == SUB2)
    {
        float_type** dense_weights = (float_type**) malloc(num_in_neurons*sizeof(float_type*));

        for(size_t j = 0; j < num_in_neurons; j++)
        {
            dense_weights[j] = (float_type*) calloc(num_in_neurons, sizeof(float_type));
        }

        if(OP == MUL)
        {
            float_type* dense_cst = (float_type*) calloc(num_in_neurons, sizeof(float_type));

            for(size_t j = 0; j < num_in_neurons; j++)
            {
                dense_weights[j][j] = cst[j];
            }

            layer_create_dense_exprs(coeffs, csts, dense_weights, dense_cst, num_in_neurons, num_in_neurons);

            free(dense_cst);
        }
        else if(OP == SUB1)
        {
            for(size_t j = 0; j < num_in_neurons; j++)
            {
                dense_weights[j][j] = -1;
            }

            layer_create_dense_exprs(coeffs, csts, dense_weights, cst, num_in_neurons, num_in_neurons);
        }
        else if(OP == SUB2)
        {
            float_type* dense_cst = (float_type*) calloc(num_in_neurons, sizeof(float_type));

            for(size_t j = 0; j < num_in_neurons; j++)
            {
                dense_weights[j][j] = 1;
            }

            for(size_t j = 0; j < num_in_neurons; j++)
            {
                dense_cst[j] = -cst[j];
            }

            layer_create_dense_exprs(coeffs, csts, dense_weights, dense_cst, num_in_neurons, num_in_neurons);

            free(dense_cst);
        }

        for(size_t j = 0; j < num_in_neurons; j++)
        {
            free(dense_weights[j]);
        }

        free(dense_weights);
    }
    else
    {
        layer_create_dense_exprs(coeffs, csts, weights, cst, num_out_neurons, num_in_neurons);
    }

    update_state_using_previous_layers(man, fp, fp->numlayers - 1, use_area_heuristic);
}


void ffn_handle_intermediate_affine_layer(elina_manager_t* man, elina_abstract0_t* element, const float_type* const * weights, const float_type* bias, const size_t num_out_neurons, const size_t num_in_neurons, size_t* predecessors, const bool use_area_heuristic)
{
    ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, NONE, true, use_area_heuristic, MATMULT);
}


void ffn_handle_intermediate_relu_layer(elina_manager_t* man, elina_abstract0_t* element, const float_type* const * weights, const float_type* bias, const size_t num_out_neurons, const size_t num_in_neurons, size_t* predecessors, const bool use_area_heuristic)
{
    ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, RELU, true, use_area_heuristic, MATMULT);
}


void ffn_handle_intermediate_sub_layer(elina_manager_t* man, elina_abstract0_t* element, const float_type* cst, const bool is_minuend, const size_t num_in_neurons, size_t* predecessors, const bool use_area_heuristic)
{
    if(is_minuend == true)
    {
        ffn_handle_intermediate_layer(man, element, nullptr, cst, num_in_neurons, num_in_neurons, predecessors, NONE, true, use_area_heuristic, SUB1);
    }
    else
    {
        ffn_handle_intermediate_layer(man, element, nullptr, cst, num_in_neurons, num_in_neurons, predecessors, NONE, true, use_area_heuristic, SUB2);
    }
}

void ffn_handle_intermediate_mul_layer(elina_manager_t* man, elina_abstract0_t* element, const float_type* bias, const size_t num_in_neurons, size_t* predecessors, const bool use_area_heuristic)
{
    ffn_handle_intermediate_layer(man, element, nullptr, bias, num_in_neurons, num_in_neurons, predecessors, NONE, true, use_area_heuristic, MUL);
}


void ffn_handle_last_layer(elina_manager_t* man, elina_abstract0_t* element, const float_type* const * weights, const float_type* bias, const size_t num_out_neurons, const size_t num_in_neurons, size_t* predecessors, const bool has_activation, const activation_type_t activation, const bool alloc, const bool use_area_heuristic)
{
    fppoly_t* fp = fppoly_of_abstract0(element);
    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    if(has_activation)
    {
        ffn_add_layer(fp, num_out_neurons, num_in_neurons, FFN, activation);
    }
    else
    {
        ffn_add_layer(fp, num_out_neurons, num_in_neurons, FFN, NONE);
    }

    float_type* coeffs = fp->layers[fp->numlayers - 1]->coeffs;
    float_type* csts = fp->layers[fp->numlayers - 1]->csts;

    fp->layers[fp->numlayers - 1]->predecessors = predecessors;

    layer_create_dense_exprs(coeffs, csts, weights, bias, num_out_neurons, num_in_neurons);

    update_state_using_previous_layers(man, fp, fp->numlayers - 1, use_area_heuristic);
}

void ffn_handle_last_relu_layer(elina_manager_t* man, elina_abstract0_t* element, const float_type* const * weights, const float_type* bias, const size_t num_out_neurons, const size_t num_in_neurons, size_t* predecessors, const bool has_relu, const bool use_area_heuristic)
{
    ffn_handle_last_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, predecessors, has_relu, RELU, true, use_area_heuristic);
}


__global__
void create_sub_expr(float_type* __restrict__ inf_coeff, float_type* __restrict__ sup_coeff, float_type* __restrict__ inf_cst, float_type* __restrict__ sup_cst, const int index, const elina_dim_t y, const elina_dim_t x)
{
    inf_cst[index] = 0;
    sup_cst[index] = 0;

    for(int i = 0; i < 10; i++)
    {
        inf_coeff[index*10 + i] = 0.;
        sup_coeff[index*10 + i] = 0.;
    }

    inf_coeff[index*10 + y] = 1.;
    sup_coeff[index*10 + y] = 1.;

    inf_coeff[index*10 + x] = -1.;
    sup_coeff[index*10 + x] = -1.;
}


__global__
void coeffs_from_previous_layer(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, float_type* __restrict__ res_inf_coeff, float_type* __restrict__ res_sup_coeff, const float_type* __restrict__ aux_coeffs, const int num_out_neurons_current_layer, const int num_in_neurons_current_layer)
{
    const int n = blockIdx.x;

    int j = threadIdx.x;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(j < num_in_neurons_current_layer)
    {
        const int b = n*num_in_neurons_current_layer + j;

        float_type inf_coeff = 0;
        float_type sup_coeff = 0;

        for(int i = 0; i < num_out_neurons_current_layer; i++)
        {
            const int a = n*num_out_neurons_current_layer + i;
            const int c = i*num_in_neurons_current_layer + j;

            const float_type prev_inf_coeff = expr_inf_coeff[a];
            const float_type prev_sup_coeff = expr_sup_coeff[a];

            if((prev_inf_coeff != 0) || (prev_sup_coeff != 0))
            {
                affine_trans_coeff(inf_coeff, sup_coeff, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, aux_coeffs[c]);
            }
        }

        res_inf_coeff[b] = inf_coeff;
        res_sup_coeff[b] = sup_coeff;

        j += blockDim.x;
    }
}


__global__
void coeffs_from_previous_layer_conv(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, float_type* __restrict__ res_inf_coeff, float_type* __restrict__ res_sup_coeff, const float_type* __restrict__ aux_coeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int input_size_x, const int input_size_y, const int input_size_z, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int n = blockIdx.x;

    const int inp_z = threadIdx.x;
    const int y_shift = threadIdx.y;
    const int x_shift = threadIdx.z;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    if(x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z < filter_size_x*filter_size_y*input_size_z)
    {
        for(int out_x = 0; out_x < output_size_x; out_x++)
        {
            for(int out_y = 0; out_y < output_size_y; out_y++)
            {
                const int x_val = out_x*stride_x + x_shift - pad_x;
                const int y_val = out_y*stride_y + y_shift - pad_y;

                if(!((y_val < 0) || (y_val >= input_size_y)))
                {
                    if(!((x_val < 0) || (x_val >= input_size_x)))
                    {
                        const int mat_in = x_val*input_size_y*input_size_z + y_val*input_size_z + inp_z;
                        const int b = n*input_size_x*input_size_y*input_size_z + mat_in;

                        float_type inf_coeff = res_inf_coeff[b];
                        float_type sup_coeff = res_sup_coeff[b];

                        for(int out_z = 0; out_z < output_size_z; out_z++)
                        {
                            const int mat_out = out_x*output_size_y*output_size_z + out_y*output_size_z + out_z;

                            const int a = n*output_size_x*output_size_y*output_size_z + mat_out;

                            const float_type prev_inf_coeff = expr_inf_coeff[a];
                            const float_type prev_sup_coeff = expr_sup_coeff[a];

                            if((prev_inf_coeff != 0) || (prev_sup_coeff != 0))
                            {
                                const int filter_index = out_z*filter_size_x*filter_size_y*input_size_z + x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z;

                                const float_type aux_coeff = aux_coeffs[filter_index];

                                affine_trans_coeff(inf_coeff, sup_coeff, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, aux_coeff);
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


__global__
void coeffs_from_previous_layer_conv_x_filter_serial(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, float_type* __restrict__ res_inf_coeff, float_type* __restrict__ res_sup_coeff, const float_type* __restrict__ aux_coeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int input_size_x, const int input_size_y, const int input_size_z, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int n = blockIdx.x;

    const int inp_z = threadIdx.x;
    const int y_shift = threadIdx.y;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    if(y_shift*input_size_z + inp_z < filter_size_y*input_size_z)
    {
        for(int out_x = 0; out_x < output_size_x; out_x++)
        {
            for(int out_y = 0; out_y < output_size_y; out_y++)
            {
                for(int x_shift = 0; x_shift < filter_size_x; x_shift++)
                {
                    const int x_val = out_x*stride_x + x_shift - pad_x;
                    const int y_val = out_y*stride_y + y_shift - pad_y;

                    if(!((y_val < 0) || (y_val >= input_size_y)))
                    {
                        if(!((x_val < 0) || (x_val >= input_size_x)))
                        {
                            const int mat_in = x_val*input_size_y*input_size_z + y_val*input_size_z + inp_z;
                            const int b = n*input_size_x*input_size_y*input_size_z + mat_in;

                            float_type inf_coeff = res_inf_coeff[b];
                            float_type sup_coeff = res_sup_coeff[b];

                            for(int out_z = 0; out_z < output_size_z; out_z++)
                            {
                                const int mat_out = out_x*output_size_y*output_size_z + out_y*output_size_z + out_z;

                                const int a = n*output_size_x*output_size_y*output_size_z + mat_out;

                                const float_type prev_inf_coeff = expr_inf_coeff[a];
                                const float_type prev_sup_coeff = expr_sup_coeff[a];

                                if((prev_inf_coeff != 0) || (prev_sup_coeff != 0))
                                {
                                    const int filter_index = out_z*filter_size_x*filter_size_y*input_size_z + x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z;

                                    const float_type aux_coeff = aux_coeffs[filter_index];

                                    affine_trans_coeff(inf_coeff, sup_coeff, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, aux_coeff);
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
}


__global__
void coeffs_from_previous_layer_conv_filter_serial(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, float_type* __restrict__ res_inf_coeff, float_type* __restrict__ res_sup_coeff, const float_type* __restrict__ aux_coeffs, const int output_size_x, const int output_size_y, const int output_size_z, const int input_size_x, const int input_size_y, const int input_size_z, const int filter_size_x, const int filter_size_y, const int stride_x, const int stride_y, const int pad_x, const int pad_y)
{
    const int n = blockIdx.x;

    int inp_z = threadIdx.x;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(inp_z < input_size_z)
    {
        for(int out_x = 0; out_x < output_size_x; out_x++)
        {
            for(int out_y = 0; out_y < output_size_y; out_y++)
            {
                for(int x_shift = 0; x_shift < filter_size_x; x_shift++)
                {
                    for(int y_shift = 0; y_shift < filter_size_y; y_shift++)
                    {
                        const int x_val = out_x*stride_x + x_shift - pad_x;
                        const int y_val = out_y*stride_y + y_shift - pad_y;

                        if(!((y_val < 0) || (y_val >= input_size_y)))
                        {
                            if(!((x_val < 0) || (x_val >= input_size_x)))
                            {
                                const int mat_in = x_val*input_size_y*input_size_z + y_val*input_size_z + inp_z;
                                const int b = n*input_size_x*input_size_y*input_size_z + mat_in;

                                float_type inf_coeff = res_inf_coeff[b];
                                float_type sup_coeff = res_sup_coeff[b];

                                for(int out_z = 0; out_z < output_size_z; out_z++)
                                {
                                    const int mat_out = out_x*output_size_y*output_size_z + out_y*output_size_z + out_z;

                                    const int a = n*output_size_x*output_size_y*output_size_z + mat_out;

                                    const float_type prev_inf_coeff = expr_inf_coeff[a];
                                    const float_type prev_sup_coeff = expr_sup_coeff[a];

                                    if((prev_inf_coeff != 0) || (prev_sup_coeff != 0))
                                    {
                                        const int filter_index = out_z*filter_size_x*filter_size_y*input_size_z + x_shift*filter_size_y*input_size_z + y_shift*input_size_z + inp_z;

                                        const float_type aux_coeff = aux_coeffs[filter_index];

                                        affine_trans_coeff(inf_coeff, sup_coeff, tmp1, tmp2, maxRes, maxMul, prev_inf_coeff, prev_sup_coeff, aux_coeff);
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

        inp_z += blockDim.x;
    }
}


__global__
void csts_from_previous_layer(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, const float_type* __restrict__ expr_inf_cst, const float_type* __restrict__ expr_sup_cst, float_type* __restrict__ res_inf_cst, float_type* __restrict__ res_sup_cst, const float_type* __restrict__ aux_csts, const int num_out_neurons_current_layer)
{
    const int n = blockIdx.x;
    int i = threadIdx.x;

    float_type inf_cst = 0;
    float_type sup_cst = 0;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(i < num_out_neurons_current_layer)
    {
        const int a = n*num_out_neurons_current_layer + i;

        affine_trans_cst(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, expr_inf_coeff[a], expr_sup_coeff[a], aux_csts[i]);

        i += blockDim.x;
    }

    block_reduce_sum_csts(inf_cst, sup_cst, blockDim.x);

    if(threadIdx.x == 0)
    {
        res_inf_cst[n] = add_rd(inf_cst, expr_inf_cst[n]);
        res_sup_cst[n] = add_ru(sup_cst, expr_sup_cst[n]);
    }
}


__global__
void csts_from_previous_layer_conv(const float_type* __restrict__ expr_inf_coeff, const float_type* __restrict__ expr_sup_coeff, const float_type* __restrict__ expr_inf_cst, const float_type* __restrict__ expr_sup_cst, float_type* __restrict__ res_inf_cst, float_type* __restrict__ res_sup_cst, const float_type* __restrict__ aux_csts, const int current_layer_out_size_x, const int current_layer_out_size_y, const int current_layer_out_size_z)
{
    const int n = blockIdx.x;
    int j = threadIdx.x;

    float_type inf_cst = 0;
    float_type sup_cst = 0;

    float_type tmp1, tmp2;
    float_type maxRes, maxMul;

    while(j < current_layer_out_size_z)
    {
        for(int i = 0; i < current_layer_out_size_x*current_layer_out_size_y; i++)
        {
            const int a = n*current_layer_out_size_x*current_layer_out_size_y*current_layer_out_size_z + i*current_layer_out_size_z + j;

            affine_trans_cst(inf_cst, sup_cst, tmp1, tmp2, maxRes, maxMul, expr_inf_coeff[a], expr_sup_coeff[a], aux_csts[j]);
        }

        j += blockDim.x;
    }

    block_reduce_sum_csts(inf_cst, sup_cst, blockDim.x);

    if(threadIdx.x == 0)
    {
        res_inf_cst[n] = add_rd(inf_cst, expr_inf_cst[n]);
        res_sup_cst[n] = add_ru(sup_cst, expr_sup_cst[n]);
    }
}


__global__
void lexpr_replace_relu_bounds(float_type* __restrict__ inf_coeff, float_type* __restrict__ sup_coeff, float_type* __restrict__ inf_cst, float_type* __restrict__ sup_cst, const float_type* __restrict__ lb_array, const float_type* __restrict__ ub_array, const int num_out_neurons_current_layer, const bool use_area_heuristic)
{
    const int n = blockIdx.x;

    int i = threadIdx.x;

    float_type res_inf_cst = 0;
    float_type res_sup_cst = 0;

    while(i < num_out_neurons_current_layer)
    {
        const int a = n*num_out_neurons_current_layer + i;

        lcoeff_replace_relu_bounds(inf_coeff[a], sup_coeff[a], res_inf_cst, res_sup_cst, lb_array[i], ub_array[i], use_area_heuristic);

        i += blockDim.x;
    }

    block_reduce_sum_rd(res_inf_cst, blockDim.x);
    block_reduce_sum_ru(res_sup_cst, blockDim.x);

    if(threadIdx.x == 0)
    {
        inf_cst[n] = add_rd(res_inf_cst, inf_cst[n]);
        sup_cst[n] = add_ru(res_sup_cst, sup_cst[n]);
    }
}


void update_state_using_predecessor_layer_lower_half(fppoly_t* fp, float_type** linf_coeff, float_type** lsup_coeff, float_type** linf_cst, float_type** lsup_cst, float_type** linf_coeff_tmp, float_type** lsup_coeff_tmp, float_type** linf_cst_tmp, float_type** lsup_cst_tmp, const size_t k, const bool use_area_heuristic, const size_t num_out_neurons_last_layer)
{
    const size_t num_out_neurons_current_layer = fp->layers[k]->num_out_neurons;
    const size_t num_in_neurons_current_layer  = fp->layers[k]->num_in_neurons;

    if(fp->layers[k]->activation == RELU)
    {
        lexpr_replace_relu_bounds<<<num_out_neurons_last_layer, num_threads>>>(*linf_coeff, *lsup_coeff, *linf_cst, *lsup_cst, fp->layers[k]->lb_array, fp->layers[k]->ub_array, num_out_neurons_current_layer, use_area_heuristic);
    }

    *linf_coeff_tmp = malloc_device<float_type>(num_out_neurons_last_layer*num_in_neurons_current_layer);
    *lsup_coeff_tmp = malloc_device<float_type>(num_out_neurons_last_layer*num_in_neurons_current_layer);

    cudaMemset(*linf_coeff_tmp, 0, num_out_neurons_last_layer*num_in_neurons_current_layer*sizeof(float_type));
    cudaMemset(*lsup_coeff_tmp, 0, num_out_neurons_last_layer*num_in_neurons_current_layer*sizeof(float_type));

    if(fp->layers[k]->type == CONV)
    {
        if(fp->layers[k]->input_size[0]*fp->layers[k]->input_size[1]*fp->layers[k]->input_size[2] > num_threads)
        {
            if(fp->layers[k]->input_size[1]*fp->layers[k]->input_size[2] > num_threads)
            {
                if(fp->layers[k]->input_size[2] > num_threads)
                {
                    coeffs_from_previous_layer_conv_filter_serial<<<num_out_neurons_last_layer, num_threads>>>(*linf_coeff, *lsup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, fp->layers[k]->filter_weights, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
                }
                else
                {
                    coeffs_from_previous_layer_conv_filter_serial<<<num_out_neurons_last_layer, fp->layers[k]->input_size[2]>>>(*linf_coeff, *lsup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, fp->layers[k]->filter_weights, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
                }
            }
            else
            {
                coeffs_from_previous_layer_conv_x_filter_serial<<<num_out_neurons_last_layer, dim3(fp->layers[k]->input_size[2], fp->layers[k]->filter_size[1], 1)>>>(*linf_coeff, *lsup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, fp->layers[k]->filter_weights, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
            }
        }
        else
        {
            coeffs_from_previous_layer_conv<<<num_out_neurons_last_layer, dim3(fp->layers[k]->input_size[2], fp->layers[k]->filter_size[1], fp->layers[k]->filter_size[0])>>>(*linf_coeff, *lsup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, fp->layers[k]->filter_weights, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->input_size[2], fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);
        }

        csts_from_previous_layer_conv<<<num_out_neurons_last_layer, fp->layers[k]->output_size[2]>>>(*linf_coeff, *lsup_coeff, *linf_cst, *lsup_cst, *linf_cst_tmp, *lsup_cst_tmp, fp->layers[k]->filter_bias, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2]);

    }
    else
    {
        if(fp->layers[k]->type == FFN)
        {
            coeffs_from_previous_layer<<<num_out_neurons_last_layer, num_threads>>>(*linf_coeff, *lsup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, fp->layers[k]->coeffs, num_out_neurons_current_layer, num_in_neurons_current_layer);

            csts_from_previous_layer<<<num_out_neurons_last_layer, num_threads>>>(*linf_coeff, *lsup_coeff, *linf_cst, *lsup_cst, *linf_cst_tmp, *lsup_cst_tmp, fp->layers[k]->csts, num_out_neurons_current_layer);
        }
        else
        {
            lcoeffs_from_previous_layer_pool<<<num_out_neurons_last_layer, fp->layers[k]->output_size[2]>>>(*linf_coeff, *lsup_coeff, *linf_coeff_tmp, *lsup_coeff_tmp, fp->layers[k]->pool_lcoeffs, fp->layers[k]->pool_ucoeffs, fp->layers[k]->output_size[0], fp->layers[k]->output_size[1], fp->layers[k]->output_size[2], fp->layers[k]->input_size[0], fp->layers[k]->input_size[1], fp->layers[k]->filter_size[0], fp->layers[k]->filter_size[1], fp->layers[k]->strides[0], fp->layers[k]->strides[1], fp->layers[k]->pad[0], fp->layers[k]->pad[1]);

            lcsts_from_input_poly<<<num_out_neurons_last_layer, 1>>>(*linf_coeff, *lsup_coeff, *linf_cst, *lsup_cst, *linf_cst_tmp, *lsup_cst_tmp, fp->layers[k]->pool_lcsts, fp->layers[k]->pool_ucsts, fp->layers[k]->lb_array, fp->layers[k]->ub_array, num_out_neurons_current_layer);
        }
    }

    std::swap(*linf_coeff, *linf_coeff_tmp);
    std::swap(*lsup_coeff, *lsup_coeff_tmp);
    std::swap(*linf_cst, *linf_cst_tmp);
    std::swap(*lsup_cst, *lsup_cst_tmp);

    cudaFree(*linf_coeff_tmp);
    cudaFree(*lsup_coeff_tmp);

    *linf_coeff_tmp = nullptr;
    *lsup_coeff_tmp = nullptr;
}


void get_lb_using_previous_layers(elina_manager_t* man, fppoly_t* const fp, const bool use_area_heuristic)
{
    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    const size_t num_out_neurons_last_layer = 90;

    const size_t num_in_neurons_0_layer = fp->layers[0]->num_in_neurons;

    float_type* lb_dev = malloc_device<float_type>(num_out_neurons_last_layer);

    float_type* linf_coeff = malloc_device<float_type>(num_out_neurons_last_layer*10);
    float_type* lsup_coeff = malloc_device<float_type>(num_out_neurons_last_layer*10);

    float_type* linf_cst = malloc_device<float_type>(num_out_neurons_last_layer);
    float_type* lsup_cst = malloc_device<float_type>(num_out_neurons_last_layer);

    size_t index = 0;

    for(elina_dim_t y = 0; y < 10; y++)
    {
        for(elina_dim_t x = 0; x < 10; x++)
        {
            if(y != x)
            {
                create_sub_expr<<<1, 1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, index, y, x);
                index++;
            }
        }
    }

    float_type* linf_coeff_tmp;
    float_type* lsup_coeff_tmp;

    float_type* linf_cst_tmp = malloc_device<float_type>(num_out_neurons_last_layer);
    float_type* lsup_cst_tmp = malloc_device<float_type>(num_out_neurons_last_layer);

    int k = fp->numlayers - 1;

    while(k >= 0)
    {
        if(fp->layers[k]->type == RESIDUAL)
        {
            if(fp->layers[k]->activation==RELU)
            {
                lexpr_replace_relu_bounds<<<num_out_neurons_last_layer, num_threads>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, fp->layers[k]->lb_array, fp->layers[k]->ub_array, fp->layers[k]->num_out_neurons, use_area_heuristic);
            }

            const size_t reslayer_size = fp->layers[k]->num_out_neurons;

            const size_t predecessor1 = fp->layers[k]->predecessors[0] - 1;
            const size_t predecessor2 = fp->layers[k]->predecessors[1] - 1;

            char* predecessor_map = (char*) calloc(k, sizeof(char));
            int iter = predecessor1;

            while(iter >= 0)
            {
                predecessor_map[iter] = 1;
                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            iter =  predecessor2;
            int common_predecessor = 0;

            while(iter >= 0)
            {
                if(predecessor_map[iter] == 1)
                {
                    common_predecessor = iter;

                    break;
                }

                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            free(predecessor_map);

            float_type* copy_linf_coeff = malloc_device<float_type>(num_out_neurons_last_layer*reslayer_size);
            float_type* copy_lsup_coeff = malloc_device<float_type>(num_out_neurons_last_layer*reslayer_size);

            float_type* copy_linf_cst = malloc_device<float_type>(num_out_neurons_last_layer);
            float_type* copy_lsup_cst = malloc_device<float_type>(num_out_neurons_last_layer);

            copy_coeffs_and_csts<<<num_out_neurons_last_layer, 1>>>(copy_linf_coeff, copy_linf_cst, linf_coeff, linf_cst, reslayer_size);
            copy_coeffs_and_csts<<<num_out_neurons_last_layer, 1>>>(copy_lsup_coeff, copy_lsup_cst, lsup_coeff, lsup_cst, reslayer_size);

            cudaMemset(copy_linf_cst, 0, num_out_neurons_last_layer*sizeof(float_type));
            cudaMemset(copy_lsup_cst, 0, num_out_neurons_last_layer*sizeof(float_type));

            float_type* copy_linf_coeff_tmp;
            float_type* copy_lsup_coeff_tmp;

            float_type* copy_linf_cst_tmp = malloc_device<float_type>(num_out_neurons_last_layer);
            float_type* copy_lsup_cst_tmp = malloc_device<float_type>(num_out_neurons_last_layer);

            iter = predecessor1;

            while(iter != common_predecessor)
            {
                update_state_using_predecessor_layer_lower_half(fp, &linf_coeff, &lsup_coeff, &linf_cst, &lsup_cst, &linf_coeff_tmp, &lsup_coeff_tmp, &linf_cst_tmp, &lsup_cst_tmp, iter, use_area_heuristic, num_out_neurons_last_layer);

                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            iter = predecessor2;

            while(iter != common_predecessor)
            {
                update_state_using_predecessor_layer_lower_half(fp, &copy_linf_coeff, &copy_lsup_coeff, &copy_linf_cst, &copy_lsup_cst, &copy_linf_coeff_tmp, &copy_lsup_coeff_tmp, &copy_linf_cst_tmp, &copy_lsup_cst_tmp, iter, use_area_heuristic, num_out_neurons_last_layer);

                iter = fp->layers[iter]->predecessors[0] - 1;
            }

            add_coeffs_and_csts<<<num_out_neurons_last_layer, 1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, copy_linf_coeff, copy_lsup_coeff, copy_linf_cst, copy_lsup_cst, fp->layers[common_predecessor]->num_out_neurons);

            cudaFree(copy_linf_coeff);
            cudaFree(copy_lsup_coeff);
            cudaFree(copy_linf_cst);
            cudaFree(copy_lsup_cst);

            cudaFree(copy_linf_cst_tmp);
            cudaFree(copy_lsup_cst_tmp);

            copy_linf_coeff = nullptr;
            copy_lsup_coeff = nullptr;
            copy_linf_cst = nullptr;
            copy_lsup_cst = nullptr;

            copy_linf_cst_tmp = nullptr;
            copy_lsup_cst_tmp = nullptr;

            k = common_predecessor;
        }
        else
        {
            update_state_using_predecessor_layer_lower_half(fp, &linf_coeff, &lsup_coeff, &linf_cst, &lsup_cst, &linf_coeff_tmp, &lsup_coeff_tmp, &linf_cst_tmp, &lsup_cst_tmp, k, use_area_heuristic, num_out_neurons_last_layer);

            k = fp->layers[k]->predecessors[0] - 1;
        }
    }

    if((fp->input_lweights != nullptr) && (fp->input_uweights != nullptr) && (fp->input_lcst != nullptr) && (fp->input_ucst != nullptr))
    {
        const size_t mu = fp->mu;

        linf_coeff_tmp = malloc_device<float_type>(num_out_neurons_last_layer*mu);
        lsup_coeff_tmp = malloc_device<float_type>(num_out_neurons_last_layer*mu);

        cudaMemset(linf_coeff_tmp, 0, num_out_neurons_last_layer*mu*sizeof(float_type));
        cudaMemset(lsup_coeff_tmp, 0, num_out_neurons_last_layer*mu*sizeof(float_type));

        lcoeffs_from_input_poly<<<num_out_neurons_last_layer, num_threads>>>(linf_coeff, lsup_coeff, linf_coeff_tmp, lsup_coeff_tmp, fp->input_lweights, fp->input_uweights, num_in_neurons_0_layer, mu);

        lcsts_from_input_poly<<<num_out_neurons_last_layer, 1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, linf_cst_tmp, lsup_cst_tmp, fp->input_lcst, fp->input_ucst, fp->input_inf, fp->input_sup, num_in_neurons_0_layer);

        std::swap(linf_coeff, linf_coeff_tmp);
        std::swap(lsup_coeff, lsup_coeff_tmp);
        std::swap(linf_cst, linf_cst_tmp);
        std::swap(lsup_cst, lsup_cst_tmp);

        cudaFree(linf_coeff_tmp);
        cudaFree(lsup_coeff_tmp);

        linf_coeff_tmp = nullptr;
        lsup_coeff_tmp = nullptr;

        compute_lb_from_expr<<<num_out_neurons_last_layer, num_threads>>>(lb_dev, linf_coeff, lsup_coeff, linf_cst, fp->input_inf + num_in_neurons_0_layer, fp->input_sup + num_in_neurons_0_layer, mu);
    }
    else
    {
        compute_lb_from_expr<<<num_out_neurons_last_layer, num_threads>>>(lb_dev, linf_coeff, lsup_coeff, linf_cst, fp->input_inf, fp->input_sup, num_in_neurons_0_layer);
    }

    cudaFree(linf_coeff);
    cudaFree(lsup_coeff);
    cudaFree(linf_cst);
    cudaFree(lsup_cst);

    cudaFree(linf_cst_tmp);
    cudaFree(lsup_cst_tmp);

    linf_coeff = nullptr;
    lsup_coeff = nullptr;
    linf_cst = nullptr;
    lsup_cst = nullptr;

    linf_cst_tmp = nullptr;
    lsup_cst_tmp = nullptr;

    float_type lb[num_out_neurons_last_layer];
    copy_to_host(lb, lb_dev, num_out_neurons_last_layer);

    cudaFree(lb_dev);

    lb_dev = nullptr;

    for(size_t i = 0; i < 10; i++)
    {
        for(size_t j = 0; j < 9; j++)
        {
            if(lb[9*i + j] > 0)
            {
                results[i][j] = true;
            }
            else
            {
                results[i][j] = false;
            }
        }
    }
}


bool is_greater(elina_manager_t* man, elina_abstract0_t* element, const elina_dim_t y, const elina_dim_t x, const bool use_area_heuristic)
{
    fppoly_t* fp = fppoly_of_abstract0(element);
    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    if(!results_calculated)
    {
        get_lb_using_previous_layers(man, fp, use_area_heuristic);
        results_calculated = true;

        if(results[0][0] == true)
        {
            output_counter_x++;
        }
        else
        {
            output_counter_y++;
        }

        return results[0][0];
    }
    else
    {
        bool result = results[output_counter_y][output_counter_x];

        if(result == true)
        {
            output_counter_x++;
        }
        else
        {
            output_counter_y++;
            output_counter_x = 0;
        }

        return result;
    }
}


void conv_add_layer(fppoly_t* const fp, const size_t num_out_neurons, const size_t num_in_neurons, const size_t num_nonzero_weights, const size_t num_biases, const layertype_t type, const activation_type_t activation)
{
    layer_t* layer = (layer_t*) malloc(sizeof(layer_t));

    layer->num_out_neurons = num_out_neurons;
    layer->num_in_neurons = num_in_neurons;

    layer->type = type;
    layer->activation = activation;

    layer->lb_array = malloc_device<float_type>(num_out_neurons);
    layer->ub_array = malloc_device<float_type>(num_out_neurons);

    layer->coeffs = nullptr;
    layer->csts = nullptr;

    layer->filter_weights = malloc_device<float_type>(num_nonzero_weights);
    layer->filter_bias = malloc_device<float_type>(num_biases);
    cudaMemset(layer->filter_bias, 0, num_biases*sizeof(float_type));

    layer->input_size = (size_t*) malloc(3*sizeof(size_t));
    layer->output_size = (size_t*) malloc(3*sizeof(size_t));
    layer->filter_size = (size_t*) malloc(2*sizeof(size_t));
    layer->strides = (size_t*) malloc(2*sizeof(size_t));
    layer->pad = (long int*) malloc(2*sizeof(long int));

    layer->pool_lcoeffs = nullptr;
    layer->pool_ucoeffs = nullptr;

    layer->pool_lcsts = nullptr;
    layer->pool_ucsts = nullptr;

    fp->layers[fp->numlayers] = layer;

    fp->numlayers++;
}


void layer_create_sparse_exprs(fppoly_t* const fp, const float_type* filter_weights, const float_type* filter_bias, const size_t* input_size, const size_t* filter_size, const size_t* strides, const size_t* output_size, const size_t pad_top, const size_t pad_left, const bool has_bias, const activation_type_t activation)
{
    const size_t num_pixels = input_size[0]*input_size[1]*input_size[2];

    const size_t num_out_neurons = output_size[0]*output_size[1]*output_size[2];
    const size_t size = filter_size[0]*filter_size[1]*input_size[2]*output_size[2];

    conv_add_layer(fp, num_out_neurons, num_pixels, size, output_size[2], CONV, activation);

    layer_t* current_layer = fp->layers[fp->numlayers - 1];

    const long int pad[2] = {(long int) pad_top, (long int) pad_left};

    float_type* filter_weights_tmp = (float_type*) malloc(size*sizeof(float_type));
    float_type* filter_bias_tmp = (float_type*) malloc(size*sizeof(float_type));

    for(size_t out_z = 0; out_z < output_size[2]; out_z++)
    {
        for(size_t inp_z = 0; inp_z < input_size[2]; inp_z++)
        {
            for(size_t x_shift = 0; x_shift < filter_size[0]; x_shift++)
            {
                for(size_t y_shift = 0; y_shift < filter_size[1]; y_shift++)
                {
                    const size_t read_index = x_shift*filter_size[1]*input_size[2]*output_size[2] + y_shift*input_size[2]*output_size[2] + inp_z*output_size[2] + out_z;
                    const size_t write_index = out_z*filter_size[0]*filter_size[1]*input_size[2] + x_shift*filter_size[1]*input_size[2] + y_shift*input_size[2] + inp_z;
                    filter_weights_tmp[write_index] = filter_weights[read_index];
                }
            }
        }

        filter_bias_tmp[out_z] = filter_bias[out_z];
    }

    copy_to_device(current_layer->filter_weights, filter_weights_tmp, size);

    if(has_bias)
    {
        copy_to_device(current_layer->filter_bias, filter_bias_tmp, output_size[2]);
    }

    memcpy(current_layer->input_size, input_size, 3*sizeof(size_t));
    memcpy(current_layer->output_size, output_size, 3*sizeof(size_t));
    memcpy(current_layer->filter_size, filter_size, 2*sizeof(size_t));
    memcpy(current_layer->strides, strides, 2*sizeof(size_t));
    memcpy(current_layer->pad, pad, 2*sizeof(long int));

    free(filter_weights_tmp);
    free(filter_bias_tmp);
}


void conv_handle_first_layer(elina_manager_t* man, elina_abstract0_t* element, const float_type* filter_weights, const float_type* filter_bias, const size_t* input_size, const size_t* filter_size, const size_t num_filters, const size_t* strides, const size_t* output_size, const size_t pad_top, const size_t pad_left, const bool has_bias, size_t* predecessors)
{
    fppoly_t* const fp = fppoly_of_abstract0(element);
    fp->layers = (layer_t**) malloc(2000*sizeof(layer_t*));

    layer_create_sparse_exprs(fp, filter_weights, filter_bias, input_size, filter_size, strides, output_size, pad_top, pad_left, has_bias, RELU);

    fp->layers[fp->numlayers - 1]->predecessors = predecessors;

    const size_t num_out_neurons_0_layer = fp->layers[0]->num_out_neurons;

    const long int offset_x = -fp->layers[0]->pad[0];
    const long int offset_y = -fp->layers[0]->pad[1];

    const long int length_x = fp->layers[0]->filter_size[0];
    const long int length_y = fp->layers[0]->filter_size[1];

    const long int shift_x = fp->layers[0]->strides[0];
    const long int shift_y = fp->layers[0]->strides[1];

    float_type* coeffs = malloc_device<float_type>(num_out_neurons_0_layer*length_x*length_y*fp->layers[0]->input_size[2]);
    float_type* csts = malloc_device<float_type>(num_out_neurons_0_layer);

    cudaMemset(coeffs, 0, num_out_neurons_0_layer*length_x*length_y*fp->layers[0]->input_size[2]*sizeof(float_type));
    cudaMemset(csts, 0, num_out_neurons_0_layer*sizeof(float_type));

    device_layer_create_sparse_exprs<<<dim3(fp->layers[0]->output_size[0], fp->layers[0]->output_size[1], fp->layers[0]->output_size[2]), 1>>>(coeffs, csts, fp->layers[0]->filter_weights, fp->layers[0]->filter_bias, 0, fp->layers[0]->input_size[0], fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], fp->layers[0]->filter_size[0], fp->layers[0]->filter_size[1], fp->layers[0]->strides[0], fp->layers[0]->strides[1], fp->layers[0]->pad[0], fp->layers[0]->pad[1]);

    if((fp->input_lweights != nullptr) && (fp->input_uweights != nullptr) && (fp->input_lcst != nullptr) && (fp->input_ucst != nullptr))
    {
        const size_t num_in_neurons_0_layer = fp->layers[0]->num_in_neurons;

        const size_t mu = fp->mu;

        float_type* linf_coeff = malloc_device<float_type>(num_out_neurons_0_layer*mu);
        float_type* lsup_coeff = malloc_device<float_type>(num_out_neurons_0_layer*mu);
        float_type* uinf_coeff = malloc_device<float_type>(num_out_neurons_0_layer*mu);
        float_type* usup_coeff = malloc_device<float_type>(num_out_neurons_0_layer*mu);

        float_type* linf_cst = malloc_device<float_type>(num_out_neurons_0_layer);
        float_type* lsup_cst = malloc_device<float_type>(num_out_neurons_0_layer);
        float_type* uinf_cst = malloc_device<float_type>(num_out_neurons_0_layer);
        float_type* usup_cst = malloc_device<float_type>(num_out_neurons_0_layer);

        cudaMemset(linf_coeff, 0, num_out_neurons_0_layer*mu*sizeof(float_type));
        cudaMemset(lsup_coeff, 0, num_out_neurons_0_layer*mu*sizeof(float_type));
        cudaMemset(uinf_coeff, 0, num_out_neurons_0_layer*mu*sizeof(float_type));
        cudaMemset(usup_coeff, 0, num_out_neurons_0_layer*mu*sizeof(float_type));

        lcoeffs_from_input_poly_sparse<<<dim3(fp->layers[0]->output_size[0], fp->layers[0]->output_size[1], fp->layers[0]->output_size[2]), num_threads>>>(coeffs, coeffs, linf_coeff, lsup_coeff, fp->input_lweights, fp->input_uweights, fp->layers[0]->input_size[0], fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], mu, offset_x, offset_y, length_x, length_y, shift_x, shift_y);
        ucoeffs_from_input_poly_sparse<<<dim3(fp->layers[0]->output_size[0], fp->layers[0]->output_size[1], fp->layers[0]->output_size[2]), num_threads>>>(coeffs, coeffs, uinf_coeff, usup_coeff, fp->input_lweights, fp->input_uweights, fp->layers[0]->input_size[0], fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], mu, offset_x, offset_y, length_x, length_y, shift_x, shift_y);

        lcsts_from_input_poly_sparse<<<dim3(fp->layers[0]->output_size[0], fp->layers[0]->output_size[1], fp->layers[0]->output_size[2]), 1>>>(coeffs, coeffs, csts, csts, linf_cst, lsup_cst, fp->input_lcst, fp->input_ucst, fp->input_inf, fp->input_sup, fp->layers[0]->input_size[0], fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y);
        ucsts_from_input_poly_sparse<<<dim3(fp->layers[0]->output_size[0], fp->layers[0]->output_size[1], fp->layers[0]->output_size[2]), 1>>>(coeffs, coeffs, csts, csts, uinf_cst, usup_cst, fp->input_lcst, fp->input_ucst, fp->input_inf, fp->input_sup, fp->layers[0]->input_size[0], fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y);

        compute_lb_from_expr<<<num_out_neurons_0_layer, num_threads>>>(fp->layers[0]->lb_array, linf_coeff, lsup_coeff, linf_cst, fp->input_inf + num_in_neurons_0_layer, fp->input_sup + num_in_neurons_0_layer, mu);
        compute_ub_from_expr<<<num_out_neurons_0_layer, num_threads>>>(fp->layers[0]->ub_array, uinf_coeff, usup_coeff, usup_cst, fp->input_inf + num_in_neurons_0_layer, fp->input_sup + num_in_neurons_0_layer, mu);

        cudaFree(linf_coeff);
        cudaFree(lsup_coeff);
        cudaFree(uinf_coeff);
        cudaFree(usup_coeff);

        linf_coeff = nullptr;
        lsup_coeff = nullptr;
        uinf_coeff = nullptr;
        usup_coeff = nullptr;
    }
    else
    {
        compute_lb_from_expr_conv_sparse<<<dim3(fp->layers[0]->output_size[0], fp->layers[0]->output_size[1], fp->layers[0]->output_size[2]), fp->layers[0]->input_size[2]>>>(fp->layers[0]->lb_array, coeffs, coeffs, csts, fp->input_inf, fp->input_sup, 1, 0, fp->layers[0]->input_size[0], fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y);
        compute_ub_from_expr_conv_sparse<<<dim3(fp->layers[0]->output_size[0], fp->layers[0]->output_size[1], fp->layers[0]->output_size[2]), fp->layers[0]->input_size[2]>>>(fp->layers[0]->ub_array, coeffs, coeffs, csts, fp->input_inf, fp->input_sup, 1, 0, fp->layers[0]->input_size[0], fp->layers[0]->input_size[1], fp->layers[0]->input_size[2], offset_x, offset_y, length_x, length_y, shift_x, shift_y);
    }

    cudaFree(coeffs);
    cudaFree(csts);

    coeffs = nullptr;
    csts = nullptr;
}


void conv_handle_intermediate_layer(elina_manager_t* man, elina_abstract0_t* element, const float_type* filter_weights, const float_type* filter_bias, const size_t* input_size, const size_t* filter_size, const size_t num_filters, const size_t* strides, const size_t* output_size, const size_t pad_top, const size_t pad_left, const bool has_bias, size_t* predecessors, const activation_type_t activation, const bool use_area_heuristic, const bool retain_training_data, const float_type* gradient)
{
    fppoly_t* const fp = fppoly_of_abstract0(element);

    layer_create_sparse_exprs(fp, filter_weights, filter_bias, input_size, filter_size, strides, output_size, pad_top, pad_left, has_bias, activation);

    fp->layers[fp->numlayers - 1]->predecessors = predecessors;

    update_state_using_previous_layers_sparse_full(man, fp, fp->numlayers - 1, use_area_heuristic, retain_training_data, gradient);
}


void conv_handle_intermediate_relu_layer(elina_manager_t* man, elina_abstract0_t* element, const float_type* filter_weights, const float_type* filter_bias, const size_t* input_size, const size_t* filter_size, const size_t num_filters, const size_t* strides, const size_t* output_size, const size_t pad_top, const size_t pad_left, const bool has_bias, size_t* predecessors, const bool use_area_heuristic, const bool retain_training_data, const float_type* gradient)
{
	conv_handle_intermediate_layer(man, element, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, output_size, pad_top, pad_left, has_bias, predecessors, RELU, use_area_heuristic, retain_training_data, gradient);
}


void conv_handle_intermediate_affine_layer(elina_manager_t* man, elina_abstract0_t* element, const float_type* filter_weights, const float_type* filter_bias, const size_t* input_size, const size_t* filter_size, const size_t num_filters, const size_t* strides, const size_t* output_size, const size_t pad_top, const size_t pad_left, const bool has_bias, size_t* predecessors, const bool use_area_heuristic, const bool retain_training_data, const float_type* gradient)
{
	conv_handle_intermediate_layer(man, element, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, output_size, pad_top, pad_left, has_bias, predecessors, NONE, use_area_heuristic, retain_training_data, gradient);
}


void pool_add_layer(fppoly_t* const fp, const size_t num_out_neurons, const size_t num_in_neurons, const layertype_t type, const activation_type_t activation)
{
    layer_t* layer = (layer_t*) malloc(sizeof(layer_t));

    layer->num_out_neurons = num_out_neurons;
    layer->num_in_neurons = num_in_neurons;

    layer->type = type;
    layer->activation = activation;

    layer->lb_array = malloc_device<float_type>(num_out_neurons);
    layer->ub_array = malloc_device<float_type>(num_out_neurons);

    layer->coeffs = nullptr;
    layer->csts = nullptr;

    layer->filter_weights = nullptr;
    layer->filter_bias = nullptr;

    layer->input_size = (size_t*) malloc(3*sizeof(size_t));
    layer->output_size = (size_t*) malloc(3*sizeof(size_t));
    layer->filter_size = (size_t*) malloc(2*sizeof(size_t));
    layer->strides = (size_t*) malloc(2*sizeof(size_t));
    layer->pad = (long int*) malloc(2*sizeof(long int));

    layer->pool_lcoeffs = nullptr;
    layer->pool_ucoeffs = nullptr;

    layer->pool_lcsts = malloc_device<float_type>(num_out_neurons);
    layer->pool_ucsts = malloc_device<float_type>(num_out_neurons);

    fp->layers[fp->numlayers] = layer;

    fp->numlayers++;
}


size_t handle_pool_layer(elina_manager_t* man, elina_abstract0_t* element, const size_t* pool_size, const size_t* input_size, const size_t* strides, const size_t dimensionality, const size_t pad_top, const size_t pad_left, const size_t* output_size, size_t* predecessors, const bool is_maxpool)
{
    assert(dimensionality == 3);
    assert(pool_size[2] == 1);

    size_t num_in_neurons = input_size[0]*input_size[1]*input_size[2];
    size_t num_out_neurons = output_size[0]*output_size[1]*output_size[2];

    fppoly_t* fp = fppoly_of_abstract0(element);
    pool_add_layer(fp, num_out_neurons, num_in_neurons, POOL, NONE);

    fp->layers[fp->numlayers - 1]->predecessors = predecessors;

    for(size_t i = 0; i < 3; i++)
    {
        fp->layers[fp->numlayers - 1]->output_size[i] = output_size[i];
        fp->layers[fp->numlayers - 1]->input_size[i] = input_size[i];
    }

    for(size_t i = 0; i < 2; i++)
    {
        fp->layers[fp->numlayers - 1]->filter_size[i] = pool_size[i];
        fp->layers[fp->numlayers - 1]->strides[i] = strides[i];
    }

    fp->layers[fp->numlayers - 1]->pad[0] = pad_top;
    fp->layers[fp->numlayers - 1]->pad[1] = pad_left;

    fp->layers[fp->numlayers - 1]->pool_lcoeffs = malloc_device<float_type>(num_out_neurons*pool_size[0]*pool_size[1]);
    fp->layers[fp->numlayers - 1]->pool_ucoeffs = malloc_device<float_type>(num_out_neurons*pool_size[0]*pool_size[1]);

    float_type* lb_in_array = (float_type*) malloc(num_in_neurons*sizeof(float_type));
    float_type* ub_in_array = (float_type*) malloc(num_in_neurons*sizeof(float_type));

    float_type* lb_out_array = (float_type*) malloc(num_out_neurons*sizeof(float_type));
    float_type* ub_out_array = (float_type*) malloc(num_out_neurons*sizeof(float_type));

    copy_to_host(lb_in_array, fp->layers[fp->numlayers - 2]->lb_array, num_in_neurons);
    copy_to_host(ub_in_array, fp->layers[fp->numlayers - 2]->ub_array, num_in_neurons);

    float_type* pool_lcoeffs_host = (float_type*) calloc(num_out_neurons*pool_size[0]*pool_size[1], sizeof(float_type));
    float_type* pool_ucoeffs_host = (float_type*) calloc(num_out_neurons*pool_size[0]*pool_size[1], sizeof(float_type));

    float_type* pool_lcsts_host = (float_type*) malloc(num_out_neurons*sizeof(float_type));
    float_type* pool_ucsts_host = (float_type*) malloc(num_out_neurons*sizeof(float_type));

    float_type* inf = (float_type*) calloc(pool_size[0]*pool_size[1], sizeof(float_type));
    float_type* sup = (float_type*) calloc(pool_size[0]*pool_size[1], sizeof(float_type));

    for(size_t out_x = 0; out_x < output_size[0]; out_x++)
    {
        for(size_t out_y = 0; out_y < output_size[1]; out_y++)
        {
            for(size_t out_z = 0; out_z < output_size[2]; out_z++)
            {
                size_t out_pos = out_x*output_size[1]*output_size[2] + out_y*output_size[2] + out_z;

                size_t inp_x = out_x*pool_size[0];
                size_t inp_y = out_y*pool_size[1];
                size_t inp_z = out_z;
                size_t inp_pos = inp_x*input_size[1]*input_size[2] + inp_y*input_size[2] + inp_z;

                size_t counter = 0;

                float_type max_u = -INFINITY;
                float_type max_l = -INFINITY;

                size_t max_l_var = 0;
                size_t max_u_var = 0;

                for(size_t x_shift = 0; x_shift < pool_size[0]; x_shift++)
                {
                    for(size_t y_shift = 0; y_shift < pool_size[1]; y_shift++)
                    {
                        const size_t l = x_shift*pool_size[1] + y_shift;

                        const int x_val = out_x*strides[0] + x_shift - pad_top;
                        const int y_val = out_y*strides[1] + y_shift - pad_left;

                        if((x_val < 0) || (x_val >= input_size[0]))
                        {
                            // Set inconsistent inf/sup so convex hull is not affected
                            inf[l] =  1.0;
                            sup[l] = -1.0;

                            continue;
                        }

                        if((y_val < 0) || (y_val >= input_size[1]))
                        {
                            // Set inconsistent inf/sup so convex hull is not affected
                            inf[l] =  1.0;
                            sup[l] = -1.0;

                            continue;
                        }

                        size_t pool_cur_dim = inp_pos + x_shift*input_size[1]*input_size[2] + y_shift*input_size[2];

                        if(pool_cur_dim >= num_in_neurons)
                        {
                            continue;
                        }

                        float_type lb = lb_in_array[pool_cur_dim];
                        float_type ub = ub_in_array[pool_cur_dim];

                        if(ub <= 0.0)
                        {
                            inf[l] = 0.0;
                            sup[l] = 0.0;
                        }
                        else if(lb > 0.0)
                        {
                            inf[l] = lb;
                            sup[l] = ub;
                        }
                        else
                        {
                            inf[l] = 0.0;
                            sup[l] = ub;
                        }

                        if(sup[l] > max_u)
                        {
                            max_u = sup[l];
                            max_u_var = l;
                        }

                        if(inf[l] > max_l)
                        {
                            max_l = inf[l];
                            max_l_var = l;
                        }

                        counter++;
                    }
                }

                if(is_maxpool)
                {
                    bool flag = false;
                    size_t var = 0;

                    for(size_t j = 0; j < pool_size[0]*pool_size[1]; j++)
                    {
                        bool is_greater = true;

                        for(size_t k = 0;  k < pool_size[0]*pool_size[1]; k++)
                        {
                            if(k == j)
                            {
                                continue;
                            }

                            if((inf[j] > sup[j]) || (inf[k] >= sup[k]))
                            {
                                continue;
                            }
                            if((inf[k] == sup[k]) && (inf[j] >= sup[k]))
                            {
                                continue;
                            }
                            else if((inf[j] == inf[k]) && (sup[j] == sup[k]) && (inf[j] == sup[j]))
                            {
                                continue;
                            }
                            else if(inf[j] <= sup[k])
                            {
                                is_greater = false;

                                break;
                            }
                        }

                        if(is_greater)
                        {
                            flag = true;
                            var = j;

                            break;
                        }
                    }

                    if(flag)
                    {
                        pool_lcoeffs_host[out_x*output_size[1]*pool_size[0]*pool_size[1]*output_size[2] + out_y*pool_size[0]*pool_size[1]*output_size[2] + var*output_size[2] + out_z] = 1.0;
                        pool_ucoeffs_host[out_x*output_size[1]*pool_size[0]*pool_size[1]*output_size[2] + out_y*pool_size[0]*pool_size[1]*output_size[2] + var*output_size[2] + out_z] = 1.0;

                        pool_lcsts_host[out_pos] = 0.0;
                        pool_ucsts_host[out_pos] = 0.0;
                    }
                    else
                    {
                        pool_lcoeffs_host[out_x*output_size[1]*pool_size[0]*pool_size[1]*output_size[2] + out_y*pool_size[0]*pool_size[1]*output_size[2] + max_l_var*output_size[2] + out_z] = 1.0;
                        pool_lcsts_host[out_pos] = 0.0;

                        dd_MatrixPtr M = nullptr;

                        if(counter <= 10)
                        {
                            M = maxpool_deeppoly_approx(inf, sup, pool_size[0]*pool_size[1]);
                        }

                        bool rel_flag = false;
                        size_t rel_index = 0;
                        float_type best_val = INFINITY;

                        if(M != NULL)
                        {
                            for(size_t i = 0; i < M->rowsize; i++)
                            {
                                const float_type Miy = dd_get_d(M->matrix[i][pool_size[0]*pool_size[1] + 1]);

                                if(Miy < 0)
                                {
                                    const float_type div = -Miy;
                                    float_type val = dd_get_d(M->matrix[i][0])/div;
                                    bool rel_cons = false;

                                    for(size_t j = 0; j < pool_size[0]*pool_size[1]; j++)
                                    {
                                        const float_type Mij = dd_get_d(M->matrix[i][j + 1])/div;

                                        if(Mij < 0)
                                        {
                                            rel_cons = true;
                                            val = val + Mij*inf[j];
                                        }
                                        else if (Mij > 0)
                                        {
                                            rel_cons = true;
                                            val = val + Mij*sup[j];
                                        }
                                    }

                                    if(rel_cons && (val < best_val))
                                    {
                                        rel_flag = true;
                                        rel_index = i;
                                        best_val = val;
                                    }
                                }
                            }
                        }
                        if((rel_flag == true) && (fabs(best_val - max_u) < 0.01))
                        {
                            const float_type div = -dd_get_d(M->matrix[rel_index][pool_size[0]*pool_size[1] + 1]);

                            for(size_t j = 0; j < pool_size[0]*pool_size[1]; j++)
                            {
                                const float_type Mij = dd_get_d(M->matrix[rel_index][j + 1]);
                                pool_ucoeffs_host[out_x*output_size[1]*pool_size[0]*pool_size[1]*output_size[2] + out_y*pool_size[0]*pool_size[1]*output_size[2] + j*output_size[2] + out_z] = Mij/div;
                            }

                            pool_ucsts_host[out_pos] = dd_get_d(M->matrix[rel_index][0])/div;
                        }
                        else
                        {
                            pool_ucoeffs_host[out_x*output_size[1]*pool_size[0]*pool_size[1]*output_size[2] + out_y*pool_size[0]*pool_size[1]*output_size[2] + max_u_var*output_size[2] + out_z] = 0.0;
                            pool_ucsts_host[out_pos] = max_u;
                        }

                        dd_FreeMatrix(M);
                    }

                    lb_out_array[out_pos] = max_l;
                    ub_out_array[out_pos] = max_u;
                }
                else
                {
                    float_type avg_l = 0.0;
                    float_type avg_u = 0.0;

                    fesetround(FE_UPWARD);
                    float_type coeff_l = (float_type)1.0/(float_type)(counter);

                    fesetround(FE_DOWNWARD);
                    float_type coeff_u = (float_type)1.0/(float_type)(counter);

                    for(size_t j = 0; j < pool_size[0]*pool_size[1]; j++)
                    {
                        pool_lcoeffs_host[out_x*output_size[1]*pool_size[0]*pool_size[1]*output_size[2] + out_y*pool_size[0]*pool_size[1]*output_size[2] + j*output_size[2] + out_z] = coeff_l;

                        if(inf[j] <= sup[j])
                        {
                            avg_l = avg_l + inf[j]*coeff_l;
                        }
                    }

                    fesetround(FE_UPWARD);
                    for(size_t j = 0; j < pool_size[0]*pool_size[1]; j++)
                    {
                        pool_ucoeffs_host[out_x*output_size[1]*pool_size[0]*pool_size[1]*output_size[2] + out_y*pool_size[0]*pool_size[1]*output_size[2] + j*output_size[2] + out_z] = coeff_u;

                        if(inf[j] <= sup[j])
                        {
                            avg_u = avg_u + sup[j]*coeff_u;
                        }
                    }

                    fesetround(FE_TONEAREST);

                    pool_lcsts_host[out_pos] = 0.0;
                    pool_ucsts_host[out_pos] = 0.0;

                    lb_out_array[out_pos] = avg_l;
                    ub_out_array[out_pos] = avg_u;
                }
            }
        }
    }

    copy_to_device(fp->layers[fp->numlayers - 1]->lb_array, lb_out_array, num_out_neurons);
    copy_to_device(fp->layers[fp->numlayers - 1]->ub_array, ub_out_array, num_out_neurons);

    copy_to_device(fp->layers[fp->numlayers - 1]->pool_lcoeffs, pool_lcoeffs_host, num_out_neurons*pool_size[0]*pool_size[1]);
    copy_to_device(fp->layers[fp->numlayers - 1]->pool_ucoeffs, pool_ucoeffs_host, num_out_neurons*pool_size[0]*pool_size[1]);

    copy_to_device(fp->layers[fp->numlayers - 1]->pool_lcsts, pool_lcsts_host, num_out_neurons);
    copy_to_device(fp->layers[fp->numlayers - 1]->pool_ucsts, pool_ucsts_host, num_out_neurons);

    free(inf);
    free(sup);

    free(lb_in_array);
    free(ub_in_array);

    free(lb_out_array);
    free(ub_out_array);

    free(pool_lcoeffs_host);
    free(pool_ucoeffs_host);

    free(pool_lcsts_host);
    free(pool_ucsts_host);

    return num_out_neurons;
}


void res_add_layer(fppoly_t* const fp, const size_t num_neurons, const layertype_t type, const activation_type_t activation)
{
    layer_t* layer = (layer_t*) malloc(sizeof(layer_t));

    layer->num_out_neurons = num_neurons;
    layer->num_in_neurons = num_neurons;

    layer->type = type;
    layer->activation = activation;

    layer->coeffs = nullptr;
    layer->csts = nullptr;

    layer->lb_array = malloc_device<float_type>(layer->num_out_neurons);
    layer->ub_array = malloc_device<float_type>(layer->num_out_neurons);

    cudaMemset(layer->lb_array, 0, layer->num_out_neurons*sizeof(float_type));
    cudaMemset(layer->ub_array, 0, layer->num_out_neurons*sizeof(float_type));

    layer->filter_weights = nullptr;
    layer->filter_bias = nullptr;

    layer->input_size = (size_t*) malloc(3*sizeof(size_t));
    layer->output_size = (size_t*) malloc(3*sizeof(size_t));
    layer->filter_size = (size_t*) malloc(2*sizeof(size_t));
    layer->strides = (size_t*) malloc(2*sizeof(size_t));
    layer->pad = (long int*) malloc(2*sizeof(long int));

    for(size_t i = 0; i < 3; i++)
    {
        layer->output_size[i] = fp->layers[fp->numlayers - 1]->output_size[i];
        layer->input_size[i] = fp->layers[fp->numlayers - 1]->output_size[i];
    }

    layer->filter_size[0] = 1;
    layer->filter_size[1] = 1;

    layer->strides[0] = 1;
    layer->strides[1] = 1;

    layer->pad[0] = 0;
    layer->pad[1] = 0;

    fp->layers[fp->numlayers] = layer;

    fp->numlayers++;
}


void handle_residual_layer(elina_manager_t* man, elina_abstract0_t* element, const size_t num_neurons, size_t* predecessors, const activation_type_t activation, const bool use_area_heuristic)
{
	fppoly_t* fp = fppoly_of_abstract0(element);
	res_add_layer(fp, num_neurons, RESIDUAL, activation);
	fp->layers[fp->numlayers - 1]->predecessors = predecessors;

    update_state_using_previous_layers_sparse_full(man, fp, fp->numlayers - 1, use_area_heuristic, false, nullptr);
}


void handle_residual_relu_layer(elina_manager_t* man, elina_abstract0_t* element, const size_t num_neurons, size_t* predecessors, const bool use_area_heuristic)
{
	handle_residual_layer(man, element, num_neurons, predecessors, RELU, use_area_heuristic);
}


void handle_residual_affine_layer(elina_manager_t* man, elina_abstract0_t* element, const size_t num_neurons, size_t* predecessors, const bool use_area_heuristic)
{
	handle_residual_layer(man, element, num_neurons, predecessors, NONE, use_area_heuristic);
}


void free_layer(layer_t* layer)
{
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

    cudaFree(layer->pool_lcoeffs);
    cudaFree(layer->pool_ucoeffs);
    cudaFree(layer->pool_lcsts);
    cudaFree(layer->pool_ucsts);

    layer->pool_lcoeffs = nullptr;
    layer->pool_ucoeffs = nullptr;
    layer->pool_lcsts = nullptr;
    layer->pool_ucsts = nullptr;

    free(layer);
    layer = nullptr;
}


void fppoly_free(elina_manager_t* man, fppoly_t* fp)
{
    for(size_t i = 0; i < fp->numlayers; i++)
    {
        free_layer(fp->layers[i]);
    }

    free(fp->layers);
    fp->layers = nullptr;

    cudaFree(fp->input_inf);
    cudaFree(fp->input_sup);
    fp->input_inf = nullptr;
    fp->input_sup = nullptr;

    cudaFree(fp->input_lweights);
    cudaFree(fp->input_uweights);
    fp->input_lweights = nullptr;
    fp->input_uweights = nullptr;

    cudaFree(fp->input_lcst);
    cudaFree(fp->input_ucst);
    fp->input_lcst = nullptr;
    fp->input_ucst = nullptr;

    free(fp);
    fp = nullptr;
}


void layer_print(const layer_t* layer)
{
    //neurons_print<<<1, 1>>>(layer->neurons, layer->num_out_neurons);
}


void fppoly_fprint(FILE* const stream, elina_manager_t* man, const fppoly_t* const fp, const char* const * name_of_dim)
{
    for(size_t i = 0; i < fp->numlayers; i++)
    {
        printf("layer: %zu\n", i);
        layer_print(fp->layers[i]);
    }
}


size_t get_num_neurons_in_layer(elina_manager_t* man, elina_abstract0_t* abs, const size_t layerno)
{
    fppoly_t* fp = fppoly_of_abstract0(abs);

    if(layerno >= fp->numlayers)
    {
	    fprintf(stdout,"the layer does not exist\n");

	    return 0;
    }

    layer_t * layer = fp->layers[layerno];

    return layer->num_out_neurons;
}


elina_interval_t** box_for_layer(elina_manager_t* man, elina_abstract0_t* abs, const size_t layerno)
{
    fppoly_t* fp = fppoly_of_abstract0(abs);

    if(layerno >= fp->numlayers)
    {
        fprintf(stdout,"the layer does not exist\n");

        return NULL;
    }

    layer_t* layer = fp->layers[layerno];
    const size_t num_neurons = layer->num_out_neurons;

    float_type* lb_array = layer->lb_array;
    float_type* ub_array = layer->ub_array;

    float_type* lb_array_host = (float_type*) malloc(num_neurons*sizeof(float_type));
    float_type* ub_array_host = (float_type*) malloc(num_neurons*sizeof(float_type));

    copy_to_host(lb_array_host, lb_array, num_neurons);
    copy_to_host(ub_array_host, ub_array, num_neurons);

    elina_interval_t** itv_arr = (elina_interval_t**) malloc(num_neurons*sizeof(elina_interval_t*));

    for(size_t i = 0; i < num_neurons; i++)
    {
        itv_arr[i] = elina_interval_alloc();
        elina_interval_set_double(itv_arr[i], lb_array_host[i], ub_array_host[i]);
    }

    free(lb_array_host);
    free(ub_array_host);

    return itv_arr;
}


__global__
void update_bounds_for_neuron_dev(float_type* lb_array, float_type* ub_array, const float_type lb, const float_type ub, const size_t neuron_no)
{
    lb_array[neuron_no] = lb;
    ub_array[neuron_no] = ub;
}


void update_bounds_for_neuron(elina_manager_t* man, elina_abstract0_t* abs, const size_t layerno, const size_t neuron_no, const float_type lb, const float_type ub)
{
    fppoly_t* fp = fppoly_of_abstract0(abs);

    if(layerno >= fp->numlayers)
    {
        fprintf(stdout,"the layer does not exist\n");

        return;
    }

    update_bounds_for_neuron_dev<<<1, 1>>>(fp->layers[layerno]->lb_array, fp->layers[layerno]->ub_array, lb, ub, neuron_no);
}


int get_num_neurons_training_layer()
{
    return num_neurons_training_layer;
}


float_type* get_adv()
{
    return adv_host;
}

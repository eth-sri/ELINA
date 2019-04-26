/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
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

#include "fppoly_gpu.h"

#include <iostream>
#include <cuda.h>
#include <chrono>


const size_t num_threads = 128;

bool results[90];
bool results_calculated;
size_t output_counter;

__constant__ const double min_denormal = 4.940656458412465441766e-324;
__constant__ const double ulp =          2.220446049250313080848e-16;


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


__device__
void elina_double_interval_mul(double* const a_inf, double* const a_sup, const double b_inf, const double b_sup, const double c_inf, const double c_sup)
{
    if(c_inf <= 0)
    {
        /* interval c is positive */
        if(b_inf <= 0)
        {
            /*interval b is positive*/
            if((b_inf == 0) || (c_inf == 0))
            {
                *a_inf = 0.0;
            }
            else
            {
                *a_inf = b_inf*-c_inf;
            }

            if((b_sup == 0) || (c_sup == 0))
            {
                *a_sup = 0.0;
            }
            else
            {
                *a_sup = b_sup*c_sup;
            }
        }
        else if(b_sup <= 0)
        {
            /* interval b is negative */
            if((c_sup == 0) || (b_inf == 0))
            {
                *a_inf = 0.0;
            }
            else
            {
                *a_inf = c_sup*b_inf;
            }

            if((c_inf == 0) || (b_sup == 0))
            {
                *a_sup = 0.0;
            }
            else
            {
                *a_sup = -c_inf*b_sup;
            }
        }
        else
        {
            /* there is 0 in between for b */
            if((c_sup == 0) || (b_inf == 0))
            {
                *a_inf = 0.0;
            }
            else
            {
                *a_inf = b_inf*c_sup;
            }

            if((c_sup == 0) || (b_sup == 0))
            {
                *a_sup = 0.0;
            }
            else
            {
                *a_sup = b_sup*c_sup;
            }
        }
    }
    else if(c_sup <= 0)
    {
        /* interval c is negative */
        if(b_inf <= 0)
        {
            /*interval b is positive*/
            if((b_sup == 0) || (c_inf == 0))
            {
                *a_inf = 0.0;
            }
            else
            {
                *a_inf = b_sup*c_inf;
            }

            if((b_inf == 0) || (c_sup == 0))
            {
                *a_sup = 0.0;
            }
            else
            {
                *a_sup = -b_inf*c_sup;
            }
        }
        else if(b_sup <= 0)
        {
            /* interval b is negative */
            if((b_sup == 0) || (c_sup == 0))
            {
                *a_inf = 0.0;
            }
            else
            {
                *a_inf = b_sup*-c_sup;
            }

            if((b_inf == 0) || (c_inf == 0))
            {
                *a_sup = 0.0;
            }
            else
            {
                *a_sup = b_inf*c_inf;
            }
        }
        else
        {
            /* there is 0 in between for b */
            if((c_inf == 0) || (b_sup == 0))
            {
                *a_inf = 0.0;
            }
            else
            {
                *a_inf = b_sup*c_inf;
            }

            if((c_inf == 0) || (b_inf == 0))
            {
                *a_sup = 0.0;
            }
            else
            {
                *a_sup = b_inf*c_inf;
            }
        }
    }
    else if(b_inf <= 0)
    {
        /* interval b is positive */
        if(c_inf <= 0)
        {
            /*interval c is positive */
            if((b_inf == 0) || (c_inf == 0))
            {
                *a_inf = 0.0;
            }
            else
            {
                *a_inf = -b_inf*c_inf;
            }

            if((b_sup == 0) || (c_sup == 0))
            {
                *a_sup = 0.0;
            }
            else
            {
                *a_sup = b_sup*c_sup;
            }
        }
        else if(c_sup <= 0)
        {
            /* interval c is negative */
            if((b_sup == 0) || (c_inf == 0))
            {
                *a_inf = 0.0;
            }
            else
            {
                *a_inf = b_sup*c_inf;
            }

            if((b_inf == 0) || (c_sup == 0))
            {
                *a_sup = 0.0;
            }
            else
            {
                *a_sup = -b_inf*c_sup;
            }
        }
        else
        {
            /* there is 0 in between for c */
            if((b_sup == 0) || (c_inf == 0))
            {
                *a_inf = 0.0;
            }
            else
            {
                *a_inf = b_sup*c_inf;
            }

            if((b_sup == 0) || (c_sup == 0))
            {
                *a_sup = 0.0;
            }
            else
            {
                *a_sup = b_sup*c_sup;
            }
        }
    }
    else if(b_sup <= 0)
    {
        /* interval b is negative */
        if(c_inf <= 0)
        {
            /* interval c is positive */
            if((b_inf == 0) || (c_sup == 0))
            {
                *a_inf = 0.0;
            }
            else
            {
                *a_inf = b_inf*c_sup;
            }

            if((b_sup == 0) || (c_inf == 0))
            {
                *a_sup = 0.0;
            }
            else
            {
                *a_sup = b_sup*-c_inf;
            }
        }
        else if(c_sup <= 0)
        {
            /* interval c is negative */
            if((b_sup == 0) || (c_sup == 0))
            {
                *a_inf = 0.0;
            }
            else
            {
                *a_inf = -b_sup*c_sup;
            }

            if((b_inf == 0) || (c_inf == 0))
            {
                *a_sup = 0.0;
            }
            else
            {
                *a_sup = b_inf*c_inf;
            }
        }
        else
        {
            /* there is 0 in between for c */
            if((b_inf == 0) || (c_sup == 0))
            {
                *a_inf = 0.0;
            }
            else
            {
                *a_inf = b_inf*c_sup;
            }

            if((b_inf == 0) || (c_inf == 0))
            {
                *a_sup = 0.0;
            }
            else
            {
                *a_sup = b_inf*c_inf;
            }
        }
    }
    else
    {
        /* there is 0 in between for both b and c */
        double tmp_inf1 = b_sup*c_inf;
        double tmp_sup1 = b_inf*c_inf;
        double tmp_inf2 = b_inf*c_sup;
        double tmp_sup2 = b_sup*c_sup;
        *a_inf = fmax(tmp_inf1, tmp_inf2);
        *a_sup = fmax(tmp_sup1, tmp_sup2);
    }
}


__device__
void elina_double_interval_div(double* const a_inf, double* const a_sup, const double b_inf, const double b_sup, const double c_inf, const double c_sup)
{
    if (c_inf < 0)
    {
        /* c is positive */
        if (b_inf <= 0)
        {
            /* b is positive */
            *a_inf = b_inf/c_sup;
            *a_sup = b_sup/-c_inf;
        }
        else if (b_sup <= 0)
        {
            /* b is negative */
            *a_inf = -b_inf/c_inf;
            *a_sup = b_sup/c_sup;
        }
        else
        {
            /* 0 is in the middle of b: one divides b by c->inf */
            *a_inf = b_inf/-c_inf;
            *a_sup = b_sup/-c_inf;
        }
    }
    else if (c_sup < 0)
    {
        /* c is negative */
        if (b_inf <= 0)
        {
            /* b is positive */
            *a_sup = b_inf/c_inf;
            *a_inf = -b_sup/c_sup;
        }
        else if (b_sup <= 0)
        {
            /* b is negative */
            *a_inf = b_sup/c_inf;
            *a_sup = -b_inf/c_sup;
        }
        else
        {
            /* 0 is in the middle of b: one cross-divide b by c->sup */
            *a_inf = b_sup/c_sup;
            *a_sup = b_inf/c_sup;
        }
    }
    else if ((b_inf == 0) && (b_sup == 0))
    {
        /* b is [0,0] */
        *a_inf = b_inf;
        *a_sup = b_sup;
    }
    else
    {
        *a_inf = INFINITY;
        *a_sup = INFINITY;
    }
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
    pr->min_denormal = ldexpl(1.0, -1074);
    pr->ulp = ldexpl(1.0, -52);

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


elina_manager_t* fppoly_manager_alloc()
{
    std::cout << "This is the GPU version of fppoly!" << std::endl;
    results_calculated = false;
    output_counter = 1;

    void** funptr;
    //fesetround(FE_UPWARD);
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


/*
__device__
void expr_print(const expr_t* const expr)
{
    if((expr->inf_coeff == nullptr) || (expr->sup_coeff == nullptr))
    {
        printf("+ [%g, %g]\n", -expr->inf_cst, expr->sup_cst);

        return;
    }

    for(size_t i = 0; i < size; i++)
    {
        if(i == 0)
        {
            printf("[%g, %g]x0 ", -expr->inf_coeff[0], expr->sup_coeff[0]);
        }
        else
        {
            printf("+ [%g, %g]x%zu ", -expr->inf_coeff[i], expr->sup_coeff[i], i);
        }
    }

    printf("+ [%g, %g]\n", -expr->inf_cst, expr->sup_cst);
}
*/


layer_t* create_layer(const size_t num_out_neurons, const size_t num_in_neurons, const layertype_t type, const activation_type_t activation)
{
    layer_t* layer = (layer_t*) malloc(sizeof(layer_t));

    layer->num_out_neurons = num_out_neurons;
    layer->num_in_neurons = num_in_neurons;

    layer->type = type;
    layer->activation = activation;

    cudaMalloc((void**) &layer->lb_array, num_out_neurons*sizeof(double));
    cudaMalloc((void**) &layer->ub_array, num_out_neurons*sizeof(double));

    cudaMalloc((void**) &layer->inf_coeff, num_out_neurons*num_in_neurons*sizeof(double));
    cudaMalloc((void**) &layer->sup_coeff, num_out_neurons*num_in_neurons*sizeof(double));

    cudaMalloc((void**) &layer->inf_cst, num_out_neurons*sizeof(double));
    cudaMalloc((void**) &layer->sup_cst, num_out_neurons*sizeof(double));

    return layer;
}


void fppoly_from_network_input_box(fppoly_t* const res, const size_t intdim, const size_t realdim, const double* inf_array, const double* sup_array)
{
    res->layers = nullptr;
    res->numlayers = 0;

    size_t num_pixels = intdim + realdim;

    double* tmp_input_inf = (double*) malloc(num_pixels*sizeof(double));
    double* tmp_input_sup = (double*) malloc(num_pixels*sizeof(double));

    for(size_t i = 0; i < num_pixels; i++)
    {
        tmp_input_inf[i] = -inf_array[i];
        tmp_input_sup[i] = sup_array[i];
    }

    cudaMalloc((void**) &(res->input_inf), num_pixels*sizeof(double));
    cudaMalloc((void**) &(res->input_sup), num_pixels*sizeof(double));

    cudaMemcpy(res->input_inf, tmp_input_inf, num_pixels*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(res->input_sup, tmp_input_sup, num_pixels*sizeof(double), cudaMemcpyHostToDevice);

    free(tmp_input_inf);
    free(tmp_input_sup);

    res->num_pixels = num_pixels;
}


elina_abstract0_t* fppoly_from_network_input(elina_manager_t* man, const size_t intdim, const size_t realdim, const double* inf_array, const double* sup_array)
{
    fppoly_t* res = (fppoly_t*) malloc(sizeof(fppoly_t));
    fppoly_from_network_input_box(res, intdim, realdim, inf_array, sup_array);

    return abstract0_of_fppoly(man, res);
}


void fppoly_add_new_layer(fppoly_t* const fp, const size_t num_out_neurons, const size_t num_in_neurons, const layertype_t type, const activation_type_t activation)
{
    const size_t numlayers = fp->numlayers;
    fp->layers[numlayers] = create_layer(num_out_neurons, num_in_neurons, type, activation);
    fp->numlayers++;
}


__device__
void elina_double_interval_add_expr_coeff(double* const res_inf, double* const res_sup, const double inf, const double sup, const double inf_expr, const double sup_expr)
{
    *res_inf = inf + inf_expr;
    *res_sup = sup + sup_expr;
    const double maxA = fmax(fabs(inf_expr), fabs(sup_expr));
    double tmp1, tmp2;
    elina_double_interval_mul(&tmp1, &tmp2, inf, sup, maxA*ulp, maxA*ulp);
    *res_inf += tmp1;
    *res_sup += tmp2;
}


__device__
void elina_double_interval_add_cst_coeff(double* const res_inf, double* const res_sup, const double inf, const double sup, const double inf_expr, const double sup_expr)
{
    elina_double_interval_add_expr_coeff(res_inf, res_sup, inf, sup, inf_expr, sup_expr);
    *res_inf += min_denormal;
    *res_sup += min_denormal;
}


__device__
void elina_double_interval_mul_expr_coeff(double* const res_inf, double* const res_sup, const double inf, const double sup, const double inf_expr, const double sup_expr)
{
    elina_double_interval_mul(res_inf, res_sup, inf, sup, inf_expr, sup_expr);
    const double maxA = fmax(fabs(inf_expr), fabs(sup_expr));
    double tmp1, tmp2;
    elina_double_interval_mul(&tmp1, &tmp2, inf, sup, maxA*ulp, maxA*ulp);
    *res_inf += tmp1;
    *res_sup += tmp2;
}


__device__
void elina_double_interval_mul_cst_coeff(double* const res_inf, double* const res_sup, const double inf, const double sup, const double inf_expr, const double sup_expr)
{
    elina_double_interval_mul_expr_coeff(res_inf, res_sup, inf, sup, inf_expr, sup_expr);
    *res_inf += min_denormal;
    *res_sup += min_denormal;
}


__global__
void compute_lb_from_expr(double* lb_array, double* inf_coeff, double* sup_coeff, const double* inf_cst, const double* input_inf, const double* input_sup, const size_t num_exprs, const size_t expr_size)
{
    const size_t n = blockIdx.x;

    double res_inf = inf_cst[n];

    double tmp1, tmp2;

    for(size_t i = 0; i < expr_size; i++)
    {
        elina_double_interval_mul(&tmp1, &tmp2, inf_coeff[n*expr_size + i], sup_coeff[n*expr_size + i], input_inf[i], input_sup[i]);
        res_inf = res_inf + tmp1;
    }

    lb_array[n] = res_inf;
}


__global__
void compute_ub_from_expr(double* ub_array, double* inf_coeff, double* sup_coeff, const double* sup_cst, const double* input_inf, const double* input_sup, const size_t num_exprs, const size_t expr_size)
{
    const size_t n = blockIdx.x;

    double res_sup = sup_cst[n];

    double tmp1, tmp2;

    for(size_t i = 0; i < expr_size; i++)
    {
        elina_double_interval_mul(&tmp1, &tmp2, inf_coeff[n*expr_size + i], sup_coeff[n*expr_size + i], input_inf[i], input_sup[i]);
        res_sup = res_sup + tmp2;
    }

   ub_array[n] = res_sup;
}


__global__
void device_layer_create_dense_expr(double* inf_coeff, double* sup_coeff, double* inf_cst, double* sup_cst, const double* weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons)
{
    const size_t i = blockIdx.x;

    const double* weight_i = weights + i*num_in_neurons;
    const double bias_i = bias[i];

    inf_cst[i] = -bias_i;
    sup_cst[i] = bias_i;

    for(size_t j = 0; j < num_in_neurons; j++)
    {
        inf_coeff[i*num_in_neurons + j] = -weight_i[j];
        sup_coeff[i*num_in_neurons + j] = weight_i[j];
    }
}


void layer_create_dense_exprs(double* inf_coeff, double* sup_coeff, double* inf_cst, double* sup_cst, const double** weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons)
{
    double* tmp_weights;
    cudaMalloc((void**) &tmp_weights, num_out_neurons*num_in_neurons*sizeof(double));

    double* tmp_bias;
    cudaMalloc((void**) &tmp_bias, num_out_neurons*sizeof(double));

    for(size_t i = 0; i < num_out_neurons; i++)
    {
        cudaMemcpy(tmp_weights + i*num_in_neurons, weights[i], num_in_neurons*sizeof(double), cudaMemcpyHostToDevice);
    }

    cudaMemcpy(tmp_bias, bias, num_out_neurons*sizeof(double), cudaMemcpyHostToDevice);

    device_layer_create_dense_expr<<<num_out_neurons, 1>>>(inf_coeff, sup_coeff, inf_cst, sup_cst, tmp_weights, tmp_bias, num_out_neurons, num_in_neurons);

    cudaFree(tmp_weights);
    cudaFree(tmp_bias);
}


__global__
void copy_expr_array(double* target_inf_coeff, double* target_sup_coeff, double* target_inf_cst, double* target_sup_cst, double* source_inf_coeff, double* source_sup_coeff, double* source_inf_cst, double* source_sup_cst, const size_t num_exprs, const size_t expr_size)
{
    const size_t i = blockIdx.x;

    for(size_t j = 0; j < expr_size; j++)
    {
        target_inf_coeff[i*expr_size + j] = source_inf_coeff[i*expr_size + j];
        target_sup_coeff[i*expr_size + j] = source_sup_coeff[i*expr_size + j];
    }

    target_inf_cst[i] = source_inf_cst[i];
    target_sup_cst[i] = source_sup_cst[i];
}


void layer_compute_bounds_from_exprs(double* inf_coeff, double* sup_coeff, double* inf_cst, double* sup_cst, double* lb_array, double* ub_array, double* input_inf, double* input_sup, const size_t num_out_neurons, const size_t num_in_neurons)
{
    compute_lb_from_expr<<<num_out_neurons, 1>>>(lb_array, inf_coeff, sup_coeff, inf_cst, input_inf, input_sup, num_out_neurons, num_in_neurons);
    compute_ub_from_expr<<<num_out_neurons, 1>>>(ub_array, inf_coeff, sup_coeff, sup_cst, input_inf, input_sup, num_out_neurons, num_in_neurons);
}


void ffn_handle_first_layer(elina_manager_t* man, elina_abstract0_t* abs, const double** weights, const double* bias, const size_t size, const size_t num_pixels, const activation_type_t activation)
{
    fppoly_t* res = fppoly_of_abstract0(abs);
    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    res->layers = (layer_t**) malloc(20*sizeof(layer_t*));
    fppoly_add_new_layer(res, size, num_pixels, FFN, activation);

    double* inf_coeff = res->layers[0]->inf_coeff;
    double* sup_coeff = res->layers[0]->sup_coeff;

    double* inf_cst = res->layers[0]->inf_cst;
    double* sup_cst = res->layers[0]->sup_cst;

    layer_create_dense_exprs(inf_coeff, sup_coeff, inf_cst, sup_cst, weights, bias, size, num_pixels);
    layer_compute_bounds_from_exprs(inf_coeff, sup_coeff, inf_cst, sup_cst, res->layers[0]->lb_array, res->layers[0]->ub_array, res->input_inf, res->input_sup, res->layers[0]->num_out_neurons, res->layers[0]->num_in_neurons);
}


void ffn_handle_first_relu_layer(elina_manager_t* man, elina_abstract0_t* abs, const double** weights, const double* bias, const size_t size, const size_t num_pixels)
{
        ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels, RELU);
}


void ffn_handle_first_sigmoid_layer(elina_manager_t* man, elina_abstract0_t* abs, const double**weights, const double* bias, const size_t size, const size_t num_pixels)
{
        //ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels, SIGMOID);
}


void ffn_handle_first_tanh_layer(elina_manager_t* man, elina_abstract0_t* abs, const double** weights, const double* bias, const size_t size, const size_t num_pixels)
{
        //ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels, TANH);
}


__global__
void lexpr_replace_relu_bounds(double* inf_coeff, double* sup_coeff, double* inf_cst, double* sup_cst, double* lb_array, double* ub_array, const size_t num_out_neurons_last_layer, const size_t num_out_neurons_current_layer)
{
    const size_t n = blockIdx.x;
    const size_t i = blockIdx.y*blockDim.x + threadIdx.x;

    if(i < num_out_neurons_current_layer)
    {
        const size_t a = n*num_out_neurons_current_layer + i;

        const double lb = lb_array[i];
        const double ub = ub_array[i];
        const double width = ub + lb;
        const double lambda_inf = -ub/width;
        const double lambda_sup = ub/width;

        const double old_inf_coeff = inf_coeff[a];
        const double old_sup_coeff = sup_coeff[a];

        if((old_sup_coeff == 0) && (old_inf_coeff == 0))
        {
            inf_coeff[a] = 0.0;
            sup_coeff[a] = 0.0;

            return;
        }
        else if(ub <= 0)
        {
            inf_coeff[a] = 0.0;
            sup_coeff[a] = 0.0;

            return;
        }
        else if(lb < 0)
        {
            inf_coeff[a] = old_inf_coeff;
            sup_coeff[a] = old_sup_coeff;
        }
        else if(old_sup_coeff < 0)
        {
            const double mu_inf = lambda_inf*lb;
            const double mu_sup = lambda_sup*lb;
            elina_double_interval_mul_expr_coeff(&inf_coeff[a], &sup_coeff[a], lambda_inf, lambda_sup, old_inf_coeff, old_sup_coeff);
            double tmp1, tmp2;
            elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, mu_inf, mu_sup, old_inf_coeff, old_sup_coeff);

            atomicAdd(&inf_cst[n], tmp1 + min_denormal);
            atomicAdd(&sup_cst[n], tmp2 + min_denormal);
        }
        else if (old_inf_coeff < 0)
        {
            const double area1 = lb*ub;
            const double area2 = 0.5*ub*width;
            const double area3 = 0.5*lb*width;

            if((area1 < area2) && (area1 < area3))
            {
                elina_double_interval_mul_expr_coeff(&inf_coeff[a], &sup_coeff[a], lambda_inf, lambda_sup, old_inf_coeff, old_sup_coeff);
            }
            else if((area2 < area1) && (area2 < area3))
            {
                inf_coeff[a] = 0.0;
                sup_coeff[a] = 0.0;
            }
            else
            {
                inf_coeff[a] = old_inf_coeff;
                sup_coeff[a] = old_sup_coeff;
            }
        }
        else
        {
            inf_coeff[a] = 0.0;
            sup_coeff[a] = 0.0;
            double tmp1, tmp2;
            elina_double_interval_mul(&tmp1, &tmp2, old_inf_coeff, old_sup_coeff, 0, ub);

            atomicAdd(&inf_cst[n], tmp1);
            atomicAdd(&sup_cst[n], -tmp1);
        }
    }
}


__global__
void uexpr_replace_relu_bounds(double* inf_coeff, double* sup_coeff, double* inf_cst, double* sup_cst, double* lb_array, double* ub_array, const size_t num_out_neurons_last_layer, const size_t num_out_neurons_current_layer)
{
    const size_t n = blockIdx.x;
    const size_t i = blockIdx.y*blockDim.x + threadIdx.x;

    if(i < num_out_neurons_current_layer)
    {
        const size_t a = n*num_out_neurons_current_layer + i;

        const double lb = lb_array[i];
        const double ub = ub_array[i];
        const double width = ub + lb;
        const double lambda_inf = -ub/width;
        const double lambda_sup = ub/width;

        const double old_inf_coeff = inf_coeff[a];
        const double old_sup_coeff = sup_coeff[a];

        if((old_sup_coeff == 0) && (old_inf_coeff == 0))
        {
            inf_coeff[a] = 0.0;
            sup_coeff[a] = 0.0;

            return;
        }
        else if(ub <= 0)
        {
            inf_coeff[a] = 0.0;
            sup_coeff[a] = 0.0;

            return;
        }
        else if(lb < 0)
        {
            inf_coeff[a] = old_inf_coeff;
            sup_coeff[a] = old_sup_coeff;
        }
        else if(old_inf_coeff < 0)
        {
            const double mu_inf = lambda_inf*lb;
            const double mu_sup = lambda_sup*lb;
            elina_double_interval_mul_expr_coeff(&inf_coeff[a], &sup_coeff[a], lambda_inf, lambda_sup, old_inf_coeff, old_sup_coeff);
            double tmp1, tmp2;
            elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, mu_inf, mu_sup, old_inf_coeff, old_sup_coeff);

            atomicAdd(&inf_cst[n], tmp1 + min_denormal);
            atomicAdd(&sup_cst[n], tmp2 + min_denormal);
        }
        else if(old_sup_coeff < 0)
        {
            const double area1 = lb*ub;
            const double area2 = 0.5*ub*width;
            const double area3 = 0.5*lb*width;

            if((area1 < area2) && (area1 < area3))
            {
                elina_double_interval_mul_expr_coeff(&inf_coeff[a], &sup_coeff[a], lambda_inf, lambda_sup, old_inf_coeff, old_sup_coeff);
            }
            else if((area2 < area1) && (area2 < area3))
            {
                inf_coeff[a] = 0.0;
                sup_coeff[a] = 0.0;
            }
            else
            {
                inf_coeff[a] = old_inf_coeff;
                sup_coeff[a] = old_sup_coeff;
            }
        }
        else
        {
            inf_coeff[a] = 0.0;
            sup_coeff[a] = 0.0;
            double tmp1, tmp2;
            elina_double_interval_mul(&tmp1, &tmp2, old_inf_coeff, old_sup_coeff, 0, ub);

            atomicAdd(&inf_cst[n], -tmp2);
            atomicAdd(&sup_cst[n], tmp2);
        }
    }
}


__global__
void coeffs_from_previous_layer(double* expr_inf_coeff, double* expr_sup_coeff, double* res_inf_coeff, double* res_sup_coeff, double* aux_inf_coeff, double* aux_sup_coeff, const size_t num_out_neurons_last_layer, const size_t num_out_neurons_current_layer, const size_t num_in_neurons_current_layer)
{
    const size_t n = blockIdx.x;
    const size_t j = blockIdx.y*blockDim.x + threadIdx.x;

    if(j < num_in_neurons_current_layer)
    {
        size_t i = 0;

        size_t a = n*num_out_neurons_current_layer + i;
        const size_t b = n*num_in_neurons_current_layer + j;
        size_t c = i*num_in_neurons_current_layer + j;

        elina_double_interval_mul_expr_coeff(&res_inf_coeff[b], &res_sup_coeff[b], expr_inf_coeff[a], expr_sup_coeff[a], aux_inf_coeff[c], aux_sup_coeff[c]);

        double tmp1, tmp2;
        double maxRes, maxMul;

        for(i = 1; i < num_out_neurons_current_layer; i++)
        {
            a++;
            c += num_in_neurons_current_layer;

            if((expr_inf_coeff[a] != 0) || (expr_sup_coeff[a] != 0))
            {
                elina_double_interval_mul_expr_coeff(&tmp1, &tmp2, expr_inf_coeff[a], expr_sup_coeff[a], aux_inf_coeff[c], aux_sup_coeff[c]);

                maxRes = fmax(fabs(res_inf_coeff[b]), fabs(res_sup_coeff[b]));
                maxMul = fmax(fabs(tmp1), fabs(tmp2));

                res_inf_coeff[b] = res_inf_coeff[b] + tmp1 + (maxRes + maxMul)*ulp;
                res_sup_coeff[b] = res_sup_coeff[b] + tmp2 + (maxRes + maxMul)*ulp;
            }
        }
    }
}


__global__
void csts_from_previous_layer(double* expr_inf_coeff, double* expr_sup_coeff, double* expr_inf_cst, double* expr_sup_cst, double* res_inf_cst, double* res_sup_cst, double* aux_inf_cst, double* aux_sup_cst, const size_t num_out_neurons_last_layer, const size_t num_out_neurons_current_layer)
{
    const size_t n = blockIdx.x;

    size_t i = 0;

    size_t a = n*num_out_neurons_current_layer + i;

    elina_double_interval_mul_cst_coeff(&res_inf_cst[n], &res_sup_cst[n], expr_inf_coeff[a], expr_sup_coeff[a], aux_inf_cst[i], aux_sup_cst[i]);

    double tmp1, tmp2;
    double maxRes, maxMul;

    for(i = 1; i < num_out_neurons_current_layer; i++)
    {
        a++;

        if((expr_inf_coeff[a] != 0) || (expr_sup_coeff[a] != 0))
        {
            elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, expr_inf_coeff[a], expr_sup_coeff[a], aux_inf_cst[i], aux_sup_cst[i]);

            maxRes = fmax(fabs(res_inf_cst[n]), fabs(res_sup_cst[n]));
            maxMul = fmax(fabs(tmp1), fabs(tmp2));

            res_inf_cst[n] += tmp1 + (maxRes + maxMul)*ulp + min_denormal;
            res_sup_cst[n] += tmp2 + (maxRes + maxMul)*ulp + min_denormal;
        }
    }

    res_inf_cst[n] = res_inf_cst[n] + expr_inf_cst[n];
    res_sup_cst[n] = res_sup_cst[n] + expr_sup_cst[n];
}


void update_state_using_previous_layers(elina_manager_t* man, fppoly_t* fp, const size_t layerno)
{
    auto start = std::chrono::system_clock::now();

    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    const size_t num_out_neurons_last_layer = fp->layers[layerno]->num_out_neurons;
    const size_t num_in_neurons_last_layer = fp->layers[layerno]->num_in_neurons;

    const size_t num_in_neurons_first_layer = fp->layers[0]->num_in_neurons;

    std::cout << "num_out_neurons_last " << num_out_neurons_last_layer << std::endl;

    double* inf_coeff = fp->layers[layerno]->inf_coeff;
    double* sup_coeff = fp->layers[layerno]->sup_coeff;

    double* inf_cst = fp->layers[layerno]->inf_cst;
    double* sup_cst = fp->layers[layerno]->sup_cst;

    double* lb_array = fp->layers[layerno]->lb_array;
    double* ub_array = fp->layers[layerno]->ub_array;

    double* linf_coeff;
    double* lsup_coeff;
    double* linf_cst;
    double* lsup_cst;

    cudaMalloc((void**) &linf_coeff, num_out_neurons_last_layer*num_in_neurons_last_layer*sizeof(double));
    cudaMalloc((void**) &lsup_coeff, num_out_neurons_last_layer*num_in_neurons_last_layer*sizeof(double));
    cudaMalloc((void**) &linf_cst, num_out_neurons_last_layer*sizeof(double));
    cudaMalloc((void**) &lsup_cst, num_out_neurons_last_layer*sizeof(double));

    double* uinf_coeff;
    double* usup_coeff;
    double* uinf_cst;
    double* usup_cst;

    cudaMalloc((void**) &uinf_coeff, num_out_neurons_last_layer*num_in_neurons_last_layer*sizeof(double));
    cudaMalloc((void**) &usup_coeff, num_out_neurons_last_layer*num_in_neurons_last_layer*sizeof(double));
    cudaMalloc((void**) &uinf_cst, num_out_neurons_last_layer*sizeof(double));
    cudaMalloc((void**) &usup_cst, num_out_neurons_last_layer*sizeof(double));

    copy_expr_array<<<num_out_neurons_last_layer, 1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, inf_coeff, sup_coeff, inf_cst, sup_cst, num_out_neurons_last_layer, num_in_neurons_last_layer);
    copy_expr_array<<<num_out_neurons_last_layer, 1>>>(uinf_coeff, usup_coeff, uinf_cst, usup_cst, inf_coeff, sup_coeff, inf_cst, sup_cst, num_out_neurons_last_layer, num_in_neurons_last_layer);

    double* linf_coeff_tmp;
    double* lsup_coeff_tmp;
    double* linf_cst_tmp;
    double* lsup_cst_tmp;

    double* uinf_coeff_tmp;
    double* usup_coeff_tmp;
    double* uinf_cst_tmp;
    double* usup_cst_tmp;

    cudaMalloc((void**) &linf_cst_tmp, num_out_neurons_last_layer*sizeof(double));
    cudaMalloc((void**) &lsup_cst_tmp, num_out_neurons_last_layer*sizeof(double));
    cudaMalloc((void**) &uinf_cst_tmp, num_out_neurons_last_layer*sizeof(double));
    cudaMalloc((void**) &usup_cst_tmp, num_out_neurons_last_layer*sizeof(double));

    for(int k = layerno - 1; k >= 0; k--)
    {
        const size_t num_out_neurons_current_layer = fp->layers[k]->num_out_neurons;
        const size_t num_in_neurons_current_layer  = fp->layers[k]->num_in_neurons;
        std::cout << "num_out_neurons_current " << num_out_neurons_current_layer << " num_in_neurons_current " << num_in_neurons_current_layer << std::endl;

        const dim3 num_blocks_relu(num_out_neurons_last_layer, num_out_neurons_current_layer/num_threads + 1, 1);
        const dim3 num_blocks_linear(num_out_neurons_last_layer, num_in_neurons_current_layer/num_threads + 1, 1);

        std::cout << "num_threads" << num_threads << " num_blocks_relu " << num_blocks_relu.y << " num_blocks_linear " << num_blocks_linear.y << std::endl;

        double* aux_inf_coeff = fp->layers[k]->inf_coeff;
        double* aux_sup_coeff = fp->layers[k]->sup_coeff;

        double* aux_inf_cst = fp->layers[k]->inf_cst;
        double* aux_sup_cst = fp->layers[k]->sup_cst;

        double* aux_lb_array = fp->layers[k]->lb_array;
        double* aux_ub_array = fp->layers[k]->ub_array;

        if(fp->layers[k]->activation == RELU)
        {
            lexpr_replace_relu_bounds<<<num_blocks_relu, num_threads>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, aux_lb_array, aux_ub_array, num_out_neurons_last_layer, num_out_neurons_current_layer);
            uexpr_replace_relu_bounds<<<num_blocks_relu, num_threads>>>(uinf_coeff, usup_coeff, uinf_cst, usup_cst, aux_lb_array, aux_ub_array, num_out_neurons_last_layer, num_out_neurons_current_layer);
        }

        cudaMalloc((void**) &linf_coeff_tmp, num_out_neurons_last_layer*num_in_neurons_current_layer*sizeof(double));
        cudaMalloc((void**) &lsup_coeff_tmp, num_out_neurons_last_layer*num_in_neurons_current_layer*sizeof(double));
        cudaMalloc((void**) &uinf_coeff_tmp, num_out_neurons_last_layer*num_in_neurons_current_layer*sizeof(double));
        cudaMalloc((void**) &usup_coeff_tmp, num_out_neurons_last_layer*num_in_neurons_current_layer*sizeof(double));

        coeffs_from_previous_layer<<<num_blocks_linear, num_threads>>>(linf_coeff, lsup_coeff, linf_coeff_tmp, lsup_coeff_tmp, aux_inf_coeff, aux_sup_coeff, num_out_neurons_last_layer, num_out_neurons_current_layer, num_in_neurons_current_layer);
        coeffs_from_previous_layer<<<num_blocks_linear, num_threads>>>(uinf_coeff, usup_coeff, uinf_coeff_tmp, usup_coeff_tmp, aux_inf_coeff, aux_sup_coeff, num_out_neurons_last_layer, num_out_neurons_current_layer, num_in_neurons_current_layer);

        csts_from_previous_layer<<<num_out_neurons_last_layer, 1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, linf_cst_tmp, lsup_cst_tmp, aux_inf_cst, aux_sup_cst, num_out_neurons_last_layer, num_out_neurons_current_layer);
        csts_from_previous_layer<<<num_out_neurons_last_layer, 1>>>(uinf_coeff, usup_coeff, uinf_cst, usup_cst, uinf_cst_tmp, usup_cst_tmp, aux_inf_cst, aux_sup_cst, num_out_neurons_last_layer, num_out_neurons_current_layer);

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

    compute_lb_from_expr<<<num_out_neurons_last_layer, 1>>>(lb_array, linf_coeff, lsup_coeff, linf_cst, fp->input_inf, fp->input_sup, num_out_neurons_last_layer, num_in_neurons_first_layer);
    compute_ub_from_expr<<<num_out_neurons_last_layer, 1>>>(ub_array, uinf_coeff, usup_coeff, usup_cst, fp->input_inf, fp->input_sup, num_out_neurons_last_layer, num_in_neurons_first_layer);

    cudaFree(linf_coeff);
    cudaFree(lsup_coeff);
    cudaFree(linf_cst);
    cudaFree(lsup_cst);

    cudaFree(uinf_coeff);
    cudaFree(usup_coeff);
    cudaFree(uinf_cst);
    cudaFree(usup_cst);

    cudaFree(linf_coeff_tmp);
    cudaFree(lsup_coeff_tmp);
    cudaFree(linf_cst_tmp);
    cudaFree(lsup_cst_tmp);

    cudaFree(uinf_coeff_tmp);
    cudaFree(usup_coeff_tmp);
    cudaFree(uinf_cst_tmp);
    cudaFree(usup_cst_tmp);

    cudaDeviceSynchronize();

    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl << std::endl;
}


void ffn_handle_intermediate_layer(elina_manager_t* man, elina_abstract0_t* element, const double** weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons, const activation_type_t activation)
{
    fppoly_t* fp = fppoly_of_abstract0(element);
    fppoly_add_new_layer(fp, num_out_neurons, num_in_neurons, FFN, activation);

    double* inf_coeff = fp->layers[fp->numlayers - 1]->inf_coeff;
    double* sup_coeff = fp->layers[fp->numlayers - 1]->sup_coeff;

    double* inf_cst = fp->layers[fp->numlayers - 1]->inf_cst;
    double* sup_cst = fp->layers[fp->numlayers - 1]->sup_cst;

    layer_create_dense_exprs(inf_coeff, sup_coeff, inf_cst, sup_cst, weights, bias, num_out_neurons, num_in_neurons);

    update_state_using_previous_layers(man, fp, fp->numlayers - 1);
}


void ffn_handle_intermediate_relu_layer(elina_manager_t* man, elina_abstract0_t* element, const double** weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons)
{
    ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, RELU);
}

void ffn_handle_intermediate_sigmoid_layer(elina_manager_t* man, elina_abstract0_t* element, const double** weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons)
{
    //ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, SIGMOID);
}

void ffn_handle_intermediate_tanh_layer(elina_manager_t* man, elina_abstract0_t* element, const double** weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons)
{
    //ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, TANH);
}


__global__
void print_bounds(double* bounds_array, const size_t num_out_neurons)
{
    for(size_t i = 0; i < num_out_neurons; i++)
    {
        printf("out inf number %i is: %g\n", i, bounds_array[i]);
    }
}


void ffn_handle_last_layer(elina_manager_t* man, elina_abstract0_t* element, const double** weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons, const bool has_activation, const activation_type_t activation)
{
    fppoly_t* fp = fppoly_of_abstract0(element);
    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    if(has_activation)
    {
        fppoly_add_new_layer(fp, num_out_neurons, num_in_neurons, FFN, activation);
    }
    else
    {
        fppoly_add_new_layer(fp, num_out_neurons, num_in_neurons, FFN, NONE);
    }

    double* inf_coeff = fp->layers[fp->numlayers - 1]->inf_coeff;
    double* sup_coeff = fp->layers[fp->numlayers - 1]->sup_coeff;

    double* inf_cst = fp->layers[fp->numlayers - 1]->inf_cst;
    double* sup_cst = fp->layers[fp->numlayers - 1]->sup_cst;

    layer_create_dense_exprs(inf_coeff, sup_coeff, inf_cst, sup_cst, weights, bias, num_out_neurons, num_in_neurons);

    update_state_using_previous_layers(man, fp, fp->numlayers - 1);

    double* lb_array = fp->layers[fp->numlayers - 1]->lb_array;
    double* ub_array = fp->layers[fp->numlayers - 1]->ub_array;

    print_bounds<<<1, 1>>>(lb_array, num_out_neurons);
    print_bounds<<<1, 1>>>(ub_array, num_out_neurons);
}

void ffn_handle_last_relu_layer(elina_manager_t* man, elina_abstract0_t* element, const double** weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons, const bool has_relu)
{
    ffn_handle_last_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, has_relu, RELU);
}

void ffn_handle_last_sigmoid_layer(elina_manager_t* man, elina_abstract0_t* element, const double** weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons, const bool has_sigmoid)
{
    //ffn_handle_last_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, has_sigmoid, SIGMOID);
}

void ffn_handle_last_tanh_layer(elina_manager_t* man, elina_abstract0_t* element, const double** weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons, const bool has_tanh)
{
    //ffn_handle_last_layer(man, element, weights, bias, num_out_neurons, num_in_neurons, has_tanh, TANH);
}


__global__
void create_sub_expr(double* inf_coeff, double* sup_coeff, double* inf_cst, double* sup_cst, const size_t index, const elina_dim_t y, const elina_dim_t x)
{
    inf_cst[index] = 0;
    sup_cst[index] = 0;

    for(size_t i = 0; i < 10; i++)
    {
        inf_coeff[index*10 + i] = 0.;
        sup_coeff[index*10 + i] = 0.;
    }

    inf_coeff[index*10 + y] = -1.;
    sup_coeff[index*10 + y] = 1.;

    inf_coeff[index*10 + x] = 1.;
    sup_coeff[index*10 + x] = -1.;
}


void get_lb_using_previous_layers(elina_manager_t* man, const fppoly_t* const fp)
{
    const size_t numlayers = fp->numlayers;
    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    const size_t num_out_neurons_last_layer = 90;

    const size_t num_in_neurons_first_layer = fp->layers[0]->num_in_neurons;

    double* lb_dev;
    cudaMalloc((void**) &lb_dev, num_out_neurons_last_layer*sizeof(double));

    double* linf_coeff;
    double* lsup_coeff;
    double* linf_cst;
    double* lsup_cst;

    cudaMalloc((void**) &linf_coeff, num_out_neurons_last_layer*10*sizeof(double*));
    cudaMalloc((void**) &lsup_coeff, num_out_neurons_last_layer*10*sizeof(double*));
    cudaMalloc((void**) &linf_cst, num_out_neurons_last_layer*10*sizeof(double));
    cudaMalloc((void**) &lsup_cst, num_out_neurons_last_layer*10*sizeof(double));

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

    double* linf_coeff_tmp;
    double* lsup_coeff_tmp;
    double* linf_cst_tmp;
    double* lsup_cst_tmp;

    cudaMalloc((void**) &linf_coeff_tmp, num_out_neurons_last_layer*sizeof(double*));
    cudaMalloc((void**) &lsup_coeff_tmp, num_out_neurons_last_layer*sizeof(double*));
    cudaMalloc((void**) &linf_cst_tmp, num_out_neurons_last_layer*sizeof(double));
    cudaMalloc((void**) &lsup_cst_tmp, num_out_neurons_last_layer*sizeof(double));

    for(int k = numlayers - 1; k >= 0; k--)
    {
        const size_t num_out_neurons_current_layer = fp->layers[k]->num_out_neurons;
        const size_t num_in_neurons_current_layer  = fp->layers[k]->num_in_neurons;

        const dim3 num_blocks_relu(num_out_neurons_last_layer, num_out_neurons_current_layer/num_threads + 1, 1);
        const dim3 num_blocks_linear(num_out_neurons_last_layer, num_in_neurons_current_layer/num_threads + 1, 1);

        double* aux_inf_coeff = fp->layers[k]->inf_coeff;
        double* aux_sup_coeff = fp->layers[k]->sup_coeff;

        double* aux_inf_cst = fp->layers[k]->inf_cst;
        double* aux_sup_cst = fp->layers[k]->sup_cst;

        double* aux_lb_array = fp->layers[k]->lb_array;
        double* aux_ub_array = fp->layers[k]->ub_array;

        if(fp->layers[k]->activation == RELU)
        {
            lexpr_replace_relu_bounds<<<num_blocks_relu, num_threads>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, aux_lb_array, aux_ub_array, num_out_neurons_last_layer, num_out_neurons_current_layer);
        }

        cudaMalloc((void**) &linf_coeff_tmp, num_out_neurons_last_layer*num_in_neurons_current_layer*sizeof(double));
        cudaMalloc((void**) &lsup_coeff_tmp, num_out_neurons_last_layer*num_in_neurons_current_layer*sizeof(double));

        coeffs_from_previous_layer<<<num_blocks_linear, num_threads>>>(linf_coeff, lsup_coeff, linf_coeff_tmp, lsup_coeff_tmp, aux_inf_coeff, aux_sup_coeff, num_out_neurons_last_layer, num_out_neurons_current_layer, num_in_neurons_current_layer);

        csts_from_previous_layer<<<num_out_neurons_last_layer, 1>>>(linf_coeff, lsup_coeff, linf_cst, lsup_cst, linf_cst_tmp, lsup_cst_tmp, aux_inf_cst, aux_sup_cst, num_out_neurons_last_layer, num_out_neurons_current_layer);

        std::swap(linf_coeff, linf_coeff_tmp);
        std::swap(lsup_coeff, lsup_coeff_tmp);
        std::swap(linf_cst, linf_cst_tmp);
        std::swap(lsup_cst, lsup_cst_tmp);

        cudaFree(linf_coeff_tmp);
        cudaFree(lsup_coeff_tmp);
    }

    compute_lb_from_expr<<<num_out_neurons_last_layer, 1>>>(lb_dev, linf_coeff, lsup_coeff, linf_cst, fp->input_inf, fp->input_sup, num_out_neurons_last_layer, num_in_neurons_first_layer);

    cudaFree(linf_coeff);
    cudaFree(lsup_coeff);
    cudaFree(linf_cst);
    cudaFree(lsup_cst);

    cudaFree(linf_coeff_tmp);
    cudaFree(lsup_coeff_tmp);
    cudaFree(linf_cst_tmp);
    cudaFree(lsup_cst_tmp);

    double lb[num_out_neurons_last_layer];
    cudaMemcpy(&lb, lb_dev, num_out_neurons_last_layer*sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(lb_dev);

    for(size_t i = 0; i < num_out_neurons_last_layer; i++)
    {
        if(lb[i] < 0)
        {
            results[i] = true;
        }
        else
        {
            results[i] = false;
        }
    }
}


bool is_greater(elina_manager_t* man, elina_abstract0_t* element, const elina_dim_t y, const elina_dim_t x)
{
    const fppoly_t* fp = fppoly_of_abstract0(element);
    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    if(!results_calculated)
    {
        get_lb_using_previous_layers(man, fp);
        results_calculated = true;

        return results[0];
    }
    else
    {
        bool result = results[output_counter];
        output_counter++;

        return result;
    }
}


void device_layer_create_sparse_exprs(double* inf_coeff, double* sup_coeff, double* inf_cst, double* sup_cst, const double* filter_weights,
                                      const double* filter_bias, const size_t* input_size, const size_t* output_size, const size_t* filter_size,
                                      const size_t* strides, const long int pad_top, const long int pad_left, const size_t num_pixels)
{
    const size_t num_out_neurons = output_size[0]*output_size[1]*output_size[2];

    double* dense_coeff = (double*) calloc(num_out_neurons*num_pixels, sizeof(double));
    double* bias = (double*) calloc(num_out_neurons, sizeof(double));

    for(size_t out_x = 0; out_x < output_size[0]; out_x++)
    {
        for(size_t out_y = 0; out_y < output_size[1]; out_y++)
        {
            for(size_t out_z = 0; out_z < output_size[2]; out_z++)
            {
                const size_t mat_x = out_x*output_size[1]*output_size[2] + out_y*output_size[2] + out_z;

                for(size_t x_shift = 0; x_shift < filter_size[0]; x_shift++)
                {
                    for(size_t y_shift = 0; y_shift < filter_size[1]; y_shift++)
                    {
                        for(size_t inp_z = 0; inp_z < input_size[2]; inp_z++)
                        {
                            const long int x_val = out_x*strides[0] + x_shift - pad_top;
                            const long int y_val = out_y*strides[1] + y_shift - pad_left;

                            if((y_val < 0) || (y_val >= (long int)input_size[1]))
                            {
                                continue;
                            }

                            if((x_val < 0) || (x_val >= (long int)input_size[0]))
                            {
                                continue;
                            }

                            const size_t mat_y = x_val*input_size[1]*input_size[2] + y_val*input_size[2] + inp_z;

                            if(mat_y >= num_pixels)
                            {
                                continue;
                            }

                            const size_t filter_index = x_shift*filter_size[1]*input_size[2]*output_size[2] + y_shift*input_size[2]*output_size[2] + inp_z*output_size[2] + out_z;
                            dense_coeff[mat_x*num_pixels + mat_y] = filter_weights[filter_index];
                        }
                    }
                }

                bias[mat_x] = filter_bias[out_z];
            }
        }
    }

    double* dense_coeff_dev;
    double* bias_dev;

    cudaMalloc((void**) &dense_coeff_dev, num_out_neurons*num_pixels*sizeof(double));
    cudaMalloc((void**) &bias_dev, num_out_neurons*sizeof(double));

    cudaMemcpy(dense_coeff_dev, dense_coeff, num_out_neurons*num_pixels*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(bias_dev, bias, num_out_neurons*sizeof(double), cudaMemcpyHostToDevice);

    device_layer_create_dense_expr<<<num_out_neurons, 1>>>(inf_coeff, sup_coeff, inf_cst, sup_cst, dense_coeff_dev, bias_dev, num_out_neurons, num_pixels);

    cudaFree(dense_coeff_dev);
    cudaFree(bias_dev);

    free(dense_coeff);
    free(bias);
}


void layer_create_sparse_exprs(fppoly_t* const fp, const double* filter_weights, const double* filter_bias,
                               const size_t* input_size, const size_t* filter_size, const size_t num_filters, const size_t* strides,
                               const bool is_valid_padding, const bool has_bias)
{
    const size_t num_pixels = input_size[0]*input_size[1]*input_size[2];

    size_t output_size[3];

    if(is_valid_padding)
    {
        output_size[0] = ceil((double)(input_size[0] - filter_size[0] + 1)/(double)strides[0]);
        output_size[1] = ceil((double)(input_size[1] - filter_size[1] + 1)/(double)strides[1]);
    }
    else
    {
        output_size[0] = ceil((double)input_size[0]/(double)strides[0]);
        output_size[1] = ceil((double)input_size[1]/(double)strides[1]);
    }

    output_size[2] = num_filters;

    const size_t num_out_neurons = output_size[0]*output_size[1]*output_size[2];
    fppoly_add_new_layer(fp, num_out_neurons, num_pixels, CONV, RELU);

    double* inf_coeff = fp->layers[fp->numlayers - 1]->inf_coeff;
    double* sup_coeff = fp->layers[fp->numlayers - 1]->sup_coeff;

    double* inf_cst = fp->layers[fp->numlayers - 1]->inf_cst;
    double* sup_cst = fp->layers[fp->numlayers - 1]->sup_cst;

    long int pad_along_height = 0;
    long int pad_along_width = 0;
    long int pad_top = 0;
    long int pad_left = 0;

    if(!is_valid_padding)
    {
        if(input_size[0]%strides[0] == 0)
        {
            const long int tmp = filter_size[0] - strides[0];
            pad_along_height = max(tmp, long(0));
        }
        else
        {
            const long int tmp = filter_size[0] - (input_size[0]%strides[0]);
            pad_along_height = max(tmp, long(0));
        }

        if(input_size[1]%strides[1] == 0)
        {
            const long int tmp = filter_size[1] - strides[1];
            pad_along_width = max(tmp, long(0));
        }
        else
        {
            const long int tmp = filter_size[1] - (input_size[1]%strides[1]);
            pad_along_width = max(tmp, long(0));
        }

        pad_top = pad_along_height/2;
        pad_left = pad_along_width/2;
    }

    const size_t size = filter_size[0]*filter_size[1]*input_size[2]*output_size[2];

    double* filter_weights_tmp = (double*) malloc(size*sizeof(double));
    double* filter_bias_tmp = (double*) calloc(output_size[2], sizeof(double));

    size_t* input_size_tmp = (size_t*) malloc(3*sizeof(size_t));
    size_t* output_size_tmp = (size_t*) malloc(3*sizeof(size_t));
    size_t* filter_size_tmp = (size_t*) malloc(2*sizeof(size_t));
    size_t* strides_tmp = (size_t*) malloc(2*sizeof(size_t));

    cudaMemcpy(filter_weights_tmp, filter_weights, size*sizeof(double), cudaMemcpyHostToHost);

    if(has_bias)
    {
        cudaMemcpy(filter_bias_tmp, filter_bias, output_size[2]*sizeof(double), cudaMemcpyHostToHost);
    }

    cudaMemcpy(input_size_tmp, input_size, 3*sizeof(size_t), cudaMemcpyHostToHost);
    cudaMemcpy(output_size_tmp, output_size, 3*sizeof(size_t), cudaMemcpyHostToHost);
    cudaMemcpy(filter_size_tmp, filter_size, 2*sizeof(size_t), cudaMemcpyHostToHost);
    cudaMemcpy(strides_tmp, strides, 2*sizeof(size_t), cudaMemcpyHostToHost);

    device_layer_create_sparse_exprs(inf_coeff, sup_coeff, inf_cst, sup_cst, filter_weights_tmp, filter_bias_tmp, input_size_tmp, output_size_tmp,
                                     filter_size_tmp, strides_tmp, pad_top, pad_left, num_pixels);

    free(filter_weights_tmp);
    free(filter_bias_tmp);

    free(input_size_tmp);
    free(output_size_tmp);
    free(filter_size_tmp);
    free(strides_tmp);
}


void conv_handle_first_layer(elina_manager_t* man, elina_abstract0_t* element, const double* filter_weights, const double* filter_bias,
                             const size_t* input_size, const size_t* filter_size, const size_t num_filters, const size_t* strides,
                             const bool is_valid_padding, const bool has_bias)
{
    fppoly_t* const fp = fppoly_of_abstract0(element);
    fp->layers = (layer_t**) malloc(20*sizeof(layer_t*));

    layer_create_sparse_exprs(fp, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, is_valid_padding, has_bias);

    double* inf_coeff = fp->layers[0]->inf_coeff;
    double* sup_coeff = fp->layers[0]->sup_coeff;

    double* inf_cst = fp->layers[0]->inf_cst;
    double* sup_cst = fp->layers[0]->sup_cst;

    layer_compute_bounds_from_exprs(inf_coeff, sup_coeff, inf_cst, sup_cst, fp->layers[0]->lb_array, fp->layers[0]->ub_array, fp->input_inf, fp->input_sup, fp->layers[0]->num_out_neurons, fp->layers[0]->num_in_neurons);
}


void conv_handle_intermediate_relu_layer(elina_manager_t* man, elina_abstract0_t* element, const double* filter_weights, const double* filter_bias,
                                         const size_t* input_size, const size_t* filter_size, const size_t num_filters, const size_t* strides,
                                         const bool is_valid_padding, const bool has_bias)
{
    fppoly_t* const fp = fppoly_of_abstract0(element);

    layer_create_sparse_exprs(fp, filter_weights, filter_bias, input_size, filter_size, num_filters, strides, is_valid_padding, has_bias);

    update_state_using_previous_layers(man, fp, fp->numlayers - 1);
}


void free_layer(layer_t* layer)
{
    cudaFree(layer->inf_coeff);
    cudaFree(layer->sup_coeff);
    cudaFree(layer->inf_cst);
    cudaFree(layer->sup_cst);

    layer->inf_coeff = nullptr;
    layer->sup_coeff = nullptr;
    layer->inf_cst = nullptr;
    layer->sup_cst = nullptr;

    cudaFree(layer->lb_array);
    cudaFree(layer->ub_array);

    layer->lb_array = nullptr;
    layer->ub_array = nullptr;

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
    fp->input_inf = nullptr;
    cudaFree(fp->input_sup);
    fp->input_sup = nullptr;

    free(fp);
    fp = nullptr;
}


void layer_print(const layer_t* layer)
{
    //neurons_print<<<1, 1>>>(layer->neurons, layer->num_out_neurons);
}


void fppoly_fprint(FILE* const stream, elina_manager_t* man, const fppoly_t* const fp, const char** name_of_dim)
{
    for(size_t i = 0; i < fp->numlayers; i++)
    {
        printf("layer: %zu\n", i);
        layer_print(fp->layers[i]);
    }
}

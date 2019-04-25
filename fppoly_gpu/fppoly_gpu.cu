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

bool initialized = false;

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
    size_t set_limit = 3*size_t(2u << 30u);
    size_t limit;

    if(!initialized)
    {
        // SET THIS ONLY ONCE BEFORE FIRST KERNEL CALL, ELSE CUDA ERROR (INVALID ARGUMENT ERROR)
        cudaDeviceSetLimit(cudaLimitMallocHeapSize, set_limit);
        initialized = true;
    }

    cudaDeviceGetLimit(&limit, cudaLimitMallocHeapSize);

    std::cout << "The Heap Limit is " << limit << std::endl;

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

    cudaMalloc((void**) &layer->expr_array, num_out_neurons*sizeof(expr_t*));

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
void compute_lb_from_expr(double* lb_array, expr_t** expr_array, double* input_inf, double* input_sup, const size_t num_exprs, const size_t expr_size)
{

    size_t n = blockIdx.x;

    if(n < num_exprs)
    {
        expr_t* expr = expr_array[n];

        double res_inf = expr->inf_cst;

        if((expr->inf_coeff == nullptr) || (expr->sup_coeff == nullptr))
        {
            lb_array[n] = 0;

            return;
        }

        double tmp1, tmp2;
        size_t k;

        for(size_t i = 0; i < expr_size; i++)
        {
            k = i;

            elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i], expr->sup_coeff[i], input_inf[k], input_sup[k]);
            res_inf = res_inf + tmp1;
        }

        lb_array[n] = res_inf;
    }
}


__global__
void compute_ub_from_expr(double* ub_array, expr_t** expr_array, double* input_inf, double* input_sup, const size_t num_exprs, const size_t expr_size)
{
    size_t n = blockIdx.x;

    if(n < num_exprs)
    {
        expr_t* expr = expr_array[n];

        double res_sup = expr->sup_cst;

        if((expr->inf_coeff == nullptr) || (expr->sup_coeff == nullptr))
        {
            ub_array[n] = 0;

            return;
        }

        double tmp1, tmp2;
        size_t k;

        for(size_t i = 0; i < expr_size; i++)
        {
            k = i;

            elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i], expr->sup_coeff[i], input_inf[k], input_sup[k]);
            res_sup = res_sup + tmp2;
        }

       ub_array[n] = res_sup;
    }
}


__global__
void device_layer_create_dense_expr(expr_t** expr_array, const double* weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons)
{
    size_t i = blockIdx.x;

    if(i < num_out_neurons)
    {
        const double* weight_i = weights + i*num_in_neurons;
        const double bias_i = bias[i];

        expr_t* expr = (expr_t*) malloc(sizeof(expr_t));

        expr->inf_coeff = (double*) malloc(num_in_neurons*sizeof(double));
        expr->sup_coeff = (double*) malloc(num_in_neurons*sizeof(double));

        expr->inf_cst = -bias_i;
        expr->sup_cst = bias_i;

        for(size_t j = 0; j < num_in_neurons; j++)
        {
            expr->inf_coeff[j] = -weight_i[j];
            expr->sup_coeff[j] = weight_i[j];
        }

        expr_array[i] = expr;
    }
}


void layer_create_dense_exprs(expr_t** expr_array, const double** weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons)
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

    device_layer_create_dense_expr<<<num_out_neurons, 1>>>(expr_array, tmp_weights, tmp_bias, num_out_neurons, num_in_neurons);

    cudaFree(tmp_weights);
    cudaFree(tmp_bias);
}


__global__
void copy_expr_array(expr_t** target_expr_array, expr_t** source_expr_array, const size_t num_exprs, const size_t expr_size)
{
    size_t i = blockIdx.x;

    if(i < num_exprs)
    {
        target_expr_array[i] = (expr_t*) malloc(sizeof(expr_t));

        target_expr_array[i]->inf_coeff = (double*) malloc(expr_size*sizeof(double));
        target_expr_array[i]->sup_coeff = (double*) malloc(expr_size*sizeof(double));

        target_expr_array[i]->inf_cst = source_expr_array[i]->inf_cst;
        target_expr_array[i]->sup_cst = source_expr_array[i]->sup_cst;

        for(size_t j = 0; j < expr_size; j++)
        {
            target_expr_array[i]->inf_coeff[j] = source_expr_array[i]->inf_coeff[j];
            target_expr_array[i]->sup_coeff[j] = source_expr_array[i]->sup_coeff[j];
        }
    }
}


__global__
void free_expr_array(expr_t** expr_array, const size_t size)
{
    size_t i = blockIdx.x;

    if(i < size)
    {
        if(expr_array[i]->inf_coeff)
        {
            free(expr_array[i]->inf_coeff);
            expr_array[i]->inf_coeff = nullptr;
        }

        if(expr_array[i]->sup_coeff)
        {
            free(expr_array[i]->sup_coeff);
            expr_array[i]->sup_coeff = nullptr;
        }

        free(expr_array[i]);
        expr_array[i] = nullptr;
    }
}


void layer_compute_bounds_from_exprs(expr_t** expr_array, double* lb_array, double* ub_array, double* input_inf, double* input_sup, const size_t num_out_neurons, const size_t num_in_neurons)
{
    // allocate
    expr_t** lexpr_array;
    expr_t** uexpr_array;

    cudaMalloc((void**) &lexpr_array, num_out_neurons*sizeof(expr_t*));
    cudaMalloc((void**) &uexpr_array, num_out_neurons*sizeof(expr_t*));

    copy_expr_array<<<num_out_neurons, 1>>>(lexpr_array, expr_array, num_out_neurons, num_in_neurons);
    copy_expr_array<<<num_out_neurons, 1>>>(uexpr_array, expr_array, num_out_neurons, num_in_neurons);

    compute_lb_from_expr<<<num_out_neurons, 1>>>(lb_array, lexpr_array, input_inf, input_sup, num_out_neurons, num_in_neurons);
    compute_ub_from_expr<<<num_out_neurons, 1>>>(ub_array, uexpr_array, input_inf, input_sup, num_out_neurons, num_in_neurons);

    // free
    free_expr_array<<<num_out_neurons, 1>>>(lexpr_array, num_out_neurons);
    free_expr_array<<<num_out_neurons, 1>>>(uexpr_array, num_out_neurons);

    cudaFree(lexpr_array);
    cudaFree(uexpr_array);
}


void ffn_handle_first_layer(elina_manager_t* man, elina_abstract0_t* abs, const double** weights, const double* bias, const size_t size, const size_t num_pixels, const activation_type_t activation)
{
    fppoly_t* res = fppoly_of_abstract0(abs);
    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    res->layers = (layer_t**) malloc(20*sizeof(layer_t*));
    fppoly_add_new_layer(res, size, num_pixels, FFN, activation);

    expr_t** expr_array = res->layers[0]->expr_array;

    layer_create_dense_exprs(expr_array, weights, bias, size, num_pixels);
    layer_compute_bounds_from_exprs(expr_array, res->layers[0]->lb_array, res->layers[0]->ub_array, res->input_inf, res->input_sup, res->layers[0]->num_out_neurons, res->layers[0]->num_in_neurons);
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
void lexpr_replace_relu_bounds(expr_t** expr_array, double* lb_array, double* ub_array, const size_t num_out_neurons_last_layer, const size_t num_out_neurons_current_layer)
{
    size_t n = blockIdx.x;
    size_t i = blockIdx.y*blockDim.x + threadIdx.x;

    if(n < num_out_neurons_last_layer)
    {
        expr_t* expr = expr_array[n];

        if(i < num_out_neurons_current_layer)
        {
            const double lb = lb_array[i];
            const double ub = ub_array[i];
            const double width = ub + lb;
            const double lambda_inf = -ub/width;
            const double lambda_sup = ub/width;

            const double old_inf_coeff = expr->inf_coeff[i];
            const double old_sup_coeff = expr->sup_coeff[i];

            if((old_sup_coeff == 0) && (old_inf_coeff == 0))
            {
                expr->inf_coeff[i] = 0.0;
                expr->sup_coeff[i] = 0.0;

                return;
            }
            else if(ub <= 0)
            {
                expr->inf_coeff[i] = 0.0;
                expr->sup_coeff[i] = 0.0;

                return;
            }
            else if(lb < 0)
            {
                expr->inf_coeff[i] = old_inf_coeff;
                expr->sup_coeff[i] = old_sup_coeff;
            }
            else if(old_sup_coeff < 0)
            {
                const double mu_inf = lambda_inf*lb;
                const double mu_sup = lambda_sup*lb;
                elina_double_interval_mul_expr_coeff(&expr->inf_coeff[i], &expr->sup_coeff[i], lambda_inf, lambda_sup, old_inf_coeff, old_sup_coeff);
                double tmp1, tmp2;
                elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, mu_inf, mu_sup, old_inf_coeff, old_sup_coeff);
                atomicAdd(&expr->inf_cst, tmp1 + min_denormal);
                atomicAdd(&expr->sup_cst, tmp2 + min_denormal);
            }
            else if (old_inf_coeff < 0)
            {
                const double area1 = lb*ub;
                const double area2 = 0.5*ub*width;
                const double area3 = 0.5*lb*width;

                if((area1 < area2) && (area1 < area3))
                {
                    elina_double_interval_mul_expr_coeff(&expr->inf_coeff[i], &expr->sup_coeff[i], lambda_inf, lambda_sup, old_inf_coeff, old_sup_coeff);
                }
                else if((area2 < area1) && (area2 < area3))
                {
                    expr->inf_coeff[i] = 0.0;
                    expr->sup_coeff[i] = 0.0;
                }
                else
                {
                    expr->inf_coeff[i] = old_inf_coeff;
                    expr->sup_coeff[i] = old_sup_coeff;
                }
            }
            else
            {
                expr->inf_coeff[i] = 0.0;
                expr->sup_coeff[i] = 0.0;
                double tmp1, tmp2;
                elina_double_interval_mul(&tmp1, &tmp2, old_inf_coeff, old_sup_coeff, 0, ub);
                atomicAdd(&expr->inf_cst, tmp1);
                atomicAdd(&expr->sup_cst, -tmp1);
            }
        }
    }
}


__global__
void uexpr_replace_relu_bounds(expr_t** expr_array, double* lb_array, double* ub_array, const size_t num_out_neurons_last_layer, const size_t num_out_neurons_current_layer)
{
    size_t n = blockIdx.x;
    size_t i = blockIdx.y*blockDim.x + threadIdx.x;

    if(n < num_out_neurons_last_layer)
    {
        expr_t* expr = expr_array[n];

        if(i < num_out_neurons_current_layer)
        {
            const double lb = lb_array[i];
            const double ub = ub_array[i];
            const double width = ub + lb;
            const double lambda_inf = -ub/width;
            const double lambda_sup = ub/width;

            const double old_inf_coeff = expr->inf_coeff[i];
            const double old_sup_coeff = expr->sup_coeff[i];

            if((old_sup_coeff == 0) && (old_inf_coeff == 0))
            {
                expr->inf_coeff[i] = 0.0;
                expr->sup_coeff[i] = 0.0;

                return;
            }
            else if(ub <= 0)
            {
                expr->inf_coeff[i] = 0.0;
                expr->sup_coeff[i] = 0.0;

                return;
            }
            else if(lb < 0)
            {
                expr->inf_coeff[i] = old_inf_coeff;
                expr->sup_coeff[i] = old_sup_coeff;
            }
            else if(old_inf_coeff < 0)
            {
                const double mu_inf = lambda_inf*lb;
                const double mu_sup = lambda_sup*lb;
                elina_double_interval_mul_expr_coeff(&expr->inf_coeff[i], &expr->sup_coeff[i], lambda_inf, lambda_sup, old_inf_coeff, old_sup_coeff);
                double tmp1, tmp2;
                elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, mu_inf, mu_sup, old_inf_coeff, old_sup_coeff);
                atomicAdd(&expr->inf_cst, tmp1 + min_denormal);
                atomicAdd(&expr->sup_cst, tmp2 + min_denormal);
            }
            else if(old_sup_coeff < 0)
            {
                const double area1 = lb*ub;
                const double area2 = 0.5*ub*width;
                const double area3 = 0.5*lb*width;

                if((area1 < area2) && (area1 < area3))
                {
                    elina_double_interval_mul_expr_coeff(&expr->inf_coeff[i], &expr->sup_coeff[i], lambda_inf, lambda_sup, old_inf_coeff, old_sup_coeff);
                }
                else if((area2 < area1) && (area2 < area3))
                {
                    expr->inf_coeff[i] = 0.0;
                    expr->sup_coeff[i] = 0.0;
                }
                else
                {
                    expr->inf_coeff[i] = old_inf_coeff;
                    expr->sup_coeff[i] = old_sup_coeff;
                }
            }
            else
            {
                expr->inf_coeff[i] = 0.0;
                expr->sup_coeff[i] = 0.0;
                double tmp1, tmp2;
                elina_double_interval_mul(&tmp1, &tmp2, old_inf_coeff, old_sup_coeff, 0, ub);
                atomicAdd(&expr->inf_cst, -tmp2);
                atomicAdd(&expr->sup_cst, tmp2);
            }
        }
    }
}


__global__
void expr_from_previous_layer(expr_t** expr_array, expr_t** res_array, expr_t** aux_expr_array, const size_t num_out_neurons_last_layer, const size_t num_out_neurons_current_layer, const size_t num_in_neurons_current_layer)
{
    size_t n = blockIdx.x;
    size_t j = blockIdx.y*blockDim.x + threadIdx.x;

    if(n < num_out_neurons_last_layer)
    {
        expr_t* expr = expr_array[n];
        expr_t* res = res_array[n];

        size_t i;

        if(j < num_in_neurons_current_layer)
        {
            if(j == 0)
            {
                i = 0;

                expr_t* aux_expr = aux_expr_array[i];

                elina_double_interval_mul_cst_coeff(&res->inf_cst, &res->sup_cst, expr->inf_coeff[0], expr->sup_coeff[0], aux_expr->inf_cst, aux_expr->sup_cst);

                double tmp1, tmp2;
                double maxRes, maxMul;

                for(i = 1; i < num_out_neurons_current_layer; i++)
                {
                    if((expr->inf_coeff[i] != 0) || (expr->sup_coeff[i] != 0))
                    {
                        aux_expr = aux_expr_array[i];

                        elina_double_interval_mul_cst_coeff(&tmp1, &tmp2, expr->inf_coeff[i], expr->sup_coeff[i], aux_expr->inf_cst, aux_expr->sup_cst);

                        maxRes = fmax(fabs(res->inf_cst), fabs(res->sup_cst));
                        maxMul = fmax(fabs(tmp1), fabs(tmp2));

                        res->inf_cst += tmp1 + (maxRes + maxMul)*ulp + min_denormal;
                        res->sup_cst += tmp2 + (maxRes + maxMul)*ulp + min_denormal;
                    }

                }

                res->inf_cst = res->inf_cst + expr->inf_cst;
                res->sup_cst = res->sup_cst + expr->sup_cst;
            }

            i = 0;

            expr_t* aux_expr = aux_expr_array[i];

            elina_double_interval_mul_expr_coeff(&res->inf_coeff[j], &res->sup_coeff[j], expr->inf_coeff[0], expr->sup_coeff[0], aux_expr->inf_coeff[j], aux_expr->sup_coeff[j]);

            double tmp1, tmp2;
            double maxRes, maxMul;

            for(i = 1; i < num_out_neurons_current_layer; i++)
            {
                if((expr->inf_coeff[i] != 0) || (expr->sup_coeff[i] != 0))
                {
                    aux_expr = aux_expr_array[i];

                    elina_double_interval_mul_expr_coeff(&tmp1, &tmp2, expr->inf_coeff[i], expr->sup_coeff[i], aux_expr->inf_coeff[j], aux_expr->sup_coeff[j]);

                    maxRes = fmax(fabs(res->inf_coeff[j]), fabs(res->sup_coeff[j]));
                    maxMul = fmax(fabs(tmp1), fabs(tmp2));

                    res->inf_coeff[j] = res->inf_coeff[j] + tmp1 + (maxRes + maxMul)*ulp;
                    res->sup_coeff[j] = res->sup_coeff[j] + tmp2 + (maxRes + maxMul)*ulp;
                }
            }
        }
    }
}


__global__
void layer_allocate_exprs(expr_t** expr_array, const size_t array_size, const size_t expr_size)
{
    size_t i = blockIdx.x;

    if(i < array_size)
    {
        expr_t* expr = (expr_t*) malloc(sizeof(expr_t));

        expr->inf_coeff = (double*) malloc(expr_size*sizeof(double));
        expr->sup_coeff = (double*) malloc(expr_size*sizeof(double));

        expr->inf_cst = 0;
        expr->sup_cst = 0;

        expr_array[i] = expr;
    }
}


void update_state_using_previous_layers(elina_manager_t* man, fppoly_t* fp, const size_t layerno)
{
    auto start = std::chrono::system_clock::now();

    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    const size_t num_out_neurons_last_layer = fp->layers[layerno]->num_out_neurons;
    const size_t num_in_neurons_last_layer = fp->layers[layerno]->num_in_neurons;

    const size_t num_in_neurons_first_layer = fp->layers[0]->num_in_neurons;

    std::cout << "num_out_neurons_last " << num_out_neurons_last_layer << std::endl;

    expr_t** expr_array = fp->layers[layerno]->expr_array;

    double* lb_array = fp->layers[layerno]->lb_array;
    double* ub_array = fp->layers[layerno]->ub_array;

    expr_t** lexpr_array;
    expr_t** uexpr_array;

    cudaMalloc((void**) &lexpr_array, num_out_neurons_last_layer*sizeof(expr_t*));
    cudaMalloc((void**) &uexpr_array, num_out_neurons_last_layer*sizeof(expr_t*));

    expr_t** lexpr_array_tmp;
    expr_t** uexpr_array_tmp;

    cudaMalloc((void**) &lexpr_array_tmp, num_out_neurons_last_layer*sizeof(expr_t*));
    cudaMalloc((void**) &uexpr_array_tmp, num_out_neurons_last_layer*sizeof(expr_t*));

    copy_expr_array<<<num_out_neurons_last_layer, 1>>>(lexpr_array, expr_array, num_out_neurons_last_layer, num_in_neurons_last_layer);
    copy_expr_array<<<num_out_neurons_last_layer, 1>>>(uexpr_array, expr_array, num_out_neurons_last_layer, num_in_neurons_last_layer);

    for(int k = layerno - 1; k >= 0; k--)
    {
        const size_t num_out_neurons_current_layer = fp->layers[k]->num_out_neurons;
        const size_t num_in_neurons_current_layer  = fp->layers[k]->num_in_neurons;
        std::cout << "num_out_neurons_current " << num_out_neurons_current_layer << " num_in_neurons_current " << num_in_neurons_current_layer << std::endl;

        const dim3 num_blocks_relu(num_out_neurons_last_layer, num_out_neurons_current_layer/num_threads + 1, 1);
        const dim3 num_blocks_linear(num_out_neurons_last_layer, num_in_neurons_current_layer/num_threads + 1, 1);

        std::cout << "num_threads" << num_threads << " num_blocks_relu " << num_blocks_relu.y << " num_blocks_linear " << num_blocks_linear.y << std::endl;

        expr_t** aux_expr_array = fp->layers[k]->expr_array;

        double* aux_lb_array = fp->layers[k]->lb_array;
        double* aux_ub_array = fp->layers[k]->ub_array;

        if(fp->layers[k]->activation == RELU)
        {
            lexpr_replace_relu_bounds<<<num_blocks_relu, num_threads>>>(lexpr_array, aux_lb_array, aux_ub_array, num_out_neurons_last_layer, num_out_neurons_current_layer);
            uexpr_replace_relu_bounds<<<num_blocks_relu, num_threads>>>(uexpr_array, aux_lb_array, aux_ub_array, num_out_neurons_last_layer, num_out_neurons_current_layer);
        }

        layer_allocate_exprs<<<num_out_neurons_last_layer, 1>>>(lexpr_array_tmp, num_out_neurons_last_layer, num_in_neurons_current_layer);
        layer_allocate_exprs<<<num_out_neurons_last_layer, 1>>>(uexpr_array_tmp, num_out_neurons_last_layer, num_in_neurons_current_layer);

        expr_from_previous_layer<<<num_blocks_linear, num_threads>>>(lexpr_array, lexpr_array_tmp, aux_expr_array, num_out_neurons_last_layer, num_out_neurons_current_layer, num_in_neurons_current_layer);
        expr_from_previous_layer<<<num_blocks_linear, num_threads>>>(uexpr_array, uexpr_array_tmp, aux_expr_array, num_out_neurons_last_layer, num_out_neurons_current_layer, num_in_neurons_current_layer);

        std::swap(lexpr_array, lexpr_array_tmp);
        std::swap(uexpr_array, uexpr_array_tmp);

        free_expr_array<<<num_out_neurons_last_layer, 1>>>(lexpr_array_tmp, num_out_neurons_current_layer);
        free_expr_array<<<num_out_neurons_last_layer, 1>>>(uexpr_array_tmp, num_out_neurons_current_layer);
    }

    compute_lb_from_expr<<<num_out_neurons_last_layer, 1>>>(lb_array, lexpr_array, fp->input_inf, fp->input_sup, num_out_neurons_last_layer, num_in_neurons_first_layer);
    compute_ub_from_expr<<<num_out_neurons_last_layer, 1>>>(ub_array, uexpr_array, fp->input_inf, fp->input_sup, num_out_neurons_last_layer, num_in_neurons_first_layer);

    free_expr_array<<<num_out_neurons_last_layer, 1>>>(lexpr_array, num_out_neurons_last_layer);
    free_expr_array<<<num_out_neurons_last_layer, 1>>>(uexpr_array, num_out_neurons_last_layer);

    cudaFree(lexpr_array);
    cudaFree(uexpr_array);

    cudaFree(lexpr_array_tmp);
    cudaFree(uexpr_array_tmp);

    cudaDeviceSynchronize();

    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl << std::endl;
}


void ffn_handle_intermediate_layer(elina_manager_t* man, elina_abstract0_t* element, const double** weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons, const activation_type_t activation)
{
    fppoly_t* fp = fppoly_of_abstract0(element);
    fppoly_add_new_layer(fp, num_out_neurons, num_in_neurons, FFN, activation);
    expr_t** expr_array = fp->layers[fp->numlayers - 1]->expr_array;

    layer_create_dense_exprs(expr_array, weights, bias, num_out_neurons, num_in_neurons);

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
void print_bounds(double* lb_array, double* ub_array, const size_t num_out_neurons)
{
    size_t i = blockIdx.x;

    if(i < num_out_neurons)
    {
        printf("out inf number %i is: %g\n", i, lb_array[i]);
        printf("out sup number %i is: %g\n", i, ub_array[i]);
    }
}


void ffn_handle_last_layer(elina_manager_t* man, elina_abstract0_t* element, const double** weights, const double* bias, const size_t num_out_neurons, const size_t num_in_neurons, const bool has_activation, const activation_type_t activation)
{
    fppoly_t* fp = fppoly_of_abstract0(element);

    if(has_activation)
    {
        fppoly_add_new_layer(fp, num_out_neurons, num_in_neurons, FFN, activation);
    }
    else
    {
        fppoly_add_new_layer(fp, num_out_neurons, num_in_neurons, FFN, NONE);
    }

    expr_t** expr_array = fp->layers[fp->numlayers - 1]->expr_array;

    double* lb_array = fp->layers[fp->numlayers - 1]->lb_array;
    double* ub_array = fp->layers[fp->numlayers - 1]->ub_array;

    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    layer_create_dense_exprs(expr_array, weights, bias, num_out_neurons, num_in_neurons);

    update_state_using_previous_layers(man, fp, fp->numlayers - 1);

    print_bounds<<<num_out_neurons, 1>>>(lb_array, ub_array, num_out_neurons);
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
void create_sub_expr(expr_t** sub, const size_t index, const elina_dim_t y, const elina_dim_t x)
{

    expr_t* expr = (expr_t*) malloc(sizeof(expr_t));

    expr->inf_cst = 0;
    expr->sup_cst = 0;

    expr->inf_coeff = (double*) malloc(10*sizeof(double));
    expr->sup_coeff = (double*) malloc(10*sizeof(double));

    for(size_t i = 0; i < 10; i++)
    {
        expr->inf_coeff[i] = 0.;
        expr->sup_coeff[i] = 0.;
    }

    expr->inf_coeff[y] = -1.;
    expr->sup_coeff[y] = 1.;

    expr->inf_coeff[x] = 1.;
    expr->sup_coeff[x] = -1.;

    sub[index] = expr;
}


void get_lb_using_previous_layers(elina_manager_t* man, const fppoly_t* const fp)
{
    const size_t numlayers = fp->numlayers;
    fppoly_internal_t* pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

    const size_t num_out_neurons_last_layer = 90;

    const size_t num_in_neurons_first_layer = fp->layers[0]->num_in_neurons;

    double* lb_dev;
    cudaMalloc((void**) &lb_dev, num_out_neurons_last_layer*sizeof(double));

    expr_t** lexpr_array;
    cudaMalloc((void**) &lexpr_array, num_out_neurons_last_layer*sizeof(expr_t*));

    size_t index = 0;

    for(elina_dim_t y = 0; y < 10; y++)
    {
        for(elina_dim_t x = 0; x < 10; x++)
        {
            if(y != x)
            {
                create_sub_expr<<<1, 1>>>(lexpr_array, index, y, x);
                index++;
            }
        }
    }

    expr_t** lexpr_array_tmp;
    cudaMalloc((void**) &lexpr_array_tmp, num_out_neurons_last_layer*sizeof(expr_t*));

    for(int k = numlayers - 1; k >= 0; k--)
    {
        const size_t num_out_neurons_current_layer = fp->layers[k]->num_out_neurons;
        const size_t num_in_neurons_current_layer  = fp->layers[k]->num_in_neurons;

        const dim3 num_blocks_relu(num_out_neurons_last_layer, num_out_neurons_current_layer/num_threads + 1, 1);
        const dim3 num_blocks_linear(num_out_neurons_last_layer, num_in_neurons_current_layer/num_threads + 1, 1);

        expr_t** aux_expr_array = fp->layers[k]->expr_array;

        double* aux_lb_array = fp->layers[k]->lb_array;
        double* aux_ub_array = fp->layers[k]->ub_array;

        if(fp->layers[k]->activation == RELU)
        {
            lexpr_replace_relu_bounds<<<num_blocks_relu, num_threads>>>(lexpr_array, aux_lb_array, aux_ub_array, num_out_neurons_last_layer, num_out_neurons_current_layer);
        }

        layer_allocate_exprs<<<num_out_neurons_last_layer, 1>>>(lexpr_array_tmp, num_out_neurons_last_layer, num_in_neurons_current_layer);

        expr_from_previous_layer<<<num_blocks_linear, num_threads>>>(lexpr_array, lexpr_array_tmp, aux_expr_array, num_out_neurons_last_layer, num_out_neurons_current_layer, num_in_neurons_current_layer);

        std::swap(lexpr_array, lexpr_array_tmp);

        free_expr_array<<<num_out_neurons_last_layer, 1>>>(lexpr_array_tmp, num_out_neurons_current_layer);
    }

    compute_lb_from_expr<<<num_out_neurons_last_layer, 1>>>(lb_dev, lexpr_array, fp->input_inf, fp->input_sup, num_out_neurons_last_layer, num_in_neurons_first_layer);

    double lb[num_out_neurons_last_layer];
    cudaMemcpy(&lb, lb_dev, num_out_neurons_last_layer*sizeof(double), cudaMemcpyDeviceToHost);

    free_expr_array<<<num_out_neurons_last_layer, 1>>>(lexpr_array, num_out_neurons_last_layer);
    cudaFree(lexpr_array);
    cudaFree(lexpr_array_tmp);
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


void device_layer_create_sparse_exprs(expr_t** expr_array, const double* filter_weights, const double* filter_bias,
                               const size_t* input_size, const size_t* output_size, const size_t* filter_size, const size_t* strides,
                               const long int pad_top, const long int pad_left, const size_t num_pixels)
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
                const size_t num_coeff = input_size[2]*filter_size[0]*filter_size[1];

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

    device_layer_create_dense_expr<<<num_out_neurons, 1>>>(expr_array, dense_coeff_dev, bias_dev, num_out_neurons, num_pixels);

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
    expr_t** expr_array = fp->layers[fp->numlayers - 1]->expr_array;

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

    device_layer_create_sparse_exprs(expr_array, filter_weights_tmp, filter_bias_tmp, input_size_tmp, output_size_tmp,
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

    const size_t num_out_neurons = fp->layers[fp->numlayers - 1]->num_out_neurons;
    const size_t num_in_neurons = fp->layers[fp->numlayers - 1]->num_in_neurons;

    expr_t** expr_array = fp->layers[fp->numlayers - 1]->expr_array;

    layer_compute_bounds_from_exprs(expr_array, fp->layers[0]->lb_array, fp->layers[0]->ub_array, fp->input_inf, fp->input_sup, num_out_neurons, num_in_neurons);
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
    free_expr_array<<<layer->num_out_neurons, 1>>>(layer->expr_array, layer->num_out_neurons);

    cudaFree(layer->expr_array);
    layer->expr_array = nullptr;

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

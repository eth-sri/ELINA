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



#ifndef __ZONOML_H_INCLUDED__
#define __ZONOML_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#if defined (HAS_APRON)
#include "apron_wrapper.h"
#else
#include "elina_coeff.h"
#include "elina_dimension.h"
#include "elina_manager.h"
#endif

#include "elina_generic.h"

elina_manager_t* zonoml_manager_alloc(void);

elina_abstract0_t *relu_zono(elina_manager_t* man, bool destructive, elina_abstract0_t * abs, elina_dim_t x);

elina_abstract0_t *relu_zono_refined(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t x, double new_inf, double new_sup);

elina_abstract0_t* maxpool_zono(elina_manager_t *man, bool destructive, elina_abstract0_t *abs, 
			   size_t *pool_size, size_t *input_size, size_t src_offset, size_t* strides, 
			   size_t dimensionality, size_t dst_offset, size_t pad_top, size_t pad_left, size_t *output_size);

elina_abstract0_t *maxpool_zono_refined(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t x, double new_inf, double new_sup);

elina_abstract0_t *handle_gather_layer(elina_manager_t* man, bool destructive, elina_abstract0_t * abs, size_t *indexes);

elina_abstract0_t* ffn_matmult_zono(elina_manager_t * man, bool destructive, elina_abstract0_t* element, elina_dim_t start_offset,
			       double **weights, double * bias,  size_t num_var, size_t expr_offset, size_t expr_size);

elina_abstract0_t* conv_matmult_zono(elina_manager_t* man, bool destructive, elina_abstract0_t* element, elina_dim_t start_offset, double *filter_weights, double * filter_bias,  
				      size_t * input_size, size_t expr_offset, size_t *filter_size, size_t num_filters, size_t *strides, size_t* output_size, size_t pad_top, size_t pad_left, bool has_bias);

bool is_greater(elina_manager_t *man, elina_abstract0_t *elem, elina_dim_t y, elina_dim_t x);

elina_abstract0_t* zonotope_from_network_input(elina_manager_t* man, size_t intdim, size_t realdim, double* inf_array, double * sup_array);

elina_abstract0_t * relu_zono_layerwise(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t start_offset, elina_dim_t num_dim, bool create_new_noise_symbol);

elina_abstract0_t * sigmoid_zono_layerwise(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t start_offset, elina_dim_t num_dim);

elina_abstract0_t * tanh_zono_layerwise(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t start_offset, elina_dim_t num_dim);

elina_abstract0_t * elina_abstract0_from_zonotope(elina_manager_t *man, size_t intdim, size_t realdim, size_t num_error_terms, double **coeffs);

void zono_add(elina_manager_t *man, elina_abstract0_t *elem, size_t dst_offset, size_t src_offset, size_t num_var);

void zono_copy_section(elina_manager_t *man, elina_abstract0_t *elem, size_t dst_offset, size_t src_offset, size_t num_var);

double get_interval_width_var_zono(elina_manager_t *man, elina_abstract0_t *elem, size_t i);

elina_abstract0_t* ffn_matmult_without_bias_zono(elina_manager_t * man, bool destructive, elina_abstract0_t* element, elina_dim_t start_offset,
			       double **weights, size_t num_var, size_t expr_offset, size_t expr_size);

elina_abstract0_t* ffn_add_bias_zono(elina_manager_t * man, bool destructive, elina_abstract0_t* element, elina_dim_t start_offset,
			        double * bias, size_t num_var);

elina_abstract0_t* ffn_sub_bias_zono(elina_manager_t * man, bool destructive, elina_abstract0_t* element, elina_dim_t start_offset,
			        double * bias, bool is_minuend, size_t num_var);

elina_abstract0_t* ffn_mul_bias_zono(elina_manager_t * man, bool destructive, elina_abstract0_t* element, elina_dim_t start_offset,
			        double * bias, size_t num_var);


bool affine_form_is_box(elina_manager_t* man, elina_abstract0_t *abs, elina_dim_t x);

double * get_affine_form_for_dim(elina_manager_t* man, elina_abstract0_t *abs, size_t dim);

static inline long int max(long int a, long int b){
	return a > b ? a : b;
}

#ifdef __cplusplus
 }
#endif

#endif

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



#ifndef __FPPOLY_H_INCLUDED__
#define __FPPOLY_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#if defined (HAS_APRON)
#include "apron_wrapper.h"
#else
#include "elina_coeff.h"
#include "elina_dimension.h"
#include "elina_linexpr0.h"
#include "elina_manager.h"
#endif

#include <fenv.h>
#include <pthread.h>
#include <unistd.h>
//#include <sys/sysinfo.h>
#include "elina_generic.h"
#include "elina_box_meetjoin.h"



typedef struct fppoly_internal_t{
  /* Name of function */
  elina_funid_t funid;

  /* local parameters for the function */
  elina_funopt_t* funopt;
  /* raised when a conversion from/to a user type resulted in an
     overapproximation
  */
  bool conv;
  double min_denormal;
  double ulp;
  /* back pointer to elina_manager*/
  elina_manager_t* man;
}fppoly_internal_t;



typedef enum layertype_t {
  FFN,  /* FFN layer */
  CONV,    /* CONV layer */
  MAXPOOL,   /* MAXPOOL layer */
  LSTM, /* LSTM layer */
  RESIDUAL, /* RESIDUAL layer */
} layertype_t;

typedef enum activation_type_t{
    RELU,
    SIGMOID,
    TANH,
    PARABOLA,  /* Parabolic assignments */
    LOG,       /* Logrithmic assignments */
    NONE,
}activation_type_t;
    
typedef enum exprtype_t{
 DENSE,
 SPARSE,
}exprtype_t;

typedef struct expr_t{
	double *inf_coeff;
	double *sup_coeff;
	double inf_cst;
	double sup_cst;
	exprtype_t type;
	size_t * dim;
    size_t size;
}expr_t;

typedef struct neuron_t{
	double lb;
	double ub;
	expr_t * expr;
	expr_t * lexpr;
	expr_t * uexpr;
}neuron_t;

typedef struct layer_t{
	size_t dims;
	layertype_t type;
        activation_type_t activation;
	neuron_t **neurons;
	double * h_t_inf;
	double * h_t_sup;
	double * c_t_inf;
	double * c_t_sup;
	size_t *predecessors;
}layer_t;

typedef struct output_abstract_t{
	double * output_inf;
	double * output_sup;
	expr_t ** lexpr;
	expr_t ** uexpr;
}output_abstract_t;

typedef struct fppoly_t{
	layer_t ** layers;
	size_t numlayers;
	double *input_inf;
	double * input_sup;
	expr_t ** input_lexpr;
	expr_t ** input_uexpr;
	size_t size;
	size_t num_pixels;
	size_t lstm_index;
	output_abstract_t * out;
}fppoly_t;


typedef struct nn_thread_t{
	size_t start;
	size_t end;
	elina_manager_t *man;
	fppoly_t *fp;
	size_t layerno;
	bool use_area_heuristic;
}nn_thread_t;


elina_manager_t* fppoly_manager_alloc(void);

elina_abstract0_t* fppoly_from_network_input(elina_manager_t *man, size_t intdim, size_t realdim, double *inf_array, double *sup_array);

void fppoly_set_network_input_box(elina_manager_t *man, elina_abstract0_t* element, size_t intdim, size_t realdim, double *inf_array, double * sup_array);
    
elina_abstract0_t* fppoly_from_network_input_poly(elina_manager_t *man, size_t intdim, size_t realdim, double *inf_array, double *sup_array, 
                                                  double * lexpr_weights, double * lexpr_cst, size_t * lexpr_dim, double * uexpr_weights,
						  double * uexpr_cst, size_t * uexpr_dim, size_t expr_size);

void ffn_handle_first_relu_layer(elina_manager_t* man, elina_abstract0_t * abs, double **weights, double *bias,  size_t size, size_t num_pixels, size_t *predecessors);

void ffn_handle_first_sigmoid_layer(elina_manager_t* man, elina_abstract0_t * abs, double **weights, double *bias,  size_t size, size_t num_pixels, size_t *predecessors);
        
void ffn_handle_first_tanh_layer(elina_manager_t* man, elina_abstract0_t * abs, double **weights, double *bias,  size_t size, size_t num_pixels, size_t *predecessors);

void ffn_handle_first_parabola_layer(elina_manager_t* man, elina_abstract0_t * abs, double **weights, double *bias, size_t size, size_t num_pixels, size_t *predecessors);

void ffn_handle_first_log_layer(elina_manager_t* man, elina_abstract0_t * abs, double **weights, double *bias,  size_t size, size_t num_pixels, size_t *predecessors);
    
void ffn_handle_first_relu_layer_no_alloc(elina_manager_t* man, elina_abstract0_t * abs, double **weights, double *bias,  size_t size, size_t num_pixels, size_t *predecessors);
    
void ffn_handle_first_sigmoid_layer_no_alloc(elina_manager_t* man, elina_abstract0_t * abs, double **weights, double *bias,  size_t size, size_t num_pixels, size_t *predecessors);
    
void ffn_handle_first_tanh_layer_no_alloc(elina_manager_t* man, elina_abstract0_t * abs, double **weights, double *bias,  size_t size, size_t num_pixels, size_t *predecessors);
    
void ffn_handle_first_parabola_layer_no_alloc(elina_manager_t* man, elina_abstract0_t * abs, double **weights, double *bias, size_t size, size_t num_pixels, size_t *predecessors);
    
void ffn_handle_first_log_layer_no_alloc(elina_manager_t* man, elina_abstract0_t * abs, double **weights, double *bias,  size_t size, size_t num_pixels, size_t *predecessors);

void ffn_handle_intermediate_affine_layer(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool use_area_heuristic);

void ffn_handle_intermediate_relu_layer(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool use_area_heuristic);
    
void ffn_handle_intermediate_sigmoid_layer(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool use_area_heuristic);
    
void ffn_handle_intermediate_tanh_layer(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool use_area_heuristic);

void ffn_handle_intermediate_parabola_layer(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool use_area_heuristic);

void ffn_handle_intermediate_log_layer(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool use_area_heuristic);
    
void ffn_handle_intermediate_affine_layer_no_alloc(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool use_area_heuristic);
    
void ffn_handle_intermediate_relu_layer_no_alloc(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool use_area_heuristic);
    
void ffn_handle_intermediate_sigmoid_layer_no_alloc(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool use_area_heuristic);
    
void ffn_handle_intermediate_tanh_layer_no_alloc(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool use_area_heuristic);
    
void ffn_handle_intermediate_parabola_layer_no_alloc(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool use_area_heuristic);
    
void ffn_handle_intermediate_log_layer_no_alloc(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool use_area_heuristic);

void fppoly_fprint(FILE* stream, elina_manager_t* man, fppoly_t* fp, char** name_of_dim);

void ffn_handle_last_relu_layer(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias,  size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool has_relu, bool use_area_heuristic);
    
void ffn_handle_last_sigmoid_layer(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias,  size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool has_sigmoid, bool use_area_heuristic);
    
void ffn_handle_last_tanh_layer(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias,  size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool has_tanh, bool use_area_heuristic);

void ffn_handle_last_parabola_layer(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, 
				    size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool has_parabola, bool use_area_heuristic);

void ffn_handle_last_log_layer(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias,  size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool has_log, bool use_area_heuristic);
    
    
void ffn_handle_last_relu_layer_no_alloc(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias,  size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool has_relu, bool use_area_heuristic);
    
void ffn_handle_last_sigmoid_layer_no_alloc(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias,  size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool has_sigmoid, bool use_area_heuristic);
    
void ffn_handle_last_tanh_layer_no_alloc(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias,  size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool has_tanh, bool use_area_heuristic);
    
void ffn_handle_last_parabola_layer_no_alloc(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool has_parabola, bool use_area_heuristic);
    
void ffn_handle_last_log_layer_no_alloc(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias,  size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, bool has_log, bool use_area_heuristic);

void fppoly_free(elina_manager_t *man, fppoly_t *fp);

bool is_greater(elina_manager_t* man, elina_abstract0_t* element, elina_dim_t y, elina_dim_t x, bool use_area_heuristic);

void conv_handle_first_layer(elina_manager_t *man, elina_abstract0_t * element, double *filter_weights, double *filter_bias,  
					  size_t *input_size, size_t *filter_size, size_t num_filters, size_t *strides, bool is_valid_padding, bool has_bias, size_t *predecessors);

void conv_handle_intermediate_relu_layer(elina_manager_t* man, elina_abstract0_t* element, double *filter_weights, double * filter_bias,  
				         size_t * input_size, size_t *filter_size, size_t num_filters, size_t *strides, bool is_valid_padding, bool has_bias, size_t *predecessors, bool use_area_heuristic);

void conv_handle_intermediate_affine_layer(elina_manager_t* man, elina_abstract0_t* element, double *filter_weights, double * filter_bias,  
				         size_t * input_size, size_t *filter_size, size_t num_filters, size_t *strides, bool is_valid_padding, bool has_bias, size_t *predecessors, bool use_area_heuristic);

size_t handle_maxpool_layer(elina_manager_t *man, elina_abstract0_t *abs, 
			   size_t *pool_size, size_t *input_size, size_t *predecessors);

void create_lstm_layer(elina_manager_t *man, elina_abstract0_t *abs, size_t h, size_t *predecessors);

void handle_lstm_layer(elina_manager_t *man, elina_abstract0_t *abs, double **weights,  double *bias, size_t d, size_t h, size_t *predecessors, bool use_area_heuristic);

void fppoly_alloc_first_layer(fppoly_t *fp, size_t size,  layertype_t type, activation_type_t activation);

elina_linexpr0_t * get_lexpr_for_output_neuron(elina_manager_t *man, elina_abstract0_t *abs, size_t i);

elina_linexpr0_t * get_uexpr_for_output_neuron(elina_manager_t *man, elina_abstract0_t *abs, size_t i);

elina_interval_t * box_for_neuron(elina_manager_t* man, elina_abstract0_t * abs, size_t layerno, size_t neuron_no);

elina_interval_t ** box_for_layer(elina_manager_t* man, elina_abstract0_t * abs, size_t layerno);

size_t get_num_neurons_in_layer(elina_manager_t* man, elina_abstract0_t * abs, size_t layerno);

void free_neuron(neuron_t *neuron);

void free_non_lstm_layer_expr(elina_manager_t *man, elina_abstract0_t *abs, size_t layerno);
    
void update_bounds_for_neuron(elina_manager_t *man, elina_abstract0_t *abs, size_t layerno, size_t neuron_no, double lb, double ub);

elina_interval_t * get_bounds_for_linexpr(elina_manager_t *man, elina_abstract0_t *element, elina_linexpr0_t *linexpr0, size_t layerno);

void handle_residual_relu_layer(elina_manager_t *man, elina_abstract0_t *element, size_t num_neurons, size_t *predecessors, bool use_area_heuristic);
void handle_residual_affine_layer(elina_manager_t *man, elina_abstract0_t *element, size_t num_neurons, size_t *predecessors, bool use_area_heuristic);

#ifdef __cplusplus
 }
#endif

#endif


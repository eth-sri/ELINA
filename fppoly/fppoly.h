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



typedef enum fnn_op{
	MATMULT,
        SUB1,
	SUB2,
        MUL,
}fnn_op;


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
	expr_t * lexpr;
	expr_t * uexpr;
}neuron_t;

typedef struct layer_t{
	size_t dims;
	neuron_t **neurons;
	double * h_t_inf;
	double * h_t_sup;
	double * c_t_inf;
	double * c_t_sup;
	size_t *predecessors;
	size_t num_predecessors;
	bool is_activation;
	bool is_concat;
	size_t *C;
	size_t num_channels;
}layer_t;


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
    size_t *spatial_indices;
    size_t *spatial_neighbors;
    size_t spatial_size;
    double spatial_gamma;
}fppoly_t;


typedef struct nn_thread_t{
	size_t start;
	size_t end;
	elina_manager_t *man;
	fppoly_t *fp;
	size_t layerno;
	elina_linexpr0_t ** linexpr0;
	double *res;
}nn_thread_t;


elina_manager_t* fppoly_manager_alloc(void);

elina_abstract0_t* fppoly_from_network_input(elina_manager_t *man, size_t intdim, size_t realdim, double *inf_array, double *sup_array);

void fppoly_set_network_input_box(elina_manager_t *man, elina_abstract0_t* element, size_t intdim, size_t realdim, double *inf_array, double * sup_array);
    
elina_abstract0_t* fppoly_from_network_input_poly(elina_manager_t *man,
        size_t intdim, size_t realdim, double *inf_array, double *sup_array,
        double * lexpr_weights, double * lexpr_cst, size_t * lexpr_dim,
        double * uexpr_weights, double * uexpr_cst, size_t * uexpr_dim,
        size_t expr_size, size_t * spatial_indices, size_t * spatial_neighbors,
        size_t spatial_size, double spatial_gamma);

fppoly_internal_t* fppoly_init_from_manager(elina_manager_t* man, elina_funid_t funid);


void handle_fully_connected_layer(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, size_t num_predecessors);

void handle_fully_connected_layer_no_alloc(elina_manager_t* man, elina_abstract0_t* element, double **weights, double * bias, size_t num_out_neurons, size_t num_in_neurons, size_t *predecessors, size_t num_predecessors);

void handle_sub_layer(elina_manager_t* man, elina_abstract0_t* element,  double * cst, bool is_minuend, size_t num_in_neurons, size_t *predecessors, size_t num_predecessors);

void handle_mul_layer(elina_manager_t* man, elina_abstract0_t* element,  double * bias,  size_t num_in_neurons, size_t *predecessors, size_t num_predecessors);

void handle_relu_layer(elina_manager_t *man, elina_abstract0_t* element, size_t num_neurons, size_t *predecessors, size_t num_predecessors, bool use_default_heuristics);

void handle_sign_layer(elina_manager_t *man, elina_abstract0_t* element, size_t num_neurons, size_t *predecessors, size_t num_predecessors);

void handle_sigmoid_layer(elina_manager_t *man, elina_abstract0_t* element, size_t num_neurons, size_t *predecessors, size_t num_predecessors);

void handle_tanh_layer(elina_manager_t *man, elina_abstract0_t* element, size_t num_neurons, size_t *predecessors, size_t num_predecessors);

void handle_parabola_layer(elina_manager_t *man, elina_abstract0_t* element, size_t num_neurons, size_t *predecessors, size_t num_predecessors);

void handle_log_layer(elina_manager_t *man, elina_abstract0_t* element, size_t num_neurons, size_t *predecessors, size_t num_predecessors);



void fppoly_fprint(FILE* stream, elina_manager_t* man, fppoly_t* fp, char** name_of_dim);


void fppoly_free(elina_manager_t *man, fppoly_t *fp);

bool is_greater(elina_manager_t* man, elina_abstract0_t* element, elina_dim_t y, elina_dim_t x);



void handle_convolutional_layer(elina_manager_t* man, elina_abstract0_t* element, double *filter_weights, double * filter_bias,  
				         size_t * input_size, size_t *filter_size, size_t num_filters, size_t *strides, size_t *output_size, size_t pad_top, size_t pad_left, bool has_bias, size_t *predecessors, size_t num_predecessors);

/*void conv_handle_intermediate_affine_layer(elina_manager_t* man, elina_abstract0_t* element, double *filter_weights, double * filter_bias,  
				         size_t * input_size, size_t *filter_size, size_t num_filters, size_t *strides, size_t *output_size, size_t pad_top, size_t pad_left, bool has_bias, size_t *predecessors, bool use_area_heuristic);*/

size_t handle_pool_layer(elina_manager_t *man, elina_abstract0_t *element,
                           size_t *pool_size, size_t *input_size, size_t *strides, size_t pad_top, size_t pad_left, size_t * output_size, size_t *predecessors, size_t num_predecessors, bool is_maxpool);

void create_lstm_layer(elina_manager_t *man, elina_abstract0_t *abs, size_t h, size_t *predecessors, size_t num_predecessors);

void handle_lstm_layer(elina_manager_t *man, elina_abstract0_t *abs, double **weights,  double *bias, size_t d, size_t h, size_t *predecessors, size_t num_predecessors);


neuron_t *neuron_alloc(void);

elina_linexpr0_t * get_lexpr_for_output_neuron(elina_manager_t *man, elina_abstract0_t *abs, size_t i);

elina_linexpr0_t * get_uexpr_for_output_neuron(elina_manager_t *man, elina_abstract0_t *abs, size_t i);

elina_interval_t * box_for_neuron(elina_manager_t* man, elina_abstract0_t * abs, size_t layerno, size_t neuron_no);

elina_interval_t ** box_for_layer(elina_manager_t* man, elina_abstract0_t * abs, size_t layerno);

size_t get_num_neurons_in_layer(elina_manager_t* man, elina_abstract0_t * abs, size_t layerno);

void free_neuron(neuron_t *neuron);

void free_non_lstm_layer_expr(elina_manager_t *man, elina_abstract0_t *abs, size_t layerno);
    
void update_bounds_for_neuron(elina_manager_t *man, elina_abstract0_t *abs, size_t layerno, size_t neuron_no, double lb, double ub);

double* get_upper_bound_for_linexpr(elina_manager_t *man, elina_abstract0_t *element, elina_linexpr0_t **linexpr0, size_t size, size_t layerno);

void handle_residual_layer(elina_manager_t *man, elina_abstract0_t *element, size_t num_neurons, size_t *predecessors, size_t num_predecessors);

void handle_concatenation_layer(elina_manager_t* man, elina_abstract0_t* element, size_t * predecessors, size_t num_predecessors, size_t *C);

void handle_tiling_layer(elina_manager_t* man, elina_abstract0_t* element, size_t * predecessors, size_t num_predecessors, size_t repeat);
//void handle_residual_affine_layer(elina_manager_t *man, elina_abstract0_t *element, size_t num_neurons, size_t *predecessors, bool use_area_heuristic);

fppoly_t* fppoly_of_abstract0(elina_abstract0_t* a);

void fppoly_add_new_layer(fppoly_t *fp, size_t size, size_t *predecessors, size_t num_predecessors, bool is_activation);

void update_activation_upper_bound_for_neuron(elina_manager_t *man, elina_abstract0_t *abs, size_t layerno, size_t neuron_no, double* coeff, size_t *dim, size_t size);

void update_activation_lower_bound_for_neuron(elina_manager_t *man, elina_abstract0_t *abs, size_t layerno, size_t neuron_no, double* coeff, size_t *dim, size_t size);

elina_linexpr0_t *get_output_lexpr_defined_over_previous_layers(elina_manager_t *man, elina_abstract0_t *element, int neuron_no, int prev_layer);

elina_linexpr0_t *get_output_uexpr_defined_over_previous_layers(elina_manager_t *man, elina_abstract0_t *element, int neuron_no, int prev_layer);

#ifdef __cplusplus
 }
#endif

#endif


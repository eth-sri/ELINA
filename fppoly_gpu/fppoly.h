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

#ifndef __FPPOLY_H_INCLUDED__
#define __FPPOLY_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#if defined(HAS_APRON)
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
#include "elina_box_meetjoin.h"
#include "elina_generic.h"

typedef struct fppoly_internal_t {
  /* Name of function */
  elina_funid_t funid;

  /* local parameters for the function */
  elina_funopt_t *funopt;
  /* raised when a conversion from/to a user type resulted in an
     overapproximation
  */
  bool conv;
  double min_denormal;
  double ulp;
  /* back pointer to elina_manager*/
  elina_manager_t *man;
} fppoly_internal_t;

typedef enum layertype_t {
  FFN,     /* FFN layer */
  CONV,    /* CONV layer */
  MAXPOOL, /* MAXPOOL layer */
} layertype_t;

typedef enum activation_type_t {
  RELU,
  SIGMOID,
  TANH,
  NONE,
} activation_type_t;

typedef enum exprtype_t {
  DENSE,
  SPARSE,
} exprtype_t;

typedef struct expr_t {
  double *inf_coeff;
  double *sup_coeff;
  double inf_cst;
  double sup_cst;
  exprtype_t type;
  size_t *dim;
  size_t size;
} expr_t;

typedef struct neuron_t {
  double lb;
  double ub;
  expr_t *expr;
  expr_t *maxpool_lexpr;
  expr_t *maxpool_uexpr;
} neuron_t;

typedef struct layer_t {
  size_t dims;
  layertype_t type;
  activation_type_t activation;
  neuron_t **neurons;
} layer_t;

typedef struct output_abstract_t {
  double *output_inf;
  double *output_sup;
  expr_t **lexpr;
  expr_t **uexpr;
} output_abstract_t;

typedef struct fppoly_t {
  layer_t **layers;
  size_t numlayers;
  double *input_inf;
  double *input_sup;
  expr_t **input_lexpr;
  expr_t **input_uexpr;
  size_t size;
  size_t num_pixels;
  output_abstract_t *out;
} fppoly_t;

typedef struct nn_thread_t {
  size_t start;
  size_t end;
  elina_manager_t *man;
  fppoly_t *fp;
  size_t layerno;
} nn_thread_t;

elina_manager_t *fppoly_manager_alloc(void);

elina_abstract0_t *fppoly_from_network_input(elina_manager_t *man,
                                             size_t intdim, size_t realdim,
                                             double *inf_array,
                                             double *sup_array);

elina_abstract0_t *fppoly_from_network_input_poly(
    elina_manager_t *man, size_t intdim, size_t realdim, double *inf_array,
    double *sup_array, double *lexpr_weights, double *lexpr_cst,
    size_t *lexpr_dim, double *uexpr_weights, double *uexpr_cst,
    size_t *uexpr_dim, size_t expr_size);

void ffn_handle_first_relu_layer(elina_manager_t *man, elina_abstract0_t *abs,
                                 double **weights, double *bias, size_t size,
                                 size_t num_pixels);

void ffn_handle_first_sigmoid_layer(elina_manager_t *man,
                                    elina_abstract0_t *abs, double **weights,
                                    double *bias, size_t size,
                                    size_t num_pixels);

void ffn_handle_first_tanh_layer(elina_manager_t *man, elina_abstract0_t *abs,
                                 double **weights, double *bias, size_t size,
                                 size_t num_pixels);

void ffn_handle_intermediate_relu_layer(elina_manager_t *man,
                                        elina_abstract0_t *element,
                                        double **weights, double *bias,
                                        size_t num_out_neurons,
                                        size_t num_in_neurons);

void ffn_handle_intermediate_sigmoid_layer(elina_manager_t *man,
                                           elina_abstract0_t *element,
                                           double **weights, double *bias,
                                           size_t num_out_neurons,
                                           size_t num_in_neurons);

void ffn_handle_intermediate_tanh_layer(elina_manager_t *man,
                                        elina_abstract0_t *element,
                                        double **weights, double *bias,
                                        size_t num_out_neurons,
                                        size_t num_in_neurons);

void fppoly_fprint(FILE *stream, elina_manager_t *man, fppoly_t *fp,
                   char **name_of_dim);

void ffn_handle_last_relu_layer(elina_manager_t *man,
                                elina_abstract0_t *element, double **weights,
                                double *bias, size_t num_out_neurons,
                                size_t num_in_neurons, bool has_relu);

void ffn_handle_last_sigmoid_layer(elina_manager_t *man,
                                   elina_abstract0_t *element, double **weights,
                                   double *bias, size_t num_out_neurons,
                                   size_t num_in_neurons, bool has_sigmoid);

void ffn_handle_last_tanh_layer(elina_manager_t *man,
                                elina_abstract0_t *element, double **weights,
                                double *bias, size_t num_out_neurons,
                                size_t num_in_neurons, bool has_tanh);

void fppoly_free(elina_manager_t *man, fppoly_t *fp);

bool is_greater(elina_manager_t *man, elina_abstract0_t *element, elina_dim_t y,
                elina_dim_t x);

void conv_handle_first_layer(elina_manager_t *man, elina_abstract0_t *element,
                             double *filter_weights, double *filter_bias,
                             size_t *input_size, size_t *filter_size,
                             size_t num_filters, size_t *strides,
                             bool is_valid_padding, bool has_bias);

void conv_handle_intermediate_relu_layer(
    elina_manager_t *man, elina_abstract0_t *element, double *filter_weights,
    double *filter_bias, size_t *input_size, size_t *filter_size,
    size_t num_filters, size_t *strides, bool is_valid_padding, bool has_bias);

size_t handle_maxpool_layer(elina_manager_t *man, elina_abstract0_t *abs,
                            size_t *pool_size, size_t *input_size);

void fppoly_alloc_first_layer(fppoly_t *fp, size_t size, size_t num_pixels,
                              layertype_t type, activation_type_t activation);

elina_linexpr0_t *get_lexpr_for_output_neuron(elina_manager_t *man,
                                              elina_abstract0_t *abs, size_t i);

elina_linexpr0_t *get_uexpr_for_output_neuron(elina_manager_t *man,
                                              elina_abstract0_t *abs, size_t i);

#ifdef __cplusplus
}
#endif

#endif

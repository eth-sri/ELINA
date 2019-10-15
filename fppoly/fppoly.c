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

#include "backsubstitute.h"



fppoly_t* fppoly_of_abstract0(elina_abstract0_t* a)
{
  return (fppoly_t*)a->value;
}

elina_abstract0_t* abstract0_of_fppoly(elina_manager_t* man, fppoly_t* fp)
{
  elina_abstract0_t* r = malloc(sizeof(elina_abstract0_t));
  assert(r);
  r->value = fp;
  r->man = elina_manager_copy(man);
  return r;
}


static inline void fppoly_internal_free(fppoly_internal_t* pr)
{
    if (pr) {
	pr->funid = ELINA_FUNID_UNKNOWN;
	free(pr);
	pr = NULL;
    }
}


static inline fppoly_internal_t* fppoly_internal_alloc(void)
{
    fppoly_internal_t* pr = (fppoly_internal_t*)malloc(sizeof(fppoly_internal_t));
    pr->funid = ELINA_FUNID_UNKNOWN;
    pr->man = NULL;
    pr->funopt = NULL; 
    pr->min_denormal = ldexpl(1.0,-1074);
    pr->ulp = ldexpl(1.0,-52);
    return pr;
}


/* back pointer to our internal structure from the manager */
fppoly_internal_t* fppoly_init_from_manager(elina_manager_t* man, elina_funid_t funid)
{
	
    fppoly_internal_t* pr = (fppoly_internal_t*)man->internal;
    pr->funid = funid;
	
    if (!(pr->man)) pr->man = man;
	
    return pr;
}


elina_manager_t * fppoly_manager_alloc(void){
	void** funptr;
	fesetround(FE_UPWARD);
	fppoly_internal_t *pr = fppoly_internal_alloc();
	elina_manager_t *man = elina_manager_alloc("fppoly",/* Library name */
			"1.0", /* version */
			pr, /* internal structure */
			(void (*)(void*))fppoly_internal_free /* free function for internal */
			);
	funptr = man->funptr;
	funptr[ELINA_FUNID_FREE] = &fppoly_free;
	/* 3.Printing */
	funptr[ELINA_FUNID_FPRINT] = &fppoly_fprint;
	return man;
}


neuron_t *neuron_alloc(void){
	neuron_t *res =  (neuron_t *)malloc(sizeof(neuron_t));
        res->expr = NULL;
        res->lb = -INFINITY;
	res->ub = INFINITY;
	res->lexpr = NULL;
	res->uexpr = NULL;
	return res;
}

layer_t *create_layer(size_t size, layertype_t type,
                      activation_type_t activation) {
  layer_t *layer = (layer_t *)malloc(sizeof(layer_t));
  layer->dims = size;
  layer->type = type;
  layer->activation = activation;
  layer->neurons = (neuron_t **)malloc(size * sizeof(neuron_t *));
  size_t i;
  for (i = 0; i < size; i++) {
    layer->neurons[i] = neuron_alloc();
  }
  layer->h_t_inf = NULL;
  layer->h_t_sup = NULL;
  layer->c_t_inf = NULL;
  layer->c_t_sup = NULL;
  return layer;
}

void fppoly_from_network_input_box(fppoly_t *res, size_t intdim, size_t realdim, double *inf_array, double *sup_array){
	
	res->layers = NULL;
	res->numlayers = 0;
	res->lstm_index = 0;
	size_t num_pixels = intdim + realdim;
	res->input_inf = (double *)malloc(num_pixels*sizeof(double));
	res->input_sup = (double *)malloc(num_pixels*sizeof(double));
	res->input_lexpr = NULL;
	res->input_uexpr = NULL;
	size_t i;
	for(i=0; i < num_pixels; i++){
		res->input_inf[i] = -inf_array[i];
		res->input_sup[i] = sup_array[i];
	}
	res->num_pixels = num_pixels;
        res->out = NULL;
}


elina_abstract0_t * fppoly_from_network_input(elina_manager_t *man, size_t intdim, size_t realdim, double *inf_array, double *sup_array){
	fppoly_t * res = (fppoly_t *)malloc(sizeof(fppoly_t));
	fppoly_from_network_input_box(res, intdim, realdim, inf_array, sup_array);
	return abstract0_of_fppoly(man,res);
}

void fppoly_set_network_input_box(elina_manager_t *man, elina_abstract0_t* element, size_t intdim, size_t realdim, double *inf_array, double * sup_array){
    fppoly_t * res = fppoly_of_abstract0(element);
    size_t num_pixels = intdim + realdim;
    res->numlayers = 0;
    size_t i;
    for(i=0; i < num_pixels; i++){
        res->input_inf[i] = -inf_array[i];
        res->input_sup[i] = sup_array[i];
    }
}

elina_abstract0_t *fppoly_from_network_input_poly(
    elina_manager_t *man, size_t intdim, size_t realdim, double *inf_array,
    double *sup_array, double *lexpr_weights, double *lexpr_cst,
    size_t *lexpr_dim, double *uexpr_weights, double *uexpr_cst,
    size_t *uexpr_dim, size_t expr_size) {
  fppoly_t *res = (fppoly_t *)malloc(sizeof(fppoly_t));

  fppoly_from_network_input_box(res, intdim, realdim, inf_array, sup_array);
  size_t num_pixels = intdim + realdim;
  res->input_lexpr = (expr_t **)malloc(num_pixels * sizeof(expr_t *));
  res->input_uexpr = (expr_t **)malloc(num_pixels * sizeof(expr_t *));

  size_t i;
  double *tmp_weights = (double *)malloc(expr_size * sizeof(double));
  size_t *tmp_dim = (size_t *)malloc(expr_size * sizeof(size_t));

  for (i = 0; i < num_pixels; i++) {

    size_t j;
    for (j = 0; j < expr_size; j++) {
      tmp_weights[j] = lexpr_weights[i * expr_size + j];
      tmp_dim[j] = lexpr_dim[i * expr_size + j];
    }
    res->input_lexpr[i] =
        create_sparse_expr(tmp_weights, lexpr_cst[i], tmp_dim, expr_size);
    sort_sparse_expr(res->input_lexpr[i]);
    // printf("w: %p %g %g %g cst: %g dim: %p %zu %zu
    // %zu\n",lexpr_weights[i],lexpr_weights[i][0],lexpr_weights[i][1],
    // lexpr_weights[i][2],lexpr_cst[i],lexpr_dim[i],lexpr_dim[i][0],lexpr_dim[i][1],
    // lexpr_dim[i][2]); expr_print(res->input_lexpr[i]); fflush(stdout);
    for (j = 0; j < expr_size; j++) {
      tmp_weights[j] = uexpr_weights[i * expr_size + j];
      tmp_dim[j] = uexpr_dim[i * expr_size + j];
    }
    res->input_uexpr[i] =
        create_sparse_expr(tmp_weights, uexpr_cst[i], tmp_dim, expr_size);
    sort_sparse_expr(res->input_uexpr[i]);
    //	expr_print(res->input_uexpr[i]);
    //	fflush(stdout);
  }
  free(tmp_weights);
  free(tmp_dim);
  return abstract0_of_fppoly(man, res);
}

void fppoly_alloc_first_layer(fppoly_t *fp, size_t size, layertype_t type,
                              activation_type_t activation) {
  layer_t *layer = create_layer(size, type, activation);
  fp->layers = (layer_t **)malloc(2000 * sizeof(layer_t *));
  fp->layers[0] = layer;
  fp->numlayers = 1;
  return;
}

void fppoly_add_new_layer(fppoly_t *fp, size_t size, layertype_t type,
                          activation_type_t activation) {
  size_t numlayers = fp->numlayers;
  fp->layers[numlayers] = create_layer(size, type, activation);
  fp->numlayers++;
  return;
}

void ffn_handle_first_layer(elina_manager_t *man, elina_abstract0_t *abs,
                            double **weights, double *bias, size_t size,
                            size_t num_pixels, size_t *predecessors,
                            activation_type_t activation, bool alloc) {
  // printf("start \n");
  // fflush(stdout);
  fppoly_t *res = fppoly_of_abstract0(abs);
  // printf("coming here\n");
  // fflush(stdout);

  if (alloc) {
    fppoly_alloc_first_layer(res, size, FFN, activation);
  }
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  size_t i, j;

  // for(i=0; i < num_pixels; i++){
  //	elina_interval_print(itv[i]);
  //	printf("\n");
  //}
  // fflush(stdout);
  neuron_t **neurons = res->layers[0]->neurons;
  res->layers[0]->predecessors = predecessors;

  for (i = 0; i < size; i++) {
    neuron_t *neuron = neurons[i];
    double *weight_i = weights[i];
    double bias_i = bias[i];
    neuron->expr = create_dense_expr(weight_i, bias_i, num_pixels);
    neuron->lb = compute_lb_from_expr(pr, neuron->expr, res, -1);
    neuron->ub = compute_ub_from_expr(pr, neuron->expr, res, -1);
  }

  // printf("return here\n");
  // fppoly_fprint(stdout,man,res,NULL);
  // fflush(stdout);
  return;
}

void ffn_handle_first_relu_layer(elina_manager_t *man, elina_abstract0_t *abs,
                                 double **weights, double *bias, size_t size,
                                 size_t num_pixels, size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, RELU, true);
}

void ffn_handle_first_sigmoid_layer(elina_manager_t *man,
                                    elina_abstract0_t *abs, double **weights,
                                    double *bias, size_t size,
                                    size_t num_pixels, size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, SIGMOID, true);
}

void ffn_handle_first_tanh_layer(elina_manager_t *man, elina_abstract0_t *abs,
                                 double **weights, double *bias, size_t size,
                                 size_t num_pixels, size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, TANH, true);
}

void ffn_handle_first_parabola_layer(elina_manager_t *man,
                                     elina_abstract0_t *abs, double **weights,
                                     double *bias, size_t size,
                                     size_t num_pixels, size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, PARABOLA, true);
}

void ffn_handle_first_log_layer(elina_manager_t *man, elina_abstract0_t *abs,
                                double **weights, double *bias, size_t size,
                                size_t num_pixels, size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, LOG, true);
}

void ffn_handle_first_relu_layer_no_alloc(elina_manager_t *man,
                                          elina_abstract0_t *abs,
                                          double **weights, double *bias,
                                          size_t size, size_t num_pixels,
                                          size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, RELU, false);
}

void ffn_handle_first_sigmoid_layer_no_alloc(elina_manager_t *man,
                                             elina_abstract0_t *abs,
                                             double **weights, double *bias,
                                             size_t size, size_t num_pixels,
                                             size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, SIGMOID, false);
}

void ffn_handle_first_tanh_layer_no_alloc(elina_manager_t *man,
                                          elina_abstract0_t *abs,
                                          double **weights, double *bias,
                                          size_t size, size_t num_pixels,
                                          size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, TANH, false);
}

void ffn_handle_first_parabola_layer_no_alloc(elina_manager_t *man,
                                              elina_abstract0_t *abs,
                                              double **weights, double *bias,
                                              size_t size, size_t num_pixels,
                                              size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, PARABOLA, false);
}

void ffn_handle_first_log_layer_no_alloc(elina_manager_t *man,
                                         elina_abstract0_t *abs,
                                         double **weights, double *bias,
                                         size_t size, size_t num_pixels,
                                         size_t *predecessors) {
  ffn_handle_first_layer(man, abs, weights, bias, size, num_pixels,
                         predecessors, LOG, false);
}

void ffn_handle_intermediate_layer(elina_manager_t *man,
                                   elina_abstract0_t *element, double **weights,
                                   double *bias, size_t num_out_neurons,
                                   size_t num_in_neurons, size_t *predecessors,
                                   activation_type_t activation, bool alloc,
                                   bool use_area_heuristic) {
  // printf("ReLU start here %zu %zu\n",num_in_neurons,num_out_neurons);
  // fflush(stdout);
  fppoly_t *fp = fppoly_of_abstract0(element);
  size_t numlayers = fp->numlayers;
  if (alloc) {
    fppoly_add_new_layer(fp, num_out_neurons, FFN, activation);
  }
  neuron_t **out_neurons = fp->layers[numlayers]->neurons;
  fp->layers[numlayers]->predecessors = predecessors;

  size_t i;
  for (i = 0; i < num_out_neurons; i++) {
    double *weight_i = weights[i];
    double bias_i = bias[i];

    out_neurons[i]->expr = create_dense_expr(weight_i, bias_i, num_in_neurons);
  }
  update_state_using_previous_layers_parallel(man, fp, numlayers,
                                              use_area_heuristic);

  // printf("return here2\n");
  // fppoly_fprint(stdout,man,fp,NULL);
  // fflush(stdout);
  return;
}

void ffn_handle_intermediate_affine_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, NONE, true,
                                use_area_heuristic);
}

void ffn_handle_intermediate_relu_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, RELU, true,
                                use_area_heuristic);
}

void ffn_handle_intermediate_sigmoid_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, SIGMOID, true,
                                use_area_heuristic);
}

void ffn_handle_intermediate_tanh_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, TANH, true,
                                use_area_heuristic);
}

void ffn_handle_intermediate_parabola_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {

  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, PARABOLA, true,
                                use_area_heuristic);
}

void ffn_handle_intermediate_log_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, LOG, true,
                                use_area_heuristic);
}

void ffn_handle_intermediate_affine_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, NONE, false,
                                use_area_heuristic);
}

void ffn_handle_intermediate_relu_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, RELU, false,
                                use_area_heuristic);
}

void ffn_handle_intermediate_sigmoid_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, SIGMOID, false,
                                use_area_heuristic);
}

void ffn_handle_intermediate_tanh_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, TANH, false,
                                use_area_heuristic);
}

void ffn_handle_intermediate_parabola_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {

  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, PARABOLA, false,
                                use_area_heuristic);
}

void ffn_handle_intermediate_log_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool use_area_heuristic) {
  ffn_handle_intermediate_layer(man, element, weights, bias, num_out_neurons,
                                num_in_neurons, predecessors, LOG, false,
                                use_area_heuristic);
}

void handle_final_relu_layer(fppoly_internal_t *pr, output_abstract_t *out,
                             neuron_t **neurons, size_t size, bool has_relu) {
  size_t i;
  if (has_relu) {
    for (i = 0; i < size; i++) {
      out->output_inf[i] = apply_relu_lexpr(pr, &out->lexpr[i], neurons[i]);
      out->output_sup[i] = apply_relu_uexpr(pr, &out->uexpr[i], neurons[i]);
    }
  } else {
    for (i = 0; i < size; i++) {
      out->output_inf[i] = neurons[i]->lb;
      out->output_sup[i] = neurons[i]->ub;
    }
  }
}

void handle_final_sigmoid_layer(fppoly_internal_t *pr, output_abstract_t *out,
                                neuron_t **neurons, size_t size,
                                bool has_sigmoid) {
  size_t i;
  if (has_sigmoid) {
    for (i = 0; i < size; i++) {
      out->output_inf[i] = apply_sigmoid_lexpr(pr, &out->lexpr[i], neurons[i]);
      out->output_sup[i] = apply_sigmoid_uexpr(pr, &out->uexpr[i], neurons[i]);
    }
  } else {
    for (i = 0; i < size; i++) {
      out->output_inf[i] = neurons[i]->lb;
      out->output_sup[i] = neurons[i]->ub;
    }
  }
}

void handle_final_tanh_layer(fppoly_internal_t *pr, output_abstract_t *out,
                             neuron_t **neurons, size_t size, bool has_tanh) {
  size_t i;
  if (has_tanh) {
    for (i = 0; i < size; i++) {
      out->output_inf[i] = apply_tanh_lexpr(pr, &out->lexpr[i], neurons[i]);
      out->output_sup[i] = apply_tanh_uexpr(pr, &out->uexpr[i], neurons[i]);
    }
  } else {
    for (i = 0; i < size; i++) {
      out->output_inf[i] = neurons[i]->lb;
      out->output_sup[i] = neurons[i]->ub;
    }
  }
}

void neuron_fprint(FILE * stream, neuron_t *neuron, char ** name_of_dim){
	//expr_fprint(stream,neuron->expr);
	fprintf(stream,"[%g, %g]\n",-neuron->lb,neuron->ub);
}

void layer_fprint(FILE * stream, layer_t * layer, char** name_of_dim){
	size_t dims = layer->dims;
	size_t i;
	for(i = 0; i < dims; i++){
		fprintf(stream,"neuron: %zu ", i);
		neuron_fprint(stream, layer->neurons[i], name_of_dim);
	}
}

void ffn_handle_last_layer(elina_manager_t *man, elina_abstract0_t *element,
                           double **weights, double *bias,
                           size_t num_out_neurons, size_t num_in_neurons,
                           size_t *predecessors, bool has_activation,
                           activation_type_t activation, bool alloc,
                           bool use_area_heuristic) {
  // printf("last\n");
  // fflush(stdout);
  fppoly_t *fp = fppoly_of_abstract0(element);
  size_t numlayers = fp->numlayers;
  if (alloc) {
    if (has_activation) {
      fppoly_add_new_layer(fp, num_out_neurons, FFN, activation);
    } else {
      fppoly_add_new_layer(fp, num_out_neurons, FFN, NONE);
    }
  }
  output_abstract_t *out =
      (output_abstract_t *)malloc(sizeof(output_abstract_t));
  out->output_inf = (double *)malloc(num_out_neurons * sizeof(double));
  out->output_sup = (double *)malloc(num_out_neurons * sizeof(double));
  out->lexpr = (expr_t **)malloc(num_out_neurons * sizeof(expr_t *));
  out->uexpr = (expr_t **)malloc(num_out_neurons * sizeof(expr_t *));
  fp->out = out;
  neuron_t **out_neurons = fp->layers[numlayers]->neurons;
  fp->layers[numlayers]->predecessors = predecessors;
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  size_t i;
  for (i = 0; i < num_out_neurons; i++) {
    double *weight_i = weights[i];
    double bias_i = bias[i];
    out_neurons[i]->expr = create_dense_expr(weight_i, bias_i, num_in_neurons);
  }
  update_state_using_previous_layers_parallel(man, fp, numlayers,
                                              use_area_heuristic);
  if (activation == RELU) {
    handle_final_relu_layer(pr, fp->out, out_neurons, num_out_neurons,
                            has_activation);
  } else if (activation == SIGMOID) {
    handle_final_sigmoid_layer(pr, fp->out, out_neurons, num_out_neurons,
                               has_activation);
  } else if (activation == TANH) {
    handle_final_tanh_layer(pr, fp->out, out_neurons, num_out_neurons,
                            has_activation);
  }

  else {
    for (i = 0; i < num_out_neurons; i++) {
      out->output_inf[i] = out_neurons[i]->lb;
      out->output_sup[i] = out_neurons[i]->ub;
    }
  }
  // printf("finish\n");
  // fppoly_fprint(stdout,man,fp,NULL);
  // fflush(stdout);
  // layer_fprint(stdout,fp->layers[1],NULL);
  // layer_fprint(stdout,fp->layers[fp->numlayers-1],NULL);
  return;
}

void ffn_handle_last_relu_layer(elina_manager_t *man,
                                elina_abstract0_t *element, double **weights,
                                double *bias, size_t num_out_neurons,
                                size_t num_in_neurons, size_t *predecessors,
                                bool has_relu, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_relu, RELU, true,
                        use_area_heuristic);
}

void ffn_handle_last_sigmoid_layer(elina_manager_t *man,
                                   elina_abstract0_t *element, double **weights,
                                   double *bias, size_t num_out_neurons,
                                   size_t num_in_neurons, size_t *predecessors,
                                   bool has_sigmoid, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_sigmoid, SIGMOID,
                        true, use_area_heuristic);
}

void ffn_handle_last_tanh_layer(elina_manager_t *man,
                                elina_abstract0_t *element, double **weights,
                                double *bias, size_t num_out_neurons,
                                size_t num_in_neurons, size_t *predecessors,
                                bool has_tanh, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_tanh, TANH, true,
                        use_area_heuristic);
}

void ffn_handle_last_parabola_layer(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool has_parabola, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_parabola, PARABOLA,
                        true, use_area_heuristic);
}

void ffn_handle_last_log_layer(elina_manager_t *man, elina_abstract0_t *element,
                               double **weights, double *bias,
                               size_t num_out_neurons, size_t num_in_neurons,
                               size_t *predecessors, bool has_log,
                               bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_log, LOG, true,
                        use_area_heuristic);
}

void ffn_handle_last_relu_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool has_relu, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_relu, RELU, false,
                        use_area_heuristic);
}

void ffn_handle_last_sigmoid_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool has_sigmoid, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_sigmoid, SIGMOID,
                        false, use_area_heuristic);
}

void ffn_handle_last_tanh_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool has_tanh, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_tanh, TANH, false,
                        use_area_heuristic);
}

void ffn_handle_last_parabola_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool has_parabola, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_parabola, PARABOLA,
                        false, use_area_heuristic);
}

void ffn_handle_last_log_layer_no_alloc(
    elina_manager_t *man, elina_abstract0_t *element, double **weights,
    double *bias, size_t num_out_neurons, size_t num_in_neurons,
    size_t *predecessors, bool has_log, bool use_area_heuristic) {
  ffn_handle_last_layer(man, element, weights, bias, num_out_neurons,
                        num_in_neurons, predecessors, has_log, LOG, false,
                        use_area_heuristic);
}

void coeff_to_interval(elina_coeff_t *coeff, double *inf, double *sup){
	double d;
	if(coeff->discr==ELINA_COEFF_SCALAR){
		elina_scalar_t * scalar = coeff->val.scalar;
		d = scalar->val.dbl;
		*inf = -d;
		*sup = d;
	}
	else{
		elina_interval_t *interval = coeff->val.interval;
		d = interval->inf->val.dbl;
		*inf = -d;
		d = interval->sup->val.dbl;
		*sup = d;	
	}
		
}

expr_t * elina_linexpr0_to_expr(elina_linexpr0_t *linexpr0){
	size_t size = linexpr0->size;
	size_t i;
	expr_t *res = (expr_t*)malloc(sizeof(expr_t));
	res->inf_coeff = (double*)malloc(size*sizeof(double));
	res->sup_coeff = (double*)malloc(size*sizeof(double));
	res->size = size;
	if(linexpr0->discr==ELINA_LINEXPR_SPARSE){
		res->type = SPARSE;
		res->dim = (size_t *)malloc(size*sizeof(size_t));
	}
	else{
		res->type = DENSE;
		res->dim = NULL;
	}
	size_t k;
	for(i=0; i< size; i++){
		elina_coeff_t *coeff;
		if(res->type==SPARSE){
			k = linexpr0->p.linterm[i].dim;
			res->dim[i] = k;
			coeff = &linexpr0->p.linterm[i].coeff;
			coeff_to_interval(coeff,&res->inf_coeff[i],&res->sup_coeff[i]);
		}
		else{
		 	k = i;
			coeff = &linexpr0->p.coeff[k];	
			coeff_to_interval(coeff,&res->inf_coeff[k],&res->sup_coeff[k]);
		}
		
	}
	elina_coeff_t *cst = &linexpr0->cst;
	coeff_to_interval(cst,&res->inf_cst,&res->sup_cst);
	return res;
}

elina_interval_t *get_bounds_for_linexpr0(elina_manager_t *man,
                                          elina_abstract0_t *element,
                                          elina_linexpr0_t *linexpr0,
                                          size_t layerno) {

  elina_interval_t *res = elina_interval_alloc();
  fppoly_t *fp = fppoly_of_abstract0(element);
  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  // printf("start %p %lu\n",fp->layers[layerno],layerno);
  // fflush(stdout);
  expr_t *tmp = elina_linexpr0_to_expr(linexpr0);
  // printf("coming here\n");
  // fflush(stdout);
  expr_t *expr = expr_from_previous_layer(pr, tmp, fp->layers[layerno]);
  // printf("end\n");
  // fflush(stdout);
  expr_t *expr2 = copy_expr(expr);

  double lb = get_lb_using_previous_layers(man, fp, expr, layerno, true);

  double ub = get_ub_using_previous_layers(man, fp, expr2, layerno, true);

  elina_interval_set_double(res, -lb, ub);
  free_expr(expr);
  free_expr(expr2);
  free_expr(tmp);
  return res;
}

bool is_greater(elina_manager_t *man, elina_abstract0_t *element, elina_dim_t y,
                elina_dim_t x, bool use_area_heuristic) {
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

    // layer_fprint(stdout,fp->layers[3],NULL);
    double lb = get_lb_using_previous_layers(man, fp, sub, fp->numlayers,
                                             use_area_heuristic);

    // free_expr(sub);

    if (lb < 0) {
      return true;
    } else {
      return false;
    }
  }
  output_abstract_t *out = fp->out;
  expr_t *exprA = out->lexpr[y];
  expr_t *exprB = out->uexpr[x];
  if (exprA == NULL) {
    return false;
  } else {
    if (exprB == NULL) {
      if (out->output_inf[y] < 0) {
        return true;
      } else {
        return false;
      }

    } else {
      // printf("before access %zu %zu\n", exprA->size,exprB->size);
      // fflush(stdout);
      size_t sizeA = exprA->size;
      size_t sizeB = exprB->size;
      // printf("after access\n");
      // fflush(stdout);
      size_t i, k;
      expr_t *sub = (expr_t *)malloc(sizeof(expr_t));
      //
      // sub->size = size;
      sub->inf_cst = exprA->inf_cst + exprB->sup_cst;
      sub->sup_cst = exprA->sup_cst + exprB->inf_cst;
      // printf("getting here\n");
      // expr_print(exprA);
      // expr_print(exprB);
      // fflush(stdout);
      if (exprA->type == DENSE) {
        sub->inf_coeff = (double *)malloc(sizeA * sizeof(double));
        sub->sup_coeff = (double *)malloc(sizeA * sizeof(double));
        sub->dim = NULL;
        sub->size = sizeA;
        sub->type = DENSE;
        if (exprB->type == DENSE) {
          for (i = 0; i < sizeA; i++) {
            sub->inf_coeff[i] = exprA->inf_coeff[i] + exprB->sup_coeff[i];
            sub->sup_coeff[i] = exprA->sup_coeff[i] + exprB->inf_coeff[i];
          }
        } else {
          k = 0;
          for (i = 0; i < sizeA; i++) {
            if (k < sizeB && exprB->dim[k] == i) {
              sub->inf_coeff[i] = exprA->inf_coeff[i] + exprB->sup_coeff[k];
              sub->sup_coeff[i] = exprA->sup_coeff[i] + exprB->inf_coeff[k];
              k++;
            } else {
              sub->inf_coeff[i] = exprA->inf_coeff[i];
              sub->sup_coeff[i] = exprA->sup_coeff[i];
            }
          }
        }

      } else {
        if (exprB->type == DENSE) {
          sub->inf_coeff = (double *)malloc(sizeB * sizeof(double));
          sub->sup_coeff = (double *)malloc(sizeB * sizeof(double));
          sub->dim = NULL;
          sub->size = sizeB;
          sub->type = DENSE;
          i = 0;
          for (k = 0; k < sizeB; k++) {
            if (i < sizeA && exprA->dim[i] == k) {
              sub->inf_coeff[k] = exprA->inf_coeff[i] + exprB->sup_coeff[k];
              sub->sup_coeff[k] = exprA->sup_coeff[i] + exprB->inf_coeff[k];
              i++;
            } else {
              sub->inf_coeff[i] = exprB->sup_coeff[k];
              sub->sup_coeff[i] = exprB->inf_coeff[k];
            }
          }
        } else {
          sub->inf_coeff = (double *)malloc((sizeA + sizeB) * sizeof(double));
          sub->sup_coeff = (double *)malloc((sizeA + sizeB) * sizeof(double));
          sub->dim = NULL;

          sub->type = SPARSE;
          size_t l = 0;
          i = 0;
          k = 0;
          sub->dim = (size_t *)malloc((sizeA + sizeB) * sizeof(size_t));
          while (i < sizeA && k < sizeB) {
            if (exprA->dim[i] < exprB->dim[k]) {
              sub->inf_coeff[l] = exprA->inf_coeff[i];
              sub->sup_coeff[l] = exprA->sup_coeff[i];
              sub->dim[l] = exprA->dim[i];
              i++;

            } else if (exprB->dim[k] < exprA->dim[i]) {
              sub->inf_coeff[l] = exprB->sup_coeff[k];
              sub->sup_coeff[l] = exprB->inf_coeff[k];
              sub->dim[l] = exprB->dim[k];
              k++;
            } else {
              sub->inf_coeff[l] = exprA->inf_coeff[i] + exprB->sup_coeff[k];
              sub->sup_coeff[l] = exprA->sup_coeff[i] + exprB->inf_coeff[k];
              sub->dim[l] = exprA->dim[i];
              i++;
              k++;
            }
            l++;
          }
          while (i < sizeA) {
            sub->inf_coeff[l] = exprA->inf_coeff[i];
            sub->sup_coeff[l] = exprA->sup_coeff[i];
            sub->dim[l] = exprA->dim[i];
            i++;
            l++;
          }
          while (k < sizeB) {
            sub->inf_coeff[l] = exprB->inf_coeff[k];
            sub->sup_coeff[l] = exprB->sup_coeff[k];
            sub->dim[l] = exprB->dim[k];
            k++;
            l++;
          }
          sub->size = l;
          sub->inf_coeff =
              (double *)realloc(sub->inf_coeff, l * sizeof(double));
          sub->sup_coeff =
              (double *)realloc(sub->sup_coeff, l * sizeof(double));
          sub->dim = (size_t *)realloc(sub->dim, l * sizeof(size_t));
        }
      }

      // expr_print(sub);
      // fflush(stdout);
      double lb = compute_lb_from_expr(pr, sub, fp, -1);
      // printf("y: %zu x: %zu lb: %g\n",y,x,lb);
      // fflush(stdout);
      free_expr(sub);
      // double lb = -out->output_inf[y] - out->output_sup[x];
      if (lb < 0) {
        return true;
      } else {
        return false;
      }
    }
  }
}

long int max(long int a, long int b){
	return a> b? a : b;

}

void conv_handle_first_layer(elina_manager_t *man, elina_abstract0_t *abs,
                             double *filter_weights, double *filter_bias,
                             size_t *input_size, size_t *filter_size,
                             size_t num_filters, size_t *strides,
                             size_t *output_size, size_t pad_top,
                             size_t pad_left, bool has_bias,
                             size_t *predecessors) {

  size_t i, j;
  size_t num_pixels = input_size[0] * input_size[1] * input_size[2];

  size_t size = output_size[0] * output_size[1] * output_size[2];

  fppoly_internal_t *pr =
      fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  fppoly_t *res = fppoly_of_abstract0(abs);
  fppoly_alloc_first_layer(res, size, CONV, RELU);

  neuron_t **neurons = res->layers[0]->neurons;
  res->layers[0]->predecessors = predecessors;
  size_t out_x, out_y, out_z;
  size_t inp_x, inp_y, inp_z;
  size_t x_shift, y_shift;

  for (out_x = 0; out_x < output_size[0]; out_x++) {
    for (out_y = 0; out_y < output_size[1]; out_y++) {
      for (out_z = 0; out_z < output_size[2]; out_z++) {
        size_t mat_x = out_x * output_size[1] * output_size[2] +
                       out_y * output_size[2] + out_z;
        size_t num_coeff = input_size[2] * filter_size[0] * filter_size[1];
        size_t actual_coeff = 0;
        double *coeff = (double *)malloc(num_coeff * sizeof(double));
        size_t *dim = (size_t *)malloc(num_coeff * sizeof(double));
        i = 0;
        for (inp_z = 0; inp_z < input_size[2]; inp_z++) {
          for (x_shift = 0; x_shift < filter_size[0]; x_shift++) {
            for (y_shift = 0; y_shift < filter_size[1]; y_shift++) {
              long int x_val = out_x * strides[0] + x_shift - pad_top;
              long int y_val = out_y * strides[1] + y_shift - pad_left;
              if (y_val < 0 || y_val >= (long int)input_size[1]) {
                continue;
              }

              if (x_val < 0 || x_val >= (long int)input_size[0]) {
                continue;
              }
              size_t mat_y = x_val * input_size[1] * input_size[2] +
                             y_val * input_size[2] + inp_z;
              if (mat_y >= num_pixels) {
                continue;
              }
              size_t filter_index =
                  x_shift * filter_size[1] * input_size[2] * output_size[2] +
                  y_shift * input_size[2] * output_size[2] +
                  inp_z * output_size[2] + out_z;
              coeff[i] = filter_weights[filter_index];
              dim[i] = mat_y;
              i++;
              actual_coeff++;
            }
          }
        }
        double cst = has_bias ? filter_bias[out_z] : 0;
        neurons[mat_x]->expr =
            create_sparse_expr(coeff, cst, dim, actual_coeff);
        sort_sparse_expr(neurons[mat_x]->expr);

        neurons[mat_x]->lb =
            compute_lb_from_expr(pr, neurons[mat_x]->expr, res, -1);
        neurons[mat_x]->ub =
            compute_ub_from_expr(pr, neurons[mat_x]->expr, res, -1);
        free(coeff);
        free(dim);
      }
    }
  }

  // printf("return here\n");
  // fppoly_fprint(stdout,man,res,NULL);
  // fflush(stdout);
  return;
}

void conv_handle_intermediate_layer(
    elina_manager_t *man, elina_abstract0_t *element, double *filter_weights,
    double *filter_bias, size_t *input_size, size_t *filter_size,
    size_t num_filters, size_t *strides, size_t *output_size, size_t pad_top,
    size_t pad_left, bool has_bias, size_t *predecessors,
    activation_type_t activation, bool use_area_heuristic) {
  // printf("conv intermediate starts here\n");
  // fflush(stdout);
  fppoly_t *fp = fppoly_of_abstract0(element);
  size_t numlayers = fp->numlayers;
  size_t i, j;
  size_t num_pixels = input_size[0] * input_size[1] * input_size[2];

  output_size[2] = num_filters;
  size_t num_out_neurons = output_size[0] * output_size[1] * output_size[2];
  // printf("num_out_neurons: %zu %zu\n",num_out_neurons,num_pixels);
  // fflush(stdout);
  fppoly_add_new_layer(fp, num_out_neurons, CONV, activation);
  neuron_t **out_neurons = fp->layers[numlayers]->neurons;
  fp->layers[numlayers]->predecessors = predecessors;
  size_t out_x, out_y, out_z;
  size_t inp_x, inp_y, inp_z;
  size_t x_shift, y_shift;

  for (out_x = 0; out_x < output_size[0]; out_x++) {
    for (out_y = 0; out_y < output_size[1]; out_y++) {
      for (out_z = 0; out_z < output_size[2]; out_z++) {
        size_t mat_x = out_x * output_size[1] * output_size[2] +
                       out_y * output_size[2] + out_z;
        size_t num_coeff = input_size[2] * filter_size[0] * filter_size[1];
        size_t actual_coeff = 0;
        double *coeff = (double *)malloc(num_coeff * sizeof(double));
        size_t *dim = (size_t *)malloc(num_coeff * sizeof(double));
        i = 0;
        for (inp_z = 0; inp_z < input_size[2]; inp_z++) {
          for (x_shift = 0; x_shift < filter_size[0]; x_shift++) {
            for (y_shift = 0; y_shift < filter_size[1]; y_shift++) {
              long int x_val = out_x * strides[0] + x_shift - pad_top;
              long int y_val = out_y * strides[1] + y_shift - pad_left;
              if (y_val < 0 || y_val >= (long int)input_size[1]) {
                continue;
              }

              if (x_val < 0 || x_val >= (long int)input_size[0]) {
                continue;
              }
              size_t mat_y = x_val * input_size[1] * input_size[2] +
                             y_val * input_size[2] + inp_z;
              if (mat_y >= num_pixels) {
                continue;
              }
              size_t filter_index =
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
        double cst = has_bias ? filter_bias[out_z] : 0;
        out_neurons[mat_x]->expr =
            create_sparse_expr(coeff, cst, dim, actual_coeff);
        sort_sparse_expr(out_neurons[mat_x]->expr);
        free(coeff);
        free(dim);
      }
    }
  }

  update_state_using_previous_layers_parallel(man, fp, numlayers,
                                              use_area_heuristic);

  // printf("return here2\n");
  // fppoly_fprint(stdout,man,fp,NULL);
  // fflush(stdout);
  return;
}

void conv_handle_intermediate_relu_layer(
    elina_manager_t *man, elina_abstract0_t *element, double *filter_weights,
    double *filter_bias, size_t *input_size, size_t *filter_size,
    size_t num_filters, size_t *strides, size_t *output_size, size_t pad_top,
    size_t pad_left, bool has_bias, size_t *predecessors,
    bool use_area_heuristic) {

  conv_handle_intermediate_layer(man, element, filter_weights, filter_bias,
                                 input_size, filter_size, num_filters, strides,
                                 output_size, pad_top, pad_left, has_bias,
                                 predecessors, RELU, use_area_heuristic);
}

void conv_handle_intermediate_affine_layer(
    elina_manager_t *man, elina_abstract0_t *element, double *filter_weights,
    double *filter_bias, size_t *input_size, size_t *filter_size,
    size_t num_filters, size_t *strides, size_t *output_size, size_t pad_top,
    size_t pad_left, bool has_bias, size_t *predecessors,
    bool use_area_heuristic) {

  conv_handle_intermediate_layer(man, element, filter_weights, filter_bias,
                                 input_size, filter_size, num_filters, strides,
                                 output_size, pad_top, pad_left, has_bias,
                                 predecessors, NONE, use_area_heuristic);
}

void handle_residual_layer(elina_manager_t *man, elina_abstract0_t *element,
                           size_t num_neurons, size_t *predecessors,
                           activation_type_t activation,
                           bool use_area_heuristic) {
  fppoly_t *fp = fppoly_of_abstract0(element);
  size_t numlayers = fp->numlayers;
  fppoly_add_new_layer(fp, num_neurons, RESIDUAL, activation);
  fp->layers[numlayers]->predecessors = predecessors;
  size_t i;
  neuron_t **neurons = fp->layers[numlayers]->neurons;
  for (i = 0; i < num_neurons; i++) {
    double *coeff = (double *)malloc(sizeof(double));
    coeff[0] = 1;
    size_t *dim = (size_t *)malloc(sizeof(size_t));
    dim[0] = i;
    neurons[i]->expr = create_sparse_expr(coeff, 0, dim, 1);
  }
  update_state_using_previous_layers_parallel(man, fp, numlayers,
                                              use_area_heuristic);
}

void handle_residual_relu_layer(elina_manager_t *man,
                                elina_abstract0_t *element, size_t num_neurons,
                                size_t *predecessors, bool use_area_heuristic) {
  handle_residual_layer(man, element, num_neurons, predecessors, RELU,
                        use_area_heuristic);
}

void handle_residual_affine_layer(elina_manager_t *man,
                                  elina_abstract0_t *element,
                                  size_t num_neurons, size_t *predecessors,
                                  bool use_area_heuristic) {
  handle_residual_layer(man, element, num_neurons, predecessors, NONE,
                        use_area_heuristic);
}

void free_neuron(neuron_t *neuron){
  if (neuron->expr) {
    free_expr(neuron->expr);
  }
        if(neuron->lexpr){
		free_expr(neuron->lexpr);
	}
        if (neuron->uexpr) {
          free_expr(neuron->uexpr);
        }
        free(neuron);
}

void free_non_lstm_layer_expr(elina_manager_t *man, elina_abstract0_t *abs, size_t layerno){
    fppoly_t *fp = fppoly_of_abstract0(abs);
    if(layerno >= fp->numlayers){
        fprintf(stdout,"the layer does not exist\n");
        return;
    }
    layer_t * layer = fp->layers[layerno];
    size_t dims = layer->dims;
    size_t i;
    for(i=0; i < dims; i++){
        neuron_t *neuron = layer->neurons[i];
        if (neuron->expr) {
          free_expr(neuron->expr);
        }
        if(neuron->lexpr){
            free_expr(neuron->lexpr);
        }
        if (neuron->uexpr) {
          free_expr(neuron->uexpr);
        }
    }
}

void layer_free(layer_t * layer){
	size_t dims = layer->dims;
	size_t i;
	for(i=0; i < dims; i++){
		free_neuron(layer->neurons[i]);
	}
	free(layer->neurons);
	layer->neurons = NULL;
	if(layer->h_t_inf!=NULL){
		free(layer->h_t_inf);
		layer->h_t_inf = NULL;
	}

	if(layer->h_t_sup!=NULL){
		free(layer->h_t_sup);
		layer->h_t_sup = NULL;
	}
	
	if(layer->c_t_inf!=NULL){
		free(layer->c_t_inf);
		layer->c_t_inf = NULL;
	}

	if(layer->c_t_sup!=NULL){
		free(layer->c_t_sup);
		layer->c_t_sup = NULL;
	}

        free(layer);
	layer = NULL;
}



void fppoly_free(elina_manager_t *man, fppoly_t *fp){
	size_t i;
	size_t output_size = fp->layers[fp->numlayers-1]->dims;
	for(i=0; i < fp->numlayers; i++){
		layer_free(fp->layers[i]);
	}

        for (i = 0; i < output_size; i++) {
          if (fp->out->lexpr[i]) {
            free_expr(fp->out->lexpr[i]);
          }
          if (fp->out->uexpr[i]) {
            free_expr(fp->out->uexpr[i]);
          }
        }

        free(fp->layers);
	fp->layers = NULL;
	free(fp->input_inf);
	fp->input_inf = NULL;
        if(fp->input_lexpr!=NULL && fp->input_uexpr!=NULL){
		for(i=0; i < fp->num_pixels; i++){
			free(fp->input_lexpr[i]);
			free(fp->input_uexpr[i]);
		}
	
		free(fp->input_lexpr);
		fp->input_lexpr = NULL;
		free(fp->input_uexpr);
		fp->input_uexpr = NULL;
        }
	free(fp->input_sup);
	fp->input_sup = NULL;
        free(fp->out->output_inf);
        fp->out->output_inf = NULL;
        free(fp->out->output_sup);
        fp->out->output_sup = NULL;
        free(fp->out->lexpr);
        fp->out->lexpr = NULL;
        free(fp->out->uexpr);
        fp->out->uexpr = NULL;
        free(fp->out);
        fp->out = NULL;
        free(fp);
	fp = NULL;
}






void fppoly_fprint(FILE* stream, elina_manager_t* man, fppoly_t* fp, char** name_of_dim){
	size_t i;
	for(i = 0; i < fp->numlayers; i++){
		fprintf(stream,"layer: %zu\n", i);
		layer_fprint(stream, fp->layers[i], name_of_dim);
	}
	size_t output_size = fp->layers[fp->numlayers-1]->dims;
        if (fp->out != NULL) {
          fprintf(stream, "OUTPUT bounds: \n");
          for (i = 0; i < output_size; i++) {
            fprintf(stream, "%zu: [%g,%g] \n", i, -fp->out->output_inf[i],
                    fp->out->output_sup[i]);
          }
        }
}


elina_interval_t * box_for_neuron(elina_manager_t* man, elina_abstract0_t * abs, size_t layerno, size_t neuron_no){
	fppoly_t *fp = fppoly_of_abstract0(abs);
	if(layerno >= fp->numlayers){
		fprintf(stdout,"the layer does not exist\n");
		return NULL;
	}
	layer_t * layer = fp->layers[layerno];
	size_t dims = layer->dims;
	if(neuron_no >= dims){
		fprintf(stdout,"the neuron does not exist\n");
		return NULL;
	}
	neuron_t * neuron = layer->neurons[neuron_no];
	elina_interval_t * res = elina_interval_alloc();
	elina_interval_set_double(res,-neuron->lb,neuron->ub);
	return res;
}

elina_interval_t ** box_for_layer(elina_manager_t* man, elina_abstract0_t * abs, size_t layerno){
	fppoly_t *fp = fppoly_of_abstract0(abs);
	if(layerno >= fp->numlayers){
		fprintf(stdout,"the layer does not exist\n");
		return NULL;
	}
	layer_t * layer = fp->layers[layerno];
	size_t dims = layer->dims;
	elina_interval_t ** itv_arr = (elina_interval_t **)malloc(dims*sizeof(elina_interval_t *));
	size_t i;
	for(i=0; i< dims; i++){
		itv_arr[i] = box_for_neuron(man, abs, layerno, i);
		
	}
	return itv_arr;
}


size_t get_num_neurons_in_layer(elina_manager_t* man, elina_abstract0_t * abs, size_t layerno){
	fppoly_t *fp = fppoly_of_abstract0(abs);
	if(layerno >= fp->numlayers){
		fprintf(stdout,"the layer does not exist\n");
		return 0;
	}
	layer_t * layer = fp->layers[layerno];
	size_t dims = layer->dims;
	
	return dims;
}

elina_linexpr0_t * get_expr_for_output_neuron(elina_manager_t *man, elina_abstract0_t *abs, size_t i, bool is_lower){
	fppoly_internal_t *pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
	fppoly_t *fp = fppoly_of_abstract0(abs);
	
	size_t output_size = fp->layers[fp->numlayers-1]->dims;
	if(i >= output_size){
		return NULL;
	}
	size_t num_pixels = fp->num_pixels;
	expr_t * expr = NULL;
	if(is_lower){
          expr = fp->out->lexpr[i];
        }
	else{
          expr = fp->out->uexpr[i];
        }
	elina_linexpr0_t * res = NULL;
	size_t j,k;
	if((fp->input_lexpr!=NULL) && (fp->input_uexpr!=NULL)){
		if(is_lower){
			expr =  replace_input_poly_cons_in_lexpr(pr, expr, fp);
		}
		else{
			expr =  replace_input_poly_cons_in_uexpr(pr, expr, fp);
		}
	}
	size_t expr_size = expr->size;
	if(expr->type==SPARSE){
		sort_sparse_expr(expr);
		res = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,expr_size);
	}
	else{
		res = elina_linexpr0_alloc(ELINA_LINEXPR_DENSE,expr_size);
	}
	elina_linexpr0_set_cst_interval_double(res,-expr->inf_cst,expr->sup_cst);
	
	for(j=0;j < expr_size; j++){
		if(expr->type==DENSE){
			k = j;
		}
		else{
			k = expr->dim[j];
		}
		elina_linexpr0_set_coeff_interval_double(res,k,-expr->inf_coeff[j],expr->sup_coeff[j]);
	}
	if((fp->input_lexpr!=NULL) && (fp->input_uexpr!=NULL)){
		free_expr(expr);
	}
	return res;
}

elina_linexpr0_t * get_lexpr_for_output_neuron(elina_manager_t *man,elina_abstract0_t *abs, size_t i){
	return get_expr_for_output_neuron(man,abs,i, true);
}

elina_linexpr0_t * get_uexpr_for_output_neuron(elina_manager_t *man,elina_abstract0_t *abs, size_t i){
	return get_expr_for_output_neuron(man,abs,i, false);
}

void update_bounds_for_neuron(elina_manager_t *man, elina_abstract0_t *abs, size_t layerno, size_t neuron_no, double lb, double ub){
	fppoly_t *fp = fppoly_of_abstract0(abs);
	if(layerno >= fp->numlayers){
		fprintf(stdout,"the layer does not exist\n");
		return;
	}
	layer_t * layer = fp->layers[layerno];
	neuron_t * neuron = layer->neurons[neuron_no];
	neuron->lb = -lb;
	neuron->ub = ub;
}

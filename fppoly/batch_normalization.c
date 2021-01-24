#include "batch_normalization.h"

expr_t * create_batch_normalization_expr(double weight, double bias, size_t i){
	expr_t * res = alloc_expr();  
	res->inf_coeff = (double *)malloc(sizeof(double));
	res->sup_coeff = (double *)malloc(sizeof(double));
	res->type = SPARSE;
	res->size = 1;  
	res->dim = (size_t*)malloc(sizeof(size_t));
	res->dim[0] = i;
	res->inf_coeff[0] = -weight;
	res->sup_coeff[0] = weight;
	res->sup_cst = -bias;
	res->inf_cst = bias;
	
	return res;
}


void handle_batch_normalization_layer(elina_manager_t *man, elina_abstract0_t* element, double * weights, double *bias, size_t num_neurons, size_t *predecessors, size_t num_predecessors){
	
	assert(num_predecessors==1);
	fppoly_t *fp = fppoly_of_abstract0(element);
	size_t numlayers = fp->numlayers;
	fppoly_add_new_layer(fp, num_neurons, predecessors, num_predecessors,true);
	neuron_t **out_neurons = fp->layers[numlayers]->neurons;
	int k = predecessors[0]-1;
	neuron_t **in_neurons = fp->layers[k]->neurons;
	size_t i;
	
	for(i=0; i < num_neurons; i++){
		out_neurons[i]->lb = -(-in_neurons[i]->lb*weights[i] + bias[i]);
		out_neurons[i]->ub = in_neurons[i]->ub*weights[i] + bias[i];
		out_neurons[i]->lexpr = create_batch_normalization_expr(weights[i], bias[i], i);
		out_neurons[i]->uexpr = out_neurons[i]->lexpr;
		
	}
}


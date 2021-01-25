#include "round_approx.h"

expr_t * create_round_expr(neuron_t *out_neuron, neuron_t *in_neuron, size_t i, bool use_default_heuristics, bool is_lower){
	expr_t * res = alloc_expr();  
	res->inf_coeff = (double *)malloc(sizeof(double));
	res->sup_coeff = (double *)malloc(sizeof(double));
	res->type = SPARSE;
	res->size = 1;  
	res->dim = (size_t*)malloc(sizeof(size_t));
	res->dim[0] = i;
	double lb = in_neuron->lb;
	double ub = in_neuron->ub;
	//res->inf_coeff[0] = 1.0;
	//res->sup_coeff[0] = 1.0;
	res->inf_coeff[0] = 0.0;
	res->sup_coeff[0] = 0.0;
	res->inf_cst = 0.0;
	res->sup_cst = 0.0;
	if(is_lower){
		res->inf_cst = -round(-lb);
		res->sup_cst = -res->inf_cst;
	}
	else{
		res->sup_cst = round(ub);
		res->inf_cst = -res->sup_cst;
	}
	
	return res;
}


void handle_round_layer(elina_manager_t *man, elina_abstract0_t* element, size_t num_neurons, size_t *predecessors, size_t num_predecessors, bool use_default_heuristics){
	
	assert(num_predecessors==1);
	fppoly_t *fp = fppoly_of_abstract0(element);
	size_t numlayers = fp->numlayers;
	fppoly_add_new_layer(fp, num_neurons, predecessors, num_predecessors,true);
	neuron_t **out_neurons = fp->layers[numlayers]->neurons;
	int k = predecessors[0]-1;
	neuron_t **in_neurons = fp->layers[k]->neurons;
	size_t i;
	
	for(i=0; i < num_neurons; i++){
		out_neurons[i]->lb = -round(-in_neurons[i]->lb);
		out_neurons[i]->ub = round(in_neurons[i]->ub);
		out_neurons[i]->lexpr = create_round_expr(out_neurons[i], in_neurons[i], i, use_default_heuristics, true);
		out_neurons[i]->uexpr = create_round_expr(out_neurons[i], in_neurons[i], i, use_default_heuristics, false);
		
	}
}


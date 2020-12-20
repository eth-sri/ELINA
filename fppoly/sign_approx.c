#include "sign_approx.h"

expr_t * create_sign_expr(neuron_t *out_neuron, neuron_t *in_neuron, size_t i, bool is_lower){
	expr_t * res = alloc_expr();  
	res->inf_coeff = (double *)malloc(sizeof(double));
	res->sup_coeff = (double *)malloc(sizeof(double));
	res->inf_cst = 0.0;
	res->sup_cst = 0.0;
	res->type = SPARSE;
	res->size = 1;  
	res->dim = (size_t*)malloc(sizeof(size_t));
	res->dim[0] = i;
	double lb = round(in_neuron->lb);
	double ub = round(in_neuron->ub);
	
	if(ub<0 || lb <= 0){
		res->inf_coeff[0] = 0.0;
		res->sup_coeff[0] = 0.0;
		if(lb<=0){
			res->inf_cst = -1.0;
			res->sup_cst = 1.0;
		}
	}
		
	else if(is_lower){
		res->inf_coeff[0] = -1/ub;
		res->sup_coeff[0] = 1/ub;
		//res->inf_coeff[0] = 0.0;
		//res->sup_coeff[0] = 0.0;
	}
	else{
		res->inf_coeff[0] = 1/lb;
		res->sup_coeff[0] = -1/lb;
		//res->inf_coeff[0] = 0.0;
		//res->sup_coeff[0] = 0.0;
		//res->inf_cst = -1.0;
		//res->sup_cst = 1.0;
	}
	return res;
}


void handle_sign_layer(elina_manager_t *man, elina_abstract0_t* element, size_t num_neurons, size_t *predecessors, size_t num_predecessors){
	
	assert(num_predecessors==1);
	fppoly_t *fp = fppoly_of_abstract0(element);
	size_t numlayers = fp->numlayers;
	fppoly_add_new_layer(fp, num_neurons, predecessors, num_predecessors,true);
	neuron_t **out_neurons = fp->layers[numlayers]->neurons;
	int k = predecessors[0]-1;
	neuron_t **in_neurons = fp->layers[k]->neurons;
	size_t i;
	
	for(i=0; i < num_neurons; i++){
		double lb = round(in_neurons[i]->lb);
		double ub = round(in_neurons[i]->ub);
		if(lb <= 0){
			out_neurons[i]->lb = -1.0;
			out_neurons[i]->ub = 1.0;
		}
		else if (ub<0){
			out_neurons[i]->lb = 0.0;
			out_neurons[i]->ub = 0.0;
		}
		else{
			out_neurons[i]->lb = 0.0;
			out_neurons[i]->ub = 1.0;
		}
		
		out_neurons[i]->lexpr = create_sign_expr(out_neurons[i], in_neurons[i], i, true);
		out_neurons[i]->uexpr = create_sign_expr(out_neurons[i], in_neurons[i], i, false);
		
	}
	
}


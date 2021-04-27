#include "leakyrelu_approx.h"

expr_t * create_leakyrelu_expr(neuron_t *out_neuron, neuron_t *in_neuron, size_t i, double alpha, bool use_default_heuristics, bool is_lower){
	expr_t * res = alloc_expr();  
	res->inf_coeff = (double *)malloc(sizeof(double));
	res->sup_coeff = (double *)malloc(sizeof(double));
	res->inf_cst = 0.0;
	res->sup_cst = 0.0;
	res->type = SPARSE;
	res->size = 1;  
	res->dim = (size_t*)malloc(sizeof(size_t));
	res->dim[0] = i;
	double lb = in_neuron->lb;
	double ub = in_neuron->ub;
	double width = ub + lb;
	double num = ub + alpha*lb;
	double lambda_inf = -num/width;
	double lambda_sup = num/width;
	
	if(ub<=0){
		res->inf_coeff[0] = -alpha;
		res->sup_coeff[0] = alpha;
	}
	else if(lb<0){
		res->inf_coeff[0] = -1.0;
		res->sup_coeff[0] = 1.0;
	}
		
	else if(is_lower){
		double area1 = 0.5*ub*width;
		double area2 = 0.5*lb*width;
		if(use_default_heuristics){
			if(area1 < area2){
				res->inf_coeff[0] = -alpha;
				res->sup_coeff[0] = alpha;
			}
			else{
				res->inf_coeff[0] = -1.0;
				res->sup_coeff[0] = 1.0;
			}
		}
		else{
				res->inf_coeff[0] = -alpha;
				res->sup_coeff[0] = alpha;
		}
	}
	else{
		double mu_num = 0.99*ub*lb;
		double mu_inf = -mu_num/width;
		double mu_sup = mu_num/width;
		res->inf_coeff[0] = lambda_inf;
		res->sup_coeff[0] = lambda_sup;
		res->inf_cst = mu_inf;
		res->sup_cst = mu_sup;
	}
	return res;
}


void handle_leakyrelu_layer(elina_manager_t *man, elina_abstract0_t* element, size_t num_neurons, size_t *predecessors, size_t num_predecessors, double alpha, bool use_default_heuristics){
	
	assert(num_predecessors==1);
	fppoly_t *fp = fppoly_of_abstract0(element);
	size_t numlayers = fp->numlayers;
	fppoly_add_new_layer(fp, num_neurons, predecessors, num_predecessors,true);
	neuron_t **out_neurons = fp->layers[numlayers]->neurons;
	int k = predecessors[0]-1;
	neuron_t **in_neurons = fp->layers[k]->neurons;
	size_t i;
	
	for(i=0; i < num_neurons; i++){
		out_neurons[i]->lb = -fmax(alpha*-in_neurons[i]->lb, -in_neurons[i]->lb);
		out_neurons[i]->ub = fmax(alpha*in_neurons[i]->ub,in_neurons[i]->ub);
		out_neurons[i]->lexpr = create_leakyrelu_expr(out_neurons[i], in_neurons[i], i, alpha, use_default_heuristics, true);
		out_neurons[i]->uexpr = create_leakyrelu_expr(out_neurons[i], in_neurons[i], i, alpha, use_default_heuristics, false);
		
	}
	
}


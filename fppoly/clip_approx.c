#include "clip_approx.h"

expr_t * create_clip_expr(neuron_t *out_neuron, double lb, double ub, size_t i, double min_input, double max_input, bool use_default_heuristics, bool is_lower){
	expr_t * res = alloc_expr();  
	res->inf_coeff = (double *)malloc(sizeof(double));
	res->sup_coeff = (double *)malloc(sizeof(double));
	res->inf_cst = 0.0;
	res->sup_cst = 0.0;
	res->type = SPARSE;
	res->size = 1;  
	res->dim = (size_t*)malloc(sizeof(size_t));
	res->dim[0] = i;

	
	if(ub<=min_input){
		res->inf_coeff[0] = 0.0;
		res->sup_coeff[0] = 0.0;
		res->inf_cst = -min_input;
		res->sup_cst = min_input;
	}
	else if(-lb>=max_input){
		res->inf_coeff[0] = 0.0;
		res->sup_coeff[0] = 0.0;
		res->inf_cst = -max_input;
		res->sup_cst = max_input;
	}
	else if(-lb < min_input){
		if(ub<=max_input){
			
			if(is_lower){
				if(use_default_heuristics){
					if(ub < -lb){
						res->inf_coeff[0] = 0.0;
						res->sup_coeff[0] = 0.0;
					}
					else{
						res->inf_coeff[0] = -1.0;
						res->sup_coeff[0] = 1.0;
					}
				}	
				else{
					res->inf_coeff[0] = 0.0;
					res->sup_coeff[0] = 0.0;
				}
			}
			else{
				double width = ub + lb;
				double lambda_inf = -ub/width;
				double lambda_sup = ub/width;
				double mu_inf = lambda_inf*lb;
				double mu_sup = lambda_sup*lb;
				res->inf_coeff[0] = lambda_inf;
				res->sup_coeff[0] = lambda_sup;
				res->inf_cst = mu_inf;
				res->sup_cst = mu_sup;
			}
		}
		else{
			res->inf_coeff[0] = 0.0;
			res->sup_coeff[0] = 0.0;
			if(is_lower){
				res->inf_cst = -min_input;
				res->sup_cst = min_input;
			}
			else{
				res->inf_cst = -max_input;
				res->sup_cst = max_input;
			}			
		}
	}
	else {
		if(ub<=max_input){
			res->inf_coeff[0] = -1.0;
			res->sup_coeff[0] = 1.0;
		}
		else{
			if(is_lower){
				double width = ub + lb;
				double num = max_input + lb;
				double lambda_inf = -num/width;
				double lambda_sup = num/width;
				num = -lb*(ub - max_input);
				double mu_inf = -num/width;
				double mu_sup = num/width;
				res->inf_coeff[0] = lambda_inf;
				res->sup_coeff[0] = lambda_sup;
				res->inf_cst = mu_inf;
				res->sup_cst = mu_sup;
			}
			else{
				res->inf_coeff[0] = 0.0;
				res->sup_coeff[0] = 0.0;
				res->inf_cst = -max_input;
				res->sup_cst = max_input;
			}
		}
	}	
	
	
	return res;
}


void handle_clip_layer(elina_manager_t *man, elina_abstract0_t* element, double min_input, double max_input, size_t num_neurons, size_t *predecessors, size_t num_predecessors, bool use_default_heuristics){
	assert(num_predecessors==1);

	fppoly_t *fp = fppoly_of_abstract0(element);
	
	size_t numlayers = fp->numlayers;
	fppoly_add_new_layer(fp, num_neurons, predecessors, num_predecessors,true);
	neuron_t **out_neurons = fp->layers[numlayers]->neurons;
	int k = predecessors[0]-1;
	neuron_t **in_neurons = NULL;
	if(k!=-1){
		in_neurons = fp->layers[k]->neurons;
	}
	size_t i;
	for(i=0; i < num_neurons; i++){
		double lb, ub;
		
		if(k==-1){
			lb = -fp->input_inf[i];
			ub = fp->input_sup[i];
		}
		else{
			
			
			lb = -in_neurons[i]->lb;
			ub = in_neurons[i]->ub;
		}
		if(ub<=min_input){
			out_neurons[i]->lb = -min_input;
			out_neurons[i]->ub = min_input;
		}
		else if(lb>=max_input){
			out_neurons[i]->lb = -max_input;
			out_neurons[i]->ub = max_input;
		}
		else if(lb < min_input){
			out_neurons[i]->lb = -min_input;
			out_neurons[i]->ub = fmin(ub, max_input);
			
		}
		else{
			out_neurons[i]->lb = -lb;
			out_neurons[i]->ub = fmin(ub, max_input);
			
		}
		
		out_neurons[i]->lexpr = create_clip_expr(out_neurons[i], -lb, ub,  i, min_input, max_input, use_default_heuristics, true);
		out_neurons[i]->uexpr = create_clip_expr(out_neurons[i], -lb, ub, i, min_input, max_input, use_default_heuristics, false);
		
	}
	
}


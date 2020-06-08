#include "log_approx.h"

expr_t * create_log_expr(fppoly_internal_t *pr, neuron_t * out_neuron, neuron_t * in_neuron, size_t i, bool is_lower){
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
	double u_plus_l_inf = -(lb+ub);
	double u_plus_l_sup = -u_plus_l_inf;
		
	//= (u_plus_l)*(u_plus_l)
	if(!is_lower){
		double one_inf = 1;
		double one_sup = -1;
		double u_plus_l_sup = ub-lb;
		double u_plus_l_inf = -(ub-lb);
		out_neuron->ub = log(ub);
		double lambda_sup = -2/u_plus_l_inf;
		double lambda_inf = -2/u_plus_l_sup;
		res->inf_coeff[0] = lambda_inf;
		res->sup_coeff[0] = lambda_sup;
				
		fesetround(FE_DOWNWARD);
		double mu_inf = -log(-u_plus_l_inf/2);
		fesetround(FE_UPWARD);
		double mu_sup = log(u_plus_l_sup/2);
		elina_double_interval_add_cst_coeff(pr,&mu_inf,&mu_sup,one_inf, one_sup, mu_inf, mu_sup);
		res->inf_cst = mu_inf;
		res->sup_cst = mu_sup;
	}
	else{
		double u_minus_l_sup = ub +lb;
		double u_minus_l_inf = -(ub+lb);
		fesetround(FE_DOWNWARD);
		out_neuron->lb = -log(-lb);
		fesetround(FE_UPWARD);
		if (u_minus_l_sup < 1e-9) {
		    double log_lb = out_neuron->lb, log_ub = -out_neuron->lb;
		    if (log_lb > 0) {
			log_lb = -log_lb;
			log_ub = -log_ub;
		    }
		    res->inf_cst = out_neuron->lb;
		    res->sup_cst = -out_neuron->lb;
	       } 
	       else {
		    double inv_u_by_l_sup = -1/u_minus_l_inf;
		    double inv_u_by_l_inf = -1/u_minus_l_sup;

		    double u_by_l_sup = -ub/lb;
		    double u_by_l_inf = ub/lb;

		    fesetround(FE_DOWNWARD);
		    double log_u_by_l_inf = -log(-u_by_l_inf);
		    double log_l_inf = -log(-lb);

		    fesetround(FE_UPWARD);
		    double log_u_by_l_sup = log(u_by_l_sup);
		    double log_l_sup = log(-lb);

		    double lambda_inf, lambda_sup;
		    elina_double_interval_mul_cst_coeff(pr,&lambda_inf,&lambda_sup,log_u_by_l_inf,log_u_by_l_sup,inv_u_by_l_inf,inv_u_by_l_sup);
		    res->inf_coeff[0] = lambda_inf;
		    res->sup_coeff[0] = lambda_sup;

		    double mu_inf, mu_sup;
		    elina_double_interval_mul_cst_coeff(pr,&mu_inf,&mu_sup,-lb,lb,lambda_inf,lambda_sup);
		    elina_double_interval_add_cst_coeff(pr,&mu_inf,&mu_sup,log_l_inf, log_l_sup, mu_inf, mu_sup);

		    
		    res->inf_cst = mu_inf;
		    res->sup_cst = mu_sup;
		}
	}
	return res;
}


void handle_log_layer(elina_manager_t *man, elina_abstract0_t* element, size_t num_neurons, size_t *predecessors, size_t num_predecessors){
	assert(num_predecessors==1);
	fppoly_t *fp = fppoly_of_abstract0(element);
	fppoly_internal_t *pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
	size_t numlayers = fp->numlayers;
	fppoly_add_new_layer(fp, num_neurons, predecessors, num_predecessors, true);
	neuron_t **out_neurons = fp->layers[numlayers]->neurons;
	int k = predecessors[0]-1;
	neuron_t **in_neurons = fp->layers[k]->neurons;
	size_t i;
	for(i=0; i < num_neurons; i++){
		out_neurons[i]->lexpr = create_log_expr(pr, out_neurons[i], in_neurons[i], i, true);
		out_neurons[i]->uexpr = create_log_expr(pr, out_neurons[i], in_neurons[i], i, false);
	}
}


#include "parabola_approx.h"

expr_t * create_parabola_expr(fppoly_internal_t *pr, neuron_t * out_neuron, neuron_t *in_neuron, size_t i, bool is_lower){
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
	double u_plus_l_sup = (ub-lb);
        double u_plus_l_inf = -(ub-lb);
        //= (u_plus_l)*(u_plus_l)
	if(!is_lower){
            double lu_sup;
            double lu_inf;
            out_neuron->ub = ub*ub;
            elina_double_interval_mul_cst_coeff(pr,&lu_inf,&lu_sup,lb,-lb,-ub,ub);
            res->inf_coeff[0] = u_plus_l_inf;
            res->sup_coeff[0] = u_plus_l_sup;
	    res->inf_cst = lu_inf;
	    res->sup_cst = lu_sup;
	}
	else{
	    res->inf_coeff[0] = 4*lb;
	    res->sup_coeff[0] = -4*lb;
	    res->inf_cst = 4*lb*lb;
	    res->sup_cst = -4*lb*lb;
	    if(lb>0){
	    	out_neuron->lb = 0;
	    }
	    else{
	    	out_neuron->lb = -lb*lb;
	    }
	}

	return res;
}

void handle_parabola_layer(elina_manager_t *man, elina_abstract0_t* element, size_t num_neurons, size_t *predecessors, size_t num_predecessors){
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
		out_neurons[i]->lexpr = create_parabola_expr(pr, out_neurons[i], in_neurons[i], i, true);
		out_neurons[i]->uexpr = create_parabola_expr(pr, out_neurons[i], in_neurons[i], i, false);
	}
}


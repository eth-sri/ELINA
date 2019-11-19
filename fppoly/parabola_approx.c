#include "parabola_approx.h"

expr_t * lexpr_replace_parabola_bounds(fppoly_internal_t * pr, expr_t * expr, neuron_t ** neurons){
	size_t num_neurons = expr->size;
	size_t i,k;
	expr_t * res = alloc_expr();  
	res->inf_coeff = (double *)malloc(num_neurons*sizeof(double));
	res->sup_coeff = (double *)malloc(num_neurons*sizeof(double));
	res->inf_cst = expr->inf_cst;
	res->sup_cst = expr->sup_cst;
	res->type = expr->type;
	res->size = num_neurons;  

	for(i = 0; i < num_neurons; i++){
		if(expr->type==DENSE){
			k = i;
		}
		else{
			k = expr->dim[i];
		}
		neuron_t *neuron_k = neurons[k];
		double lb = neurons[k]->lb;
		double ub = neurons[k]->ub;
		res->inf_coeff[i] = 0.0;
		res->sup_coeff[i] = 0.0;
		if((expr->sup_coeff[i]==0) && (expr->inf_coeff[i]==0)){
			continue;
		}
		double u_plus_l_sup = (ub-lb);
		double u_plus_l_inf = -(ub-lb);
		
		
                //= (u_plus_l)*(u_plus_l)
		if(expr->sup_coeff[i]<0){
			
            double lu_sup;
            double lu_inf;
            elina_double_interval_mul_cst_coeff(pr,&lu_inf,&lu_sup,lb,-lb,-ub,ub);
			//res->coeff[i] = lambda*expr->coeff[i];
			//res->cst = res->cst + expr->coeff[i]*mu;
			elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],u_plus_l_inf,u_plus_l_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
			double tmp1, tmp2;
			elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,lu_inf,lu_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
			res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
			res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
		}
		else if (expr->inf_coeff[i]<0){
			res->inf_coeff[i] = 0;
			res->sup_coeff[i] = 0;
			//double u_plus_l_sq_inf, u_plus_l_sq_sup; 
			//u_plus_l_sq_inf = u_plus_l_inf/2;
			//u_plus_l_sq_sup = u_plus_l_sup/2;
                        //elina_double_interval_mul_cst_coeff(pr,&u_plus_l_sq_inf,&u_plus_l_sq_sup,u_plus_l_inf/2,u_plus_l_sup/2,u_plus_l_inf/2,u_plus_l_sup/2);
			//elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],u_plus_l_inf,u_plus_l_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
			//double tmp1, tmp2;
			//elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,-u_plus_l_sq_inf,-u_plus_l_sq_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
			double tmp1, tmp2;
			elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i], expr->sup_coeff[i], -lb*lb, lb*lb);
			res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
			res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
		}
		else{
			res->inf_coeff[i] = 0.0;
			res->sup_coeff[i] = 0.0;
            
			double tmp1, tmp2;
            if(lb*lb < ub*ub){
			elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i], -lb*lb, ub*ub);
            }
            else{
                elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i], -ub*ub, lb*lb);
            }
			res->inf_cst = res->inf_cst + tmp1;
			res->sup_cst = res->sup_cst + tmp2;
		}
	}
	if(expr->type==SPARSE){
		res->dim = (size_t*)malloc(num_neurons*sizeof(size_t));
		for(i=0; i < num_neurons; i++){
			res->dim[i] = expr->dim[i];
		}
	}
	return res;
}

expr_t * uexpr_replace_parabola_bounds(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
	size_t num_neurons = expr->size;
	size_t i, k;
	expr_t * res = alloc_expr();
	res->inf_coeff = (double *)malloc(num_neurons*sizeof(double));
	res->sup_coeff = (double *)malloc(num_neurons*sizeof(double));
	res->inf_cst = expr->inf_cst;
	res->sup_cst = expr->sup_cst;
	res->type = expr->type;
	res->size = num_neurons;  
	for(i = 0; i < num_neurons; i++){
		if(expr->type==DENSE){
			k = i;
		}
		else{
			k = expr->dim[i];
		}
		neuron_t *neuron_k = neurons[k];
		double lb = neurons[k]->lb;
		double ub = neurons[k]->ub;
		res->inf_coeff[i] = 0.0;
		res->sup_coeff[i] = 0.0;
		if((expr->sup_coeff[i]==0) && (expr->inf_coeff[i]==0)){
			continue;
		}
		double u_plus_l_sup = (ub-lb);
                double u_plus_l_inf = -(ub-lb);
		
		
                //= (u_plus_l)*(u_plus_l)
		if(expr->inf_coeff[i]<0){
			
            double lu_sup;
            double lu_inf;
            elina_double_interval_mul_cst_coeff(pr,&lu_inf,&lu_sup,lb,-lb,-ub,ub);
			//res->coeff[i] = lambda*expr->coeff[i];
			//res->cst = res->cst + expr->coeff[i]*mu;
			elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],u_plus_l_inf,u_plus_l_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
			double tmp1, tmp2;
			elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,lu_inf,lu_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
			res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
			res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
		}
		else if (expr->sup_coeff[i]<0){
			res->inf_coeff[i] = 0;
			res->sup_coeff[i] = 0;
			//double u_plus_l_sq_inf, u_plus_l_sq_sup; 
			//u_plus_l_sq_inf = u_plus_l_inf/2;
			//u_plus_l_sq_sup = u_plus_l_sup/2;
                        //elina_double_interval_mul_cst_coeff(pr,&u_plus_l_sq_inf,&u_plus_l_sq_sup,u_plus_l_inf/2,u_plus_l_sup/2,u_plus_l_inf/2,u_plus_l_sup/2);
			//elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],u_plus_l_inf,u_plus_l_sup,expr->inf_coeff[i],expr->sup_coeff[i]);

			double tmp1, tmp2;
			elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,-lb*lb,lb*lb,expr->inf_coeff[i],expr->sup_coeff[i]);
			//double tmp =lb*lb;
			res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
			res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
		}
		else{
			res->inf_coeff[i] = 0.0;
			res->sup_coeff[i] = 0.0;
			double tmp1, tmp2;
            if(lb*lb < ub*ub){
                elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i], -lb*lb, ub*ub);
            }
            else{
                elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i], -ub*ub, lb*lb);
            }
			res->inf_cst = res->inf_cst + tmp2;
			res->sup_cst = res->sup_cst + tmp2;
		}
		
	}

	
	if(expr->type==SPARSE){
		res->dim = (size_t*)malloc(num_neurons*sizeof(size_t));
		for(i=0; i < num_neurons; i++){
			res->dim[i] = expr->dim[i];
		}
	}
	return res;
}

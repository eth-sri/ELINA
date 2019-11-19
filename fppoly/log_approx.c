#include "log_approx.h"



expr_t * lexpr_replace_log_bounds(fppoly_internal_t * pr, expr_t * expr, neuron_t ** neurons){
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
		double u_plus_l_inf = -(lb+ub);
		double u_plus_l_sup = -u_plus_l_inf;
		
		
		//= (u_plus_l)*(u_plus_l)
		if(expr->sup_coeff[i]<0){
			double one_inf = 1;
			double one_sup = -1;
			double u_plus_l_sup = ub-lb;
			double u_plus_l_inf = -(ub-lb);

			double lambda_sup = -2/u_plus_l_inf;
			double lambda_inf = -2/u_plus_l_sup;
			elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],lambda_inf,lambda_sup,expr->inf_coeff[i],expr->sup_coeff[i]);			
	
			
			fesetround(FE_DOWNWARD);
			double mu_inf = -log(-u_plus_l_inf/2);
			fesetround(FE_UPWARD);
			double mu_sup = log(u_plus_l_sup/2);			
			elina_double_interval_add_cst_coeff(pr,&mu_inf,&mu_sup,one_inf, one_sup, mu_inf, mu_sup);
			double tmp1, tmp2;
			elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,mu_inf,mu_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
			res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
			res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
		}
		else if (expr->inf_coeff[i]<0){
			double u_minus_l_sup = ub +lb;
			double u_minus_l_inf = -(ub+lb);

			if (u_minus_l_sup < 1e-9) {
			  double tmp1, tmp2;
			  double log_lb = -log(-lb), log_ub = log(-lb);
			  if (log_lb > 0) {
				log_lb = -log_lb;
				log_ub = -log_ub;
			  }
			  elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i], expr->sup_coeff[i], -log(-lb), log(-lb));
			  res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
			  res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
			} else {
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
			  elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],lambda_inf,lambda_sup,expr->inf_coeff[i],expr->sup_coeff[i]);

			  double mu_inf, mu_sup;
			  elina_double_interval_mul_cst_coeff(pr,&mu_inf,&mu_sup,-lb,lb,lambda_inf,lambda_sup);
			  elina_double_interval_add_cst_coeff(pr,&mu_inf,&mu_sup,log_l_inf, log_l_sup, mu_inf, mu_sup);

			  double tmp1, tmp2;
			  elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,mu_inf,mu_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
			  res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
			  res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
			}
		}
		else{
			
			res->inf_coeff[i] = 0.0;
			res->sup_coeff[i] = 0.0;
            double log_lb = -log(-lb);
            double log_ub = log(ub);
            double tmp1, tmp2;
            elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],log_lb,log_ub);
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

expr_t * uexpr_replace_log_bounds(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
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
		double u_plus_l_inf = -(lb+ub);
		double u_plus_l_sup = -u_plus_l_inf;
		
		//= (u_plus_l)*(u_plus_l)
		if(expr->inf_coeff[i]<0){
			double one_inf = 1;
			double one_sup = -1;
			double u_plus_l_sup = ub-lb;
			double u_plus_l_inf = -(ub-lb);

			double lambda_sup = -2/u_plus_l_inf;
			double lambda_inf = -2/u_plus_l_sup;
			elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],lambda_inf,lambda_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
	
			
			fesetround(FE_DOWNWARD);
			double mu_inf = -log(-u_plus_l_inf/2);
			fesetround(FE_UPWARD);
			double mu_sup = log(u_plus_l_sup/2);
			elina_double_interval_add_cst_coeff(pr,&mu_inf,&mu_sup,one_inf, one_sup, mu_inf, mu_sup);
			double tmp1, tmp2;
			elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,mu_inf,mu_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
			res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
			res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
		}
		else if (expr->sup_coeff[i]<0){
			double u_minus_l_sup = ub +lb;
			double u_minus_l_inf = -(ub+lb);

			if (u_minus_l_sup < 1e-9) {
			  double tmp1, tmp2;
			  double log_lb = -log(-lb), log_ub = log(-lb);
			  if (log_lb > 0) {
				log_lb = -log_lb;
				log_ub = -log_ub;
			  }
			  elina_double_interval_mul(&tmp1, &tmp2, expr->inf_coeff[i], expr->sup_coeff[i], -log(-lb), log(-lb));
			  res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
			  res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
			} else {
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
			  elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],lambda_inf,lambda_sup,expr->inf_coeff[i],expr->sup_coeff[i]);

			  double mu_inf, mu_sup;
			  elina_double_interval_mul_cst_coeff(pr,&mu_inf,&mu_sup,-lb,lb,lambda_inf,lambda_sup);
			  elina_double_interval_add_cst_coeff(pr,&mu_inf,&mu_sup,log_l_inf, log_l_sup, mu_inf, mu_sup);

			  double tmp1, tmp2;
			  elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,mu_inf,mu_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
			  res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
			  res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
			}
		}
		else{
			
			res->inf_coeff[i] = 0.0;
			res->sup_coeff[i] = 0.0;
            double log_lb = -log(-lb);
            double log_ub = log(ub);
            double tmp1, tmp2;
            elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],log_lb,log_ub);
            res->inf_cst = res->inf_cst + tmp1;
            res->sup_cst = res->sup_cst + tmp2;
			//res->inf_cst = -INFINITY;
			//res->sup_cst = INFINITY;
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


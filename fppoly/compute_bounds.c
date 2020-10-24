#ifdef GUROBI
#include <stdlib.h>
#include <stdio.h>

#include "gurobi_c.h"
#endif

#include "compute_bounds.h"
#include "math.h"

expr_t * replace_input_poly_cons_in_lexpr(fppoly_internal_t *pr, expr_t * expr, fppoly_t * fp){
	size_t dims = expr->size;
	size_t i,k;
	double tmp1, tmp2;
	expr_t * res;
	if(expr->type==DENSE){
		k = 0;
	}
	else{
		k = expr->dim[0];		
	}
	expr_t * mul_expr = NULL;
			
	if(expr->sup_coeff[0] <0){
		mul_expr = fp->input_uexpr[k];
	}
	else if(expr->inf_coeff[0] < 0){
		mul_expr = fp->input_lexpr[k];
	}
		
	if(mul_expr!=NULL){
		if(mul_expr->size==0){
			res = multiply_cst_expr(pr,mul_expr,expr->inf_coeff[0],expr->sup_coeff[0]);
		}
		else{
			res = multiply_expr(pr,mul_expr,expr->inf_coeff[0],expr->sup_coeff[0]);
		}
	}
		
	else{
		elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,expr->inf_coeff[0],expr->sup_coeff[0],fp->input_inf[k],fp->input_sup[k]);
		res = create_cst_expr(tmp1, -tmp1);
	}
	for(i=1; i < dims; i++){
		if(expr->type==DENSE){
			k = i;
		}
		else{
			k = expr->dim[i];
		}
			
		expr_t * mul_expr = NULL;
		expr_t * sum_expr = NULL;
		if(expr->sup_coeff[i] <0){
			mul_expr = fp->input_uexpr[k];
		}
		else if(expr->inf_coeff[i] <0){
			mul_expr = fp->input_lexpr[k];
		}
			
		if(mul_expr!=NULL){
			if(mul_expr->size==0){
				sum_expr = multiply_cst_expr(pr,mul_expr, expr->inf_coeff[i],expr->sup_coeff[i]);
				add_cst_expr(pr,res,sum_expr);
			}	
			else if(expr->inf_coeff[i]!=0 && expr->sup_coeff[i]!=0){
				sum_expr = multiply_expr(pr,mul_expr, expr->inf_coeff[i],expr->sup_coeff[i]);
				add_expr(pr,res,sum_expr);
			}
				//free_expr(mul_expr);
			if(sum_expr!=NULL){
				free_expr(sum_expr);
			}
		}
		else{
			elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],fp->input_inf[k],fp->input_sup[k]);
			res->inf_cst = res->inf_cst + tmp1;
			res->sup_cst = res->sup_cst - tmp1;
		}
	}
		
	res->inf_cst = res->inf_cst + expr->inf_cst; 
	res->sup_cst = res->sup_cst + expr->sup_cst; 
	return res;
}


expr_t * replace_input_poly_cons_in_uexpr(fppoly_internal_t *pr, expr_t * expr, fppoly_t * fp){
	size_t dims = expr->size;
	size_t i,k;
	double tmp1, tmp2;
	expr_t * res;
	if(expr->type==DENSE){
		k = 0;
	}
	else{
		k = expr->dim[0];		
	}
	expr_t * mul_expr = NULL;
			
	if(expr->sup_coeff[0] <0){
		mul_expr = fp->input_lexpr[k];
	}
	else if(expr->inf_coeff[0] < 0){
		mul_expr = fp->input_uexpr[k];
	}
		
	if(mul_expr!=NULL){
		if(mul_expr->size==0){
			res = multiply_cst_expr(pr,mul_expr,expr->inf_coeff[0],expr->sup_coeff[0]);
		}
		else{
			res = multiply_expr(pr,mul_expr,expr->inf_coeff[0],expr->sup_coeff[0]);
		}
	}
	else{
		elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,expr->inf_coeff[0],expr->sup_coeff[0],fp->input_inf[k],fp->input_sup[k]);
		res = create_cst_expr(-tmp2, tmp2);
	}
                //printf("finish\n");
		//fflush(stdout);
	for(i=1; i < dims; i++){
		if(expr->type==DENSE){
			k = i;
		}
		else{
			k = expr->dim[i];
		}
		expr_t * mul_expr = NULL;
		expr_t * sum_expr = NULL;
		if(expr->sup_coeff[i] <0){
			mul_expr = fp->input_lexpr[k];
		}
		else if(expr->inf_coeff[i] <0){
			mul_expr = fp->input_uexpr[k];
		}
			
		if(mul_expr!=NULL){
			if(mul_expr->size==0){
				sum_expr = multiply_cst_expr(pr,mul_expr, expr->inf_coeff[i],expr->sup_coeff[i]);
				add_cst_expr(pr,res,sum_expr);
			}	
			else if(expr->inf_coeff[i]!=0 && expr->sup_coeff[i]!=0){
				sum_expr = multiply_expr(pr,mul_expr, expr->inf_coeff[i],expr->sup_coeff[i]);
				add_expr(pr,res,sum_expr);
			}
				//free_expr(mul_expr);
			if(sum_expr!=NULL){
				free_expr(sum_expr);
			}
		}
		else{
			elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],fp->input_inf[k],fp->input_sup[k]);
			res->inf_cst = res->inf_cst - tmp2;
			res->sup_cst = res->sup_cst + tmp2;
		}
	}
	res->inf_cst = res->inf_cst + expr->inf_cst; 
	res->sup_cst = res->sup_cst + expr->sup_cst; 
	return res;
}

#ifdef GUROBI
void handle_gurobi_error(int error, GRBenv *env) {
    if (error) {
        printf("Gurobi error: %s\n", GRBgeterrormsg(env));
        exit(1);
    }
}

double substitute_spatial_gurobi(expr_t *expr, fppoly_t *fp, const int opt_sense) {

    GRBenv *env = NULL;
    GRBmodel *model = NULL;

    int error;

    error = GRBemptyenv(&env);
    handle_gurobi_error(error, env);
    error = GRBsetintparam(env, "OutputFlag", 0);
    handle_gurobi_error(error, env);
    error = GRBsetintparam(env, "NumericFocus", 2);
    handle_gurobi_error(error, env);
    error = GRBstartenv(env);
    handle_gurobi_error(error, env);

    double *lb, *ub, *obj;
    const size_t dims = expr->size;
    const size_t numvars = 3 * dims;

    lb = malloc(numvars * sizeof(double));
    ub = malloc(numvars * sizeof(double));
    obj = malloc(numvars * sizeof(double));

    for (size_t i = 0; i < dims; ++i) {
        const size_t k = expr->type == DENSE ? i : expr->dim[i];
        lb[i] = -fp->input_inf[k];
        ub[i] = fp->input_sup[k];
        obj[i] = opt_sense == GRB_MINIMIZE ? -expr->inf_coeff[i] : expr->sup_coeff[i];

        for (size_t j = 0; j < 2; ++j) {
            const size_t l = fp->input_uexpr[k]->dim[j];
            lb[dims + 2 * i + j] = -fp->input_inf[l];
            ub[dims + 2 * i + j] = fp->input_sup[l];
            obj[dims + 2 * i + j] = 0;
        }
    }

    error = GRBnewmodel(env, &model, NULL, numvars, obj, lb, ub, NULL, NULL);
    handle_gurobi_error(error, env);
    error = GRBsetintattr(model, "ModelSense", opt_sense);
    handle_gurobi_error(error, env);

    for (size_t i = 0; i < dims; ++i) {
        const size_t k = expr->type == DENSE ? i : expr->dim[i];

        int ind[] = {i, dims + 2 * i, dims + 2 * i + 1};

        double lb_val[] = {
            -1, -fp->input_lexpr[k]->inf_coeff[0], -fp->input_lexpr[k]->inf_coeff[1]
        };
        error = GRBaddconstr(model, 3, ind, lb_val, GRB_LESS_EQUAL, fp->input_lexpr[k]->inf_cst, NULL);
        handle_gurobi_error(error, env);

        double ub_val[] = {
            1, -fp->input_uexpr[k]->sup_coeff[0], -fp->input_uexpr[k]->sup_coeff[1]
        };
        error = GRBaddconstr(model, 3, ind, ub_val, GRB_LESS_EQUAL, fp->input_uexpr[k]->sup_cst, NULL);
        handle_gurobi_error(error, env);
    }

    size_t idx, nbr, s_idx, s_nbr;
    const size_t num_pixels = fp->num_pixels;

    for (size_t i = 0; i < fp->spatial_size; ++i) {
        idx = fp->spatial_indices[i];
        nbr = fp->spatial_neighbors[i];

        if (expr->type == DENSE) {
            s_idx = idx;
            s_nbr = nbr;
        } else {
            s_idx = s_nbr = num_pixels;

            for (size_t j = 0; j < dims; ++j) {
                if (expr->dim[j] == idx) {
                    s_idx = j;
                }
                if (expr->dim[j] == nbr) {
                    s_nbr = j;
                }
            }

            if ((s_idx == num_pixels) || (s_nbr == num_pixels)) {
                continue;
            }
        }

        int ind_x[] = {dims + 2 * s_idx, dims + 2 * s_nbr};
        int ind_y[] = {dims + 2 * s_idx + 1, dims + 2 * s_nbr + 1};
        double val[] = {1., -1.};

        error = GRBaddconstr(model, 2, ind_x, val, GRB_LESS_EQUAL, fp->spatial_gamma, NULL);
        handle_gurobi_error(error, env);
        error = GRBaddconstr(model, 2, ind_y, val, GRB_LESS_EQUAL, fp->spatial_gamma, NULL);
        handle_gurobi_error(error, env);
        error = GRBaddconstr(model, 2, ind_x, val, GRB_GREATER_EQUAL, -fp->spatial_gamma, NULL);
        handle_gurobi_error(error, env);
        error = GRBaddconstr(model, 2, ind_y, val, GRB_GREATER_EQUAL, -fp->spatial_gamma, NULL);
        handle_gurobi_error(error, env);
    }

    int opt_status;
    double obj_val;

    error = GRBoptimize(model);
    handle_gurobi_error(error, env);
    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &opt_status);
    handle_gurobi_error(error, env);
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &obj_val);
    handle_gurobi_error(error, env);

    if (opt_status != GRB_OPTIMAL) {
        printf("Gurobi model status not optimal %i\n", opt_status);
        exit(1);
    }

    free(lb);
    free(ub);
    free(obj);

    GRBfreemodel(model);
    GRBfreeenv(env);

    return obj_val;
}
#endif

double compute_lb_from_expr(fppoly_internal_t *pr, expr_t * expr, fppoly_t * fp, int layerno){

#ifdef GUROBI
    if ((fp->input_lexpr!=NULL) && (fp->input_uexpr!=NULL) && layerno==-1 && fp->spatial_size > 0) {
        return expr->inf_cst - substitute_spatial_gurobi(expr, fp, GRB_MINIMIZE);
    }
#endif

	size_t i,k;
	double tmp1, tmp2;
        //printf("start\n");
        //fflush(stdout);
	if((fp->input_lexpr!=NULL) && (fp->input_uexpr!=NULL) && layerno==-1){
		expr =  replace_input_poly_cons_in_lexpr(pr, expr, fp);
	}
        //expr_print(expr);
	//fflush(stdout);
	size_t dims = expr->size;
	double res_inf = expr->inf_cst;
	if(expr->inf_coeff==NULL || expr->sup_coeff==NULL){
		return res_inf;
	}
	for(i=0; i < dims; i++){
		//if(expr->inf_coeff[i]<0){
		if(expr->type==DENSE){
			k = i;
		}
		else{
			k = expr->dim[i];
		}
			if(layerno==-1){
				elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],fp->input_inf[k],fp->input_sup[k]);
			}
			else{
				elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],fp->layers[layerno]->neurons[k]->lb,fp->layers[layerno]->neurons[k]->ub);
			}
			//printf("tmp1: %g\n",tmp1);
			res_inf = res_inf + tmp1;
			
	}
//	printf("inf: %g\n",-res_inf);
//	fflush(stdout);
        if(fp->input_lexpr!=NULL && fp->input_uexpr!=NULL && layerno==-1){
		free_expr(expr);
	}
        //printf("finish\n");
        //fflush(stdout);
	return res_inf;
}

double compute_ub_from_expr(fppoly_internal_t *pr, expr_t * expr, fppoly_t * fp, int layerno){

#ifdef GUROBI
    if ((fp->input_lexpr!=NULL) && (fp->input_uexpr!=NULL) && layerno==-1 && fp->spatial_size > 0) {
        return expr->sup_cst + substitute_spatial_gurobi(expr, fp, GRB_MAXIMIZE);
    }
#endif

	size_t i,k;
	double tmp1, tmp2;

	if((fp->input_lexpr!=NULL) && (fp->input_uexpr!=NULL) && layerno==-1){
		expr =  replace_input_poly_cons_in_uexpr(pr, expr, fp);
	}

	size_t dims = expr->size;
	double res_sup = expr->sup_cst;
	if(expr->inf_coeff==NULL || expr->sup_coeff==NULL){
		return res_sup;
	}
	for(i=0; i < dims; i++){
		//if(expr->inf_coeff[i]<0){
		if(expr->type==DENSE){
			k = i;
		}
		else{
			k = expr->dim[i];
		}		
		if(layerno==-1){
			elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],fp->input_inf[k],fp->input_sup[k]);
		}
		else{
			elina_double_interval_mul(&tmp1,&tmp2,expr->inf_coeff[i],expr->sup_coeff[i],fp->layers[layerno]->neurons[k]->lb,fp->layers[layerno]->neurons[k]->ub);
		}
		res_sup = res_sup + tmp2;
			
	}
	//printf("sup: %g\n",res_sup);
	//fflush(stdout);
	if(fp->input_lexpr!=NULL && fp->input_uexpr!=NULL && layerno==-1){
		free_expr(expr);
	}
	return res_sup;
}


double get_lb_using_predecessor_layer(fppoly_internal_t * pr,fppoly_t *fp, expr_t **lexpr_ptr, int k){
	expr_t * tmp_l;
	neuron_t ** aux_neurons = fp->layers[k]->neurons;
	expr_t *lexpr = *lexpr_ptr;
	double res = INFINITY;
	res = compute_lb_from_expr(pr,lexpr,fp,k);
	tmp_l = lexpr;
	*lexpr_ptr = lexpr_replace_bounds(pr,lexpr,aux_neurons, fp->layers[k]->is_activation);
	free_expr(tmp_l);
	return res;
}


double get_ub_using_predecessor_layer(fppoly_internal_t * pr,fppoly_t *fp, expr_t **uexpr_ptr, int k){
	expr_t * tmp_u;
	neuron_t ** aux_neurons = fp->layers[k]->neurons;
	expr_t *uexpr = *uexpr_ptr;
	double res = INFINITY;
	tmp_u = uexpr;
	res = compute_ub_from_expr(pr,uexpr,fp,k);
	*uexpr_ptr = uexpr_replace_bounds(pr,uexpr,aux_neurons, fp->layers[k]->is_activation);
	free_expr(tmp_u);
	return res;
}

double get_lb_using_previous_layers(elina_manager_t *man, fppoly_t *fp, expr_t *expr, size_t layerno){
	size_t i;
	int k;
	//size_t numlayers = fp->numlayers;
	//printf("COMING HERE\n");
	//fflush(stdout);
	expr_t * lexpr = copy_expr(expr);
        fppoly_internal_t * pr = fppoly_init_from_manager(man,ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
	if(fp->numlayers==layerno){
		
		k = layerno-1;
	}
	else if((fp->layers[layerno]->is_concat == true) || (fp->layers[layerno]->num_predecessors==2)){
		k = layerno;
	}
	else{
		k = fp->layers[layerno]->predecessors[0]-1;
	}	
	double res = INFINITY;
	while(k >=0){
	        if(fp->layers[k]->is_concat==true){
		//	expr_print(lexpr);
		//	fflush(stdout);
			//printf("k: %zu\n", k);
			size_t i;
			size_t *predecessors = fp->layers[k]->predecessors;
			size_t num_predecessors = fp->layers[k]->num_predecessors;
			int common_predecessor = INT_MAX;
			expr_t ** sub_expr = (expr_t**)malloc(num_predecessors*sizeof(expr_t*));
			size_t index_start = 0;
			for(i=0; i < num_predecessors; i++){
				int pred = predecessors[i]-1;
				size_t num_neurons = fp->layers[pred]->dims;
				if(pred < common_predecessor){
					common_predecessor = pred;
				}
				sub_expr[i] = extract_subexpr(lexpr,index_start, num_neurons);
			printf("index start %zu %zu %zu\n", i,index_start,num_neurons);
			fflush(stdout);
				index_start = index_start + num_neurons;
			}
			printf("common %zu %zu\n", sub_expr[0]->size, sub_expr[1]->size);
			expr_print(lexpr);
			expr_print(sub_expr[0]);
			expr_print(sub_expr[1]);
			fflush(stdout);
			for(i=0; i < num_predecessors; i++){
				int iter = predecessors[i]-1;
				if(sub_expr[i]->size>0){
					while(iter!=common_predecessor){
						get_lb_using_predecessor_layer(pr,fp, &sub_expr[i],  iter);
						printf("iter %zu %d\n",sub_expr[i]->size, iter);
						fflush(stdout);
						iter = fp->layers[iter]->predecessors[0]-1;
					}
				}
			}
			double inf_cst = lexpr->inf_cst;
			double sup_cst = lexpr->sup_cst;
			free_expr(lexpr);
			bool flag = true;
			//lexpr = copy_expr(sub_expr[0]);
			for(i=0; i < num_predecessors; i++){
				if(sub_expr[i]->size>0){
					if(flag==true){
						lexpr = copy_expr(sub_expr[i]);
						flag = false;
					}
					else{
		//				sort_sparse_expr(lexpr);
						printf("ADDING %zu %zu\n", lexpr->size, sub_expr[i]->size);
						fflush(stdout);
						add_expr(pr, lexpr, sub_expr[i]);
						printf("after adding: %zu\n",lexpr->size);
						fflush(stdout);
		//				if(lexpr->size!=sub_expr[i]->size)
					}				//printf("sizes: %zu %zu\n", lexpr->size,sub_expr[i]->size);
					
				}
		//		if(lexpr->size!=sub_expr[i]->size)
		//			printf("sizes %zu %zu\n", lexpr->size, sub_expr[i]->size);
				free_expr(sub_expr[i]);
			}	
			lexpr->inf_cst = lexpr->inf_cst + inf_cst;
			lexpr->sup_cst = lexpr->sup_cst + sup_cst;
			//free_expr(sub_expr[0]);
			free(sub_expr);
			k = common_predecessor;
			//printf("IS ACTIVATION: %d\n",fp->layers[k]->is_activation);
			//expr_print(lexpr);
			//fflush(stdout);
		}
		else if(fp->layers[k]->num_predecessors==2){
				expr_t * lexpr_copy = copy_expr(lexpr);
				lexpr_copy->inf_cst = 0;
				lexpr_copy->sup_cst = 0;
				size_t predecessor1 = fp->layers[k]->predecessors[0]-1;
				size_t predecessor2 = fp->layers[k]->predecessors[1]-1;
				
				char * predecessor_map = (char *)calloc(k,sizeof(char));
				// Assume no nested residual layers
				int iter = fp->layers[predecessor1]->predecessors[0]-1;
				while(iter>=0){
					predecessor_map[iter] = 1;
					iter = fp->layers[iter]->predecessors[0]-1;
				}
				iter =  fp->layers[predecessor2]->predecessors[0]-1;
				int common_predecessor = 0;
				while(iter>=0){
					if(predecessor_map[iter] == 1){
						common_predecessor = iter;
						break;
					}
					iter = fp->layers[iter]->predecessors[0]-1;
				}
				
				iter = predecessor1;
				while(iter!=common_predecessor){
					get_lb_using_predecessor_layer(pr,fp, &lexpr,  iter);
					iter = fp->layers[iter]->predecessors[0]-1;
				}
				iter =  predecessor2;
				while(iter!=common_predecessor){
					get_lb_using_predecessor_layer(pr,fp, &lexpr_copy,  iter);
					iter = fp->layers[iter]->predecessors[0]-1;					
				}
				free(predecessor_map);
				add_expr(pr,lexpr,lexpr_copy);
				
				free_expr(lexpr_copy);
				
				// Assume at least one non-residual layer between two residual layers
				k = common_predecessor;		
				continue;
			}
			else {
								
				 res =fmin(res,get_lb_using_predecessor_layer(pr,fp, &lexpr, k));
				 k = fp->layers[k]->predecessors[0]-1;
				
			}
			
	}
		
	res = fmin(res,compute_lb_from_expr(pr,lexpr,fp,-1)); 
        free_expr(lexpr);
        //if(fp->layers[layerno]->is_concat == true){
		printf("res: %g\n",res);
		fflush(stdout);
		//expr_print(expr);
		fflush(stdout);
	//}
	return res;
	
}


double get_ub_using_previous_layers(elina_manager_t *man, fppoly_t *fp, expr_t *expr, size_t layerno){
	size_t i;
	int k;
	//size_t numlayers = fp->numlayers;
	expr_t * uexpr = copy_expr(expr);
        fppoly_internal_t * pr = fppoly_init_from_manager(man,ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
        
	if(fp->numlayers==layerno){
		k = layerno-1;
	}
	else if((fp->layers[layerno]->is_concat == true) || (fp->layers[layerno]->num_predecessors==2)){
		k = layerno;
	}
	else{
		k = fp->layers[layerno]->predecessors[0]-1;
	}	
	double res =INFINITY;
	while(k >=0){
		if(fp->layers[k]->is_concat==true){
                        //sort_expr(lexpr);
                        size_t i;
                        size_t *predecessors = fp->layers[k]->predecessors;
                        size_t num_predecessors = fp->layers[k]->num_predecessors;
                        int common_predecessor = INT_MAX;
                        expr_t ** sub_expr = (expr_t**)malloc(num_predecessors*sizeof(expr_t*));
                        size_t index_start = 0;
                        for(i=0; i < num_predecessors; i++){
                                int pred = predecessors[i]-1;
                                size_t num_neurons = fp->layers[pred]->dims;
                                if(pred < common_predecessor){
                                        common_predecessor = pred;
                                }
                                sub_expr[i] = extract_subexpr(uexpr,index_start, num_neurons);
                                index_start = index_start + num_neurons;
                        }
                        for(i=0; i < num_predecessors; i++){
                                int iter = predecessors[i]-1;
				if(sub_expr[i]->size>0){
                                	while(iter!=common_predecessor){
                                        	get_ub_using_predecessor_layer(pr,fp, &sub_expr[i],  iter);
                                        	iter = fp->layers[iter]->predecessors[0]-1;
                                	}
				}
                        }
			double inf_cst = uexpr->inf_cst;
			double sup_cst = uexpr->sup_cst;
                        free_expr(uexpr);
			bool flag = true;
                        for(i=0; i < num_predecessors; i++){
				if(sub_expr[i]->size>0){
					if(flag==true){
						uexpr = copy_expr(sub_expr[i]);
						flag = false;
					}
					else{
		//				sort_sparse_expr(uexpr);
                                		add_expr(pr, uexpr, sub_expr[i]);
					}
				}
                                free_expr(sub_expr[i]);
                        }
                        //free_expr(sub_expr[0]);
                        free(sub_expr);
			uexpr->inf_cst = uexpr->inf_cst + inf_cst;
			uexpr->sup_cst = uexpr->sup_cst + sup_cst;
                        k = common_predecessor;
                }

		
		else if(fp->layers[k]->num_predecessors==2){
				expr_t * uexpr_copy = copy_expr(uexpr);
				uexpr_copy->inf_cst = 0;
				uexpr_copy->sup_cst = 0;
				size_t predecessor1 = fp->layers[k]->predecessors[0]-1;
				size_t predecessor2 = fp->layers[k]->predecessors[1]-1;
				
				char * predecessor_map = (char *)calloc(k,sizeof(char));
				// Assume no nested residual layers
				int iter = fp->layers[predecessor1]->predecessors[0]-1;
				while(iter>=0){
					predecessor_map[iter] = 1;
					iter = fp->layers[iter]->predecessors[0]-1;
				}
				iter =  fp->layers[predecessor2]->predecessors[0]-1;
				int common_predecessor = 0;
				while(iter>=0){
					if(predecessor_map[iter] == 1){
						common_predecessor = iter;
						break;
					}
					iter = fp->layers[iter]->predecessors[0]-1;
				}
				
				iter = predecessor1;
				while(iter!=common_predecessor){
					get_ub_using_predecessor_layer(pr,fp, &uexpr,  iter);
					iter = fp->layers[iter]->predecessors[0]-1;
				}
				iter =  predecessor2;
				while(iter!=common_predecessor){
					get_ub_using_predecessor_layer(pr,fp, &uexpr_copy,  iter);
					iter = fp->layers[iter]->predecessors[0]-1;					
				}
				free(predecessor_map);
				add_expr(pr,uexpr,uexpr_copy);
				
				free_expr(uexpr_copy);
				
				// Assume at least one non-residual layer between two residual layers
				k = common_predecessor;
				
				continue;
			}
			else {
				
				 res= fmin(res,get_ub_using_predecessor_layer(pr,fp, &uexpr, k));
				 k = fp->layers[k]->predecessors[0]-1;
				 
			}
			
	}
		
	res = fmin(res,compute_ub_from_expr(pr,uexpr,fp,-1)); 
        printf("UPPER BOUND: %g\n",res);
        fflush(stdout);
        free_expr(uexpr);
	return res;
	
}

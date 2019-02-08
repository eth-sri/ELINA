/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
 *  This software is distributed under GNU Lesser General Public License Version 3.0.
 *  For more information, see the ELINA project website at:
 *  http://elina.ethz.ch
 *
 *  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
 *  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
 *  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
 *  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 *  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY     
 *  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
 *  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
 *  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
 *  CONTRACT, TORT OR OTHERWISE).
 *
 */

#include "zonoml_reduced_product.h"

/* transfer information from zonotope to octagon */
elina_lincons0_array_t get_meet_lincons_array(elina_dim_t y, double inf_l, double inf_u, double sup_l, double sup_u){
	elina_lincons0_array_t lincons = elina_lincons0_array_make(2);
	lincons.p[0].constyp = ELINA_CONS_SUPEQ;
	elina_linexpr0_t * expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
	elina_coeff_set_interval_double(&expr->cst,-inf_u,-inf_l);
	elina_linterm_t * lterm = &expr->p.linterm[0];
	lterm->dim = y;
	elina_coeff_set_interval_double(&lterm->coeff,1,1); 
	lincons.p[0].linexpr0 = expr;
	
	lincons.p[1].constyp = ELINA_CONS_SUPEQ;
	expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
	elina_coeff_set_interval_double(&expr->cst,sup_l, sup_u);
	lterm = &expr->p.linterm[0];
	lterm->dim = y;
	elina_coeff_set_interval_double(&lterm->coeff,-1,-1); 
	lincons.p[1].linexpr0 = expr;	
	
	return lincons;
}


zonotope_aff_t * create_affine_form_for_box(zonotope_internal_t *pr, double inf, double sup, zonotope_noise_symbol_t *epsilon, char *is_used){
	double mid_inf = 0.0;
    	double mid_sup = 0.0;
    	double dev_inf = 0.0;
    	double dev_sup = 0.0;
	zonotope_aff_t * res_box = zonotope_aff_alloc_init(pr);
	elina_interval_middev(&mid_inf, &mid_sup, &dev_inf, &dev_sup, inf,sup);
	res_box->c_inf = res_box->c_inf + mid_inf;
	res_box->c_sup = res_box->c_sup + mid_sup;
	res_box->itv_inf = inf;
	res_box->itv_sup = sup;
				
   	if (dev_inf!=0 || dev_sup!=0) {
	    zonotope_aaterm_t* ptr = zonotope_aaterm_alloc_init();
	    ptr->inf = dev_inf;
	    ptr->sup = dev_sup;
	    ptr->pnsym = epsilon;
	    *is_used = 1;
            if (res_box->end) res_box->end->n = ptr;
	    else res_box->q = ptr;
	    res_box->end = ptr;
	    res_box->l++;
        }
	return res_box;
}


elina_abstract0_t *relu_zono(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t x){
	
	//zonoml_fprint(stdout,man,zo,NULL);
	//elina_lincons0_array_t arr = zonoml_to_lincons_array(man,zo);
	//elina_lincons0_array_fprint(stdout,&arr,NULL);
	//elina_lincons0_array_clear(&arr);
	zonotope_t * input = zonotope_of_abstract0(abs);
	zonotope_t *zo = destructive? input : zonotope_copy(man,input);
	if(zonotope_is_bottom(man,zo)){
		
		return abstract0_of_zonotope(man,zo);
	}
	zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
	//elina_abstract0_t * oct = zo->oct;
	elina_dimension_t dims = zonotope_dimension(man,zo);
	elina_dim_t num_dim = dims.intdim + dims.realdim;
	elina_interval_t * bound_x = zonotope_bound_dimension(man,zo,x);
	double sup = bound_x->sup->val.dbl;
	double inf = bound_x->inf->val.dbl;
	elina_interval_free(bound_x);
	if(sup<=0){
		elina_linexpr0_t * assign_expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0); 
		elina_coeff_set_scalar_double(&assign_expr->cst,0);
		//elina_linterm_t *lterm = &assign_expr->p.linterm[0];
		//lterm->dim = x;
		//elina_coeff_set_scalar_double(&lterm->coeff,1); 
		zo = zonotope_assign_linexpr_array(man,true,zo,&x,&assign_expr,1,NULL);
		return abstract0_of_zonotope(man,zo);		
	}
	else if(inf>0){
		return abstract0_of_zonotope(man,zo);
	}
	else{
		elina_dimchange_t * add_dim = elina_dimchange_alloc(0,1);
		elina_dim_t y = num_dim;
		add_dim->dim[0] = y;
		zo = zonotope_add_dimensions(man,true,zo,add_dim,false);
		
		double inf_l = -inf;
		double inf_u = inf;
		double sup_l = -sup;
		double sup_u = sup;
		
			
		double width_l = sup_l + inf_u;
		double width_u = sup_u + inf_l;
		//if(-inf>sup){
			//flag = 0;
			//count++;
			//relu_zono_neg(man,zo,y,x);
			
			//return;
		//}
		elina_lincons0_array_t lincons;
		double bound_l, bound_u;
		elina_double_interval_mul(&bound_l, &bound_u, inf_l, inf_u, sup_l, sup_u);
		double tmp = bound_l;
		bound_l = bound_u;
		bound_u = tmp;
		elina_double_interval_div(&bound_l, &bound_u, bound_l, bound_u, width_l, width_u);
		
		double lambda_l, lambda_u;
		elina_double_interval_div(&lambda_l, &lambda_u, sup_l, sup_u, width_l, width_u);
		//double bound = (-inf*sup)/width;
		
		if((-inf>sup)&&-bound_l>1){
		//if(0){
			double alpha_l, alpha_u;
			elina_double_interval_div(&alpha_l,&alpha_u,-1,1,bound_l,bound_u);
			//double alpha_inf = -1/bound;
			//double alpha_sup = 1/bound;
			double beta_l = -1 + alpha_u;
			double beta_u = 1 + alpha_l;
			elina_dimchange_t * dimchange = elina_dimchange_alloc(0,2);
			dimchange->dim[0] = num_dim+1;
			dimchange->dim[1] = num_dim+1;
			zo = zonotope_add_dimensions(man,true,zo,dimchange,false);
			lincons = get_meet_lincons_array(num_dim+1,0,0, -sup_l, sup_u);
			zo = zonotope_meet_lincons_array(man,true,zo,&lincons);
			elina_lincons0_array_clear(&lincons);
			lincons = get_meet_lincons_array(num_dim+2,0,0,-bound_l, bound_u);
			zo = zonotope_meet_lincons_array(man,true,zo,&lincons);
			elina_lincons0_array_clear(&lincons);
			
			elina_dim_t z = num_dim+2;
			elina_linterm_t * lterm;
			elina_linexpr0_t * assign_expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2); 
			elina_coeff_set_interval_double(&assign_expr->cst,0,0);
			lterm = &assign_expr->p.linterm[0];
			lterm->dim = z;
			elina_coeff_set_interval_double(&lterm->coeff,1,1); 
			lterm = &assign_expr->p.linterm[1];
			lterm->dim = x;
			elina_coeff_set_interval_double(&lterm->coeff,-lambda_l, lambda_u); 
		//elina_dim_t extra_var = num_dim +1;		
			zo = zonotope_assign_linexpr_array(man,true,zo,&z,&assign_expr,1,NULL);
			
			lterm = &assign_expr->p.linterm[0];
			lterm->dim = num_dim+1;
			elina_coeff_set_interval_double(&lterm->coeff,-beta_l, beta_u); 
			lterm = &assign_expr->p.linterm[1];
			lterm->dim = z;
			elina_coeff_set_interval_double(&lterm->coeff,-alpha_l,alpha_u); 
		//elina_dim_t extra_var = num_dim +1;		
			zo = zonotope_assign_linexpr_array(man,true,zo,&y,&assign_expr,1,NULL);
					
			elina_linexpr0_free(assign_expr);
			zo = zonotope_remove_dimensions(man,true,zo,dimchange);
			elina_dimchange_free(dimchange);
			
			//return abstract0_of_zonotope(man,zo);
		}
		//size_t z;
		
		
		//elina_lincons0_array_t lincons = elina_lincons0_array_make(4*size*num_dim - 2*size);
		//size_t k = 0;
		
		//double min_width = (sup - ((sup*inf)/(sup-inf)));
		
		else{
			
			lincons = get_meet_lincons_array(y,0,0,-bound_l,bound_u);
			zo = zonotope_meet_lincons_array(man,true,zo,&lincons);
			elina_lincons0_array_clear(&lincons);	
			elina_linterm_t * lterm;
			elina_linexpr0_t * assign_expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2); 
			elina_coeff_set_interval_double(&assign_expr->cst,0,0);
			lterm = &assign_expr->p.linterm[0];
			lterm->dim = y;
			elina_coeff_set_interval_double(&lterm->coeff,1,1); 
			lterm = &assign_expr->p.linterm[1];
			lterm->dim = x;
			//printf("div: %g %g\n", lambda_l, lambda_u);
			//fflush(stdout);
			elina_coeff_set_interval_double(&lterm->coeff,-lambda_l, lambda_u); 
			//elina_dim_t extra_var = num_dim +1;
			//elina_linexpr0_fprint(stdout,assign_expr,NULL);
			//fflush(stdout);		
			zo = zonotope_assign_linexpr_array(man,true,zo,&y,&assign_expr,1,NULL);
			elina_linexpr0_free(assign_expr);
		}
		elina_dimperm_t* dimperm = elina_dimperm_alloc(num_dim+1);
		elina_dimperm_set_id(dimperm);
		dimperm->dim[x] = num_dim;
		dimperm->dim[num_dim] = x;
		zo = zonotope_permute_dimensions(man,true,zo,dimperm);
		zo = zonotope_remove_dimensions(man,true,zo,add_dim);
		elina_dimchange_free(add_dim);
		elina_dimperm_free(dimperm);
		return abstract0_of_zonotope(man,zo);
	}
}		
	
void * handle_relu_zono_parallel(void *args){
	zonoml_relu_thread_t * data = (zonoml_relu_thread_t *)args;
	elina_manager_t * man = data->man;
	//elina_abstract0_t * abs = data->abs;
	zonotope_t * zo = data->z;
	zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);


	size_t start_offset = data->start_offset;
	size_t idx_start = data->start;
	size_t idx_end = data->end;
	
	zonotope_noise_symbol_t **epsilon_map = data->epsilon_map;
	zonotope_noise_symbol_t **epsilon_map_extra = data->epsilon_map_extra;
	size_t offset = start_offset + idx_start;
	size_t i;
	for (i=idx_start; i< idx_end; i++) {
        //printf("i: %zu\n",i);
        //fflush(stdout);
		elina_interval_t * bound_x = zonotope_bound_dimension(man,zo,offset);
		double sup = bound_x->sup->val.dbl;
		double inf = bound_x->inf->val.dbl;
		elina_interval_free(bound_x);
		if(sup<=0){
			zonotope_aff_check_free(pr, zo->paf[offset]);
			zonotope_aff_t * res = zonotope_aff_alloc_init(pr);
    			res->c_inf = 0.0;
    			res->c_sup = 0.0;
    			res->itv_inf = 0.0;
    			res->itv_sup = 0.0;
			zo->paf[offset] = res;
			zo->paf[offset]->pby++;
			zo->box_inf[offset] = 0.0;
			zo->box_sup[offset] = 0.0;
		}
		else if(inf>=0){
			
		}
		else{
			//zonotope_aff_check_free(pr, zo->paf[offset]);
			double inf_l = -inf;
			double inf_u = inf;
			double sup_l = -sup;
			double sup_u = sup;
		
			
			double width_l = sup_l + inf_u;
			double width_u = sup_u + inf_l;
		
			
			double bound_l, bound_u;
			elina_double_interval_mul(&bound_l, &bound_u, inf_l, inf_u, sup_l, sup_u);
			double tmp = bound_l;
			bound_l = bound_u;
			bound_u = tmp;
			elina_double_interval_div(&bound_l, &bound_u, bound_l, bound_u, width_l, width_u);
		
			double lambda_l, lambda_u;
			elina_double_interval_div(&lambda_l, &lambda_u, sup_l, sup_u, width_l, width_u);
			
			
			
			if((-inf>sup)&&-bound_l>1){
				
				double alpha_l, alpha_u;
				elina_double_interval_div(&alpha_l,&alpha_u,-1,1,bound_l,bound_u);
				
				double beta_l = -1 + alpha_u;
				double beta_u = 1 + alpha_l;
				
				double mid_inf = 0.0;
    				double mid_sup = 0.0;
    				double dev_inf = 0.0;
    				double dev_sup = 0.0;
				zonotope_aff_t * res_box = zonotope_aff_alloc_init(pr);
				elina_interval_middev(&mid_inf, &mid_sup, &dev_inf, &dev_sup, 0,sup_u);
				res_box->c_inf = res_box->c_inf + mid_inf;
				res_box->c_sup = res_box->c_sup + mid_sup;
				res_box->itv_inf = 0.0;
				res_box->itv_sup = sup_u;
				
   				 if (dev_inf!=0 || dev_sup!=0) {
					zonotope_aaterm_t* ptr = zonotope_aaterm_alloc_init();
					ptr->inf = dev_inf;
					ptr->sup = dev_sup;
					ptr->pnsym = epsilon_map_extra[i];
                     			if (res_box->end) res_box->end->n = ptr;
					else res_box->q = ptr;
					res_box->end = ptr;
					res_box->l++;
    				}
				
				elina_interval_t *beta = elina_interval_alloc();
				elina_interval_set_double(beta,-beta_l,beta_u);
            			zonotope_aff_t *tmp = res_box;
				res_box = zonotope_aff_mul_itv(pr, tmp, beta);
				zonotope_aff_check_free(pr, tmp);
				elina_interval_free(beta);


				elina_interval_t *lambda = elina_interval_alloc();
				elina_interval_set_double(lambda,-lambda_l,lambda_u);
				zonotope_aff_t * res_zono = zonotope_aff_mul_itv(pr, zo->paf[offset], lambda);
				zonotope_aff_check_free(pr, zo->paf[offset]);
				elina_interval_free(lambda);
				
				elina_interval_middev(&mid_inf, &mid_sup, &dev_inf, &dev_sup, 0,bound_u);
				res_zono->c_inf = res_zono->c_inf + mid_inf;
				res_zono->c_sup = res_zono->c_sup + mid_sup;
				
				
   				 if (dev_inf!=0 || dev_sup!=0) {
					zonotope_aaterm_t* ptr = zonotope_aaterm_alloc_init();
					ptr->inf = dev_inf;
					ptr->sup = dev_sup;
					ptr->pnsym = epsilon_map[i];
                     			if (res_zono->end) res_zono->end->n = ptr;
					else res_zono->q = ptr;
					res_zono->end = ptr;
					res_zono->l++;
    				}
								
				res_zono->itv_sup+= bound_u;

				elina_interval_t *alpha = elina_interval_alloc();
				elina_interval_set_double(alpha,-alpha_l,alpha_u);
            			tmp = res_zono;
				res_zono = zonotope_aff_mul_itv(pr, tmp, alpha);
				zonotope_aff_check_free(pr, tmp);
				elina_interval_free(alpha);

				zonotope_aff_t *res = zonotope_aff_add(pr, res_box, res_zono, zo);
				zo->paf[offset] = res;
				zo->box_inf[offset] = zo->paf[offset]->itv_inf;
				zo->box_sup[offset] = zo->paf[offset]->itv_sup;
				zonotope_aff_check_free(pr, res_box);
				zonotope_aff_check_free(pr, res_zono);
				
			
			}
			else{	
    			
				elina_interval_t *lambda = elina_interval_alloc();
				elina_interval_set_double(lambda,-lambda_l,lambda_u);
            		
				zonotope_aff_t * res = zonotope_aff_mul_itv(pr, zo->paf[offset], lambda);
				zonotope_aff_check_free(pr, zo->paf[offset]);
				elina_interval_free(lambda);

				double mid_inf = 0.0;
    				double mid_sup = 0.0;
    				double dev_inf = 0.0;
    				double dev_sup = 0.0;
				elina_interval_middev(&mid_inf, &mid_sup, &dev_inf, &dev_sup, 0,bound_u);
				res->c_inf = res->c_inf + mid_inf;
				res->c_sup = res->c_sup + mid_sup;
				
   				if (dev_inf!=0 || dev_sup!=0) {
					zonotope_aaterm_t* ptr = zonotope_aaterm_alloc_init();
					ptr->inf = dev_inf;
					ptr->sup = dev_sup;
					ptr->pnsym = epsilon_map[i];
                     			if (res->end) res->end->n = ptr;
					else res->q = ptr;
					res->end = ptr;
					res->l++;
    				}				           

				//res->itv_inf+= ;
				res->itv_sup+= bound_u;
				zo->paf[offset] = res;
				zo->box_inf[offset] = zo->paf[offset]->itv_inf;
				zo->box_sup[offset] = zo->paf[offset]->itv_sup;
			}
			//zo->paf[offset] = res;
			zo->paf[offset]->pby++;
			
		}
		offset++;
    	}
	return NULL; 
}

elina_abstract0_t *relu_zono_refined(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t x, double new_inf, double new_sup){
	
	//zonoml_fprint(stdout,man,zo,NULL);
	//elina_lincons0_array_t arr = zonoml_to_lincons_array(man,zo);
	//elina_lincons0_array_fprint(stdout,&arr,NULL);
	//elina_lincons0_array_clear(&arr);
	zonotope_t * input = zonotope_of_abstract0(abs);
	zonotope_t *zo = destructive? input : zonotope_copy(man,input);
	if(zonotope_is_bottom(man,zo)){
		
		return abstract0_of_zonotope(man,zo);
	}
	zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
	//elina_abstract0_t * oct = zo->oct;
	elina_dimension_t dims = zonotope_dimension(man,zo);
	elina_dim_t num_dim = dims.intdim + dims.realdim;
	
	
	elina_dimchange_t * add_dim = elina_dimchange_alloc(0,1);
	elina_dim_t y = num_dim;
	add_dim->dim[0] = y;
	zo = zonotope_add_dimensions(man,true,zo,add_dim,false);
	elina_linexpr0_t * assign_expr;
	elina_linterm_t * lterm;
	elina_lincons0_array_t lincons;
	if(new_sup<=0){
		assign_expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0); 
		elina_coeff_set_interval_double(&assign_expr->cst,0,0);
		zo = zonotope_assign_linexpr_array(man,true,zo,&y,&assign_expr,1,NULL);	
	}
	else if(new_inf>=0){
		assign_expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1); 
		elina_coeff_set_interval_double(&assign_expr->cst,0,0);
		lterm = &assign_expr->p.linterm[0];
		lterm->dim = x;
		elina_coeff_set_interval_double(&lterm->coeff,1,1); 
		zo = zonotope_assign_linexpr_array(man,true,zo,&y,&assign_expr,1,NULL);	
		lincons = elina_lincons0_array_make(2);
		lincons.p[0].constyp = ELINA_CONS_SUPEQ;
		elina_linexpr0_t * expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
		elina_coeff_set_scalar_double(&expr->cst,-new_inf);
		elina_linterm_t * lterm = &expr->p.linterm[0];
		lterm->dim = y;
		elina_coeff_set_interval_double(&lterm->coeff,1,1); 
		lincons.p[0].linexpr0 = expr;

                lincons.p[1].constyp = ELINA_CONS_SUPEQ;
		expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
		elina_coeff_set_scalar_double(&expr->cst,new_sup);
		lterm = &expr->p.linterm[0];
		lterm->dim = y;
		elina_coeff_set_interval_double(&lterm->coeff,-1,-1); 
		lincons.p[1].linexpr0 = expr;


		zo = zonotope_meet_lincons_array(man,true,zo,&lincons);
		elina_lincons0_array_clear(&lincons);	
	}
	else{
		
		lincons = elina_lincons0_array_make(2);
		lincons.p[0].constyp = ELINA_CONS_SUPEQ;
		elina_linexpr0_t * expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
		elina_coeff_set_scalar_double(&expr->cst,-new_inf);
		elina_linterm_t * lterm = &expr->p.linterm[0];
		lterm->dim = x;
		elina_coeff_set_scalar_double(&lterm->coeff,1); 
		lincons.p[0].linexpr0 = expr;


		lincons.p[1].constyp = ELINA_CONS_SUPEQ;
		expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
		elina_coeff_set_scalar_double(&expr->cst,new_sup);
		lterm = &expr->p.linterm[0];
		lterm->dim = x;
		elina_coeff_set_scalar_double(&lterm->coeff,-1); 
		lincons.p[1].linexpr0 = expr;


		zo = zonotope_meet_lincons_array(man,true,zo,&lincons);
		elina_lincons0_array_clear(&lincons);
			
		double inf_l = -new_inf;
		double inf_u = new_inf;
		double sup_l = -new_sup;
		double sup_u = new_sup;		
		double width_l = sup_l + inf_u;
		double width_u = sup_u + inf_l;
		
		double bound_l, bound_u;
		elina_double_interval_mul(&bound_l, &bound_u, inf_l, inf_u, sup_l, sup_u);
		double tmp = bound_l;
		bound_l = bound_u;
		bound_u = tmp;
		elina_double_interval_div(&bound_l, &bound_u, bound_l, bound_u, width_l, width_u);
			
		double lambda_l, lambda_u;
		elina_double_interval_div(&lambda_l, &lambda_u, sup_l, sup_u, width_l, width_u);
			
		lincons = get_meet_lincons_array(y,0,0,-bound_l,bound_u);
		zo = zonotope_meet_lincons_array(man,true,zo,&lincons);
		elina_lincons0_array_clear(&lincons);	
		
		assign_expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2); 
		elina_coeff_set_interval_double(&assign_expr->cst,0,0);
		lterm = &assign_expr->p.linterm[0];
		lterm->dim = y;
		elina_coeff_set_interval_double(&lterm->coeff,1,1); 
		lterm = &assign_expr->p.linterm[1];
		lterm->dim = x;
		elina_coeff_set_interval_double(&lterm->coeff,-lambda_l, lambda_u); 
		zo = zonotope_assign_linexpr_array(man,true,zo,&y,&assign_expr,1,NULL);	
		
	}	
	
	elina_linexpr0_free(assign_expr);
	elina_dimperm_t* dimperm = elina_dimperm_alloc(num_dim+1);
	elina_dimperm_set_id(dimperm);
	dimperm->dim[x] = num_dim;
	dimperm->dim[num_dim] = x;
	zo = zonotope_permute_dimensions(man,true,zo,dimperm);
	zo = zonotope_remove_dimensions(man,true,zo,add_dim);
	elina_dimchange_free(add_dim);
	elina_dimperm_free(dimperm);
	return abstract0_of_zonotope(man,zo);
	
}		

bool affine_form_is_box(elina_manager_t* man, elina_abstract0_t *abs, elina_dim_t x){
	zonotope_t * zo = zonotope_of_abstract0(abs);
	return (zo->paf[x]->l<=1);
}
	

elina_abstract0_t *maxpool_zono_refined(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t x, double new_inf, double new_sup){
	zonotope_t * input = zonotope_of_abstract0(abs);
	zonotope_t *zo = destructive? input : zonotope_copy(man,input);
	if(zonotope_is_bottom(man,zo)){
		
		return abstract0_of_zonotope(man,zo);
	}
	
	elina_lincons0_array_t lincons = elina_lincons0_array_make(2);

	lincons.p[0].constyp = ELINA_CONS_SUPEQ;
	elina_linexpr0_t * expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
	elina_coeff_set_scalar_double(&expr->cst,-new_inf);
	elina_linterm_t * lterm = &expr->p.linterm[0];
	lterm->dim = x;
	elina_coeff_set_interval_double(&lterm->coeff,1,1); 
	lincons.p[0].linexpr0 = expr;

        lincons.p[1].constyp = ELINA_CONS_SUPEQ;
	expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
	elina_coeff_set_scalar_double(&expr->cst,new_sup);
	lterm = &expr->p.linterm[0];
	lterm->dim = x;
	elina_coeff_set_interval_double(&lterm->coeff,-1,-1); 
	lincons.p[1].linexpr0 = expr;


	zo = zonotope_meet_lincons_array(man,true,zo,&lincons);
	elina_lincons0_array_clear(&lincons);
	return abstract0_of_zonotope(man,zo);
}

elina_abstract0_t * relu_zono_layerwise(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t start_offset, elina_dim_t num_dim){
	//elina_dim_t i;
	//elina_dim_t end = start_offset + num_dim;
	elina_dimension_t dimension = elina_abstract0_dimension(man,abs);
	
	elina_abstract0_t *res = destructive? abs : elina_abstract0_copy(man,abs);
   
	//for(i=start_offset; i < end; i++){
	//	res= relu_zono(man,true,res,i);
	//}
        zonotope_t *zo = zonotope_of_abstract0(res);
        relu_zono_parallel(man, zo, start_offset, num_dim, handle_relu_zono_parallel);
        res = abstract0_of_zonotope(man,zo);
       
    return res;
}


elina_abstract0_t * relu_zono_layerwise_split(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t start_offset, elina_dim_t num_dim, 
			elina_dim_t* split_index, int *split_path, size_t split_size){
	elina_dim_t i;
	elina_dim_t end = start_offset + num_dim;
	elina_abstract0_t *res = destructive? abs : elina_abstract0_copy(man,abs);
	zonotope_t *zo = zonotope_of_abstract0(res);
        char *map = (char*)calloc(num_dim,sizeof(char));
	for(i =0; i < split_size; i++){
		map[split_index[i]+start_offset] = 1;
	}
	for(i=start_offset; i < end; i++){
		
		if(map[i]){
			size_t j;
			for(j=0; j < split_size; j++){
				if(split_index[j]==i-start_offset){
					break;
				}
			}
			
			int sgn = split_path[j];
			if(sgn==1){
				
				elina_lincons0_array_t lincons = elina_lincons0_array_make(1);
				lincons.p[0].constyp = ELINA_CONS_SUPEQ;
				elina_linexpr0_t * expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
				elina_coeff_set_interval_double(&expr->cst,0,0);
				elina_linterm_t * lterm = &expr->p.linterm[0];
				lterm->dim = i;
				elina_coeff_set_interval_double(&lterm->coeff,1,1); 
				lincons.p[0].linexpr0 = expr;
				/*elina_interval_t * bound_i = zonotope_bound_dimension(man,zo,i);
				double ub = bound_i->sup->val.dbl;
				elina_interval_free(bound_i);
				zo = zonotope_forget_array(man,true,zo,&i,1,false);
				elina_lincons0_array_t lincons = get_meet_lincons_array(i,0,0,ub,ub);*/
				zo = zonotope_meet_lincons_array(man,true,zo,&lincons);
				
				elina_lincons0_array_clear(&lincons);
			}
			else{
				
				elina_lincons0_array_t lincons = elina_lincons0_array_make(1);
				lincons.p[0].constyp = ELINA_CONS_SUPEQ;
				elina_linexpr0_t * expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
				elina_coeff_set_interval_double(&expr->cst,0,0);
				elina_linterm_t * lterm = &expr->p.linterm[0];
				lterm->dim = i;
				elina_coeff_set_interval_double(&lterm->coeff,-1,-1); 
				lincons.p[0].linexpr0 = expr;
				zo = zonotope_meet_lincons_array(man,true,zo,&lincons);
				
				elina_lincons0_array_clear(&lincons);
				elina_linexpr0_t * assign_expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0); 
				elina_coeff_set_interval_double(&assign_expr->cst,0,0);
	
				zo = zonotope_assign_linexpr_array(man,true,zo,&i,&assign_expr,1,NULL);
				
			}
			res->value = zo;
		}
		else{
			res= relu_zono(man,true,res,i);
			zo = res->value;
		}
	}
	free(map);
	
	return res;
}



void * handle_s_curve_zono_parallel(void *args){
	zonoml_s_curve_thread_t * data = (zonoml_s_curve_thread_t *)args;
	elina_manager_t * man = data->man;
	//elina_abstract0_t * abs = data->abs;
	zonotope_t * zo = data->z;
	zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);


	size_t start_offset = data->start_offset;
	size_t idx_start = data->start;
	size_t idx_end = data->end;
	
	zonotope_noise_symbol_t **epsilon_map = data->epsilon_map;
	bool is_sigmoid = data->is_sigmoid;
	char *is_used = data->is_used;
	size_t offset = start_offset + idx_start;
	size_t i;
	for (i=idx_start; i< idx_end; i++) {
        	
		elina_interval_t * bound_x = zonotope_bound_dimension(man,zo,offset);
		double sup = bound_x->sup->val.dbl;
		double inf = bound_x->inf->val.dbl;
		elina_interval_free(bound_x);
		double inf_l = -inf;
		double inf_u = inf;
		double sup_l = -sup;
		double sup_u = sup;
		if(inf ==sup){
		   fesetround(FE_DOWNWARD);
		   double val_inf = is_sigmoid ? exp(inf) : tanh(inf);
		   if(is_sigmoid){
		      val_inf = val_inf/(1+val_inf);
		   }
		   fesetround(FE_UPWARD);
		   double val_sup = is_sigmoid ? exp(inf) : tanh(inf);
		   if(is_sigmoid){
			val_sup = val_sup/(1+val_sup);
		   }
		   zonotope_aff_check_free(pr,zo->paf[offset]);
		   zonotope_aff_t *res = create_affine_form_for_box(pr,-val_inf,val_sup,epsilon_map[i],&is_used[i]);
		   zo->paf[offset] = res;
		   zo->box_inf[offset] = -val_inf;
		   zo->box_sup[offset] = val_sup;
		   zo->paf[offset]->pby++;
	 	}
		else{
			fesetround(FE_DOWNWARD);
			double e_sup_l = is_sigmoid ? -exp(sup) : -tanh(sup);
			double e_inf_l = is_sigmoid ? -exp(inf) : -tanh(inf);
			fesetround(FE_UPWARD); 
			double e_sup_u = is_sigmoid ? exp(sup) : tanh(sup);
			double e_inf_u = is_sigmoid ? exp(inf) : tanh(inf);
			double f_sup_l, f_sup_u;
			double f_inf_l, f_inf_u;
			double den_sup_l, den_sup_u;
			double den_inf_l, den_inf_u;
			double connecting_slope_l, connecting_slope_u;
		
			if(is_sigmoid){
				den_sup_l = -1 + e_sup_l;
				den_sup_u = 1 + e_sup_u;				
				den_inf_l = -1 + e_inf_l;
				den_inf_u = 1 + e_inf_u;		
				elina_double_interval_div(&f_sup_l, &f_sup_u, e_sup_l, e_sup_u, den_sup_l, den_sup_u);
				elina_double_interval_div(&f_inf_l, &f_inf_u, e_inf_l, e_inf_u, den_inf_l, den_inf_u);
			}
			else{
				f_inf_l = e_inf_l;
				f_inf_u = e_inf_u;
				f_sup_l = e_sup_l;
				f_sup_u = e_sup_u;
				den_inf_l = e_inf_l;
				den_inf_u = e_inf_u;
				den_sup_l = e_sup_l;
				den_sup_u = e_sup_u;
			}
		
			double slope_l, slope_u;
			bool boxify = false;
			if(inf>0){
				double sq_den_sup_l, sq_den_sup_u;
				elina_double_interval_mul(&sq_den_sup_l, &sq_den_sup_u, den_sup_l, den_sup_u, den_sup_l, den_sup_u);
				if(is_sigmoid){
					elina_double_interval_div(&slope_l, &slope_u, e_sup_l, e_sup_u, sq_den_sup_l, sq_den_sup_u);
				}
				else{
					slope_l = -1 + sq_den_sup_u;
					slope_u = 1 + sq_den_sup_l; 
				}

			}
			else if(sup<0){
				double sq_den_inf_l, sq_den_inf_u;
				elina_double_interval_mul(&sq_den_inf_l, &sq_den_inf_u, den_inf_l, den_inf_u, den_inf_l, den_inf_u);	
				if(is_sigmoid){
						
					elina_double_interval_div(&slope_l, & slope_u, e_inf_l, e_inf_u, sq_den_inf_l, sq_den_inf_u);
				}
				else{
				
					slope_l = -1 + sq_den_inf_u;
					slope_u = 1 + sq_den_inf_l; 
				}
			}
			else{
				double sq_den_sup_l, sq_den_sup_u;
				double slope1_l, slope1_u;
				elina_double_interval_mul(&sq_den_sup_l, &sq_den_sup_u, den_sup_l, den_sup_u, den_sup_l, den_sup_u);
				if(is_sigmoid){
					elina_double_interval_div(&slope1_l, &slope1_u, e_sup_l, e_sup_u, sq_den_sup_l, sq_den_sup_u);
				}

				else{
					slope1_l = -1 + sq_den_sup_u;
					slope1_u = 1 + sq_den_sup_l; 
				}

				double sq_den_inf_l, sq_den_inf_u;
				double slope2_l, slope2_u;
				elina_double_interval_mul(&sq_den_inf_l, &sq_den_inf_u, den_inf_l, den_inf_u, den_inf_l, den_inf_u);
				if(is_sigmoid){
					elina_double_interval_div(&slope2_l, & slope2_u, e_inf_l, e_inf_u, sq_den_inf_l, sq_den_inf_u);
				}
				else{
					slope2_l = -1 + sq_den_inf_u;
					slope2_u = 1 + sq_den_inf_l; 
				}
				if(slope1_u < -slope2_l){
					slope_l = slope1_l;
					slope_u = slope1_u;
				}
				else if(slope2_u < -slope1_l){
					slope_l = slope2_l;
					slope_u = slope2_u;				
				}
				else{
					boxify = true;
					
				}	
			}
		
			if(boxify){
				// meet with linear constraint
				
				zonotope_aff_check_free(pr,zo->paf[offset]);
				zonotope_aff_t *res = create_affine_form_for_box(pr,f_inf_l,f_sup_u,epsilon_map[i],&is_used[i]);
		   		zo->paf[offset] = res;
		   		zo->box_inf[offset] = f_inf_l;
		   		zo->box_sup[offset] = f_sup_u;
		   		zo->paf[offset]->pby++;
				
			}
			else{
				double tmp_l, tmp_u;
				elina_double_interval_mul(&tmp_l, &tmp_u, sup_l, sup_u, slope_l, slope_u);
				double bound1_l, bound1_u;
				bound1_l = f_sup_l + tmp_u;			
				bound1_u = f_sup_u + tmp_l;			

				elina_double_interval_mul(&tmp_l, &tmp_u, inf_l, inf_u, slope_l, slope_u);
				double bound2_l, bound2_u;
				bound2_l = f_inf_l + tmp_u;
				bound2_u = f_inf_u + tmp_l;
				
				if(-bound1_l > bound2_u){	
					elina_interval_t *lambda = elina_interval_alloc();
					elina_interval_set_double(lambda,-slope_l,slope_u);
            		
					zonotope_aff_t * res = zonotope_aff_mul_itv(pr, zo->paf[offset], lambda);
					zonotope_aff_check_free(pr,zo->paf[offset]);
					elina_interval_free(lambda);

					double mid_inf = 0.0;
    					double mid_sup = 0.0;
    					double dev_inf = 0.0;
    					double dev_sup = 0.0;
					elina_interval_middev(&mid_inf, &mid_sup, &dev_inf, &dev_sup, bound2_l,bound1_u);
					res->c_inf = res->c_inf + mid_inf;
					res->c_sup = res->c_sup + mid_sup;
				
   					if (dev_inf!=0 || dev_sup!=0) {
						zonotope_aaterm_t* ptr = zonotope_aaterm_alloc_init();
						ptr->inf = dev_inf;
						ptr->sup = dev_sup;
						ptr->pnsym = epsilon_map[i];
						is_used[i] = 1;
                     				if (res->end) res->end->n = ptr;
						else res->q = ptr;
						res->end = ptr;
						res->l++;
    					}				           

					res->itv_inf+= bound2_l;
					res->itv_sup+= bound1_u;
					zo->paf[offset] = res;
					zo->box_inf[offset] = zo->paf[offset]->itv_inf;
					zo->box_sup[offset] = zo->paf[offset]->itv_sup;
					
				}	
				else{
					//printf("boxify2\n");
					//fflush(stdout);
					zonotope_aff_check_free(pr,zo->paf[offset]);
					zonotope_aff_t *res = create_affine_form_for_box(pr,f_inf_l,f_sup_u,epsilon_map[i],&is_used[i]);
		   			zo->paf[offset] = res;
		   			zo->box_inf[offset] = f_inf_l;
		   			zo->box_sup[offset] = f_sup_u;
		   			zo->paf[offset]->pby++;
				}
			}
		}
		offset++;
    	}
	return NULL; 
}


elina_abstract0_t * s_curve_zono(elina_manager_t *man, bool destructive,
				 elina_abstract0_t *abs, elina_dim_t x, bool is_sigmoid){

	zonotope_t * input = zonotope_of_abstract0(abs);
	zonotope_t *zo = destructive? input : zonotope_copy(man,input);
	if(zonotope_is_bottom(man,zo)){
		return abstract0_of_zonotope(man,zo);
	}

	zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
	
	elina_dimension_t dims = zonotope_dimension(man,zo);
	elina_dim_t num_dim = dims.intdim + dims.realdim;

	elina_dimchange_t *dimchange = elina_dimchange_alloc(0,1);
	dimchange->dim[0] = num_dim;
	zo = zonotope_add_dimensions(man,true,zo,dimchange,false);
   
	elina_interval_t * bound_x = zonotope_bound_dimension(man,zo,x);
	double sup = bound_x->sup->val.dbl;
	double inf = bound_x->inf->val.dbl;
	elina_interval_free(bound_x);
	double inf_l = -inf;
	double inf_u = inf;
	double sup_l = -sup;
	double sup_u = sup;
	if(inf ==sup){
		
		fesetround(FE_DOWNWARD);
		double val_inf = is_sigmoid ? exp(inf) : tanh(inf);
		if(is_sigmoid){
			val_inf = val_inf/(1+val_inf);
		}
		fesetround(FE_UPWARD);
		double val_sup = is_sigmoid ? exp(inf) : tanh(inf);
		if(is_sigmoid){
			val_sup = val_sup/(1+val_sup);
		}
		//printf("inf: %.300f %.300f %.300f\n",inf,val_inf,val_sup);
		elina_linexpr0_t * assign_expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0); 
		elina_coeff_set_interval_double(&assign_expr->cst,val_inf, val_sup);
		zo = zonotope_assign_linexpr_array(man,true,zo,&num_dim,&assign_expr,1,NULL);
		//return abstract0_of_zonotope(man,zo);	
	}
	
	else{
		
		fesetround(FE_DOWNWARD);
		double e_sup_l = is_sigmoid ? -exp(sup) : -tanh(sup);
		double e_inf_l = is_sigmoid ? -exp(inf) : -tanh(inf);
		fesetround(FE_UPWARD); 
		double e_sup_u = is_sigmoid ? exp(sup) : tanh(sup);
		double e_inf_u = is_sigmoid ? exp(inf) : tanh(inf);
		double f_sup_l, f_sup_u;
		double f_inf_l, f_inf_u;
		double den_sup_l, den_sup_u;
		double den_inf_l, den_inf_u;
		double connecting_slope_l, connecting_slope_u;
		
		if(is_sigmoid){
			den_sup_l = -1 + e_sup_l;
			den_sup_u = 1 + e_sup_u;
			
			//double e_inf_l = -e_inf;
			//double e_inf_u = e_inf;
							
			den_inf_l = -1 + e_inf_l;
			den_inf_u = 1 + e_inf_u;
				
						
			elina_double_interval_div(&f_sup_l, &f_sup_u, e_sup_l, e_sup_u, den_sup_l, den_sup_u);

			
			elina_double_interval_div(&f_inf_l, &f_inf_u, e_inf_l, e_inf_u, den_inf_l, den_inf_u);
		}
		else{
			f_inf_l = e_inf_l;
			f_inf_u = e_inf_u;
			f_sup_l = e_sup_l;
			f_sup_u = e_sup_u;
			den_inf_l = e_inf_l;
			den_inf_u = e_inf_u;
			den_sup_l = e_sup_l;
			den_sup_u = e_sup_u;
		}
		elina_lincons0_array_t lincons;
		double slope_l, slope_u;
		bool boxify = false;
		if(inf>0){
			double sq_den_sup_l, sq_den_sup_u;
			elina_double_interval_mul(&sq_den_sup_l, &sq_den_sup_u, den_sup_l, den_sup_u, den_sup_l, den_sup_u);
			if(is_sigmoid){
				elina_double_interval_div(&slope_l, &slope_u, e_sup_l, e_sup_u, sq_den_sup_l, sq_den_sup_u);
			}
			else{
				slope_l = -1 + sq_den_sup_u;
				slope_u = 1 + sq_den_sup_l; 
			}

		}
		else if(sup<0){
			double sq_den_inf_l, sq_den_inf_u;
			elina_double_interval_mul(&sq_den_inf_l, &sq_den_inf_u, den_inf_l, den_inf_u, den_inf_l, den_inf_u);	
			if(is_sigmoid){
						
				elina_double_interval_div(&slope_l, & slope_u, e_inf_l, e_inf_u, sq_den_inf_l, sq_den_inf_u);
			}
			else{
				
				slope_l = -1 + sq_den_inf_u;
				slope_u = 1 + sq_den_inf_l; 
			}
		}
		else{
			double sq_den_sup_l, sq_den_sup_u;
			double slope1_l, slope1_u;
			elina_double_interval_mul(&sq_den_sup_l, &sq_den_sup_u, den_sup_l, den_sup_u, den_sup_l, den_sup_u);
			if(is_sigmoid){
				
			
				elina_double_interval_div(&slope1_l, &slope1_u, e_sup_l, e_sup_u, sq_den_sup_l, sq_den_sup_u);
			}

			else{
				
				slope1_l = -1 + sq_den_sup_u;
				slope1_u = 1 + sq_den_sup_l; 
			}

			double sq_den_inf_l, sq_den_inf_u;
			double slope2_l, slope2_u;
			elina_double_interval_mul(&sq_den_inf_l, &sq_den_inf_u, den_inf_l, den_inf_u, den_inf_l, den_inf_u);
			if(is_sigmoid){
				
			
				elina_double_interval_div(&slope2_l, & slope2_u, e_inf_l, e_inf_u, sq_den_inf_l, sq_den_inf_u);
			}
			else{
				
				slope2_l = -1 + sq_den_inf_u;
				slope2_u = 1 + sq_den_inf_l; 
			}
			if(slope1_u < -slope2_l){
				slope_l = slope1_l;
				slope_u = slope1_u;
			}
			else if(slope2_u < -slope1_l){
				slope_l = slope2_l;
				slope_u = slope2_u;				
			}
			else{
				boxify = true;
			}	
		}
		
		if(boxify){
			// meet with linear constraint
			
			lincons = get_meet_lincons_array(num_dim,-f_inf_l,f_inf_u,-f_sup_l,f_sup_u); 
			zo = zonotope_meet_lincons_array(man,true,zo,&lincons);
			elina_lincons0_array_clear(&lincons);
		}
		else{
			double tmp_l, tmp_u;
			elina_double_interval_mul(&tmp_l, &tmp_u, sup_l, sup_u, slope_l, slope_u);
			double bound1_l, bound1_u;
			bound1_l = f_sup_l + tmp_u;			
			bound1_u = f_sup_u + tmp_l;			

			elina_double_interval_mul(&tmp_l, &tmp_u, inf_l, inf_u, slope_l, slope_u);
			double bound2_l, bound2_u;
			bound2_l = f_inf_l + tmp_u;
			bound2_u = f_inf_u + tmp_l;
			
			if(-bound1_l > bound2_u){	
			
				// meet with linear constraint
				lincons = get_meet_lincons_array(num_dim,-bound2_l,bound2_u,-bound1_l,bound1_u); 
				//printf("bounds: %g %g %g %g\n", bound2_l, bound2_u, bound1_l, bound1_u);
				//elina_lincons0_array_fprint(stdout,&lincons,NULL);
				zo = zonotope_meet_lincons_array(man,true,zo,&lincons);
				elina_lincons0_array_clear(&lincons);
				
				
				// assignment
				elina_linterm_t * lterm;
				elina_linexpr0_t * assign_expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2); 
				elina_coeff_set_interval_double(&assign_expr->cst,0,0);
				lterm = &assign_expr->p.linterm[0];
				lterm->dim = num_dim;
				elina_coeff_set_interval_double(&lterm->coeff,1,1); 
				lterm = &assign_expr->p.linterm[1];
				lterm->dim = x;
					
				elina_coeff_set_interval_double(&lterm->coeff,-slope_l, slope_u); 
						
				zo = zonotope_assign_linexpr_array(man,true,zo,&num_dim,&assign_expr,1,NULL);
				elina_linexpr0_free(assign_expr);
			}	
			else{
				
				//printf("WTF %.30f %.30f %.30f %.30f %.30f %.30f\n",-f_inf_l,f_inf_u,-f_sup_l,f_sup_u,inf,sup);
				lincons = get_meet_lincons_array(num_dim,-f_inf_l,f_inf_u,-f_sup_l,f_sup_u); 
				zo = zonotope_meet_lincons_array(man,true,zo,&lincons);
				elina_lincons0_array_clear(&lincons);
			}
		}
	}
	elina_dimperm_t* dimperm = elina_dimperm_alloc(num_dim+1);
	elina_dimperm_set_id(dimperm);
	dimperm->dim[x] = num_dim;
	dimperm->dim[num_dim] = x;
	zo = zonotope_permute_dimensions(man,true,zo,dimperm);
	zo = zonotope_remove_dimensions(man,true,zo,dimchange);
	elina_dimchange_free(dimchange);
	elina_dimperm_free(dimperm);
	
	return abstract0_of_zonotope(man,zo);

}

elina_abstract0_t * sigmoid_zono(elina_manager_t *man, bool destructive,
				 elina_abstract0_t *abs, elina_dim_t x){
		return s_curve_zono(man,destructive,abs,x,true);
}


elina_abstract0_t * sigmoid_zono_layerwise(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t start_offset, elina_dim_t num_dim){
	//elina_dim_t i;
	//elina_dim_t end = start_offset + num_dim;
	elina_abstract0_t *res = destructive? abs : elina_abstract0_copy(man,abs);
	
	//for(i=start_offset; i < end; i++){
	//	res= sigmoid_zono(man,true,res,i);
	//}
	zonotope_t *zo = zonotope_of_abstract0(res);
        s_curve_zono_parallel(man, zo, start_offset, num_dim, handle_s_curve_zono_parallel,true);
        res = abstract0_of_zonotope(man,zo);
	
	return res;
}

elina_abstract0_t * tanh_zono(elina_manager_t *man, bool destructive,
				 elina_abstract0_t *abs, elina_dim_t x){
		return s_curve_zono(man,destructive,abs,x,false);
}


elina_abstract0_t * tanh_zono_layerwise(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t start_offset, elina_dim_t num_dim){
	//elina_dim_t i;
	//elina_dim_t end = start_offset + num_dim;
	elina_abstract0_t *res = destructive? abs : elina_abstract0_copy(man,abs);
	//for(i=start_offset; i < end; i++){
	//	res= tanh_zono(man,true,res,i);
	//}
	
	zonotope_t *zo = zonotope_of_abstract0(res);
        s_curve_zono_parallel(man, zo, start_offset, num_dim, handle_s_curve_zono_parallel,false);
        res = abstract0_of_zonotope(man,zo);
	
	return res;
}

bool is_greater(elina_manager_t *man, elina_abstract0_t *elem, elina_dim_t y, elina_dim_t x){
    
	zonotope_t * zo = zonotope_of_abstract0(elem);
    if(zo->box_inf[y]>zo->box_sup[x]){
        return true;
    }
	elina_dimension_t dims = zonotope_dimension(man,zo);
	elina_dim_t num_dim = dims.intdim + dims.realdim;

	elina_dimchange_t *dimchange = elina_dimchange_alloc(0,1);
	dimchange->dim[0] = num_dim;
	zo = zonotope_add_dimensions(man,true,zo,dimchange,false);
	elina_linexpr0_t * assign_expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2); 
	elina_coeff_set_interval_double(&assign_expr->cst,0,0);
	elina_linterm_t *lterm = &assign_expr->p.linterm[0];
	lterm->dim = y;
	elina_coeff_set_interval_double(&lterm->coeff,1,1); 
	lterm = &assign_expr->p.linterm[1];
	lterm->dim = x;
	elina_coeff_set_interval_double(&lterm->coeff,-1, -1); 
	//elina_dim_t extra_var = num_dim +1;		
	zo = zonotope_assign_linexpr_array(man,true,zo,&num_dim,&assign_expr,1,NULL);
	elina_interval_t * bound_x = zonotope_bound_dimension(man,zo,num_dim);
        
	bool res = elina_scalar_sgn(bound_x->inf)>0;
	elina_interval_free(bound_x);
	zo = zonotope_remove_dimensions(man,true,zo,dimchange);
	elina_dimchange_free(dimchange);
	elina_linexpr0_free(assign_expr);
	return res;
}


void * handle_maxpool_zono_parallel(void *args){
	zonoml_maxpool_thread_t * data = (zonoml_maxpool_thread_t *)args;
	zonotope_internal_t * pr = data->pr;
	zonotope_t * z = data->z;
	size_t src_offset = data->src_offset;
	size_t idx_start = data->start;
	size_t idx_end = data->end;
	size_t *pool_size = data->pool_size;

	size_t dst_offset = data->dst_offset;
	size_t *input_size = data->input_size;
	
	
	size_t *strides = data->strides;
	size_t *output_size = data->output_size;		
	long int pad_top = data->pad_top;
	long int pad_left = data->pad_left;  
	
	zonotope_noise_symbol_t ** epsilon_map = data->epsilon_map;
	char * is_used = data->is_used;	

	size_t o12 = output_size[1]*output_size[2];
   	size_t i12 = input_size[1]*input_size[2];
    	size_t p01 = pool_size[0]*pool_size[1];

	
	size_t out_pos, mat_x;
	double * inf = (double *) calloc(p01,sizeof(double));
	double * sup = (double *) calloc(p01,sizeof(double));
	size_t * pool_map = (size_t *)calloc(p01,sizeof(double));
	
	
        size_t out_x, out_y, out_z;
	size_t i,j,k;
        size_t m = input_size[0]*input_size[1]*input_size[2];
	for(mat_x=idx_start; mat_x< idx_end; mat_x++){ 
	    out_x = mat_x / o12;
	    out_y = (mat_x-out_x*o12) / output_size[2];
	    out_z = mat_x-out_x*o12 - out_y*output_size[2];
		//size_t inp_x = out_x*pool_size[0];
		//size_t inp_y = out_y*pool_size[1];
	    size_t inp_z = out_z;
		//long int inp_pos = (inp_x*strides[0]-pad_top)*i12 + (inp_y*strides[1]-pad_left)*input_size[2] + inp_z;
		//size_t pool_start_dim = mat_x*pool_size[0]*pool_size[1];
		//printf("inpXYZ: %zu, %zu, %zu\n", inp_x, inp_y, inp_z);
        	//printf("outXYZ: %zu, %zu, %zu\n", out_x, out_y, out_z);
	    size_t x_shift, y_shift, l = 0;
	    double sum_u = 0.0;
	    double sum_l = 0.0;
	    double max_u = -INFINITY;
	    double max_l = -INFINITY;
	    for(x_shift = 0; x_shift < pool_size[0]; x_shift++){
		for(y_shift = 0; y_shift < pool_size[1]; y_shift++){
		    long int x_val = out_x*strides[0] + x_shift - pad_top;
		   
		    if(x_val<0 || x_val>=(long int)input_size[0]){
			continue;			
		    }
		    long int y_val = out_y*strides[1] + y_shift - pad_left;
			
		    if(y_val<0 || y_val >= (long int)input_size[1]){
			continue;
		    }				
		    size_t mat_offset = x_val*i12 + y_val*input_size[2] + inp_z;
		    if(mat_offset>=m){
			continue;
		    }				
		    size_t  pool_cur_dim = src_offset + mat_offset;
				
		    pool_map[l] = pool_cur_dim;
		    elina_interval_t * bound = zonotope_bound_dimension(pr->man,z,pool_cur_dim);
		    inf[l] = bound->inf->val.dbl;
		    sup[l] = bound->sup->val.dbl;
		    elina_interval_free(bound);
		    sum_u = sum_u + sup[l];
		    sum_l = sum_l - inf[l];
		    if(sup[l]>max_u){
			max_u = sup[l];
		    }
		    if(inf[l] > max_l){
			max_l = inf[l];
		    }
				
		    l++;
				
		}
	    }
		
	    bool flag = false;
	    elina_dim_t var = 0;
	    for(j=0; j < l; j++){
		bool g_flag = true;
		for(k = 0;  k < l; k++){
		    if(k==j)continue;
		    if((inf[k]==sup[k]) && (inf[j]>=sup[k])){
			continue;
		    }
		    else if((inf[j]==inf[k]) && (sup[j]==sup[k]) && (inf[j]==sup[j])){
			continue;
		    }
		    
		    zonotope_aff_t * x = z->paf[pool_map[j]];
		    zonotope_aff_t * y = z->paf[pool_map[k]];
		    elina_scalar_t *scalar = elina_scalar_alloc();
		    elina_scalar_set_double(scalar,-1);		
		    zonotope_aff_t *tmp = zonotope_aff_mul_scalar(pr, y, scalar);
		    zonotope_aff_t * res= zonotope_aff_add(pr, x, tmp, z);	
		    
		    g_flag = (res->itv_inf<0);
		    elina_scalar_free(scalar);
		    zonotope_aff_check_free(pr,tmp);
		    zonotope_aff_check_free(pr,res);	 
		    if(!g_flag){
			break;
		    }
				
		}
		if(g_flag){
			flag = true;
			var =(elina_dim_t)pool_map[j];
			break;
		}
	  }
	  elina_coeff_t *coeff, *cst;
	  elina_linexpr0_t * linexpr0;
	  elina_linterm_t *linterm;
	  out_pos = mat_x + dst_offset;
          zonotope_aff_check_free(pr,z->paf[out_pos]);
	 
	  if(flag){
	     //x_new = x_var
	     elina_scalar_t *scalar = elina_scalar_alloc();
	     elina_scalar_set_double(scalar,1);		
	     z->paf[out_pos] = zonotope_aff_mul_scalar(pr, z->paf[var], scalar);
		
	     z->box_inf[out_pos] = z->paf[out_pos]->itv_inf; 
	     z->box_sup[out_pos] = z->paf[out_pos]->itv_sup; 
	    
	     elina_scalar_free(scalar);
	     z->paf[out_pos]->pby++;
	  }
	  else{
	     //max_l<= x_new <= max_u
		
	     size_t noise_index = out_pos - dst_offset;
	     
	     is_used[noise_index] = 1;
	     double mid_inf = 0.0;
    	     double mid_sup = 0.0;
    	     double dev_inf = 0.0;
    	     double dev_sup = 0.0;
	     zonotope_aff_t * res_box = zonotope_aff_alloc_init(pr);
	     elina_interval_middev(&mid_inf, &mid_sup, &dev_inf, &dev_sup, -max_l,max_u);
	     res_box->c_inf = res_box->c_inf + mid_inf;
	     res_box->c_sup = res_box->c_sup + mid_sup;
	     res_box->itv_inf = -max_l;
	     res_box->itv_sup = max_u;
				
   	     if (dev_inf!=0 || dev_sup!=0) {
		 zonotope_aaterm_t* ptr = zonotope_aaterm_alloc_init();
		 ptr->inf = dev_inf;
		 ptr->sup = dev_sup;
		 ptr->pnsym = epsilon_map[noise_index];
                 if (res_box->end) res_box->end->n = ptr;
		     else res_box->q = ptr;
		     res_box->end = ptr;
		     res_box->l++;
    		 }
	     z->paf[out_pos] = res_box;
             z->box_inf[out_pos] = -max_l;
	     z->box_sup[out_pos] = max_u;
	     z->paf[out_pos]->pby++;		
	  }
		
	}
	
	free(inf);
	free(sup);
	free(pool_map);
	return NULL; 
}


elina_abstract0_t* maxpool_zono(elina_manager_t *man, bool destructive, elina_abstract0_t *abs, 
			   size_t *pool_size, size_t *input_size, size_t src_offset, size_t* strides, 
			   size_t dimensionality, size_t dst_offset, bool is_valid_padding){
	assert(dimensionality==3);
	assert(pool_size[2]==1);
	//assert(stride[0]==2 && stride[1]==2 && stride[2]==1);
	
	zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);

	size_t * output_size = (size_t *)malloc(dimensionality*sizeof(size_t));
	//for(i=0; i < dimensionality; i++){
	//	output_size[i] = input_size[i]/pool_size[i];
	//}
	
    if(is_valid_padding){
        output_size[0] = ceil((double)(input_size[0] - pool_size[0]+1) / (double)strides[0]);
        output_size[1] = ceil((double)(input_size[1] - pool_size[1]+1) / (double)strides[1]);
    }
    else{
        output_size[0] = ceil((double)input_size[0] / (double)strides[0]);
        output_size[1] = ceil((double)input_size[1] / (double)strides[1]);
    }
    output_size[2] = input_size[2]/pool_size[2];

	
	size_t num_out_neurons = output_size[0]*output_size[1]*output_size[2];

    	
	long int pad_along_height=0, pad_along_width=0;
	long int pad_top=0,  pad_left=0,  tmp=0;
	if(!is_valid_padding){
		if (input_size[0] % strides[0] == 0){
			long int tmp = pool_size[0] - strides[0];
	  		pad_along_height = max(tmp, 0);
		}
		else{
			tmp = pool_size[0] - (input_size[0] % strides[0]);
	  		pad_along_height = max(tmp, 0);
		}
		if (input_size[1] % strides[1] == 0){
			tmp = pool_size[1] - strides[1];
	  		pad_along_width = max(tmp, 0);
		}
		else{
			tmp = pool_size[1] - (input_size[1] % strides[1]);
	  		pad_along_width = max(tmp, 0);
		}
		pad_top = pad_along_height / 2;
		pad_left = pad_along_width / 2;
		
	}
        //printf("pad top: %ld %ld\n",output_size[0],output_size[1]);
	//fflush(stdout);
	zonotope_t * input = zonotope_of_abstract0(abs);

	zonotope_t *res = destructive? input : zonotope_copy(man,input);
        elina_dimension_t dims = zonotope_dimension(pr->man,res);
        elina_dim_t num_var = dims.intdim + dims.realdim;
 	//printf("start %u\n",num_var);
	
	//fflush(stdout);
        maxpool_zono_parallel(pr, res, src_offset, pool_size, num_out_neurons, dst_offset, input_size, strides, output_size, pad_top, pad_left, handle_maxpool_zono_parallel);
	//dims = zonotope_dimension(pr->man,res);
	//num_var = dims.intdim + dims.realdim;
	//printf("end %u\n",num_var);
	//fflush(stdout);
	free(output_size);
	return abstract0_of_zonotope(man,res);
}





void zono_add(elina_manager_t *man, elina_abstract0_t *elem, size_t dst_offset, size_t src_offset, size_t num_var){
	zonotope_t * zo = zonotope_of_abstract0(elem);
	elina_linexpr0_t ** expr_array = (elina_linexpr0_t **)malloc(num_var*sizeof(elina_linexpr0_t*));
	elina_dim_t *tdim = (elina_dim_t *)malloc(num_var*sizeof(elina_dim_t));
	size_t i;
	for(i=0; i < num_var; i++){
		elina_dim_t dst_ind = dst_offset + i;
		elina_dim_t src_ind = src_offset + i;
		tdim[i] = dst_ind;
		expr_array[i] = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2); 
		elina_coeff_set_interval_double(&expr_array[i]->cst,0,0);
		elina_linterm_t *lterm = &expr_array[i]->p.linterm[0];
		lterm->dim = dst_ind;
		elina_coeff_set_interval_double(&lterm->coeff,1,1); 
		lterm = &expr_array[i]->p.linterm[1];
		lterm->dim = src_ind;
		elina_coeff_set_interval_double(&lterm->coeff,1, 1);
	}
	zo = zonotope_assign_linexpr_array(man,true,zo,tdim,expr_array,num_var,NULL);
	for(i=0; i < num_var; i++){
		elina_linexpr0_free(expr_array[i]);
	}
	free(tdim);
	free(expr_array);
	elem->value = zo;
}

void zono_copy_section(elina_manager_t *man, elina_abstract0_t *elem, size_t dst_offset, size_t src_offset, size_t num_var){
	zonotope_t * zo = zonotope_of_abstract0(elem);
	elina_linexpr0_t ** expr_array = (elina_linexpr0_t **)malloc(num_var*sizeof(elina_linexpr0_t*));
	elina_dim_t *tdim = (elina_dim_t *)malloc(num_var*sizeof(elina_dim_t));
	size_t i;
	for(i=0; i < num_var; i++){
		elina_dim_t dst_ind = dst_offset + i;
		elina_dim_t src_ind = src_offset + i;
		tdim[i] = dst_ind;
		expr_array[i] = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1); 
		elina_coeff_set_interval_double(&expr_array[i]->cst,0,0);
		elina_linterm_t *lterm = &expr_array[i]->p.linterm[0];
		lterm->dim = src_ind;
		elina_coeff_set_interval_double(&lterm->coeff,1,1); 
	}
	zo = zonotope_assign_linexpr_array(man,true,zo,tdim,expr_array,num_var,NULL);
	for(i=0; i < num_var; i++){
		elina_linexpr0_free(expr_array[i]);
	}
	free(tdim);
	free(expr_array);
	elem->value=zo;
}

double get_interval_width_var_zono(elina_manager_t *man, elina_abstract0_t *elem, size_t i){
	zonotope_t * zo = zonotope_of_abstract0(elem);
	double inf = zo->box_inf[i];
	double sup = zo->box_sup[i];
	return inf + sup;
}

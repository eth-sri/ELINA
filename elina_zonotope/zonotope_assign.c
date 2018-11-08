#include "zonotope_assign.h"

zonotope_t* zonotope_assign_linexpr_array(elina_manager_t* man,
		bool destructive,
		zonotope_t* z,
		elina_dim_t* tdim, elina_linexpr0_t** lexpr,
		size_t size,
		zonotope_t* dest)
{
    start_timing();
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
    zonotope_t* res = zonotope_copy(man, z);
	//printf("assign input %d\n",z->size);
	//zonotope_fprint(stdout,man,z,NULL);
	//fflush(stdout);
    size_t i = 0;
    for (i=0; i<res->dims; i++) {
        //fflush(stdout);
	res->paf[i]->itv_inf = res->box_inf[i];
	res->paf[i]->itv_sup = res->box_sup[i];
    }
    for (i=0; i<size; i++) {
	zonotope_aff_check_free(pr, res->paf[tdim[i]]);
	//printf("statement x%d:= \n",tdim[i]);
	//elina_linexpr0_fprint(stdout,lexpr[i],NULL);
       //printf("\n %d\n",lexpr[i]->size);
	//fflush(stdout);
        //
	res->paf[tdim[i]] = zonotope_aff_from_linexpr0(pr, lexpr[i], z);
       // printf("affine expression\n");
        //zonotope_aff_fprint(pr,stdout,res->paf[tdim[i]]);
        //
	//printf("center here1 %g %g\n",res->paf[tdim[i]]->c_inf,res->paf[tdim[i]]->c_sup);
	//zonotope_aff_reduce(pr, res->paf[tdim[i]]);
        //printf("center here2 %g %g\n",res->paf[tdim[i]]->c_inf,res->paf[tdim[i]]->c_sup);
	if (zonotope_aff_is_top(pr, res->paf[tdim[i]])) {
	    zonotope_aff_check_free(pr, res->paf[tdim[i]]);
	    res->paf[tdim[i]] = pr->top;
	} else if (zonotope_aff_is_bottom(pr, res->paf[tdim[i]])) {
	    zonotope_aff_check_free(pr, res->paf[tdim[i]]);
	    res->paf[tdim[i]] = pr->bot;
	}
	res->box_inf[tdim[i]] = res->paf[tdim[i]]->itv_inf;
	res->box_sup[tdim[i]] = res->paf[tdim[i]]->itv_sup;
	res->paf[tdim[i]]->pby++;
	//printf("assign output\n");
	//zonotope_fprint(stdout,man,res,NULL);
	//fflush(stdout);
    }
    man->result.flag_best = false;
    man->result.flag_exact = false;
	//printf("assign output %d\n",res->size);
	//zonotope_fprint(stdout,man,res,NULL);
	//fflush(stdout);
    record_timing(zonotope_assign_linexpr_time);
    return res;
}


void zonotope_relu(elina_manager_t *man, elina_abstract0_t * abs, elina_dim_t y, elina_dim_t x){
	zonotope_t * zo = (zonotope_t *)abs->value;
	zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
	//printf("relu input %d\n",zo->dims);
	//zonotope_aff_fprint(pr,stdout,zo->paf[x]);
	//fflush(stdout);
	
	elina_interval_t * bound_x = zonotope_bound_dimension(man,zo,x);
	double lb = bound_x->inf->val.dbl;
	double ub = bound_x->sup->val.dbl;
	double sup = bound_x->sup->val.dbl;
	double inf = bound_x->inf->val.dbl;
	elina_interval_free(bound_x);
	double width = sup -inf;
	elina_lincons0_array_t lincons;
	double bound = (-inf*sup)/width;
	//printf("bound: %g slope: %g\n",bound,sup/width);
	//fflush(stdout);
	zonotope_aff_check_free(pr, zo->paf[y]);
	zo->paf[y] = zonotope_aff_alloc_init(pr);
	if(1){
		/*double lambda = 1/bound;
		if(lamba>=(sup/width)){
			double new_bound = lambda*-inf;
			zo->paf[y]->c_inf = 0;
        		zo->paf[y]->c_sup = new_bound;
			elina_linterm_t * lterm;
			elina_linexpr0_t * assign_expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2); 
			elina_coeff_set_scalar_double(&assign_expr->cst,0);
		lterm = &assign_expr->p.linterm[0];
		lterm->dim = z;
		elina_coeff_set_scalar_double(&lterm->coeff,1); 
		lterm = &assign_expr->p.linterm[1];
		lterm->dim = x;
		elina_coeff_set_scalar_double(&lterm->coeff,sup/width); 
	
		zo->zono = elina_abstract0_assign_linexpr_array(pr->man_zono,true,zo->zono,&z,&assign_expr,1,NULL);
		
		lterm = &assign_expr->p.linterm[0];
		lterm->dim = num_dim;
		elina_coeff_set_scalar_double(&lterm->coeff,1-lambda); 
		lterm = &assign_expr->p.linterm[1];
		lterm->dim = z;
		elina_coeff_set_scalar_double(&lterm->coeff,lambda); 
	
		zo->zono = elina_abstract0_assign_linexpr_array(pr->man_zono,true,zo->zono,&y,&assign_expr,1,NULL);
				
		elina_linexpr0_free(assign_expr);
		
		}
		else{
			zo->paf[y]->c_inf = 0;
        		zo->paf[y]->c_sup = sup;
		}*/
		zo->paf[y]->c_inf = 0;
		zo->paf[y]->c_sup = sup;
		zo->paf[y]->itv_inf = 0;
		zo->paf[y]->itv_sup = sup;
		zo->paf[y]->pby++;
		zo->box_inf[y] = 0;
		zo->box_sup[y] = sup;
		return;
	}
	
	//zo->paf[y]->c_inf = 0;
        //zo->paf[y]->c_sup = bound;
	//zo->paf[y]->itv_inf = 0;
    	//zo->paf[y]->itv_sup = bound;
	//zo->paf[y]->pby++;	
	//zo->box_inf[y] = 0;
	//zo->box_sup[y] = bound;
	//elina_linterm_t * lterm;
	//elina_linexpr0_t * assign_expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2); 
	//elina_coeff_set_scalar_double(&assign_expr->cst,0);
	//lterm = &assign_expr->p.linterm[0];
	//lterm->dim = y;
	//elina_coeff_set_scalar_double(&lterm->coeff,1); 
	//lterm = &assign_expr->p.linterm[1];
	//lterm->dim = x;
	//elina_coeff_set_scalar_double(&lterm->coeff,sup/width); 
	
	zonotope_aff_t *tmp;
	elina_interval_t *interval = elina_interval_alloc();
	elina_interval_set_double(interval,sup/width,sup/width);
	tmp = zonotope_aff_mul_itv(pr,zo->paf[x],interval);
        
        zonotope_aff_t *tmp1 = zonotope_aff_alloc_init(pr);
        tmp1->c_inf = 0;
        tmp1->c_sup = bound;
	tmp1->itv_inf = 0;
    	tmp1->itv_sup = bound;
	tmp1->pby++;
        zo->paf[y] = zonotope_aff_add(pr,tmp1,tmp,zo);
	zonotope_aff_free(pr,tmp);
	elina_interval_free(interval);
	//zo->paf[y] = zonotope_aff_from_linexpr0(pr, &assign_expr, zo);
	if (zonotope_aff_is_top(pr, zo->paf[y])) {
	    zonotope_aff_check_free(pr, zo->paf[y]);
	    zo->paf[y] = pr->top;
	} else if (zonotope_aff_is_bottom(pr, zo->paf[y])) {
	    zonotope_aff_check_free(pr, zo->paf[y]);
	    zo->paf[y] = pr->bot;
	}
	zo->box_inf[y] = zo->paf[y]->itv_inf;
	zo->box_sup[y] = zo->paf[y]->itv_sup;
	zo->paf[y]->pby++;

	//elina_linexpr0_free(assign_expr);
	//printf("relu output %d\n",zo->dims);
	//zonotope_aff_fprint(pr,stdout,zo->paf[y]);
	//printf("center: %g %g",-zo->paf[y]->c_inf,zo->paf[y]->c_sup);
	//fflush(stdout);
	return;
}

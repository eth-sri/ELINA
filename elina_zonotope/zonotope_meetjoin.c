/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
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

#include "zonotope.h"

/**********************/
/* 1. Meet Lincons    */
/**********************/
zonotope_t* zonotope_meet_lincons_array(elina_manager_t* man, bool destructive, zonotope_t* z, elina_lincons0_array_t* array)
{
    start_timing();
     zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_MEET_LINCONS_ARRAY);
 //printf("meet start lincons %d\n",z->hypercube);
   // zonotope_fprint(stdout,man,z,NULL);
   // elina_lincons0_array_fprint(stdout,array,NULL);
     //   fflush(stdout);
   
    //arg_assert(a && array, abort(););
  
    size_t i = 0;
    zonotope_t* res = destructive ? z : zonotope_copy(man, z);
    bool is_bottom = false;
    double box_inf = 0.0;
    double box_sup = 0.0;
    bool* tchange = (bool*)calloc(2*z->dims, sizeof(bool));
    
    size_t kmax = 2; /* specifies the maximum number of iterations */
    /* intervalonly is set to false which means try to improve all dimensions, not only the ones with an interval coefficient */
    //elina_interval_t **env = elina_interval_array_alloc(res->dims);
    //for(i=0; i < res->dims; i++){
	//env[i] = elina_interval_alloc();
//	elina_interval_set_double(env[i],res->box_inf[i],res->box_sup[i]);
	//printf("i: %d %g %g\n",i,-res->box_inf[i],res->box_sup[i]);
    //}
    if (elina_double_boxize_lincons0_array(res->box_inf, res->box_sup, tchange, array, res->box_inf, res->box_sup,res->intdim, kmax, false,ELINA_SCALAR_DOUBLE)) {
	 //if (elina_boxize_lincons0_array(env, tchange, array, env, res->intdim, kmax, false,ELINA_SCALAR_DOUBLE)) {
	    /* there is some inferred bounds */
	    //for(i=0; i < res->dims; i++){
		//res->box_inf[i] = env[i]->inf->val.dbl;
		//res->box_sup[i] = env[i]->sup->val.dbl;
    	    //}
		
	    for (i=0; i<res->dims; i++) {
		
		if (tchange[2*i] || tchange[2*i + 1]) {
		    if (-res->box_inf[i] > res->box_sup[i]) {
			
			is_bottom = true;
			break;
		    } else if (res->paf[i]->q == NULL) {
			
			zonotope_aff_check_free(pr, res->paf[i]);
			res->paf[i] = zonotope_aff_alloc_init(pr);
			
			if((res->box_sup[i]!=INFINITY) && (res->box_inf[i]!=INFINITY)){
				
			    zonotope_aff_add_itv(pr, res->paf[i], res->box_inf[i],res->box_sup[i], IN);
				
			} else {
				
			    res->paf[i]->c_inf = res->box_inf[i];
			    res->paf[i]->c_sup = res->box_sup[i];
			}
			
			res->paf[i]->pby++;
		    } 
		}
	    }
	} else {
		
	    /* nothing change */
	}
	//for(i=0; i < res->dims; i++){
	//	elina_interval_free(env[i]);	
	//}
	//free(env);
    //} 
	//printf("meet licnons finish1 %p %p\n",man,man->internal);
    //zonotope_fprint(stdout,man,res,NULL);
   //fflush(stdout);
    if (!is_bottom) {
	/* texpr -> aff forme */
	
	zonotope_aff_t** aff = (zonotope_aff_t **)malloc(array->size*sizeof(zonotope_aff_t*));
	for (i=0; i<array->size; i++) {
		aff[i]=NULL;
	}
	elina_linexpr0_t* linexpr0 = NULL;
	elina_lincons0_t lincons0;
	elina_interval_t cst;
	for (i=0; i<z->dims; i++){ 
		res->paf[i]->itv_inf = res->box_inf[i];
		res->paf[i]->itv_sup = res->box_sup[i];
        }
	for (i=0; i<array->size; i++) {
	    aff[i] = zonotope_aff_from_linexpr0(pr, array->p[i].linexpr0, res);
		//printf("after linexpr0: %d %d %p %d \n",i,array->size,aff[i],aff[i]->pby);
		//zonotope_aff_fprint(pr,stdout,aff[i]);
		//fflush(stdout);
	    linexpr0 = elina_linexpr0_from_zonotope(pr,aff[i],res);
        //elina_dimension_t dimension = elina_abstract0_dimension(pr->manNS,z->abs);
        //printf("dimension: %d %d\n",dimension.intdim,dimension.realdim);
        //fflush(stdout);
        
        //elina_dimension_t dimension2 = elina_abstract0_dimension(pr->manNS,res->abs);
        //printf("dimension2: %d %d\n",dimension2.intdim,dimension2.realdim);
        //fflush(stdout);
        
	    if (aff[i]->q != NULL) {
		/* only the centers are involved in this constraint, already treated while updating res->box */
		
		
		/* infer constraints on noise symbols */
		//linexpr0 = zonotope_elina_linexpr0_set_aff(pr, aff[i], res);
		
		lincons0.constyp = array->p[i].constyp;
		lincons0.linexpr0 = linexpr0;
		lincons0.scalar = array->p[i].scalar;
		elina_lincons0_array_t eps_lincons_array;
		eps_lincons_array.size = 1;
		eps_lincons_array.p = &lincons0;
            //printf("epsilon\n");
             //zonotope_aff_fprint(pr, stdout, aff[i]);
            //elina_abstract0_fprint(stdout,pr->manNS,res->abs,NULL);
            //elina_lincons0_array_fprint(stdout,&eps_lincons_array,NULL);
            //fflush(stdout);
            //printf("library before %s %s %d\n",pr->manNS->library,res->abs->man->library,destructive);
           //fflush(stdout);
		elina_abstract0_meet_lincons_array(pr->manNS, true, res->abs, &eps_lincons_array);
		//printf("output \n");
            //elina_abstract0_fprint(stdout,pr->manNS,res->abs,NULL);
           // printf("library after %s %p %d\n",pr->manNS->library,res->abs->man,destructive);
            //fflush(stdout);
		if (elina_abstract0_is_bottom(pr->manNS, res->abs)) {
			//printf("DEF\n");
			//fflush(stdout);
		    is_bottom = true;
		    break;
		}
	    }
	    elina_linexpr0_free(linexpr0);
	}
	/* update res->gamma */
	zonotope_update_noise_symbol_cons_gamma(pr, res);
	
	for (i=0; i<array->size; i++) {
	    /* update the abstract object with the new affine forms */
	    if (array->p[i].constyp == ELINA_CONS_EQ) {
		elina_interval_t *dummy = elina_interval_alloc();
		zonotope_aff_t* tmp, *tmp1;
		size_t j = 0;
		for (j=0; j<res->dims; j++) {
		    if (res->paf[j]->q) {
			
			zonotope_aff_cons_eq_lambda(pr, &dummy, res->paf[j], aff[i], res);
			tmp = zonotope_aff_mul_itv(pr, aff[i], dummy);
			res->paf[j]->itv_inf = res->box_inf[j];
			res->paf[j]->itv_sup = res->box_sup[j];
			tmp1 = zonotope_aff_add(pr, tmp, res->paf[j], res);
			zonotope_aff_check_free(pr, res->paf[j]);
			res->paf[j] = tmp1;
			/* update res->box */
			res->box_inf[j] = fmin(res->box_inf[j], res->paf[j]->itv_inf);
			res->box_sup[j] = fmin(res->box_sup[j], res->paf[j]->itv_sup);
			res->paf[j]->pby++;
			zonotope_aff_free(pr, tmp);
		    }
		}
		elina_interval_free(dummy);
	    } else {
		/* do nothing, just update res->box */
		size_t j = 0;
		//zonotope_fprint(stdout,man,res,NULL);
		for (j=0; j<res->dims; j++) {
			//zonotope_fprint(stdout,man,res,NULL);
			//zonotope_aff_fprint(pr, stdout, res->paf[j]);
		    zonotope_aff_bound(pr, &box_inf, &box_sup, res->paf[j], res);
			//printf("bounds; %g %g\n", box_inf, box_sup);
			//fflush(stdout);
		    res->box_inf[j] = fmin(res->box_inf[j], box_inf);
		    res->box_sup[j] = fmin(res->box_sup[j], box_sup);
		}
		//zonotope_fprint(stdout,man,res,NULL);
	    }
		//printf("aff: %d %d %p %d %d\n",i,array->size,aff[i],aff[i]->pby,is_bottom);
		//zonotope_aff_fprint(pr,stdout,aff[i]);
		//fflush(stdout);
	    zonotope_aff_check_free(pr, aff[i]);
	}
	free(aff);
    }
    if (is_bottom) {
	
	size_t intdim = res->intdim;
	size_t realdim = res->dims - res->intdim;
	zonotope_free(man, res);
	res = zonotope_bottom(man, intdim, realdim);
    }
    
    free(tchange);
   //printf("meet licnons finish %d\n",res->hypercube);
   // zonotope_fprint(stdout,man,res,NULL);
   //fflush(stdout);
    
    record_timing(zonotope_meet_lincons_time);
    return res;
}


/************************************************/
/* 2. Join					*/
/************************************************/

zonotope_t* zonotope_join(elina_manager_t* man, bool destructive, zonotope_t* z1, zonotope_t* z2)
    /* TODO destructive not used  */
{
    start_timing();
    size_t i = 0;
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_JOIN);
    if((z1->dims!=z2->dims) || (z1->intdim!=z2->intdim)){
	return NULL;
    }
    zonotope_t* res;
    size_t intdim = z1->intdim;
    size_t realdim = z1->dims - z1->intdim;
    if (zonotope_is_eq(man, z1, z2)) {
	if (destructive) res = z1;
	else res = zonotope_copy(man, z1);
    }
    else if (zonotope_is_top(man, z1) ||  zonotope_is_top(man, z2)) {
	if (destructive) {
	    zonotope_free(man, z1);
	    res = z1 = zonotope_top(man,intdim,realdim);
	} else {
	    res = zonotope_top(man,intdim,realdim);
	}
    } else if (zonotope_is_bottom(man, z1)) {
	if (destructive) {
	    zonotope_free(man, z1);
	    res = z1 = zonotope_copy(man,z2);
	} else {
	    res = zonotope_copy(man,z2);
	}
    } else if (zonotope_is_bottom(man, z2)) {
	if (destructive) res = z1;
	else res = zonotope_copy(man,z1);
    } else {
	/* TODO: destructive not yet supported */
	elina_interval_t *tmp = elina_interval_alloc();
	res = zonotope_alloc(man, intdim, realdim);
	/* update res->box */
	for (i=0; i<(intdim+realdim); i++){
	     res->box_inf[i] = fmax(z1->box_inf[i],z2->box_inf[i]);
	     res->box_sup[i] = fmax(z1->box_sup[i],z2->box_sup[i]);
	}
	
	if ((z1->hypercube && z2->hypercube)) {
	    for (i=0; i<(intdim+realdim); i++) {
		//printf("%d: ",i);
		if (zonotope_aff_is_bottom(pr, z1->paf[i])){
		    res->paf[i] = z2->paf[i];
		}
		else if (zonotope_aff_is_bottom(pr, z2->paf[i])){
		    res->paf[i] = z1->paf[i];
		}
		else if (zonotope_aff_is_top(pr, z1->paf[i]) || zonotope_aff_is_top(pr, z2->paf[i])){
			 res->paf[i] = pr->top;
		}		
		else {
		    if (!isfinite(z1->box_inf[i]) || !isfinite(z1->box_sup[i]) 
			|| !isfinite(z2->box_inf[i]) || !isfinite(z2->box_sup[i])) {
			/* Do nothing, the join of concretisations is already done and stored in res->box */
			res->paf[i] = zonotope_aff_alloc_init(pr);
			res->paf[i]->c_inf = res->box_inf[i];
			res->paf[i]->c_sup = res->box_sup[i];
			res->paf[i]->itv_inf = INFINITY;
			res->paf[i]->itv_sup = INFINITY;
		    } else {
			/* join two affine form expressions */
			z1->paf[i]->itv_inf = z1->box_inf[i];
			z1->paf[i]->itv_sup = z1->box_sup[i];
			z2->paf[i]->itv_inf = z2->box_inf[i];
			z2->paf[i]->itv_sup = z2->box_sup[i];
			res->paf[i] = zonotope_aff_join_constrained6(pr, z1->paf[i], z2->paf[i], z1, z2, res);
		    }
		}
		res->paf[i]->pby++;
	    }

	} else {
	    size_t k = 0;
	    elina_dim_t j = 0;
	    size_t dims1 = zonotope_noise_symbol_cons_get_dimension(pr, z1);
	    size_t dims2 = zonotope_noise_symbol_cons_get_dimension(pr, z2);
	 
	    if (dims1 && dims2) {
		
		size_t dim2 = 0;
		elina_dimchange_t* dimchange2 = elina_dimchange_alloc(0, dims1);
		//elina_dimchange_t* dimchange1 = elina_dimchange_alloc(0, dims2);
		if ((dims1+dims2) > res->size) {
		    res->nsymcons = (elina_dim_t*)realloc(res->nsymcons, (dims1+dims2)*sizeof(elina_dim_t));
		    res->gamma = (elina_interval_t**)realloc(res->gamma, (dims1+dims2)*sizeof(elina_interval_t*));
		    for (k=res->size;k<(dims1+dims2);k++) res->gamma[k] = NULL;
		    res->size = dims1+dims2;
		}
		res->nsymcons = memcpy((void *)res->nsymcons, (const void *)z1->nsymcons, dims1*sizeof(elina_dim_t));
		elina_abstract0_free(pr->manNS, res->abs);
		
		res->abs = elina_abstract0_copy(pr->manNS, z1->abs);
		for (k=0; k<dims2; k++) zonotope_insert_constrained_noise_symbol(pr, &j, z2->nsymcons[k], res);
		for (k=0; k<dims1; k++) {
		    if (!zonotope_noise_symbol_cons_get_dimpos(pr, &j, z1->nsymcons[k], z2)) {
			dimchange2->dim[dim2] = j;
			dim2++;
		    }
		}
		dimchange2->realdim = dim2;
		
		/* destructive, without projection (new dimension set to top) */
		elina_abstract0_t * tmp_z2 = elina_abstract0_add_dimensions(pr->manNS, false, z2->abs, dimchange2, false);
		//elina_abstract0_add_dimensions(pr->manNS, true, res->abs, dimchange1, false);
		elina_dimchange_add_invert(dimchange2);
		for (k=0; k<dim2; k++) {
		    zonotope_set_lincons_dim(pr, dimchange2->dim[k]);
		    elina_abstract0_meet_lincons_array(pr->manNS, true, tmp_z2, &pr->moo);
		}

			
       	        elina_abstract0_join(pr->manNS, true, res->abs, tmp_z2);
		
		
		/* update res->gamma */
		
		zonotope_update_noise_symbol_cons_gamma(pr, res);
		
		elina_abstract0_free(pr->manNS,tmp_z2);
		//elina_abstract0_remove_dimensions(pr->manNS, true, z2->abs, dimchange2);
		//dimchange2->realdim = dims2; 
		elina_dimchange_free(dimchange2);
		//elina_dimchange_free(dimchange1);
		size_t nsymcons_size = zonotope_noise_symbol_cons_get_dimension(pr, res);
		pr->dimtoremove = (elina_dim_t*)realloc(pr->dimtoremove, (nsymcons_size)*sizeof(elina_dim_t));
		memset((void *)pr->dimtoremove, (int)0, nsymcons_size*sizeof(int));
	    } else {
		/* res->abs is a hypercube */
	    }
	    for (i=0; i<(intdim+realdim); i++) {
		if (zonotope_aff_is_bottom(pr, z1->paf[i])){
		    res->paf[i] = z2->paf[i];
		}
		else if (zonotope_aff_is_bottom(pr, z2->paf[i])){
		    res->paf[i] = z1->paf[i];
	        }
		else if (zonotope_aff_is_top(pr, z1->paf[i]) || zonotope_aff_is_top(pr, z2->paf[i])){
		    res->paf[i] = pr->top;
		}
		else if (zonotope_aff_is_eq(pr, z1->paf[i], z2->paf[i])){
			
		    res->paf[i] = z1->paf[i];
		}
		else {
		    if (!isfinite(z1->box_inf[i]) || !isfinite(z1->box_sup[i])
			 || !isfinite(z2->box_inf[i]) || !isfinite(z2->box_sup[i])) {
			/* Do nothing, the join of concretisations is already done and stored in res->box */
			res->paf[i] = zonotope_aff_alloc_init(pr);
			res->paf[i]->c_inf = res->box_inf[i];
			res->paf[i]->c_sup = res->box_sup[i];
			res->paf[i]->itv_inf = INFINITY;
			res->paf[i]->itv_sup = INFINITY;
			
		    } else {
			/* join two affine form expressions */
			z1->paf[i]->itv_inf = z1->box_inf[i];
			z1->paf[i]->itv_sup = z1->box_sup[i];
			z2->paf[i]->itv_inf = z2->box_inf[i];
			z2->paf[i]->itv_sup = z2->box_sup[i];
			res->paf[i] = zonotope_aff_join_constrained6(pr, z1->paf[i], z2->paf[i], z1, z2, res);
		    }
		}
		res->paf[i]->pby++;
	    }

	    man->result.flag_best = false;
	    man->result.flag_exact = false;
	}
	
	man->result.flag_best = true;
	man->result.flag_exact = false;
	elina_interval_free(tmp);
    }
    man->result.flag_best = true;
    man->result.flag_exact = true;
   
    record_timing(zonotope_join_time);
    return res;
}



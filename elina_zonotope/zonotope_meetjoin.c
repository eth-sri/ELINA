/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
 *  This software is distributed under GNU Lesser General Public License
 * Version 3.0. For more information, see the ELINA project website at:
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





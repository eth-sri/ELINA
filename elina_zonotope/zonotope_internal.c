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

#include <stdlib.h>
#include <unistd.h>

#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>


#include "zonotope_internal.h"
//#include "zonotope_itv_utils.h"

#include "zonotope_representation.h"
#include "zonotope_constructor.h"
#include "zonotope_meetjoin.h"
#include "zonotope_assign.h"
#include "zonotope_resize.h"
#include "zonotope_otherops.h"

double zonotope_copy_time=0;
double zonotope_is_equal_time=0;
double zonotope_is_lequal_time=0;
double zonotope_permute_dimension_time=0;
double zonotope_add_dimension_time=0;
double zonotope_remove_dimension_time=0;
double zonotope_top_time=0;
double zonotope_bottom_time=0;
double zonotope_join_time=0;
double zonotope_free_time=0;
double zonotope_forget_array_time=0;
double zonotope_meet_lincons_time=0;
double zonotope_to_box_time=0;
double zonotope_of_box_time=0;
double zonotope_is_top_time=0;
double zonotope_is_bottom_time=0;
double zonotope_assign_linexpr_time=0;

elina_manager_t* zonotope_manager_alloc(void)
{
	//CALL();
	elina_manager_t* man;
	void** funptr;

	
	elina_manager_t* manNS = elina_box_manager_alloc();
	zonotope_internal_t *zt = zonotope_internal_alloc(manNS);

	man = elina_manager_alloc("Zonotope",/* Library name */
			"1.0", /* version */
			zt, /* internal structure */
			(void (*)(void*))zonotope_internal_free /* free function for internal */
			);

	funptr = man->funptr;

	/* Internal representation */
	/***************************/
	/* 1.Memory */
	funptr[ELINA_FUNID_COPY] = &zonotope_copy;
	funptr[ELINA_FUNID_FREE] = &zonotope_free;
	/*funptr[ELINA_FUNID_SIZE] = &zonotope_size;*/
	//funptr[ELINA_FUNID_ASIZE] = &zonotope_size;
	/* 2.Control of internal representation */
	//funptr[ELINA_FUNID_MINIMIZE] = &zonotope_minimize;
	//funptr[ELINA_FUNID_CANONICALIZE] = &zonotope_canonicalize;
	//funptr[ELINA_FUNID_HASH] = &zonotope_hash;
	//funptr[ELINA_FUNID_APPROXIMATE] = &zonotope_approximate;
	/* 3.Printing */
	funptr[ELINA_FUNID_FPRINT] = &zonotope_fprint;
	//funptr[ELINA_FUNID_FPRINTDIFF] = &zonotope_fprintdiff;
	//funptr[ELINA_FUNID_FDUMP] = &zonotope_fdump;
	/* 4.Serialisation */
	//funptr[ELINA_FUNID_SERIALIZE_RAW] = &zonotope_serialize_raw;
	//funptr[ELINA_FUNID_DESERIALIZE_RAW] = &zonotope_deserialize_raw;

	/* Constructors */
	/****************/
	/* 1.Basic constructors */
	funptr[ELINA_FUNID_BOTTOM] = &zonotope_bottom;
	funptr[ELINA_FUNID_TOP] = &zonotope_top;
	funptr[ELINA_FUNID_OF_BOX] = &zonotope_of_box;
	/* 2.Accessors */
	funptr[ELINA_FUNID_DIMENSION] = &zonotope_dimension;
	/* 3.Tests */
	funptr[ELINA_FUNID_IS_BOTTOM] = &zonotope_is_bottom;
	funptr[ELINA_FUNID_IS_TOP] = &zonotope_is_top;
	//funptr[ELINA_FUNID_IS_LEQ] = &zonotope_is_leq;
	funptr[ELINA_FUNID_IS_EQ] = &zonotope_is_eq;
	//funptr[ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED] = &zonotope_is_dimension_unconstrained;
	//funptr[ELINA_FUNID_SAT_TCONS] = &zonotope_sat_tcons; /*  */
	//funptr[ELINA_FUNID_SAT_INTERVAL] = &zonotope_sat_interval;
	//funptr[ELINA_FUNID_SAT_LINCONS] = &zonotope_sat_lincons;
	/* 4.Extraction of properties */
	//funptr[ELINA_FUNID_BOUND_TEXPR] = &zonotope_bound_texpr; /*  */
	funptr[ELINA_FUNID_BOUND_DIMENSION] = &zonotope_bound_dimension;
	//funptr[ELINA_FUNID_BOUND_LINEXPR] = &zonotope_bound_linexpr;
	funptr[ELINA_FUNID_TO_BOX] = &zonotope_to_box;
	//funptr[ELINA_FUNID_TO_TCONS_ARRAY] = &zonotope_to_tcons_array; /*  */
	funptr[ELINA_FUNID_TO_LINCONS_ARRAY] = &zonotope_to_lincons_array;
	//funptr[ELINA_FUNID_TO_GENERATOR_ARRAY] = &zonotope_to_generator_array;

	/* Meet and Join */
	/*****************/
	/* 1.Meet */
	//funptr[ELINA_FUNID_MEET] = &zonotope_meet; /* */
	//funptr[ELINA_FUNID_MEET_ARRAY] = &zonotope_meet_array; /*  */
	funptr[ELINA_FUNID_MEET_LINCONS_ARRAY] = &zonotope_meet_lincons_array; /*  */
	//funptr[ELINA_FUNID_MEET_TCONS_ARRAY] = &zonotope_meet_tcons_array; /*  */
	/* 2.Join */
	funptr[ELINA_FUNID_JOIN] = &zonotope_join;
	//funptr[ELINA_FUNID_JOIN_ARRAY] = &zonotope_join_array;

	//funptr[ELINA_FUNID_ADD_RAY_ARRAY] = &zonotope_add_ray_array;

	/* Assign and Substitute */
	/*************************/
	funptr[ELINA_FUNID_ASSIGN_LINEXPR_ARRAY] = &zonotope_assign_linexpr_array;
	//funptr[ELINA_FUNID_SUBSTITUTE_LINEXPR_ARRAY] = &zonotope_substitute_linexpr_array;
	//funptr[ELINA_FUNID_ASSIGN_TEXPR_ARRAY] = &zonotope_assign_texpr_array;
	//funptr[ELINA_FUNID_SUBSTITUTE_TEXPR_ARRAY] = &zonotope_substitute_texpr_array;

	/* Resize dimensions */
	/*********************/
	funptr[ELINA_FUNID_ADD_DIMENSIONS] = &zonotope_add_dimensions;
	funptr[ELINA_FUNID_REMOVE_DIMENSIONS] = &zonotope_remove_dimensions;
	funptr[ELINA_FUNID_PERMUTE_DIMENSIONS] = &zonotope_permute_dimensions;

	/* Other functions */
	/*******************/
	funptr[ELINA_FUNID_FORGET_ARRAY] = &zonotope_forget_array;
	//funptr[ELINA_FUNID_EXPAND] = &zonotope_expand;
	//funptr[ELINA_FUNID_FOLD] = &zonotope_fold;
	//funptr[ELINA_FUNID_WIDENING] = &zonotope_widening;
	//funptr[ELINA_FUNID_CLOSURE] = &zonotope_closure;

	/* ?! */
	/******/
	/*
	   ELINA_FUNID_CHANGE_ENVIRONMENT
	   ELINA_FUNID_RENAME_ARRAY
	   ELINA_FUNID_SIZE2
	 */
	man->option.abort_if_exception[ELINA_EXC_INVALID_ARGUMENT] = false;
	return man;
}

/* back pointer to our internal structure from the manager */
zonotope_internal_t* zonotope_init_from_manager(elina_manager_t* man, elina_funid_t funid)
{
	//printf("manager %p\n",man);
	//fflush(stdout);
    zonotope_internal_t* pr = (zonotope_internal_t*)man->internal;
    pr->funid = funid;
	//printf("funid: %d\n",funid);
	//fflush(stdout);
    if (!(pr->man)) pr->man = man;
	//printf("finish\n");
	//fflush(stdout);
    return pr;
}




zonotope_aff_t* zonotope_aff_add(zonotope_internal_t* pr, zonotope_aff_t* exprA, zonotope_aff_t* exprB, zonotope_t* abs)
{
    double box_inf = 0.0;
    double box_sup = 0.0;
    double tmp_inf = 0.0;
    double tmp_sup = 0.0;
    
    zonotope_aff_t* res = zonotope_aff_alloc_init(pr);
    zonotope_aaterm_t *p, *q, *ptr;
    double maxA = fmax(fabs(exprA->c_inf),fabs(exprA->c_sup));
    double maxB = fmax(fabs(exprB->c_inf),fabs(exprB->c_sup));
    //double fp_err_inf = (maxA + maxB)*pr->ulp + pr->min_denormal;
    //double fp_err_sup = (maxA + maxB)*pr->ulp + pr->min_denormal;
    res->c_inf = exprA->c_inf + exprB->c_inf + (maxA + maxB)*pr->ulp + pr->min_denormal;
    res->c_sup = exprA->c_sup + exprB->c_sup + (maxA + maxB)*pr->ulp + pr->min_denormal;
	
    box_inf = res->c_inf;
    box_sup = res->c_sup;
    //size_t count = 0;
    if (exprA->q || exprB->q) {
        
        ptr = zonotope_aaterm_alloc_init();
        for(p = exprA->q, q = exprB->q; p || q;) {
            //count++;
            if (p && q) {
                
                if (p->pnsym->index == q->pnsym->index) {
                    maxA = fmax(fabs(p->inf),fabs(p->sup));
                    maxB = fmax(fabs(q->inf),fabs(q->sup));
                    
                    ptr->inf = p->inf + q->inf + (maxA + maxB)*pr->ulp;
                    ptr->sup = p->sup + q->sup + (maxA + maxB)*pr->ulp;
                    //fp_err_inf += (maxA + maxB)*pr->ulp;
                    //fp_err_sup += (maxA + maxB)*pr->ulp;
                    ptr->pnsym = p->pnsym;
                    p = p->n ;
                    q = q->n ;
			
                } else if (p->pnsym->index < q->pnsym->index) {
		    ptr->inf = p->inf;
		    ptr->sup = p->sup;
                    ptr->pnsym = p->pnsym;
                    p = p->n ;
                } else {
			
		    ptr->inf = q->inf;
		    ptr->sup = q->sup;
                    ptr->pnsym = q->pnsym;
                    q = q->n ;
                }
                
            } else if (p) {
		ptr->inf = p->inf;
		ptr->sup = p->sup;
                ptr->pnsym = p->pnsym;
                p = p->n ;
            } else {
		ptr->inf = q->inf;
		ptr->sup = q->sup;
                ptr->pnsym = q->pnsym;
                q = q->n ;
            }
            
            if (!ptr->inf && !ptr->sup) {
                if (!(p||q)) {
                    /* the last iteration */
                    zonotope_aaterm_free(pr, ptr);
                    if (res->end) res->end->n = NULL;
                }
            } else {
                /* keep this term */
               
                if (!res->q) res->q = ptr;
                res->end = ptr;
                res->l++;
               
                zonotope_noise_symbol_cons_get_gamma(pr, &tmp_inf, &tmp_sup, ptr->pnsym->index, abs);
		
		double inf = 0.0;
		double sup = 0.0;
                elina_double_interval_mul(&inf,&sup,tmp_inf,tmp_sup,ptr->inf,ptr->sup);
		
                box_inf =  box_inf + inf;
                box_sup =  box_sup + sup;
		//printf("%g %g %g %g %g %g %d\n",-inf,sup,-tmp_inf ,tmp_sup,-box_inf,box_sup,abs->hypercube);
                 //start_timing();
                if (p||q) {
                    /* continuing */
                    ptr->n = zonotope_aaterm_alloc_init();
                    ptr=ptr->n;
                }
                 //record_timing(zonotope_assign_linexpr_time);
            }
            
        }
        
    }
    
    //elina_interval_add(box, box, tmp, ELINA_SCALAR_DOUBLE);
    //elina_interval_add(res->itv, exprA->itv, exprB->itv, ELINA_SCALAR_DOUBLE);
    res->itv_inf =  exprA->itv_inf + exprB->itv_inf;// + fp_err_inf;
    res->itv_sup =  exprA->itv_sup + exprB->itv_sup;// + fp_err_sup;
	
	//printf("Box %g %g\n",-box_inf,box_sup);
    res->itv_inf = fmin(res->itv_inf,box_inf);
    res->itv_sup = fmin(res->itv_sup,box_sup);

    //printf("Output %g %g\n",res->itv_inf,res->itv_sup);
   //zonotope_aff_fprint(pr,stdout,res);
	
    return res;
}

zonotope_aff_t * zonotope_aff_from_linexpr0(zonotope_internal_t* pr, elina_linexpr0_t * expr, zonotope_t *z){
    size_t i;
    elina_dim_t dim;
    elina_coeff_t *coeff;
    zonotope_aff_t *res = zonotope_aff_alloc_init(pr);
    elina_coeff_t * cst = &(expr->cst);

    if(cst->discr==ELINA_COEFF_SCALAR) {
        if (cst->val.scalar->discr == ELINA_SCALAR_DOUBLE) {
            res->c_inf = -cst->val.scalar->val.dbl;
            res->c_sup = cst->val.scalar->val.dbl;
            res->itv_inf = -cst->val.scalar->val.dbl;
            res->itv_sup = cst->val.scalar->val.dbl;
        } else {
            double inf,sup;
            elina_double_set_scalar(&inf,cst->val.scalar,GMP_RNDD);
            elina_double_set_scalar(&sup,cst->val.scalar,GMP_RNDU);
            res->c_inf = -inf;
            res->c_sup = sup;
            res->itv_inf = -inf;
            res->itv_sup = sup;
        }
    } else {
        if (cst->val.interval->inf->discr == ELINA_SCALAR_DOUBLE) {
            res->c_inf = -cst->val.interval->inf->val.dbl;
            res->itv_inf = -cst->val.interval->inf->val.dbl;
        } else {
            double inf;
            elina_double_set_scalar(&inf,cst->val.scalar,GMP_RNDD);
            res->c_inf = -inf;
            res->itv_inf = -inf;
        }
        if (cst->val.interval->sup->discr == ELINA_SCALAR_DOUBLE) {
            res->c_sup = cst->val.interval->sup->val.dbl;
            res->itv_sup = cst->val.interval->sup->val.dbl;
        } else {
            double sup;
            elina_double_set_scalar(&sup,cst->val.scalar,GMP_RNDU);
            res->c_sup = sup;
            res->itv_sup = sup;
        }
    }

	//printf("size: %zu\n",expr->size);
	//fflush(stdout);
    elina_linexpr0_ForeachLinterm(expr,i,dim,coeff) {
        zonotope_aff_t *aff = z->paf[dim];
        zonotope_aff_t *tmp;
        if(coeff->discr==ELINA_COEFF_SCALAR){
            
            elina_scalar_t * scalar = elina_scalar_alloc();
            //tmp = zonotope_aff_mul_scalar(pr,aff,scalar);
            elina_scalar_set(scalar,coeff->val.scalar);
            elina_coeff_reinit(coeff,ELINA_COEFF_INTERVAL,ELINA_SCALAR_DOUBLE);
            elina_coeff_set_interval_scalar(coeff,scalar,scalar);
            elina_scalar_free(scalar);
        }
        //else{
          
        //fflush(stdout);
        elina_interval_t *interval = coeff->val.interval;
        //start_timing();
	//elina_interval_fprint(stdout,interval);
	//fflush(stdout);
        tmp = zonotope_aff_mul_itv(pr,aff,interval);
        //}
         //
        //printf("mul itv\n");
        //zonotope_aff_fprint(pr,stdout,tmp);
        //elina_interval_fprint(stdout,interval);
        
        //printf("i: %zu k: %u %g %g\n",i,dim,-tmp->itv_inf,tmp->itv_sup);
        
        
        zonotope_aff_t *tmp1 = res;
       
	//if(tmp->q==NULL){
		 //printf("add input %p\n",tmp->q);
        	//zonotope_aff_fprint(pr,stdout,tmp1);
		//zonotope_aff_fprint(pr,stdout,tmp);
        	//fflush(stdout);
	//}
        res = zonotope_aff_add(pr,tmp1,tmp,z);
	//if(tmp->q==NULL){
	//printf("add output\n");
       // zonotope_aff_fprint(pr,stdout,res);
	//fflush(stdout);
	//}
        //record_timing(zonotope_assign_linexpr_time);
        zonotope_aff_free(pr,tmp);
        zonotope_aff_free(pr,tmp1);
        //printf("i: %d dim: %d\n",i,dim);
	//printf("Output %g %g\n",-res->itv_inf,res->itv_sup);
        //zonotope_aff_fprint(pr,stdout,res);
    }
	
    return res;
}

zonotope_aff_t* zonotope_aff_mul_itv(zonotope_internal_t* pr, zonotope_aff_t* src, elina_interval_t *lambda)
{
   
    if ((!elina_scalar_sgn(lambda->inf) && !elina_scalar_sgn(lambda->sup) )|| zonotope_aff_is_known_to_be_zero(pr, src)) {
        return zonotope_aff_alloc_init(pr);
    } else if (zonotope_aff_is_bottom(pr, src) || elina_interval_is_bottom(lambda)) {
        return zonotope_aff_bottom(pr);
    } else if (zonotope_aff_is_top(pr, src) || elina_interval_is_top(lambda)) {
        return zonotope_aff_top(pr);
    } else {
         //start_timing();
        zonotope_aff_t* dst = NULL;
        zonotope_aaterm_t *p,*q;
        q = NULL;
        dst = zonotope_aff_alloc_init(pr);
         //
        elina_double_interval_mul(&dst->c_inf,&dst->c_sup, -lambda->inf->val.dbl, lambda->sup->val.dbl, src->c_inf,src->c_sup);
        //printf("coming here %g %g %g %g %g %g\n",dst->c_inf,dst->c_sup,-lambda->inf->val.dbl, lambda->sup->val.dbl, src->c_inf,src->c_sup);
	//fflush(stdout);
        if (src->q) {
            
            dst->q = q = zonotope_aaterm_alloc_init();
            for (p=src->q; p; p=p->n) {
		
                double tmp_inf = 0.0;
		double tmp_sup = 0.0;
                elina_double_interval_mul(&tmp_inf, &tmp_sup, -lambda->inf->val.dbl, lambda->sup->val.dbl, p->inf, p->sup);
		q->inf = tmp_inf;
		q->sup = tmp_sup;
                q->pnsym = p->pnsym;
                if (p->n) {
                    /* continue */
                    q->n = zonotope_aaterm_alloc_init();
                    q = q->n;
                } else {
                    /* the last iteration */
                    dst->end = q; 
                }
            }
            
        }
        
        dst->l = src->l;
	//printf("inputs: %g %g %g %g\n",-lambda->inf->val.dbl,lambda->sup->val.dbl,src->itv_inf, src->itv_sup);
        elina_double_interval_mul(&dst->itv_inf, &dst->itv_sup,  -lambda->inf->val.dbl, lambda->sup->val.dbl,src->itv_inf, src->itv_sup);
	//printf("outputs %g %g\n",dst->itv_inf,dst->itv_sup);
       // record_timing(zonotope_assign_linexpr_time);
        return dst;
    }
}

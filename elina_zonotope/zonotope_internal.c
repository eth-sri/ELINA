/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
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
	//funptr[ELINA_FUNID_BOUND_DIMENSION] = &zonotope_bound_dimension;
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
	//printf("funid: %d %d\n",funid,funid==ELINA_FUNID_UNKNOWN);
	//fflush(stdout);
    zonotope_internal_t* pr = (zonotope_internal_t*)man->internal;
    pr->funid = funid;
    if (!(pr->man)) pr->man = man;
    return pr;
}

zonotope_aff_t* zonotope_aff_add(zonotope_internal_t* pr, zonotope_aff_t* exprA, zonotope_aff_t* exprB, zonotope_t* abs)
{
    elina_interval_t *box = elina_interval_alloc();
    elina_interval_t *tmp = elina_interval_alloc();
    
    zonotope_aff_t* res = zonotope_aff_alloc_init(pr);
    zonotope_aaterm_t *p, *q, *ptr;
    elina_interval_add(res->c, exprA->c, exprB->c, ELINA_SCALAR_DOUBLE);
    elina_interval_set(box, res->c);
    
    if (exprA->q || exprB->q) {
        
        ptr = zonotope_aaterm_alloc_init();
        for(p = exprA->q, q = exprB->q; p || q;) {
            
            if (p && q) {
                
                if (p->pnsym->index == q->pnsym->index) {
                    elina_interval_add(ptr->coeff, p->coeff, q->coeff, ELINA_SCALAR_DOUBLE);
                    ptr->pnsym = p->pnsym;
                    p = p->n ;
                    q = q->n ;
                } else if (p->pnsym->index < q->pnsym->index) {
                    elina_interval_set(ptr->coeff, p->coeff);
                    ptr->pnsym = p->pnsym;
                    p = p->n ;
                } else {
                    elina_interval_set(ptr->coeff, q->coeff);
                    ptr->pnsym = q->pnsym;
                    q = q->n ;
                }
                
            } else if (p) {
                elina_interval_set(ptr->coeff, p->coeff);
                ptr->pnsym = p->pnsym;
                p = p->n ;
            } else {
                elina_interval_set(ptr->coeff, q->coeff);
                ptr->pnsym = q->pnsym;
                q = q->n ;
            }
            
            if (!ptr->coeff->inf->val.dbl && !ptr->coeff->sup->val.dbl) {
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
                ;
                zonotope_noise_symbol_cons_get_gamma(pr, tmp, ptr->pnsym->index, abs);
               
                //start_timing();
                //elina_interval_fprint(stdout,tmp);
                //elina_interval_fprint(stdout,ptr->coeff);
                //printf("\n");
                //if(elina_scalar_equal(ptr->coeff->inf,ptr->coeff->sup)) {
                
               //
                if((tmp->inf->val.dbl==pr->muu->inf->val.dbl) && (tmp->sup->val.dbl==pr->muu->sup->val.dbl)){
                //if(elina_interval_equal(tmp,pr->muu)){
                    if(ptr->coeff->sup->val.dbl>=0){
                        tmp->inf->val.dbl =- ptr->coeff->sup->val.dbl;
                        tmp->sup->val.dbl = ptr->coeff->sup->val.dbl;
                        //elina_scalar_neg(tmp->inf,ptr->coeff->sup);
                        //elina_scalar_set(tmp->sup,ptr->coeff->sup);
                    }
                    else{
                        //elina_scalar_t * add = elina_scalar_alloc_set(tmp->sup);
                        //elina_scalar_neg(tmp->sup,ptr->coeff->sup);
                        //elina_scalar_set(tmp->inf,ptr->coeff->sup);
                        //elina_scalar_free(add);
                        tmp->sup->val.dbl = -ptr->coeff->sup->val.dbl;
                        tmp->inf->val.dbl = ptr->coeff->sup->val.dbl;
                    }
                }
                else{
                    if(elina_scalar_sgn(ptr->coeff->sup)>=0){
                        elina_scalar_mul(tmp->inf,tmp->inf,ptr->coeff->sup,ELINA_SCALAR_DOUBLE);
                        elina_scalar_mul(tmp->sup,tmp->sup,ptr->coeff->sup,ELINA_SCALAR_DOUBLE);
                    }
                    else{
                        elina_scalar_t * add = elina_scalar_alloc_set(tmp->sup);
                        elina_scalar_mul(tmp->sup,tmp->inf,ptr->coeff->sup,ELINA_SCALAR_DOUBLE);
                        elina_scalar_mul(tmp->inf,add,ptr->coeff->sup,ELINA_SCALAR_DOUBLE);
                        elina_scalar_free(add);
                    }
                }
                
                
                //}
                //else{
                //  elina_interval_mul(tmp, tmp, ptr->coeff, ELINA_SCALAR_DOUBLE);
                // }
                //printf("SCALAR: %d %d %d\n",box->inf->discr==ELINA_SCALAR_MPQ,tmp->inf->discr==ELINA_SCALAR_MPQ, tmp->sup->discr==ELINA_SCALAR_MPQ);
                //elina_interval_add(box, box, tmp, ELINA_SCALAR_DOUBLE);
                box->inf->val.dbl =  box->inf->val.dbl + tmp->inf->val.dbl;
                box->sup->val.dbl =  box->sup->val.dbl + tmp->sup->val.dbl;
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
    res->itv->inf->val.dbl =  exprA->itv->inf->val.dbl + exprB->itv->inf->val.dbl;
    res->itv->sup->val.dbl =  exprA->itv->sup->val.dbl + exprB->itv->sup->val.dbl;
    elina_scalar_max(res->itv->inf, res->itv->inf, box->inf);
    elina_scalar_min(res->itv->sup, res->itv->sup, box->sup);
    elina_interval_free(box);
    elina_interval_free(tmp);
   
    return res;
}

zonotope_aff_t * zonotope_aff_from_linexpr0(zonotope_internal_t* pr, elina_linexpr0_t * expr, zonotope_t *z){
    size_t i;
    elina_dim_t dim;
    elina_coeff_t *coeff;
    zonotope_aff_t *res = zonotope_aff_alloc_init(pr);
    elina_coeff_t * cst = &(expr->cst);
    if(cst->discr==ELINA_COEFF_SCALAR){
        elina_interval_set_scalar(res->c, cst->val.scalar,cst->val.scalar);
        elina_interval_set_scalar(res->itv, cst->val.scalar,cst->val.scalar);
    }
    else{
        elina_interval_set(res->c, cst->val.interval);
        elina_interval_set(res->itv, cst->val.interval);
    }
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
        //  printf("interval\n");
        //fflush(stdout);
        elina_interval_t *interval = coeff->val.interval;
        //start_timing();
        tmp = zonotope_aff_mul_itv(pr,aff,interval);
        //}
         //
        //printf("mul itv\n");
        //zonotope_aff_fprint(pr,stdout,aff);
        //elina_interval_fprint(stdout,interval);
        
        
        //printf("result\n");
        //zonotope_aff_fprint(pr,stdout,tmp);
        //fflush(stdout);
        
        zonotope_aff_t *tmp1 = res;
        
        res = zonotope_aff_add(pr,tmp1,tmp,z);
        //record_timing(zonotope_assign_linexpr_time);
        zonotope_aff_free(pr,tmp);
        zonotope_aff_free(pr,tmp1);
        
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
        elina_interval_mul(dst->c, lambda, src->c,ELINA_SCALAR_DOUBLE);
       
        if (src->q) {
            
            dst->q = q = zonotope_aaterm_alloc_init();
            for (p=src->q; p; p=p->n) {
                
                //printf("SCALAR: %d %d %d\n",q->coeff->inf->discr,lambda->inf->discr,p->coeff->inf->discr);
                elina_interval_mul(q->coeff, lambda, p->coeff,ELINA_SCALAR_DOUBLE);
                
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
        elina_interval_mul(dst->itv, src->itv, lambda,ELINA_SCALAR_DOUBLE);
        //record_timing(zonotope_assign_linexpr_time);
        return dst;
    }
}

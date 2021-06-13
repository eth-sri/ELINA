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

#include "elina_linexpr0_arith.h"

void elina_linexpr0_reinit(elina_linexpr0_t* expr, size_t size){
  if(size==expr->size){
	
	return;
  }
  size_t i;
  switch(expr->discr){
	case ELINA_LINEXPR_DENSE:
		for  (i=size; i<expr->size; i++){
			elina_coeff_clear(&expr->p.coeff[i]);
		}
		expr->p.coeff = realloc(expr->p.coeff,size*sizeof(elina_coeff_t));
  		for (i=expr->size;i<size;i++){
    			elina_coeff_init(&expr->p.coeff[i],ELINA_COEFF_SCALAR);
  		}
		break;
	case ELINA_LINEXPR_SPARSE:
		for  (i=size; i<expr->size; i++){
			elina_coeff_clear(&expr->p.linterm[i].coeff);
		}
		expr->p.linterm = realloc(expr->p.linterm,size*sizeof(elina_linterm_t));
  		for (i=expr->size;i<size;i++){
    			elina_coeff_init(&expr->p.linterm[i].coeff,ELINA_COEFF_SCALAR);
			expr->p.linterm[i].dim = ELINA_DIM_MAX;
  		}
		break;
  }
  
  expr->size = size;
  return;
}


void elina_linexpr0_init(elina_linexpr0_t* expr, size_t size)
{
  size_t i;
  expr->size = size;
  expr->discr = ELINA_LINEXPR_SPARSE;
  elina_coeff_init(&expr->cst,ELINA_COEFF_SCALAR);
  expr->p.linterm = malloc(size*sizeof(elina_linterm_t));
  for (i=0;i<size;i++){
    	elina_coeff_init(&expr->p.linterm[i].coeff,ELINA_COEFF_SCALAR);
	expr->p.linterm[i].dim = ELINA_DIM_MAX;
  }
}


void elina_linexpr0_clear(elina_linexpr0_t* e)
{
  size_t i;

  if(e!=NULL){
    elina_coeff_clear(&e->cst);
    switch (e->discr){
    case ELINA_LINEXPR_DENSE:
      if (e->p.coeff!=NULL){
	for (i=0; i<e->size; i++) elina_coeff_clear(&e->p.coeff[i]);
	free(e->p.coeff);
      }
      break;
    case ELINA_LINEXPR_SPARSE:
      if (e->p.linterm!=NULL){
	for (i=0; i<e->size; i++) elina_coeff_clear(&e->p.linterm[i].coeff);
	free(e->p.linterm);
      }
      break;
    }
  }
}

void elina_linexpr0_neg(elina_linexpr0_t* expr)
{
  size_t i;
  elina_dim_t dim;
  elina_coeff_t *coeff;

  elina_coeff_neg(&expr->cst,&expr->cst);
  elina_linexpr0_ForeachLinterm(expr,i,dim,coeff){
    elina_coeff_neg(coeff,coeff);
  }
  return;
}




void elina_linexpr0_scale(elina_linexpr0_t* expr, elina_interval_t *interval, elina_scalar_discr_t discr)
{
  size_t i;
  elina_dim_t dim;
  elina_coeff_t *coeff, *cst;
  cst = &expr->cst;
  if (!elina_scalar_sgn(interval->inf) && !elina_scalar_sgn(interval->sup)){
    elina_coeff_reinit(cst,ELINA_COEFF_SCALAR,discr);
    elina_scalar_set_to_int(cst->val.scalar,0,discr);
    elina_linexpr0_reinit(expr,0);
    return;
  }
  elina_coeff_mul_interval(cst,cst,interval,discr);
  
  if (cst->discr==ELINA_COEFF_INTERVAL && elina_interval_is_top(cst->val.interval)){
    elina_linexpr0_reinit(expr,0);
    return;
  }
  else {
    elina_linexpr0_ForeachLinterm(expr,i,dim,coeff){
      elina_coeff_mul_interval(coeff,coeff,interval,discr);
    }
  }
  return;
}

void elina_linexpr0_add(elina_linexpr0_t **linres, elina_linexpr0_t** linexprA, elina_linexpr0_t** linexprB, elina_scalar_discr_t discr)
{
  elina_linexpr0_t * res = *linres;
  elina_linexpr0_t * exprA = *linexprA;
  elina_linexpr0_t * exprB = *linexprB;
  size_t i,j,k;
  elina_linexpr0_t *expr;
  if (res==exprA || res==exprB){
    //elina_linexpr0_init(&expr,exprA->size+exprB->size);
    expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,exprA->size+exprB->size);
  }
  else {
    expr = res;
    elina_linexpr0_reinit(expr,exprA->size+exprB->size);
  }
  bool endA,endB;
  elina_coeff_t *cst;
  elina_linterm_t *term, *termA, *termB;
  i = j = k = 0;
  endA = endB = false;
  cst = &expr->cst;
  elina_coeff_add(cst,&exprA->cst,&exprB->cst,discr);
  //expr.equality = exprA->equality && exprB->equality && itv_is_point(intern,expr.cst);
	
  if (cst->discr==ELINA_COEFF_INTERVAL && elina_interval_is_top(cst->val.interval))
    goto _elina_linexpr0_add_return;
  while (true){
    termA = &exprA->p.linterm[i];
    termB = &exprB->p.linterm[j];
    endA = endA || (i==exprA->size) || termA->dim == ELINA_DIM_MAX;
    endB = endB || (j==exprB->size) || termB->dim == ELINA_DIM_MAX;
    if (endA && endB)
      break;
    term = &expr->p.linterm[k];
    if (endA || (!endB && termB->dim < termA->dim)){
      elina_coeff_set(&term->coeff, &termB->coeff);
      term->dim = termB->dim; 
      k++; j++;
    }
    else if (endB || (!endA && termA->dim < termB->dim)){
      elina_coeff_set(&term->coeff,&termA->coeff);
      term->dim = termA->dim;
      k++; i++;
    }
    else {
      elina_coeff_add(&term->coeff, &termA->coeff, &termB->coeff,discr);
      term->dim = termA->dim;
      if (!elina_coeff_zero(&term->coeff)){
	k++;
      }
      i++; j++;
    }
  }
 _elina_linexpr0_add_return:
  elina_linexpr0_reinit(expr,k);
   if (res==exprA || res==exprB){
    elina_linexpr0_clear(res);
    //free(res);
  }
  //res = elina_linexpr0_copy(&expr);
  **linres = *expr;
  return;
}

void elina_linexpr0_sub(elina_linexpr0_t ** linres, elina_linexpr0_t** linexprA, elina_linexpr0_t** linexprB, elina_scalar_discr_t discr)
{
  elina_linexpr0_t * res = *linres;
  elina_linexpr0_t * exprA = *linexprA;
  elina_linexpr0_t * exprB = *linexprB;
  if (exprA==exprB){
    elina_linexpr0_t *expr = elina_linexpr0_copy(exprB);
    elina_linexpr0_neg(expr);
    elina_linexpr0_add(linres,linexprA,&expr,discr);
    elina_linexpr0_free(expr);
  }
  else {
    elina_linexpr0_neg(exprB);
    elina_linexpr0_add(linres,linexprA,linexprB,discr);
    if(res!=exprB){
    	elina_linexpr0_neg(exprB);
    }
  }
  return;
}

void elina_linexpr0_div(elina_linexpr0_t* expr, elina_interval_t *interval, elina_scalar_discr_t discr)
{
  size_t i;
  elina_dim_t dim;
  elina_coeff_t *coeff;
  elina_interval_t * tmp = elina_interval_alloc();
  elina_coeff_t *cst = &expr->cst;
  if(cst->discr==ELINA_COEFF_SCALAR){
	elina_interval_set_scalar(tmp,cst->val.scalar,cst->val.scalar);
  }
  else{
	elina_interval_set(tmp,cst->val.interval);
  }
  elina_interval_div(tmp,tmp,interval,discr);
  if(elina_scalar_equal(tmp->inf,tmp->sup)){
	elina_coeff_set_scalar(cst,tmp->inf);
  }
  else{
	elina_coeff_set_interval(cst,tmp);
  }
  elina_linexpr0_ForeachLinterm(expr,i,dim,coeff){
    if(coeff->discr==ELINA_COEFF_SCALAR){
	elina_interval_set_scalar(tmp,coeff->val.scalar,coeff->val.scalar);
    }
    else{
	elina_interval_set(tmp,coeff->val.interval);
    }
    elina_interval_div(tmp,tmp,interval,discr);
    if(elina_scalar_equal(tmp->inf,tmp->sup)){
	elina_coeff_set_scalar(coeff,tmp->inf);
    }
    else{
	elina_coeff_set_interval(coeff,tmp);
    }
    
  }
  elina_interval_free(tmp);
  return;
}



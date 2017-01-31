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


/* ************************************************************************* */
/* elina_linearize: generic functions for (quasi)linearisation ofinterval expressions */
/* ************************************************************************* */


#include "elina_linearize.h"


/* ********************************************************************** */
/* I. Evaluation of interval linear expressions */
/* ********************************************************************** */

/* Evaluate an ELINA interval linear expression */
bool elina_interval_eval_elina_linexpr0(elina_interval_t * itv, elina_linexpr0_t* expr, elina_interval_t** env, elina_scalar_discr_t discr)
{
  size_t i;
  elina_dim_t dim;
  elina_coeff_t* coeff;
  assert(env);
  elina_scalar_t * scalar;
  elina_interval_t * interval;
  elina_coeff_t * cst = &expr->cst;
  elina_interval_set_elina_coeff(itv,&expr->cst);
  elina_interval_t * tmp = elina_interval_alloc();
  elina_linexpr0_ForeachLinterm(expr,i,dim,coeff){
    bool eq = (coeff->discr==ELINA_COEFF_SCALAR);
    if (eq){
      scalar = coeff->val.scalar;
      if (elina_scalar_sgn(scalar)!=0){
	elina_interval_mul_scalar(tmp,
		      env[dim],
		      scalar,discr);
	elina_interval_add(itv, itv, tmp,discr);
      }
    }
    else {
      interval = coeff->val.interval;
      elina_interval_mul(tmp,env[dim],interval,discr);
      elina_interval_add(itv, itv, tmp,discr);
    }
    if (elina_interval_is_top(itv)){
      break;
    }
  }
  elina_interval_free(tmp);
  return true;
}




elina_interval_t*
eval_elina_linexpr0(elina_manager_t* man,
		 elina_abstract0_t* abs,
		 elina_linexpr0_t* expr,
		 elina_scalar_discr_t discr,
		 bool* pexact)
{
  bool exact;
  elina_dimension_t dim;
  elina_interval_t** aenv;
  elina_interval_t* r = elina_interval_alloc();
  if (pexact) *pexact = true;
  aenv = elina_abstract0_to_box(man,abs);
  if (!aenv) {
    elina_interval_set_bottom(r);
    return r;
  }
  dim = elina_abstract0_dimension(man,abs);
  exact = elina_interval_eval_elina_linexpr0(r,expr,aenv,discr);
  if (pexact) *pexact = exact;
  elina_interval_array_free(aenv,dim.intdim+dim.realdim);
  return r;
}



/*****************************
	Linearization
******************************/

/* Evaluate a constraint, composed of a constant (interval) expression */
char eval_elina_cstlincons0(elina_lincons0_t* cons)
{
  char res;
  elina_coeff_t * cst = &(cons->linexpr0->cst);
  bool equality = (cst->discr==ELINA_COEFF_SCALAR);
  elina_scalar_t *sup, *inf;
  assert (cons->linexpr0->size==0);
  if (!equality && elina_interval_is_bottom(cst->val.interval)){
    return 0;
  }
  switch (cons->constyp){
  case ELINA_CONS_EQ:
    if (equality){
      int sgn = elina_scalar_sgn(cst->val.scalar);
      res = (sgn==0 ? 1 : 0);
    }
    else {
      sup = cst->val.interval->sup;
      inf = cst->val.interval->inf;
      if (elina_scalar_sgn(sup)<0 || elina_scalar_sgn(inf)>0){
	res = 0;
      }
      else{
	res = 2;
      }
    }
    break;
  case ELINA_CONS_DISEQ:
    if(equality){
	int sgn = elina_scalar_sgn(cst->val.scalar);
	res = (sgn==0 ? 2 : 1);
    }
    else{
      sup = cst->val.interval->sup;
      inf = cst->val.interval->inf;
      if(elina_scalar_sgn(sup) < 0 || elina_scalar_sgn(inf)>0){
	 res = 1;
      }
      else{
	 res = 2;
      }
    }
    break;
  case ELINA_CONS_SUPEQ:
    if(equality){
	int sgn = elina_scalar_sgn(cst->val.scalar);
	if(sgn>=0){
		res = 1;
	}
	else{
		res = 0;
	}
    }
    else{
	sup = cst->val.interval->sup;
        inf = cst->val.interval->inf;
	if (elina_scalar_sgn(inf)>=0){
	      res = 1;
	}
	else if (elina_scalar_sgn(sup)<0){
	      res = 0;
	}
	else{
	      res = 2;
	}
    }
    break;
  case ELINA_CONS_SUP:
    if(equality){
	int sgn = elina_scalar_sgn(cst->val.scalar);
	if(sgn>0){
		res = 1;
	}
	else{
		res = 0;
	}
    }
    else{
	sup = cst->val.interval->sup;
        inf = cst->val.interval->inf;
    	if (elina_scalar_sgn(inf)>0){
      	    res = 1;
	}
        else if (elina_scalar_sgn(sup)<=0){
      	    res = 0;
	}
        else{
      	    res = 2;
	}
    }
    break;
  case ELINA_CONS_EQMOD:
    if (equality){
      elina_scalar_t * scalar = cst->val.scalar;
      if (elina_scalar_sgn(scalar)==0){
	res = 1;
      }
      else if (elina_scalar_cmp_int(cons->scalar,0)){
	res = 2;
      }
      else {
	elina_rat_t* rat = (elina_rat_t *)malloc(3*sizeof(elina_rat_t));
	elina_rat_set_elina_scalar(rat+1,scalar);
	elina_rat_set_elina_scalar(rat+2, cons->scalar);
	elina_rat_div(rat,rat+1,rat+2);
	if (rat->d==1){
	  res = 1;
	}
	else {
	  res = 2;
	}
	free(rat);
      }
    }
    else {
      res = 2;
    }
    break;
  default:
    abort();
  }
  return res;
}

static inline void elina_lincons0_swap(elina_lincons0_t* a, elina_lincons0_t* b){
	if (a!=b){ 
		elina_lincons0_t t=*a; 
		*a=*b; 
		*b=t; 
	} 
}

void elina_lincons0_array_reinit(elina_lincons0_array_t* array, size_t size)
{
  
  size_t i;
  if (size == array->size) return;
  if (size < array->size){
    for (i=size; i<array->size; i++){
      elina_lincons0_clear(&array->p[i]);
    }
    array->p = realloc(array->p,size*sizeof(elina_lincons0_t));
  }
  else { /* size > array->size */
    array->p = realloc(array->p,size*sizeof(elina_lincons0_t));
    for (i=array->size; i<size; i++){
      elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
      elina_scalar_t *scalar = elina_scalar_alloc();
      array->p[i].linexpr0 = linexpr0;
      array->p[i].scalar = scalar;
    }
  }
  array->size = size;
  return;
}

void elina_lincons0_set_bool(elina_lincons0_t* cons, bool value, elina_scalar_discr_t discr)
{
  /* constraint 0=0 if value, 1=0 otherwise */
  elina_linexpr0_realloc(cons->linexpr0,0);
  elina_coeff_t *cst = &cons->linexpr0->cst;
  elina_coeff_reinit(cst,ELINA_COEFF_SCALAR,discr);
  elina_scalar_set_to_int(cst->val.scalar,value ? 0 : 1, discr);
  cons->constyp = ELINA_CONS_EQ;
}

static bool 
elina_lincons0_is_useless_for_meet(elina_lincons0_t* cons)
{
  bool res = false;
  elina_coeff_t * cst = &(cons->linexpr0->cst);
  bool equality = (cst->discr==ELINA_COEFF_SCALAR);
  if (cons->linexpr0->size==0){
    if (eval_elina_cstlincons0(cons)!=0){
      res = true;
    }
  }
  else {
    if (!equality){
      bool sup = elina_scalar_infty(cst->val.interval->sup);
      switch (cons->constyp){
      case ELINA_CONS_EQ:
      case ELINA_CONS_DISEQ:
      case ELINA_CONS_EQMOD:
	{
	  bool inf = elina_scalar_infty(cst->val.interval->inf);
	  res = inf && sup;
	}
	break;
      case ELINA_CONS_SUPEQ:
      case ELINA_CONS_SUP:
	res = sup;
	break;
      default:
	break;
      }
    }
  }
  return res;
}

bool sat_elina_lincons0_is_false(elina_lincons0_t* cons)
{
  bool res = false;

  elina_coeff_t * cst = &(cons->linexpr0->cst);
  bool equality = (cst->discr==ELINA_COEFF_SCALAR);
  bool inf= false;
  if(equality){
	if(elina_scalar_infty(cst->val.scalar)){
		inf = true;
	}
  }
  else{
	if(elina_scalar_infty(cst->val.interval->inf)){
		inf = true;
	}
  }
  
  switch (cons->constyp){
  case ELINA_CONS_EQ:
  case ELINA_CONS_EQMOD:
    res = !equality;
    break;
  case ELINA_CONS_DISEQ:
    if(equality){
	res = inf;
    }
    else{
	res = inf && elina_scalar_infty(cst->val.interval->sup);
    }
    
    break;
  case ELINA_CONS_SUPEQ:
  case ELINA_CONS_SUP:
    res = inf;
    break;
  default:
    break;
  }
  if (!res && 
      cons->linexpr0->size==0 && 
      !eval_elina_cstlincons0(cons)){
    res = true;
  }
  return res;
}

char elina_lincons0_array_reduce(elina_lincons0_array_t* array, bool meet, elina_scalar_discr_t discr)
{
  char res;
  size_t i,size;

  res = 2;
  i = 0;
  size = array->size;
  while (i<size){
    elina_lincons0_t * cons = &array->p[i];
    if (cons->linexpr0->size==0){
      char sat = eval_elina_cstlincons0(cons);
      if (sat==1){
      elina_lincons0_array_reduce_remove:
	size --;
	elina_lincons0_swap(cons,&array->p[size]);
	continue;
      }
      else if (sat==0){
      elina_lincons0_array_reduce_false:
	elina_lincons0_array_reinit(array,1);
	elina_lincons0_set_bool(&array->p[0],false, discr);
	return 0;
      }
    }
    if (meet && elina_lincons0_is_useless_for_meet(&array->p[i]))
      goto elina_lincons0_array_reduce_remove;
    else if (!meet && sat_elina_lincons0_is_false(&array->p[i]))
      goto elina_lincons0_array_reduce_false;
    else {
      i++;
    }
  }
  elina_lincons0_array_reinit(array,size);
  if (size==0) 
    res = 1;
  else if (size==1 && array->p[0].linexpr0->size==0)
    res = eval_elina_cstlincons0(&array->p[0]);
  return res;
}

/* Transform sets of quasilinear constraint as follows:
   e.x + [a,b] >= 0 ==> e.x + b >= 0
   e.x + [a,b] > 0  ==> e.x + b > 0
   e.x + [a,b] = 0  ==> e.x + b >= 0 and e.x + a <= 0 added at the end
   e.x + [a,b] = 0 mod k ==> unchanged

   Also remove (some) trivially true constraints e.x + oo >= 0

*/
static void elina_lincons0_select_sup(elina_lincons0_t* cons)
{
  elina_coeff_t * cst = &(cons->linexpr0->cst);
  assert(cst->discr==ELINA_COEFF_INTERVAL);
  elina_scalar_t * sup = cst->val.interval->sup;
  elina_scalar_t *tmp = elina_scalar_alloc();
  elina_scalar_set(tmp,sup);
  elina_coeff_set_scalar(cst,tmp);
  elina_scalar_free(tmp);
}


static void elina_lincons0_select_inf(elina_lincons0_t* cons)
{
  elina_linexpr0_t * expr = cons->linexpr0;
  elina_coeff_t * cst = &(expr->cst);
  assert(cst->discr==ELINA_COEFF_INTERVAL);
  elina_scalar_t * inf = cst->val.interval->inf;
  elina_scalar_t *tmp = elina_scalar_alloc();
  elina_scalar_set(tmp,inf);
  elina_coeff_set_scalar(cst,tmp);
  elina_scalar_free(tmp);
  elina_linexpr0_neg(expr); 
}


void elina_lincons0_set(elina_lincons0_t * dst, elina_lincons0_t * src){
	if (dst!=src){ 
		if(dst->linexpr0){
			elina_linexpr0_free(dst->linexpr0);
		}
		//elina_linexpr0_clear(dst);
		dst->linexpr0 = elina_linexpr0_copy(src->linexpr0); 
		if(src->scalar){
			elina_scalar_set(dst->scalar,src->scalar);
		}
		else if(dst->scalar){
			elina_scalar_free(dst->scalar);
		} 
		dst->constyp = src->constyp; 
	}
}

void linearize_elina_lincons0_array(elina_lincons0_array_t* array, bool meet, elina_scalar_discr_t discr)
{
	
  size_t index,size,sizeorg;

  char res = elina_lincons0_array_reduce(array,meet, discr);
  if (res!=2) return;
 
  /* One now remove intervals when we can */
  sizeorg = array->size;
  size = sizeorg;
  for (index=0; index<sizeorg; index++){
     
    elina_lincons0_t* cons = &array->p[index];
    elina_linexpr0_t* expr = cons->linexpr0;
    elina_coeff_t * cst = &expr->cst;
    bool equality = (cst->discr==ELINA_COEFF_SCALAR);
    if (!equality){
      elina_scalar_t * isup = cst->val.interval->sup;
      elina_scalar_t * iinf = cst->val.interval->inf;
      bool sup = !elina_scalar_infty(isup);
	
      switch (cons->constyp){
      case ELINA_CONS_EQ:
	assert (meet); /* otherwise, already removed */
	{
	  bool inf = !elina_scalar_infty(iinf);
	  assert (inf || sup); /* otherwise, already removed */
	  if (inf && sup){
	    if (size>=array->size){
	      elina_lincons0_array_reinit(array,1+5*array->size/4);
	    }
	    /* Be cautious: cons and cst may be invalid now */
	    elina_lincons0_set(&array->p[size],&array->p[index]);
	    array->p[index].constyp = ELINA_CONS_SUPEQ;
	    array->p[size].constyp  = ELINA_CONS_SUPEQ;
	    elina_lincons0_select_sup(&array->p[index]);
	    elina_lincons0_select_inf(&array->p[size]);
	    size++;
	  }
	  else if (inf){
	    array->p[index].constyp = ELINA_CONS_SUPEQ;
	    elina_lincons0_select_inf(&array->p[index]);
	  }
	  else if (sup){
	    array->p[index].constyp = ELINA_CONS_SUPEQ;
	    elina_lincons0_select_sup(&array->p[index]);
	  }
	  else
	    assert(false);
	}
	break;
      case ELINA_CONS_SUPEQ:
      case ELINA_CONS_SUP:
	
	if (meet){
	  assert(sup);
	  elina_lincons0_select_sup(&array->p[index]);
	}
	else {
	  assert(!elina_scalar_infty(iinf));
	  elina_coeff_t * cst = &cons->linexpr0->cst;
	  elina_scalar_t * tinf = cst->val.interval->inf;
	   
	  elina_coeff_set_scalar(cst,tinf);
	}
	break;
      default:
	break;
      }
    }
  }
  elina_lincons0_array_reinit(array,size);
  
}

/* ********************************************************************** */
/* III. (Quasi)linearisation of interval linear expressions */
/* ********************************************************************** */

/* These functions quasilinearize in-place expressions and constraints.  They
   optimize (sets of) constraints when the parameter meet is true, by deducing
   things. If constraints are quasilinearized for testing satisfaction, meet
   should be set to false.
*/


/* Choose the middle of interval coefficient coeff for quasilinearisation */
/* Applies the following choice:

  if coeff=[-oo,+oo], choose 0
  if coeff=[-oo,x], choose x
  if coeff=[x,+oo], choose x
  if coeff = [inf,sup]
    if for_meet_inequality,
      (* always choose in favour of a finite sup bound in the constant
	 of the quasilinear expression *)
      if var=[-oo,a]
	choose inf,
	because it gives inf.var + [0,sup-inf].[-oo,a] = inf.var + [-oo,b]
      if var=[a,+oo]
	choose sup,
	because it gives sup.var + [inf-sup,0].[a,+oo] = sup.var + [-oo,b]
      if var=[a,b], choose middle
    else
      (* always choose in favour of at least a finite bound in the evaluation
	 of the quasilinear expression *)
      if var=[-oo,a]
	if inf >= 0, choose inf,
	  because inf.var + [0,sup-inf][-oo,a] evaluates to a finite sup bound
	if sup<=0, choose sup
	  because sup.var + [inf-sup,0][-oo,a] evaluates to a finite inf bound
	otherwise arbitrary choice (middle)
      if var=[a,+oo]
	if inf >= 0, choose inf,
	  because inf.var + [0,sup-inf][a,+oo] evaluates to a finite inf bound
	if sup <= 0, choose sup,
	  because sup.var + [inf-sup,0][a,+oo] evaluates to a finite sup bound
	otherwise arbitrary choice (middle)
      if var=[a,b], choose middle
*/

static void
elina_quasilinearize_choose_middle(elina_scalar_t * middle, /* the result */
				 elina_coeff_t *coeff,    /* the coefficient in which
						    middle is to be picked */
				 elina_interval_t *var,      /* the variable interval */
				 bool for_meet_inequality, /* is it for the linearisation of an inequality ? */
				elina_scalar_discr_t discr			     
							     
				 )
{
	
	
  bool equality = (coeff->discr==ELINA_COEFF_SCALAR);
	
  elina_scalar_t * inf, *sup, *scalar;
  if (!equality){
     inf = coeff->val.interval->inf;
     sup = coeff->val.interval->sup;
	
     if(elina_scalar_infty(inf)){
    	 if (elina_scalar_infty(sup)){
      	     elina_scalar_set_to_int(middle,0,discr);
         }
    	 else{
      	     elina_scalar_set(middle,sup);
	 }
     }
     else if(elina_scalar_infty(sup)){
	 elina_scalar_set(middle,inf);
     }
     else {
     	/* if coeff = [inf,sup] */
	if (for_meet_inequality){
	    if (elina_scalar_infty(var->inf))
		elina_scalar_set(middle,inf);
	    else if (elina_scalar_infty(var->sup))
		elina_scalar_set(middle,sup);
	    else /* Arbitrary choice: we take the middle */
		goto elina_quasilinearize_choose_middle_default;
	}
     	else {
         if (elina_scalar_infty(var->inf) ?
	     !elina_scalar_infty(var->sup) :
	     elina_scalar_infty(var->sup)){
		if (elina_scalar_sgn(inf)>=0)
		  elina_scalar_set(middle,inf);
		else if (elina_scalar_sgn(sup)<=0)
		  elina_scalar_set(middle,sup);
		else /* Arbitrary choice: we take the middle */
		  goto elina_quasilinearize_choose_middle_default;
      	  }
      	  else {
      		elina_quasilinearize_choose_middle_default:
			elina_scalar_add(middle,sup,inf,discr);
			elina_scalar_div_2(middle,middle);
      	  }
      }
    }
  }
  else{
	scalar = coeff->val.scalar;
	elina_scalar_set(middle,scalar);
  }
}
 


bool quasilinearize_elina_linexpr0(elina_linexpr0_t* expr, elina_interval_t** env, bool for_meet_inequality, elina_scalar_discr_t discr)
{
  
  size_t i,k,size;
  elina_dim_t dim=0;
  elina_coeff_t * coeff=NULL;
  bool peq;
  k = 0;
  i = 0;
  size = expr->size;
  while (i<size){
    elina_linterm_t * term;
    switch(expr->discr){
	case ELINA_LINEXPR_DENSE:
		dim = i;
		coeff = &expr->p.coeff[i];
		break;
	case ELINA_LINEXPR_SPARSE:
		term = &expr->p.linterm[i];
		dim = term->dim;
		coeff = &term->coeff;	
		break;
    }
	
    if(dim==ELINA_DIM_MAX){
	i++;
	continue;
    }
    peq = (coeff->discr==ELINA_COEFF_SCALAR);
    elina_interval_t * interval = env[dim];
    elina_scalar_t * inf = interval->inf;
    elina_scalar_t * sup = interval->sup;
    if (!elina_scalar_infty(inf) && elina_scalar_equal(inf,sup)){
      /* If a variable has a constant value, simplification */
      elina_coeff_mul_scalar(coeff,coeff,sup,discr);
      elina_coeff_add(&expr->cst,&expr->cst,coeff,discr);
      elina_coeff_t  tcoeff;
      elina_linterm_t tmp;
      switch(expr->discr){
	      case ELINA_LINEXPR_DENSE:
			tcoeff = expr->p.coeff[i];
			elina_coeff_reinit(&tcoeff,ELINA_COEFF_SCALAR,discr);
			elina_scalar_set_to_int(tcoeff.val.scalar,0,discr); 
			i++;
			break;
	      case ELINA_LINEXPR_SPARSE:
			size--;
			tmp = expr->p.linterm[i];
			expr->p.linterm[i] = expr->p.linterm[size];
			expr->p.linterm[size] = tmp;			
			break;
      }
    }
    else {
      if (peq == false){
	/* Compute the middle of the interval */
	elina_scalar_t * middle = elina_scalar_alloc();
	
	elina_quasilinearize_choose_middle(middle,
					 coeff,env[dim],for_meet_inequality,discr);
	
	/* Residue (interval-middle) */
	elina_coeff_t * tcoeff = elina_coeff_alloc(ELINA_COEFF_INTERVAL);
	elina_coeff_t * tcoeff2 = elina_coeff_alloc(ELINA_COEFF_INTERVAL);
	elina_coeff_sub_num(tcoeff,coeff,middle,discr);
	
	/* Multiplication of residue by variable range */
	elina_coeff_mul_interval(tcoeff2,tcoeff,env[dim],discr);
	/* Addition to the constant coefficient */
	elina_coeff_add(&expr->cst,&expr->cst,tcoeff2,discr);
        elina_coeff_t *cst = &expr->cst;
	if (cst->discr==ELINA_COEFF_INTERVAL && elina_scalar_infty(cst->val.interval->inf) && elina_scalar_infty(cst->val.interval->sup)){
	  k = 0;
	  break;
	}
	
	/* Addition of the linear term */
	if (elina_scalar_sgn(middle)!=0){
	  switch(expr->discr){
		case ELINA_LINEXPR_DENSE:
			elina_coeff_set_scalar(&expr->p.coeff[i],middle);
			break;
		case ELINA_LINEXPR_SPARSE:
			elina_coeff_set_scalar(&expr->p.linterm[k].coeff,middle);
			expr->p.linterm[k].dim = dim;
			break;
	  }
	  k++;
	}
      }
      else {
	if (k!=i){
	  switch(expr->discr){
		case ELINA_LINEXPR_DENSE:
			//elina_coeff_set(&expr->p.coeff[k],coeff);
			break;
		case ELINA_LINEXPR_SPARSE:
			elina_coeff_set(&expr->p.linterm[k].coeff,coeff);
			expr->p.linterm[k].dim = dim;
			break;
	  }
	}
	k++;
      }
      i++;
    }
  }
  
  if(expr->discr==ELINA_LINEXPR_SPARSE){
  	elina_linexpr0_reinit(expr,k);
  }
  
#if defined(NUM_FLOAT) || defined(NUM_DOUBLE) || defined(NUM_LONGDOUBLE) || defined(NUM_MPFR) || defined (NUM_NUMINT)
  return false;
#else
  return true;
#endif
  
}

bool quasilinearize_elina_lincons0(elina_lincons0_t* cons, elina_interval_t** env, bool meet, elina_scalar_discr_t discr)
{
  if (cons->linexpr0->size==0){
    /* constant expression */
    char sat = eval_elina_cstlincons0(cons);
    if (sat!=2){
      elina_lincons0_set_bool(cons,sat==1, discr);
    }
    return true;
  }
  else {
    return quasilinearize_elina_linexpr0(cons->linexpr0,env,
				      meet &&
				      (cons->constyp==ELINA_CONS_SUP ||
				       cons->constyp==ELINA_CONS_SUPEQ),discr);
  }
}


bool quasilinearize_elina_lincons0_array(elina_lincons0_array_t* array, elina_interval_t ** env, bool meet, elina_scalar_discr_t discr)
{
  size_t i,j,size;
  bool exact = true;
  
  elina_lincons0_array_reduce(array,meet, discr);
   
  size = array->size;
  for (i=0; i<size; i++){
    
    if (meet &&
	array->p[i].constyp == ELINA_CONS_EQ &&
	!elina_linexpr0_is_quasilinear((&array->p[i])->linexpr0)){
      /* Split an equality constraint into two inequalities if it is really
	 interval linear ==> better precision because quasilinearisation
	 choices differ between expr>=0 and expr<=0 */
      if (size>=array->size){
	elina_lincons0_array_reinit(array,1+(5*array->size)/4);
      }
	
      array->p[i].constyp = ELINA_CONS_SUPEQ;
      elina_lincons0_set(&array->p[size],&array->p[i]);
      elina_linexpr0_neg((&array->p[size])->linexpr0);
	
      size++;
    }
	
    exact = quasilinearize_elina_lincons0(&array->p[i],env,meet,discr) && exact;
	
    if ((&array->p[i])->linexpr0->size==0 &&
	eval_elina_cstlincons0(&array->p[i]) == 0){
      elina_lincons0_array_reinit(array,1);
      elina_lincons0_set_bool(&array->p[0],false, discr);
      return true;
    }
  }
  elina_lincons0_array_reinit(array,size);
 
  return exact;
}

void elina_lincons0_reduce_integer(elina_lincons0_t* cons, size_t intdim, elina_scalar_discr_t discr)
{
  elina_linexpr0_t* expr;
  size_t i;
  elina_coeff_t * coeff;
  elina_dim_t dim;
  bool* peq;
  bool integer;
  switch (cons->constyp){
  case ELINA_CONS_EQ:
  case ELINA_CONS_SUPEQ:
  case ELINA_CONS_SUP:
    break;
  default:
    return;
  }
  expr = cons->linexpr0;
  /* Tests if only integer variables are involved */
  if (!elina_linexpr0_is_integer(expr,intdim))
    return;
  /* Check that there are only scalar coefficients for dimensions */
  elina_linexpr0_ForeachLinterm(expr,i,dim,coeff) {
    if (coeff->discr!=ELINA_COEFF_SCALAR)
      return;
  }
    elina_rat_t * rat = (elina_rat_t *)malloc(sizeof(elina_rat_t));
    elina_rat_t * tmp = (elina_rat_t *)malloc(sizeof(elina_rat_t));
    /* compute lcm of denominators and gcd of numerators */
    rat->n = 0;
    rat->d = 1;
    elina_linexpr0_ForeachLinterm(expr,i,dim,coeff) {
      elina_rat_set_elina_scalar(tmp,coeff->val.scalar);
      rat->d = elina_int_lcm(rat->d,tmp->d);
      rat->n = elina_int_gcd(rat->n,tmp->n);
    }
    if (elina_int_sgn(rat->n)==0)
      return;
    elina_linexpr0_ForeachLinterm(expr,i,dim,coeff) {
      elina_rat_set_elina_scalar(tmp,coeff->val.scalar);
      tmp->n = tmp->n/rat->n;
      tmp->n = tmp->n*rat->d;
      tmp->n = tmp->n/tmp->d;
      tmp->d = 1;
	
      elina_scalar_set_elina_rat(coeff->val.scalar,tmp);
    }
    free(tmp);
    
    elina_rat_inv(rat,rat);
    elina_scalar_t * scalar = elina_scalar_alloc();
    elina_scalar_set_elina_rat(scalar,rat);
    elina_coeff_mul_scalar(&expr->cst,&expr->cst,scalar,discr);
    
    elina_scalar_free(scalar);

  /* Constrain bounds */
  elina_coeff_t *cst = &expr->cst;
  scalar = cst->discr==ELINA_COEFF_SCALAR ? cst->val.scalar : cst->val.interval->sup;
  elina_rat_set_elina_scalar(rat,scalar);
  if (cst->discr==ELINA_COEFF_SCALAR ||  !elina_scalar_infty(scalar)){
    if (cons->constyp==ELINA_CONS_SUP){
      if (rat->d==1){
	elina_rat_sub_uint(rat,rat,1);
      }
      else {
	rat->n = elina_int_fdiv_q(rat->n,rat->d);
        rat->d = 1;
      }
      cons->constyp = ELINA_CONS_SUPEQ;
    }
    else {
      rat->n = elina_int_fdiv_q(rat->n,rat->d);
      rat->d = 1;
    }
  }
  elina_scalar_set_elina_rat(scalar,rat);
  if (cons->constyp == ELINA_CONS_EQ){
    if(cst->discr==ELINA_COEFF_INTERVAL){
	    scalar = cst->val.interval->inf;
	    elina_rat_set_elina_scalar(rat,scalar);
	    if (!elina_scalar_infty(scalar)){
		rat->n = elina_int_fdiv_q(rat->n,rat->d);
        	rat->d = 1;
	    }
	    elina_scalar_set_elina_rat(scalar,rat);
	    if (elina_interval_is_bottom(cst->val.interval)){
	      elina_lincons0_set_bool(cons,false, discr);
	    }
    }
  }
  else {
    if (!elina_scalar_infty(scalar)){
     elina_coeff_reinit(&expr->cst,ELINA_COEFF_SCALAR,scalar->discr);
    }
  }
  free(rat);
}


char elina_lincons0_array_reduce_integer(elina_lincons0_array_t* array, size_t intdim, elina_scalar_discr_t discr)
{
  size_t i;
  for (i=0; i<array->size; i++){
    elina_lincons0_reduce_integer(&array->p[i],intdim,discr);
  }
  return elina_lincons0_array_reduce(array,true, discr);
}



/* ********************************************************************** */
/* II. Quasilinearization of interval linear expressions */
/* ********************************************************************** */

static
bool elina_quasilinearize_alloc(elina_manager_t* man, elina_abstract0_t* abs,
			  elina_interval_t** penv, elina_dimension_t* pdim)
{
  bool exact;
  assert(!elina_abstract0_is_bottom(man,abs));
  penv = elina_abstract0_to_box(man,abs);
  exact = man->result.flag_exact;
  *pdim = elina_abstract0_dimension(man,abs);
  return exact;
}

static inline
void elina_quasilinearize_free(elina_interval_t** env, elina_dimension_t dim)
{
  if(env)elina_interval_array_free(env,dim.intdim+dim.realdim);
}

/* Evaluate a interval linear expression on the abstract
   value such as to transform it into a quasilinear expression.

   discr allows to choose the type of scalars used for computations and for the
   result.  pexact is a pointer to a Boolean, which is set to true if all
   the conversions and computations were exact.
*/

elina_linexpr0_t* elina_quasilinearize_linexpr0(elina_manager_t* man,
				   void * abs,
				   elina_linexpr0_t* linexpr0,
				   bool* pexact, elina_scalar_discr_t discr)
{
  elina_dimension_t dim;
  elina_linexpr0_t* rlinexpr0;
  elina_interval_t* env;
  bool exact,exact2;

  exact = elina_quasilinearize_alloc(man,abs,&env,&dim);
  rlinexpr0 = elina_linexpr0_copy(linexpr0);
  exact = quasilinearize_elina_linexpr0(rlinexpr0,&env,false,discr)
    && exact;
  elina_quasilinearize_free(&env,dim);
  *pexact = exact;
  return rlinexpr0;
}

/* Same for elina_lincons0_t */
elina_lincons0_t elina_quasilinearize_lincons0(elina_manager_t* man,
				   void* abs, elina_lincons0_t* lincons0,
				   bool* pexact,  elina_scalar_discr_t discr, bool meet)
{
  elina_dimension_t dim;
  elina_lincons0_t rlincons0;
  elina_interval_t* env;
  bool exact;

  exact = elina_quasilinearize_alloc(man,abs,&env,&dim);
  rlincons0 = elina_lincons0_copy(lincons0);
  exact = quasilinearize_elina_lincons0(&rlincons0,&env,meet,discr)
    && exact;
  elina_quasilinearize_free(&env,dim);
  *pexact = exact;
  return rlincons0;
}

/* Same for arrays of elina_linexpr0_t */
elina_linexpr0_t** elina_quasilinearize_linexpr0_array(elina_manager_t* man,
					 void * abs, elina_linexpr0_t** texpr, size_t size,
					 bool* pexact,elina_scalar_discr_t discr)
{
  elina_dimension_t dim;
  elina_linexpr0_t** tab;
  elina_interval_t* env;
  bool exact,exact2;
  size_t i;

  exact = elina_quasilinearize_alloc(man,abs,&env,&dim);
  tab = (elina_linexpr0_t**)malloc(size*sizeof(elina_linexpr0_t*));
  for (i=0; i<size; i++){
    tab[i] = elina_linexpr0_copy(texpr[i]);
    exact = quasilinearize_elina_linexpr0(tab[i],&env,false,discr)
      && exact;
  }
  elina_quasilinearize_free(&env,dim);
  *pexact = exact;
  return tab;
}

/* Same for elina_lincons0_array_t */
elina_lincons0_array_t elina_quasilinearize_lincons0_array(elina_manager_t* man,
					 void * abs, elina_lincons0_array_t* array,
					 bool* pexact, elina_scalar_discr_t discr, bool linearize, bool meet)
{
  elina_interval_t* env;
  elina_dimension_t dim;
  bool exact;
  elina_lincons0_array_t res;
  size_t i;
  size_t size = array->size;
  res = elina_lincons0_array_make(size);
  for(i=0; i < size; i++){
	res.p[i] = elina_lincons0_copy(&array->p[i]);
  }
  exact = elina_quasilinearize_alloc(man,abs,&env,&dim);
  exact = quasilinearize_elina_lincons0_array(&res,&env,meet,discr);
  if (linearize) 
    linearize_elina_lincons0_array(&res,meet, discr);
  
  elina_quasilinearize_free(&env,dim);
  return res;
}

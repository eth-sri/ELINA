/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
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


#include <string.h>
#include <stdio.h>

#include "elina_box_internal.h"
#include "elina_box_representation.h"
#include "elina_box_constructor.h"
#include "elina_box_meetjoin.h"
#include "elina_scalar_arith.h"

/* ============================================================ */
/* Meet and Join */
/* ============================================================ */
elina_box_t* elina_box_meet(elina_manager_t* man, bool destructive, elina_box_t* a1, elina_box_t* a2)
{
    size_t i;
    size_t nbdims;
    elina_box_t* res;
    man->result.flag_best = true;
    man->result.flag_exact = false;
    res = destructive ? a1 : elina_box_alloc(a1->intdim,a1->realdim);
    if ((a1->inf==NULL && a1->sup==NULL) || (a2->inf==NULL && a2->sup==NULL)){
        elina_box_set_bottom(res);
        return res;
    }
    if (!destructive){
        elina_box_init(res);
    }
    
    nbdims = a1->intdim + a2->realdim;
    for (i=0; i<nbdims; i++){
        res->inf[i] = fmin(a1->inf[i],a2->inf[i]);
        res->sup[i] = fmin(a1->sup[i],a2->sup[i]);
    }
    return res;
}



elina_box_t* elina_box_join(elina_manager_t* man, bool destructive, elina_box_t* a1, elina_box_t* a2)
{
  size_t i;
  size_t nbdims;
  elina_box_t* res;
  man->result.flag_best = true;
  man->result.flag_exact = false;
  res = destructive ? a1 : elina_box_alloc(a1->intdim,a1->realdim);
  if ((a1->inf==NULL) && (a1->sup==NULL)){
    if ((a2->inf!=NULL) && (a2->sup!=NULL)){
      man->result.flag_exact = true;
      elina_box_set(res,a2);
    }
    return res;
  }
  else if (a2->inf==NULL && a2->sup==NULL){
    man->result.flag_exact = true;
    if (!destructive) elina_box_set(res,a1);
    return res;
  }
  man->result.flag_exact = false;
  if (!destructive){
    elina_box_init(res);
  }
  
  nbdims = a1->intdim + a2->realdim;
  for (i=0; i<nbdims; i++){
    res->inf[i] = fmax(a1->inf[i],a2->inf[i]);
    res->sup[i] = fmax(a1->sup[i],a2->sup[i]);
  }
   
  return res;
}

bool elina_double_interval_canonicalize(double *inf, double *sup, bool integer, elina_scalar_discr_t discr)
{
  bool exc;

  if (integer){
    *inf = floor(*inf);
    *sup = floor(*sup); 
  }
  if ((*inf==INFINITY) || *sup==INFINITY) return false;

  /* Check that it is not bottom */
  exc = false;
  if (*sup< -*inf )
    exc = true;
  return exc;
}



void elina_double_interval_mul(double *a_inf, double *a_sup, double b_inf, double b_sup, double c_inf, double c_sup){
	if(c_inf<=0){
		/* interval c is positive */
		if(b_inf<=0){
			/*interval b is positive*/
			if((b_inf==0) || (c_inf==0)){
				*a_inf = 0.0;
			}
			else{
				*a_inf = b_inf * -c_inf;
			}
			if((b_sup==0) || (c_sup==0)){
				*a_sup = 0.0;
			}
			else{
				*a_sup = b_sup * c_sup;
			}
		}
		else if(b_sup<=0){
			/* interval b is negative */
			if((c_sup==0) || (b_inf==0)){
				*a_inf = 0.0;
			}
			else{
			 	*a_inf = c_sup*b_inf;
			}
			if((c_inf==0) || (b_sup==0)){
				*a_sup = 0.0;
			}
			else{
			 	*a_sup = -c_inf*b_sup;
			}
		}
		else{
			/* there is 0 in between for b */
			if((c_sup==0) || (b_inf==0)){
				*a_inf = 0.0;
			}
			else{
				*a_inf = b_inf * c_sup;
			}
			if((c_sup==0) || (b_sup==0)){
				*a_sup = 0.0;
			}
			else{
				*a_sup = b_sup * c_sup;
			}
		}
	}
	else if(c_sup<=0){
		/* interval c is negative */
		if(b_inf<=0){
			/*interval b is positive*/
			if((b_sup==0) || (c_inf==0)){
				*a_inf = 0.0;
			}
			else{
				*a_inf = b_sup*c_inf;
			}
			if((b_inf==0) || (c_sup==0)){
				*a_sup = 0.0;
			}
			else{
				*a_sup = -b_inf*c_sup;
			}
		}
		else if(b_sup<=0){
			/* interval b is negative */
			if((b_sup==0) || (c_sup==0)){
				*a_inf = 0.0;
			}
			else{
				*a_inf = b_sup * -c_sup;
			}
			if((b_inf==0) || (c_inf == 0 )){
				*a_sup = 0.0;
			}
			else{
				*a_sup = b_inf * c_inf;
			} 
		}
		else{
			/* there is 0 in between for b */
			if((c_inf==0) || (b_sup==0)){
				*a_inf = 0.0;
			}
			else{
				*a_inf = b_sup*c_inf;
			}
			if((c_inf==0) || (b_inf==0)){
				*a_sup = 0.0;
			}
			else{
				*a_sup = b_inf*c_inf;
			} 
		}
	}
	else if(b_inf<=0){
		/* interval b is positive */
		if(c_inf<=0){
			/*interval c is positive */
			if((b_inf==0) || (c_inf==0)){
				*a_inf = 0.0;
			}
			else{
				*a_inf = -b_inf * c_inf;
			}
			if((b_sup==0) || (c_sup==0)){
				*a_sup = 0.0;
			}
			else{
				*a_sup = b_sup * c_sup;
			}
		}
		else if(c_sup<=0){
			/* interval c is negative */
			if((b_sup==0) || (c_inf == 0)){
				*a_inf = 0.0;
			}
			else{
				*a_inf = b_sup*c_inf;
			}
			if((b_inf==0) || (c_sup==0)){
				*a_sup = 0.0;
			}
			else{
				*a_sup = -b_inf*c_sup;
			}
		}
		else{
			/* there is 0 in between for c */
			if((b_sup==0) || (c_inf==0)){
				*a_inf = 0.0;
			}
			else{
				*a_inf = b_sup * c_inf;
			}
			if((b_sup==0) || (c_sup==0)){
				*a_sup = 0.0;
			}
			else{
				*a_sup = b_sup * c_sup;
			}
		}
	}
	else if(b_sup<=0){
		/* interval b is negative */
		if(c_inf <= 0){
			/* interval c is positive */
			if((b_inf==0) || (c_sup==0)){
				*a_inf = 0.0;
			}
			else{
				*a_inf = b_inf * c_sup;
			}
			if((b_sup==0) || (c_inf==0)){
				*a_sup = 0.0;
			}
			else{
				*a_sup = b_sup * -c_inf;
			}
		}
		else if(c_sup<=0){
			/* interval c is negative */
			if((b_sup==0) || (c_sup==0)){
				*a_inf = 0.0;
			}
			else{
				*a_inf = -b_sup * c_sup;
			}
			if((b_inf==0) || (c_inf==0)){
				*a_sup = 0.0;
			}
			else{
				*a_sup = b_inf * c_inf;
			} 
		}
		else{
			/* there is 0 in between for c */
			if((b_inf == 0) || (c_sup==0)){
				*a_inf = 0.0;
			}
			else{
				*a_inf = b_inf * c_sup;
			}
			if((b_inf==0) || (c_inf==0)){
				*a_sup = 0.0;
			}
			else{
				*a_sup = b_inf * c_inf;
			}
		}
	}
	else{
		/* there is 0 in between for both b and c */
		double tmp_inf1 = b_sup*c_inf;
		double tmp_sup1 = b_inf*c_inf;
		double tmp_inf2 = b_inf*c_sup;
		double tmp_sup2 = b_sup*c_sup;
		*a_inf = fmax(tmp_inf1, tmp_inf2);
		*a_sup = fmax(tmp_sup1, tmp_sup2);
	}
}




/* ====================================================================== */
/* Division */
/* ====================================================================== */


void elina_double_interval_div(double *a_inf, double *a_sup, double b_inf, double b_sup, double c_inf, double c_sup)
{
  if (c_inf<0){
    /* c is positive */
    if (b_inf<=0){
         /* b is positive */
         *a_inf = b_inf/c_sup;
	 *a_sup = b_sup/-c_inf;
    }
    else if (b_sup<=0){
       /* b is negative */
       *a_inf = -b_inf/c_inf;
       *a_sup = b_sup/c_sup;
    }
    else {
        /* 0 is in the middle of b: one divides b by c->inf */
        *a_inf = b_inf/-c_inf;
        *a_sup = b_sup/-c_inf;
    }
  }
  else if (c_sup<0){
    /* c is negative */
    if (b_inf<=0){
        /* b is positive */
	*a_sup = b_inf/c_inf;
        *a_inf = -b_sup/c_sup;
  	
     }
     else if (b_sup<=0){
       /* b is negative */
       *a_inf = b_sup/c_inf;
       *a_sup = -b_inf/c_sup;
     }
     else {
        /* 0 is in the middle of b: one cross-divide b by c->sup */
        *a_inf = b_sup/c_sup;
        *a_sup = b_inf/c_sup;
     }
  }
  else if ((b_inf==0) && (b_sup==0)){
    /* b is [0,0] */
    *a_inf = b_inf;
    *a_sup = b_sup;
  }
  else {
    *a_inf = INFINITY;
    *a_sup = INFINITY;
  }
}



bool elina_double_interval_eval_elina_linexpr0(double * itv_inf, double *itv_sup,  elina_linexpr0_t* expr, double* env_inf, double * env_sup, elina_scalar_discr_t discr)
{
  size_t i;
  elina_dim_t dim;
  elina_coeff_t* coeff;
  assert(env_inf&&env_sup);
  elina_scalar_t * scalar;
  elina_coeff_t * cst = &expr->cst;
  if(cst->discr==ELINA_COEFF_SCALAR){
	*itv_inf = -cst->val.scalar->val.dbl;
	*itv_sup = cst->val.scalar->val.dbl;
  }
  else{
	*itv_inf = -cst->val.interval->inf->val.dbl;
	*itv_sup = cst->val.interval->sup->val.dbl;
  }
  
  double tmp_inf = 0.0;
  double tmp_sup = 0.0;
  elina_linexpr0_ForeachLinterm(expr,i,dim,coeff){
    bool eq = (coeff->discr==ELINA_COEFF_SCALAR);
    if (eq){
      scalar = coeff->val.scalar;
      if (elina_scalar_sgn(scalar)>0){
	
	tmp_inf = env_inf[dim]*scalar->val.dbl;
	tmp_sup = env_sup[dim]*scalar->val.dbl;
	*itv_inf = *itv_inf + tmp_inf;
	*itv_sup = *itv_sup + tmp_sup;
      }
      else if(elina_scalar_sgn(scalar)<0){
	tmp_inf = env_sup[dim]*-scalar->val.dbl;
	tmp_sup = env_inf[dim]*-scalar->val.dbl;
        *itv_inf = *itv_inf + tmp_inf;
        *itv_sup = *itv_sup + tmp_sup;
      }
    }
    else {
      double inf = -coeff->val.interval->inf->val.dbl;
      double sup = coeff->val.interval->sup->val.dbl;
	//if((inf==sup) && (inf==0) ){
	//	tmp_inf = 0.0;
	//	tmp_sup = 0.0;
	//}
	//else{
      		elina_double_interval_mul(&tmp_inf,&tmp_sup,env_inf[dim],env_sup[dim],inf,sup);
	//}
      *itv_inf = *itv_inf + tmp_inf;
      *itv_sup = *itv_sup + tmp_sup;
     
    }
    if (*itv_inf==INFINITY && *itv_sup==INFINITY){
	
      break;
    }
  }
 
  return true;
}

static inline bool elina_double_scalar_sgn(double d){
	if(d > 0.0){
		return 1;
	}
	else if(d < 0.0){
		return -1;
	}
	else{
		return 0;
	}
	
}

static bool elina_double_boxize_lincons0(double* res_inf, double * res_sup,
			       bool* tchange,
			       elina_lincons0_t* cons,
			       double* env_inf, double * env_sup,
			       size_t intdim,
			       bool intervalonly, elina_scalar_discr_t discr)
{
  size_t i;
  elina_linexpr0_t* expr;
  bool change,globalchange;
  bool exc;
  assert(cons->constyp == ELINA_CONS_EQ ||
	 cons->constyp == ELINA_CONS_SUPEQ ||
	 cons->constyp == ELINA_CONS_SUP);

  expr = cons->linexpr0;
  globalchange = false;

  /* Iterates on coefficients */
  double itv_inf = 0.0;
  double itv_sup = 0.0;
  double val = 0.0;
  elina_scalar_t * scalar2 = elina_scalar_alloc();
  for (i=0; i<expr->size; i++){
    elina_coeff_t *coeff;
    elina_dim_t dim;
    if(expr->discr==ELINA_LINEXPR_SPARSE){
	elina_linterm_t *term = &expr->p.linterm[i];
        dim = term->dim;
        coeff = &term->coeff;
    }
    else{
	dim = i;
	coeff = &expr->p.coeff[i];
    }
    
    elina_coeff_t * tmp = elina_coeff_alloc(coeff->discr);
    /* 1. We decompose the expression e = ax+e' */
    elina_coeff_swap(tmp,coeff);
    double inf = 0.0;
    double sup = 0.0;
    double d;

    if (tmp->discr==ELINA_COEFF_SCALAR) {
        if (tmp->val.scalar->discr == ELINA_SCALAR_DOUBLE) {
            inf = -tmp->val.scalar->val.dbl;
            sup = tmp->val.scalar->val.dbl;
        } else {
            elina_double_set_scalar(&d,tmp->val.scalar,GMP_RNDD);
            inf = -d;
            elina_double_set_scalar(&d,tmp->val.scalar,GMP_RNDU);
            sup = d;
        }
    } else {
        if (tmp->val.interval->inf->discr == ELINA_SCALAR_DOUBLE) {
            inf = -tmp->val.interval->inf->val.dbl;
        } else {
            elina_double_set_scalar(&d,tmp->val.interval->inf,GMP_RNDD);
            inf = -d;
        }
        if (tmp->val.interval->sup->discr == ELINA_SCALAR_DOUBLE) {
            sup = tmp->val.interval->sup->val.dbl;
        } else {
            elina_double_set_scalar(&d,tmp->val.interval->inf,GMP_RNDU);
            sup = d;
        }
    }
    /* 2. evaluate e' */
    elina_double_interval_eval_elina_linexpr0(&itv_inf, & itv_sup, expr,env_inf, env_sup, discr);
    /* 3. Perform deduction from ax+e' = [-m,M]x + [-e,E] >= 0
	  we can deduce [-m,M]x + E >= 0
	  (for equality, [-m,M]x - e <= 0)
    */
    bool equality = tmp->discr==ELINA_COEFF_SCALAR;
    change = false;
    if (itv_inf!=INFINITY || itv_sup!=INFINITY){
      if (equality && !intervalonly){
	/* [-m,M]=[a,a] */
	
	int sgn = elina_double_scalar_sgn(sup);
	if (sgn!=0){
	  /* From ax+E >= 0, we deduce
	     (1) If a>0, then x>=-E/a
	     (2) If a<0, then x<=-E/a
	     From ax-e <= 0, we deduce
	     (3) If a>0, then x<=e/a
	     (4) If a<0, then x>=e/a
	  */
	  if (sgn>0 || cons->constyp == ELINA_CONS_EQ){
	    /* 1,4: inf bound
	       If we have a>0 (1), we compute E/a
	       If we have a<0 (4), we compute -e/a
	    */
		
	    if (sgn>0){
		val = itv_sup/sup;
	    }
	    else {
		val = itv_inf/inf;
	    }
	    if (dim<intdim && isfinite(val)){
	      if ((cons->constyp==ELINA_CONS_SUP) && (ceil(val) == val)){
		  val = val - 1;
	      }
	      else {
		  val = floor(val);
	      }
	    }
	    /* We update the interval */
	    if (val < res_inf[dim]){
	      change = true;
	      if (tchange) tchange[2*dim] = true;
	      res_inf[dim] = val;
	    }
	  }
	  if (sgn<0 || cons->constyp == ELINA_CONS_EQ){
	    /* 2,3: sup bound
	       If we have a<0 (2), we compute -E/a
	       If we have a>0 (3), we compute e/a
	    */
	    if (sgn<0){
	       val = itv_sup/inf;
	    }
	    else {
	       val = itv_inf/sup;
	    }
	    if (dim<intdim && isfinite(val)){
	      if ((cons->constyp==ELINA_CONS_SUP) &&
		  (ceil(val)==val)){
		   val = val - 1;
	      }
	      else {
		 val = floor(val);
	      }
	    }
	    /* We update the interval */
	    if (val < res_sup[dim]){
	      change = true;
	      if (tchange) tchange[2*dim+1] = true;
	      res_sup[dim] = val;
	    }
	  }
	}
      }
      else if (!equality){
	/* We have a real interval [-m,M] */
	/* Case [-m,M]x+E >= 0
	  (1) If -m>0, we rewrite [-m,M]x >= -E, and we have -1/m > 1/M > 0
	  (1A) If E<=0 <=> -E>=0 then x >= -E/M
	  (1B) If E>0  <=> -E<0  then x >= E/m
	  (2) If M<0, we rewrite [-M,m]x <= E, and we have 0 < 1/m < -1/M
	  (2A) If E<=0           then x <= E/m
	  (2B) If E>0            then x <= -E/M
	  Case [-m,M]x-e <= 0
	  (3) If -m>0, we rewrite [-m,M]x <= e, and we have 0 < 1/M < -1/m
	  (3A) If e<=0           then x <= e/M
	  (3B) If e>0            then x <= -e/m
	  (4) If M<0, we rewrite [-M,m]x >= -e, and we have -1/M > 1/m > 0
	  (4A) If e<=0 <=> -e>=0 then x >= -e/m
	  (4B) If e>0  <=> -e<0  then x >= e/M
	*/
	int sgnitv = inf<0 ? 1 : sup<0 ? -1 : 0;
	if (sgnitv != 0){
	  int sgne = elina_double_scalar_sgn(itv_inf);
	  int sgnE = elina_double_scalar_sgn(itv_sup);
	  if (sgnitv>0 || (cons->constyp==ELINA_CONS_EQ && sgnitv<0)){
	    /* 1,4: inf bound */
	    if (sgnitv>0){ /* 1 */
	      if (sgnE<=0){ /* 1A */
		/* We compute E/M */
		 val = itv_sup/sup;
		
	      } else { /* 1B */
		/* We compute -E/m */
		val = -(itv_sup/inf);
	      }
	    }
	    else { /* 4 */
	      if (sgne>=0){ /* 4A */
		/* We compute e/m */
		val = itv_inf/inf;
	      } else { /* 4B */
		/* We compute -e/M */
		val = -(itv_inf/sup);
	      }
	    }
	    if (dim<intdim && isfinite(val)){
	      if ((cons->constyp==ELINA_CONS_SUP) && (ceil(val)==val)){
		val = val - 1;
	      }
	      else {
		val = floor(val);
	      }
	    }
	    /* We update the interval */
	    if (val < res_inf[dim]){
	      change = true;
	      if (tchange) tchange[2*dim] = true;
	      res_inf[dim] = val;
	    }
	  }
	  if (sgnitv<0 || (cons->constyp==ELINA_CONS_EQ && sgnitv>0)){
	    /* 2,3: sup bound */
	    if (sgnitv<0){ /* 2 */
	      if (sgnE>=0){ /* 2B */
		/* We compute -E/M */
		val = -(itv_sup/sup);
	      } else { /* 2A */
		/* We compute E/m */
		val = itv_sup/inf;
	      }
	    }
	    else { /* 3 */
	      if (sgne<=0){ /* 3B */
		/* We compute -e/m */
		val = -(itv_inf/inf);
	      }
	      else { /* 3A */
		/* We compute e/M */
		val = itv_inf/sup;
	      }
	    }
	    if (dim<intdim && isfinite(val)){
	      if ((cons->constyp==ELINA_CONS_SUP) && (ceil(val)==val)){
		val = val - 1;
	      }
	      else {
		val = floor(val);
	      }
	    }
	    /* We update the interval */
	    if (val < res_sup[dim]){
	      change = true;
	      if (tchange) tchange[2*dim+1] = true;
	      res_sup[dim] = val;
	    }
	  }
	}
      }
    }
    elina_coeff_swap(tmp,coeff);
    elina_coeff_free(tmp);
    if (change){
      globalchange = true;
      exc = elina_double_interval_canonicalize(res_inf + dim,res_sup + dim,dim<intdim,discr);
      if (exc){
	res_inf[0] = -1;
	res_sup[0] = -1;
	return true;
      }
    }
  }
  if (expr->size==0 &&
      eval_elina_cstlincons0(cons)==0){
	res_inf[0] = -1;
	res_sup[0] = -1;
    globalchange = true;
  }
  return globalchange;
}



bool elina_double_boxize_lincons0_array(double * res_inf, double * res_sup, bool* tchange,
				      elina_lincons0_array_t* array,
				      double* env_inf, double * env_sup, size_t intdim,
				      size_t kmax,
				      bool intervalonly, elina_scalar_discr_t discr)
{
  size_t i,k;
  bool change,globalchange;

  if (kmax<1) kmax=1;
  if ((res_inf!=env_inf) && (res_sup!=env_sup)) kmax=1;

  globalchange = false;
  /* we possibly perform kmax passes */
  for (k=0;k<(size_t)kmax;k++){
    change = false;
    for (i=0; i<array->size; i++){
      if (array->p[i].constyp==ELINA_CONS_EQ ||
	  array->p[i].constyp==ELINA_CONS_SUPEQ ||
	  array->p[i].constyp==ELINA_CONS_SUP){
	
	change =
	  elina_double_boxize_lincons0(res_inf, res_sup, tchange,&array->p[i],env_inf,env_sup,intdim,intervalonly,discr)
	  ||
	  change
	  ;
	globalchange = globalchange || change;
	if (elina_double_interval_canonicalize(res_inf,res_sup,false,ELINA_SCALAR_DOUBLE)){
	  return true;
	}
      }
    }
    if (!change) break;
  }
  return globalchange;
}




/* ============================================================ */
/* Meet_lincons */
/* ============================================================ */

elina_box_t* elina_box_meet_lincons_array(elina_manager_t* man,
			      bool destructive,
			      elina_box_t* a,
			      elina_lincons0_array_t* array)
{
    //printf("start %p\n",man);
    //fflush(stdout);
  elina_box_t* res;
  size_t kmax;
  elina_lincons0_array_t tlincons;
  elina_box_internal_t* intern = (elina_box_internal_t*)man->internal;

  res = destructive ? a : elina_box_copy(man,a);
  if (a->inf==NULL && a->sup==NULL){
    man->result.flag_best = true;
    man->result.flag_exact = true;
  }
  else {
   
    man->result.flag_best = array->size==1;
    man->result.flag_exact = false;
    kmax = man->option.funopt[ELINA_FUNID_MEET_LINCONS_ARRAY].algorithm;
    if (kmax<1) kmax=2;
    tlincons = elina_lincons0_array_make(array->size);
    for(size_t i =0; i < array->size; i++){
	tlincons.p[i] = elina_lincons0_copy(&array->p[i]);
    }
    char tb = elina_lincons0_array_reduce_integer(&tlincons,a->intdim,ELINA_SCALAR_DOUBLE);
    if (tb==0){
      goto _elina_box_meet_lincons_array_bottom;
    }
    
    elina_double_boxize_lincons0_array(res->inf,res->sup,NULL,
			     &tlincons,res->inf,res->sup,a->intdim,kmax,false,ELINA_SCALAR_DOUBLE);
    //if (res->inf[0]<res->sup[0]){
    if(elina_double_interval_canonicalize(res->inf,res->sup,false,ELINA_SCALAR_DOUBLE)){
    _elina_box_meet_lincons_array_bottom:
      elina_box_set_bottom(res);
    }
  
    elina_lincons0_array_clear(&tlincons);
    
  }
    //printf("finish %p\n",man);
    //fflush(stdout);
  return res;
}


/* ============================================================ */
/* Widening */
/* ============================================================ */
elina_box_t* elina_box_widening(elina_manager_t* man,
                    elina_box_t* a1, elina_box_t* a2)
{
    size_t i;
    size_t nbdims;
    elina_box_t* res;
    
    man->result.flag_best = true;
    man->result.flag_exact = true;
    nbdims = a1->intdim+a1->realdim;
    if ((a1->inf==NULL) && (a1->sup==NULL)){
        return elina_box_copy(man,a2);
    }
    
    res = elina_box_copy(man,a1);
    for (i=0; i<nbdims; i++){
        if(a1->sup[i] < a2->sup[i]){
            res->sup[i] = INFINITY;
        }
        else{
            res->sup[i] = a1->sup[i];
        }
        if(a1->inf[i] < a2->inf[i]){
            res->inf[i] = INFINITY;
        }
        else{
            res->inf[i] = a1->inf[i];
        }
    }
    return res;
}



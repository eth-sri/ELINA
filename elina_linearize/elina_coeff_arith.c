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

#include "elina_coeff_arith.h"

void elina_coeff_mul_scalar(elina_coeff_t * dst, elina_coeff_t * src, elina_scalar_t * mul, elina_scalar_discr_t discr){
 if(dst!=src){
 	elina_coeff_reinit(dst,src->discr,discr);
 }
 elina_scalar_t *dscalar;
 elina_scalar_t *sscalar;
 switch(src->discr){
	case ELINA_COEFF_SCALAR:
		 dscalar = dst->val.scalar;
		 sscalar = src->val.scalar;
		 elina_scalar_mul(dscalar,sscalar,mul,discr);
		break;
	case ELINA_COEFF_INTERVAL:
		 elina_interval_mul_scalar(dst->val.interval,src->val.interval,mul,discr);
		 break;
  }
}



void elina_coeff_add(elina_coeff_t * op, elina_coeff_t * op1, elina_coeff_t * op2, elina_scalar_discr_t discr){
	elina_scalar_t *sup, *inf, *scalar;
 	elina_scalar_t *sup1, *inf1, *scalar1;
	elina_scalar_t *sup2, *inf2, *scalar2;
	
	if((op1->discr==ELINA_COEFF_INTERVAL) || (op2->discr==ELINA_COEFF_INTERVAL)){
		if(op!=op1){
			
			elina_coeff_reinit(op,ELINA_COEFF_INTERVAL,discr);
		}
		else{
			if(op->discr==ELINA_COEFF_SCALAR){
				elina_scalar_t * tmp = elina_scalar_alloc();
				elina_scalar_set(tmp,op->val.scalar);
				elina_coeff_reinit(op,ELINA_COEFF_INTERVAL,discr);
				elina_scalar_set(op->val.interval->inf,tmp);
				elina_scalar_set(op->val.interval->sup,tmp);
				elina_scalar_free(tmp);
			}
		}
	}
	else{
		if(op!=op1){
			elina_coeff_reinit(op,ELINA_COEFF_SCALAR,discr);
		}
		
	}
	if(op1->discr==ELINA_COEFF_INTERVAL){
		sup = op->val.interval->sup;
		sup1 = op1->val.interval->sup;
		inf = op->val.interval->inf;
		inf1 = op1->val.interval->inf;
		if(op2->discr==ELINA_COEFF_INTERVAL){
			sup2 = op2->val.interval->sup;
			inf2 = op2->val.interval->inf;			
		}
		else{
			sup2 = op2->val.scalar;
			inf2 = op2->val.scalar;
		}
		elina_scalar_add(sup,sup1,sup2,discr);
		elina_scalar_add(inf,inf1,inf2,discr);
	}
	else{
		sup1 = op1->val.scalar;
		inf1 = op1->val.scalar;
		if(op2->discr==ELINA_COEFF_INTERVAL){
			sup = op->val.interval->sup;
			inf = op->val.interval->inf;
			sup2 = op2->val.interval->sup;
			inf2 = op2->val.interval->inf;
			elina_scalar_add(sup,sup1,sup2,discr);
			elina_scalar_add(inf,inf1,inf2,discr);
		}
		else{
			
			scalar = op->val.scalar;
			scalar1 = op1->val.scalar;
			scalar2 = op2->val.scalar;
			elina_scalar_add(scalar,scalar1,scalar2,discr);
			
		}
	}
}


void elina_coeff_sub_num(elina_coeff_t * dst, elina_coeff_t * src, elina_scalar_t * sub, elina_scalar_discr_t discr){
	if(dst!=src){
		elina_coeff_reinit(dst,src->discr,discr);
	}
	elina_scalar_t *tmp;
	tmp = elina_scalar_alloc();
	elina_scalar_neg(tmp,sub);
	switch(src->discr){
		case ELINA_COEFF_SCALAR:
			elina_scalar_add(dst->val.scalar,src->val.scalar,tmp,discr);
			break;
		case ELINA_COEFF_INTERVAL:
			elina_scalar_add(dst->val.interval->inf,src->val.interval->inf,tmp,discr);
			elina_scalar_add(dst->val.interval->sup,src->val.interval->sup,tmp,discr);
			break;
	}
	elina_scalar_free(tmp);
}




void elina_coeff_mul_interval(elina_coeff_t * dst, elina_coeff_t *src, elina_interval_t * interval, elina_scalar_discr_t discr){
	if(dst->discr!=ELINA_COEFF_INTERVAL){
		if(dst!=src){
			elina_coeff_reinit(dst,ELINA_COEFF_INTERVAL,discr);
		}
		else{
			elina_scalar_t *tmp = elina_scalar_alloc();
			elina_scalar_set(tmp,dst->val.scalar);
			elina_coeff_reinit(dst,ELINA_COEFF_INTERVAL,discr);
			elina_interval_set_scalar(dst->val.interval,tmp,tmp);
			elina_scalar_free(tmp);
		}
	}
	elina_scalar_t *inf, *sup;
	switch(src->discr){
		case ELINA_COEFF_SCALAR:
			elina_interval_mul_scalar(dst->val.interval,interval,src->val.scalar,discr);
			break;
		case ELINA_COEFF_INTERVAL:
			elina_interval_mul(dst->val.interval,src->val.interval,interval,discr);
			break;
	}
	inf = dst->val.interval->inf;
	sup = dst->val.interval->sup;
	if(!elina_scalar_infty(inf) && elina_scalar_equal(inf,sup)){
		elina_scalar_t *tmp = elina_scalar_alloc();
		elina_scalar_set(tmp,sup);
		elina_coeff_reinit(dst,ELINA_COEFF_SCALAR,discr);
		elina_scalar_set(src->val.scalar,tmp);
		elina_scalar_free(tmp);
	}
}

void elina_interval_set_elina_coeff(elina_interval_t *interval, elina_coeff_t *coeff){
	elina_scalar_t * scalar;
	switch(coeff->discr){
		case ELINA_COEFF_SCALAR:
			scalar = coeff->val.scalar;
			elina_interval_set_scalar(interval,scalar,scalar);
			break;
		case ELINA_COEFF_INTERVAL:
			elina_scalar_set(interval->inf,coeff->val.interval->inf);
			elina_scalar_set(interval->sup,coeff->val.interval->sup);
			break;
	}
	
}


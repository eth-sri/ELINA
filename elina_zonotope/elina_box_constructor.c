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

#include <string.h>
#include <stdio.h>

#include "elina_box_constructor.h"

/* ********************************************************************** */
/* 1. Basic constructors */
/* ********************************************************************** */

/* We assume that dimensions [0..intdim-1] correspond to integer variables, and
   dimensions [intdim..intdim+realdim-1] to real variables */

/* Create a bottom (empty) value */
elina_box_t* elina_box_bottom(elina_manager_t* man, size_t intdim, size_t realdim)
{
  
  man->result.flag_best = true;
  man->result.flag_exact = true;
  elina_box_t *a = elina_box_alloc(intdim,realdim);
  return a;
}

/* Create a top (universe) value */
elina_box_t* elina_box_top(elina_manager_t* man, size_t intdim, size_t realdim)
{
  size_t i;
  
  elina_box_t* a = elina_box_alloc(intdim,realdim);
  elina_box_init(a);
  for(i=0;i<a->intdim+a->realdim; i++){
     a->inf[i] = INFINITY;
     a->sup[i] = INFINITY;
  }
  man->result.flag_best = true;
  man->result.flag_exact = true;
  return a;
}

/* Abstract an hypercube defined by the array of intervals
   of size intdim+realdim */
elina_box_t* elina_box_of_box(elina_manager_t* man,
		  size_t intdim, size_t realdim,
		  elina_interval_t** tinterval)
{
  size_t i;
  elina_box_internal_t* intern = elina_box_init_from_manager(man,ELINA_FUNID_OF_BOX);

  elina_box_t* a = elina_box_alloc(intdim,realdim);
  if (intdim+realdim!=0){
    elina_box_init(a);
    for(i=0;i<intdim+realdim; i++){
      //assert(tinterval[i]->sup->val.discr==ELINA_SCALAR_DOUBLE && );
      a->inf[i] = -tinterval[i]->inf->val.dbl;
      a->sup[i] = tinterval[i]->sup->val.dbl;
      
    }
  }
  man->result.flag_best = true;
  man->result.flag_exact = true;
  return a;
}

/* ********************************************************************** */
/* 2. Accessors */
/* ********************************************************************** */

elina_dimension_t elina_box_dimension(elina_manager_t* man, elina_box_t* a)
{ 
  elina_dimension_t res;
  res.intdim = a->intdim;
  res.realdim = a->realdim;
   // printf("intdim: %d realdim: %d\n",res.intdim,res.realdim);
  return res;
}

/* ********************************************************************** */
/* 3. Tests */
/* ********************************************************************** */

bool elina_box_is_bottom(elina_manager_t* man, elina_box_t* a)
{
  man->result.flag_best = true;
  man->result.flag_exact = true;
  return (a->inf==NULL) && (a->sup==NULL);
}

bool elina_box_is_top(elina_manager_t* man, elina_box_t* a)
{
  size_t i;
  bool res;
  size_t nbdims = a->intdim + a->realdim;

  man->result.flag_best = true;
  man->result.flag_exact = true;
  if ((a->inf==NULL) && (a->sup==NULL))
    return false;

  res = true;
  for (i=0;i<nbdims;i++){
    if (a->inf[i]!=INFINITY || a->sup[i]!=INFINITY){
      res = false;
      break;
    }
  }
  return res;
}

/* inclusion check */
bool elina_box_is_leq(elina_manager_t* man, elina_box_t* a, elina_box_t* b)
{
  size_t i;
  bool res;
  size_t nbdims;

  man->result.flag_best = true;
  man->result.flag_exact = true;
  nbdims = a->intdim + a->realdim;
  if ((a->inf==NULL) && (a->sup==NULL))
    return true;
  else if ((b->inf==NULL) && (b->sup==NULL))
    return false;

  res = true;
  for (i=0;i<nbdims;i++){
    if ((a->inf[i] > b->inf[i]) || (a->sup[i] > b->sup[i])){
      res = false;
      break;
    }
  }
  return res;
}

/* equality check */
bool elina_box_is_eq(elina_manager_t* man, elina_box_t* a, elina_box_t* b)
{
  size_t i;
  bool res;
  size_t nbdims;

  man->result.flag_best = true;
  man->result.flag_exact = true;
  nbdims = a->intdim + a->realdim;
  if ((a->inf==NULL) &&(a->sup==NULL))
    return (b->inf==NULL) && (b->sup==NULL);
  else if ((b->inf==NULL)&&(b->sup==NULL))
    return false;

  res = true;
  for (i=0;i<nbdims;i++){
    if ((a->inf[i]!=b->inf[i]) || (a->sup[i]!=b->sup[i])){
      res = false;
      break;
    }
  }
  return res;
}




/* ********************************************************************** */
/* 4 Extraction of properties */
/* ********************************************************************** */

elina_interval_t* elina_box_bound_dimension(elina_manager_t* man,
				   elina_box_t* a, elina_dim_t dim)
{
  bool exact;
  elina_box_internal_t* intern = elina_box_init_from_manager(man,ELINA_FUNID_BOUND_DIMENSION);
  elina_interval_t* interval = elina_interval_alloc();
  if ((a->inf==NULL) && (a->sup==NULL)){
    elina_interval_set_bottom(interval);
    exact = true;
  }
  else {
    elina_interval_set_double(interval,-a->inf[dim],a->sup[dim]);
    
  }
  man->result.flag_best = true;
  man->result.flag_exact = true;
  return interval;
}

/* Returns the interval taken by a linear expression
   over the abstract value */
elina_interval_t* elina_box_bound_linexpr(elina_manager_t* man,
				 elina_box_t* a, elina_linexpr0_t* expr)
{
  bool exact;
  elina_interval_t* interval = elina_interval_alloc();
  elina_box_internal_t* intern =  elina_box_init_from_manager(man,ELINA_FUNID_BOUND_LINEXPR);

  if ((a->inf==NULL) && (a->sup==NULL)){
    elina_interval_set_bottom(interval);
    exact = true;
  }
  else {
    //elina_interval_t **env = (elina_interval_t **)malloc((a->intdim+a->realdim)*sizeof(elina_interval_t*));
    size_t i;
    //for(i=0; i < a->intdim+a->realdim; i++){
	//env[i] = elina_interval_alloc();
	//elina_interval_set_double(env[i],a->inf[i],a->sup[i]);
    //}
    double itv_inf = 0;
    double itv_sup = 0;	
    exact = elina_double_interval_eval_elina_linexpr0(&itv_inf, &itv_sup,expr,a->inf,a->sup,ELINA_SCALAR_DOUBLE);
    elina_interval_set_double(interval,-itv_inf,itv_sup);
    //for(i=0; i < a->intdim+a->realdim; i++){
	//elina_interval_free(env[i]);
    //}
    //free(env);			
	
  }
  man->result.flag_best = true;
  man->result.flag_exact = exact;
  return interval;
}

/* convert to lincons array */
elina_lincons0_array_t elina_box_to_lincons_array(elina_manager_t* man, elina_box_t* a)
{
    size_t i;
    elina_lincons0_array_t array;
    
    size_t nbdims = a->intdim + a->realdim;
    
    man->result.flag_best = true;
    man->result.flag_exact = true;
    if ((a->inf==NULL) && (a->sup==NULL)){
        array = elina_lincons0_array_make(1);
        array.p[0] = elina_lincons0_make_unsat();
    }
    else if (nbdims==0){
        array = elina_lincons0_array_make(0);
    }
    else {
        size_t size;
        elina_linexpr0_t* expr;
        elina_scalar_t* scalar;
        
        size = 0;
        for (i=0;i<nbdims;i++){
            if (a->inf[i]!=INFINITY) size++;
            bool point = (a->inf[i]!=INFINITY) && (a->inf[i] == a->sup[i]);
            if (!point && (a->sup[i]!=INFINITY)) size++;
        }
        array = elina_lincons0_array_make(size);
        size = 0;
        for (i=0;i<nbdims;i++){
            bool point = false;
            if ((a->inf[i]!=INFINITY)){
                expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
                elina_coeff_set_scalar_int(&expr->p.linterm[0].coeff, 1);
                expr->p.linterm[0].dim = i;
                
                elina_coeff_reinit(&expr->cst,ELINA_COEFF_SCALAR,ELINA_SCALAR_DOUBLE);
                scalar = expr->cst.val.scalar;
                elina_scalar_set_double(scalar,-a->inf[i]);
                
                point = (a->inf[i]!=INFINITY) && (a->inf[i] == a->sup[i]);
                array.p[size].constyp = point ? ELINA_CONS_EQ : ELINA_CONS_SUPEQ;
                array.p[size].linexpr0 = expr;
                size++;
            }
            if (!point && (a->sup[i]!=INFINITY)){
                expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
                elina_coeff_set_scalar_int(&expr->p.linterm[0].coeff, -1);
                expr->p.linterm[0].dim = i;
                
                elina_coeff_reinit(&expr->cst,ELINA_COEFF_SCALAR,ELINA_SCALAR_DOUBLE);
                elina_scalar_set_double(expr->cst.val.scalar,a->sup[i]);
                
                array.p[size].constyp = ELINA_CONS_SUPEQ;
                array.p[size].linexpr0 = expr;
                size++;
            }
        }
    }
    return array;
}

elina_interval_t** elina_box_to_box(elina_manager_t* man, elina_box_t* a)
{
    size_t i;
    elina_interval_t** interval;
    size_t nbdims;
    elina_box_internal_t* intern = (elina_box_internal_t*)man->internal;
    
    man->result.flag_best = true;
    man->result.flag_exact = true;
    nbdims = a->intdim+a->realdim;
    if (nbdims==0){
        interval = NULL;
    }
    else {
        interval = elina_interval_array_alloc(nbdims);
        for (i=0; i<nbdims; i++){
            if ((a->inf==NULL) && (a->sup==NULL)){
                elina_interval_set_bottom(interval[i]);
            } else {
                elina_interval_set_double(interval[i],-a->inf[i],a->sup[i]);
            }
        }
    }
    return interval;
}


elina_box_t* elina_box_of_abstract0(elina_abstract0_t* a)
{
  return (elina_box_t*)a->value;
}

static inline elina_abstract0_t* abstract0_of_elina_box(elina_manager_t* man, elina_box_t* z)
{
  elina_abstract0_t* r = malloc(sizeof(elina_abstract0_t));
  assert(r);
  r->value = z;
  r->man = elina_manager_copy(man);
  return r;
}

elina_abstract0_t *relu_box(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t x){
	elina_box_t * input = elina_box_of_abstract0(abs);
	elina_box_t *a = destructive? input : elina_box_copy(man,input);
	a->sup[x] = fmax(0,a->sup[x]);
        a->inf[x] = fmin(0,a->inf[x]);
	return abstract0_of_elina_box(man,a);
}

elina_abstract0_t * relu_box_layerwise(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t start_offset, elina_dim_t num_dim){
	elina_dim_t i;
	//printf("ReLU input\n");
	//elina_abstract0_fprint(stdout,man,abs,NULL);
	//fflush(stdout);
	elina_dim_t end = start_offset + num_dim;
	elina_abstract0_t *res = destructive? abs : elina_abstract0_copy(man,abs);
	for(i=start_offset; i < end; i++){
		res= relu_box(man,true,res,i);
	}
	//printf("ReLU output\n");
	//elina_abstract0_fprint(stdout,man,res,NULL);
	//fflush(stdout);
	return res;
}

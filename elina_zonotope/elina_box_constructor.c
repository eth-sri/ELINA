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
    elina_interval_set_top(a->p[i]);
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

  elina_box_t* a = elina_box_alloc(intdim,realdim);
  if (intdim+realdim!=0){
    elina_box_init(a);
    for(i=0;i<intdim+realdim; i++){
      elina_interval_set(a->p[i], tinterval[i]);
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
  return a->p == NULL;
}

bool elina_box_is_top(elina_manager_t* man, elina_box_t* a)
{
  size_t i;
  bool res;
  size_t nbdims = a->intdim + a->realdim;

  man->result.flag_best = true;
  man->result.flag_exact = true;
  if (a->p == NULL)
    return false;

  res = true;
  for (i=0;i<nbdims;i++){
    if (!elina_interval_is_top(a->p[i])) {
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
  if (a->p == NULL)
    return true;
  else if (b->p == NULL)
    return false;

  res = true;
  for (i=0;i<nbdims;i++){
    if (!elina_interval_is_leq(a->p[i], b->p[i])) {
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
  if (a->p == NULL)
    return b->p == NULL;
  else if (b->p == NULL)
    return false;

  res = true;
  for (i=0;i<nbdims;i++){
    if (!elina_interval_equal(a->p[i], b->p[i])) {
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
  elina_interval_t* interval = elina_interval_alloc();
  if (a->p == NULL) {
    elina_interval_set_bottom(interval);
    // exact = true;
  } else {
    elina_interval_set(interval, a->p[dim]);
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

  if (a->p == NULL) {
    elina_interval_set_bottom(interval);
    exact = true;
  } else {
    exact = elina_interval_eval_elina_linexpr0(interval, expr, a->p,
                                               ELINA_SCALAR_DOUBLE);
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
    if (a->p == NULL) {
      array = elina_lincons0_array_make(1);
      array.p[0] = elina_lincons0_make_unsat();
    } else if (nbdims == 0) {
      array = elina_lincons0_array_make(0);
    } else {
      size_t size;
      elina_linexpr0_t *expr;
      elina_scalar_t *scalar;

      size = 0;
      for (i = 0; i < nbdims; i++) {
        if (!elina_scalar_infty(a->p[i]->inf))
          size++;
        bool point = !elina_scalar_infty(a->p[i]->inf) &&
                     elina_scalar_equal(a->p[i]->inf, a->p[i]->sup);
        if (!point && !elina_scalar_infty(a->p[i]->sup))
          size++;
      }
      array = elina_lincons0_array_make(size);
      size = 0;
      for (i = 0; i < nbdims; i++) {
        bool point = false;
        if (!elina_scalar_infty(a->p[i]->inf)) {
          expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 1);
          elina_coeff_set_scalar_int(&expr->p.linterm[0].coeff, 1);
          expr->p.linterm[0].dim = i;

          elina_coeff_reinit(&expr->cst, ELINA_COEFF_SCALAR,
                             ELINA_SCALAR_DOUBLE);
          scalar = expr->cst.val.scalar;
          elina_scalar_set(scalar, a->p[i]->inf);

          point = !elina_scalar_infty(a->p[i]->inf) &&
                  elina_scalar_equal(a->p[i]->inf, a->p[i]->sup);
          array.p[size].constyp = point ? ELINA_CONS_EQ : ELINA_CONS_SUPEQ;
          array.p[size].linexpr0 = expr;
          size++;
        }
        if (!point && !elina_scalar_infty(a->p[i]->sup)) {
          expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 1);
          elina_coeff_set_scalar_int(&expr->p.linterm[0].coeff, -1);
          expr->p.linterm[0].dim = i;

          elina_coeff_reinit(&expr->cst, ELINA_COEFF_SCALAR,
                             ELINA_SCALAR_DOUBLE);
          elina_scalar_set(expr->cst.val.scalar, a->p[i]->sup);

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
    
    man->result.flag_best = true;
    man->result.flag_exact = true;
    nbdims = a->intdim+a->realdim;
    if (nbdims==0){
        interval = NULL;
    }
    else {
        interval = elina_interval_array_alloc(nbdims);
        for (i=0; i<nbdims; i++){
          if (a->p == NULL) {
            elina_interval_set_bottom(interval[i]);
          } else {
            elina_interval_set(interval[i], a->p[i]);
          }
        }
    }
    return interval;
}


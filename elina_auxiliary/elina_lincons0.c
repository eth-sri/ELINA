/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
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

/* ************************************************************************* */
/* elina_lincons0.c: linear constraints and arrays */
/* ************************************************************************* */

#include "elina_lincons0.h"

/* ********************************************************************** */
/* I. Linear constraints */
/* ********************************************************************** */

void elina_lincons0_print(elina_lincons0_t* cons,
		       char** name_of_dim)
{ elina_lincons0_fprint(stdout,cons,name_of_dim); }
void elina_lincons0_fprint(FILE* stream, elina_lincons0_t* cons, char** name_of_dim)
{
  elina_linexpr0_fprint(stream,cons->linexpr0,name_of_dim);
  fprintf(stream,
          cons->constyp == ELINA_CONS_EQ || cons->constyp == ELINA_CONS_EQMOD
              ?
              //" = 0" :
              " 0 "
              : (cons->constyp == ELINA_CONS_SUPEQ
                     ?
                     // " >= 0" :
                     " 1 "
                     : (cons->constyp == ELINA_CONS_SUP
                            ?
                            //" > 0" :
                            " 2 "
                            : (cons->constyp == ELINA_CONS_DISEQ
                                   ?
                                   //" != 0" :
                                   " 3 "
                                   : "\"ERROR in elina_lincons0_fprint\""))));
  if (cons->constyp == ELINA_CONS_EQMOD){
    assert(cons->scalar!=NULL);
    fprintf(stream," mod ");
    elina_scalar_fprint(stream,cons->scalar);
  }
}

elina_lincons0_t elina_lincons0_make_unsat()
{
  elina_linexpr0_t* expr;

  expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
  elina_coeff_set_scalar_double(&expr->cst,-1.0);
  return elina_lincons0_make(ELINA_CONS_SUPEQ,expr,NULL);
}

bool elina_lincons0_is_unsat(elina_lincons0_t* cons)
{
  size_t i,nbcoeffs;
  elina_dim_t dim;
  elina_coeff_t* coeff;
  int sgn;
  elina_linexpr0_t* expr = cons->linexpr0;

  nbcoeffs = 0;
  elina_linexpr0_ForeachLinterm(expr,i,dim,coeff){
    if (!elina_coeff_zero(coeff)){
      nbcoeffs++;
      if (nbcoeffs>0) break;
    }
  }
  if (nbcoeffs==0){
    switch (expr->cst.discr){
    case ELINA_COEFF_SCALAR:
      sgn = elina_scalar_sgn(expr->cst.val.scalar);
      switch(cons->constyp){
      case ELINA_CONS_EQ:
      case ELINA_CONS_EQMOD:
	return (sgn!=0);
      case ELINA_CONS_DISEQ:
	return (sgn==0);
      case ELINA_CONS_SUPEQ:
	return (sgn<0);
      case ELINA_CONS_SUP:
	return (sgn<=0);
      }
    case ELINA_COEFF_INTERVAL:
      sgn = elina_scalar_sgn(expr->cst.val.interval->sup);
      switch(cons->constyp){
      case ELINA_CONS_EQ:
      case ELINA_CONS_EQMOD:
	return
	  sgn < 0 ||
	  elina_scalar_sgn(expr->cst.val.interval->inf)>0;
      case ELINA_CONS_DISEQ:
	return
	  sgn>=0 &&
	  elina_scalar_sgn(expr->cst.val.interval->inf)<=0;
      case ELINA_CONS_SUPEQ:
	return sgn<0;
      case ELINA_CONS_SUP:
	return (sgn<=0);
      }
    default:
      abort();
    }
  }
  else
    return false;
}

bool elina_lincons0_is_sat(elina_lincons0_t* cons)
{
  size_t i,nbcoeffs;
  elina_dim_t dim;
  elina_coeff_t* coeff;
  elina_linexpr0_t* expr = cons->linexpr0;

  nbcoeffs = 0;
  elina_linexpr0_ForeachLinterm(expr,i,dim,coeff){
    if (!elina_coeff_zero(coeff)){
      nbcoeffs++;
      if (nbcoeffs>0) break;
    }
  }
  if (nbcoeffs==0){
    switch (expr->cst.discr){
    case ELINA_COEFF_SCALAR:
      switch(cons->constyp){
      case ELINA_CONS_EQ:
      case ELINA_CONS_EQMOD:
	return elina_scalar_sgn(expr->cst.val.scalar)==0;
      case ELINA_CONS_DISEQ:
	return elina_scalar_sgn(expr->cst.val.scalar)!=0;
      case ELINA_CONS_SUPEQ:
	return elina_scalar_sgn(expr->cst.val.scalar)>=0;
      case ELINA_CONS_SUP:
	return elina_scalar_sgn(expr->cst.val.scalar)>0;
      }
    case ELINA_COEFF_INTERVAL:
      if (elina_interval_is_bottom(expr->cst.val.interval)) 
        return true;
      switch(cons->constyp){
      case ELINA_CONS_EQ:
      case ELINA_CONS_EQMOD:
	return
          elina_scalar_sgn(expr->cst.val.interval->inf)==0 &&
          elina_scalar_sgn(expr->cst.val.interval->sup)==0;
      case ELINA_CONS_DISEQ:
        return
          elina_scalar_sgn(expr->cst.val.interval->inf)<0 ||
          elina_scalar_sgn(expr->cst.val.interval->sup)>0;
      case ELINA_CONS_SUPEQ:
	return 
          elina_scalar_sgn(expr->cst.val.interval->inf)>=0;
      case ELINA_CONS_SUP:
	return 
          elina_scalar_sgn(expr->cst.val.interval->inf)>0;
      }
    default:
      abort();
    }
  }
  else
    return false;
}

/* ********************************************************************** */
/* II. Array of linear constraints */
/* ********************************************************************** */

elina_lincons0_array_t elina_lincons0_array_make(size_t size)
{
  elina_lincons0_array_t array;
  size_t i;
  array.size = size;
  array.p = (size==0) ? NULL : (elina_lincons0_t*)malloc(size*sizeof(elina_lincons0_t));
  for (i=0; i<size; i++){
    array.p[i].linexpr0 = NULL;
    array.p[i].scalar = NULL;
  }
  return array;
}
void elina_lincons0_array_resize(elina_lincons0_array_t* array, size_t size)
{
  size_t i;
  for (i=size; i<array->size; i++){
    elina_lincons0_clear(&array->p[i]);
  }
  array->p = (elina_lincons0_t*)realloc(array->p,size*sizeof(elina_lincons0_t));
   for (i=array->size; i<size; i++){
    array->p[i].linexpr0 = NULL;
    array->p[i].scalar = NULL;
  }
   array->size = size;
}

void elina_lincons0_array_clear(elina_lincons0_array_t* array)
{
  size_t i;

  if (array->p!=NULL){
    for (i=0; i<array->size; i++)
      elina_lincons0_clear(&array->p[i]);
    free(array->p);
    array->p=NULL;
  }
}

void elina_lincons0_array_print(elina_lincons0_array_t* array,
			     char** name_of_dim)
{ elina_lincons0_array_fprint(stdout,array,name_of_dim); }
void elina_lincons0_array_fprint(FILE* stream,
			      elina_lincons0_array_t* array,
			      char** name_of_dim)
{
  size_t i;
  fprintf(stream,"%lu\n",
	    (unsigned long)array->size);
  // if (array->size==0){
  // fprintf(stream,"empty array of constraints\n");
  //} else {
  // fprintf(stream,"array of constraints of size %lu\n",
  //    (unsigned long)array->size);
  for (i = 0; i < array->size; i++) {
    // fprintf(stream,"%2lu: ",(unsigned long)i);
    elina_lincons0_fprint(stream, &array->p[i], name_of_dim);
    fprintf(stream, "\n");
    }
    // }
}

elina_linexpr_type_t elina_lincons0_array_type(elina_lincons0_array_t* array)
{
  size_t i;
  elina_linexpr_type_t res = ELINA_LINEXPR_LINEAR;
  for (i=0; i<array->size; i++){
    elina_linexpr_type_t type = elina_linexpr0_type(array->p[i].linexpr0);
    if (type < res)
      type = res;
    if (res==ELINA_LINEXPR_INTLINEAR)
      break;
  }
  return res;
}


bool elina_lincons0_array_is_quasilinear(elina_lincons0_array_t* array)
{
  size_t i;
  bool res = true;
  for (i=0; i<array->size; i++){
    if (!elina_linexpr0_is_quasilinear(array->p[i].linexpr0)){
      res = false;
      break;
    }
  }
  return res;
}
bool elina_lincons0_array_is_linear(elina_lincons0_array_t* array)
{
  size_t i;
  bool res = true;
  for (i=0; i<array->size; i++){
    if (!elina_linexpr0_is_linear(array->p[i].linexpr0)){
      res = false;
      break;
    }
  }
  return res;
}

/* ====================================================================== */
/* II.1 Change of dimensions and permutations */
/* ====================================================================== */
void elina_lincons0_array_add_dimensions_with(elina_lincons0_array_t* array,
					   elina_dimchange_t* dimchange)
{
  size_t i;
  for(i=0; i<array->size; i++){
    elina_linexpr0_t* expr = array->p[i].linexpr0;
    if (expr) elina_linexpr0_add_dimensions_with(expr,dimchange);
  }
}
elina_lincons0_array_t elina_lincons0_array_add_dimensions(elina_lincons0_array_t* array,
						     elina_dimchange_t* dimchange)
{
  size_t i;
  elina_lincons0_array_t narray;

  narray = elina_lincons0_array_make(array->size);
  for(i=0; i<array->size; i++){
    narray.p[i] =  elina_lincons0_add_dimensions(&array->p[i],dimchange);
  }
  return narray;
}

void elina_lincons0_array_permute_dimensions_with(elina_lincons0_array_t* array,
					       elina_dimperm_t* perm)
{
  size_t i;
  for(i=0; i<array->size; i++){
    elina_linexpr0_t* expr = array->p[i].linexpr0;
    if (expr) elina_linexpr0_permute_dimensions_with(expr,perm);
  }
}
elina_lincons0_array_t elina_lincons0_array_permute_dimensions(elina_lincons0_array_t* array,
							 elina_dimperm_t* perm)
{
  size_t i;
  elina_lincons0_array_t narray;

  narray = elina_lincons0_array_make(array->size);
  for(i=0; i<array->size; i++){
    narray.p[i] =  elina_lincons0_permute_dimensions(&array->p[i],perm);
  }
  return narray;
}


/* ********************************************************************** */
/* III. Auxilliary functions */
/* ********************************************************************** */

 elina_lincons0_t elina_lincons0_make(elina_constyp_t constyp, elina_linexpr0_t* linexpr, elina_scalar_t* scalar)
{
  elina_lincons0_t cons;
  cons.constyp = constyp;
  cons.linexpr0 = linexpr;
  cons.scalar = scalar;
  return cons;
}
 elina_lincons0_t elina_lincons0_copy(elina_lincons0_t* cons)
{
  return elina_lincons0_make(cons->constyp, 
			  cons->linexpr0 ? elina_linexpr0_copy(cons->linexpr0) : NULL,
			  cons->scalar ? elina_scalar_alloc_set(cons->scalar) : NULL);
}
 void elina_lincons0_clear(elina_lincons0_t* lincons)
{
  if (lincons->linexpr0){
    elina_linexpr0_free(lincons->linexpr0);
  }
  lincons->linexpr0 = NULL;
  if (lincons->scalar){
        elina_scalar_free(lincons->scalar);
  }
  lincons->scalar = NULL;
}


void elina_lincons0_add_dimensions_with(elina_lincons0_t* cons,
				     elina_dimchange_t* dimchange)
{ elina_linexpr0_add_dimensions_with(cons->linexpr0,dimchange); }

elina_lincons0_t elina_lincons0_add_dimensions(elina_lincons0_t* cons,
					 elina_dimchange_t* dimchange)
{
  return elina_lincons0_make(cons->constyp,
			  elina_linexpr0_add_dimensions(cons->linexpr0,dimchange),
			  cons->scalar ? elina_scalar_alloc_set(cons->scalar) : NULL);
}

void elina_lincons0_permute_dimensions_with(elina_lincons0_t* cons,
					 elina_dimperm_t* perm)
{ elina_linexpr0_permute_dimensions_with(cons->linexpr0,perm); }

elina_lincons0_t elina_lincons0_permute_dimensions(elina_lincons0_t* cons,
					     elina_dimperm_t* perm)
{
  return elina_lincons0_make(cons->constyp,
			  elina_linexpr0_permute_dimensions(cons->linexpr0,perm),
			  cons->scalar ? elina_scalar_alloc_set(cons->scalar) : NULL);
}

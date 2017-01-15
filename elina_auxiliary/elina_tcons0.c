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
/* elina_tcons0.c: tree constraints and arrays */
/* ************************************************************************* */


#include "elina_tcons0.h"

/* ********************************************************************** */
/* I. Linear constraints */
/* ********************************************************************** */

void elina_tcons0_fprint(FILE* stream, elina_tcons0_t* cons, char** name_of_dim)
{
  elina_texpr0_fprint(stream,cons->texpr0,name_of_dim);
  fprintf(stream,
	  cons->constyp == ELINA_CONS_EQ || cons->constyp == ELINA_CONS_EQMOD ?
	  " = 0" :
	  ( cons->constyp == ELINA_CONS_SUPEQ ?
	    " >= 0" :
	    (cons->constyp == ELINA_CONS_SUP ?
	     " > 0" :
	     (cons->constyp == ELINA_CONS_DISEQ ?
	      " != 0" :
	      "\"ERROR in elina_tcons0_fprint\""))));
  if (cons->constyp == ELINA_CONS_EQMOD){
    assert(cons->scalar!=NULL);
    fprintf(stream," mod ");
    elina_scalar_fprint(stream,cons->scalar);
  }
}

elina_tcons0_t elina_tcons0_make_unsat()
{
  elina_texpr0_t* expr;

  expr = elina_texpr0_cst_scalar_int(-1);
  return elina_tcons0_make(ELINA_CONS_SUPEQ,expr,NULL);
}

/* ********************************************************************** */
/* II. Array of linear constraints */
/* ********************************************************************** */

elina_tcons0_array_t elina_tcons0_array_make(size_t size)
{
  elina_tcons0_array_t array;
  size_t i;
  array.size = size;
  array.p = (size==0) ? NULL : (elina_tcons0_t*)malloc(size*sizeof(elina_tcons0_t));
  for (i=0; i<size; i++){
    array.p[i].texpr0 = NULL;
    array.p[i].scalar = NULL;
  }
  return array;
}
void elina_tcons0_array_resize(elina_tcons0_array_t* array, size_t size)
{
  size_t i;
  for (i=size; i<array->size; i++){
    elina_tcons0_clear(&array->p[i]);
  }
  array->p = (elina_tcons0_t*)realloc(array->p,size*sizeof(elina_tcons0_t));
   for (i=array->size; i<size; i++){
    array->p[i].texpr0 = NULL;
    array->p[i].scalar = NULL;
  }
   array->size = size;
}

void elina_tcons0_array_clear(elina_tcons0_array_t* array)
{
  size_t i;

  if (array->p!=NULL){
    for (i=0; i<array->size; i++)
      elina_tcons0_clear(&array->p[i]);
    free(array->p);
    array->p=NULL;
  }
}

void elina_tcons0_array_fprint(FILE* stream,
			      elina_tcons0_array_t* array,
			      char** name_of_dim)
{
  size_t i;

  if (array->size==0){
    fprintf(stream,"empty array of constraints\n");
  } else {
    fprintf(stream,"array of constraints of size %lu\n",
	    (unsigned long)array->size);
    for (i=0; i<array->size; i++){
      fprintf(stream,"%2lu: ",(unsigned long)i);
      elina_tcons0_fprint(stream,&array->p[i],name_of_dim);
      fprintf(stream,"\n");
    }
  }
}

bool elina_tcons0_array_is_interval_linear(elina_tcons0_array_t* array)
{
  size_t i;
  bool res = true;
  for (i=0; i<array->size; i++){
    if (!elina_texpr0_is_interval_linear(array->p[i].texpr0)){
      res = false;
      break;
    }
  }
  return res;
}

/* ====================================================================== */
/* II.1 Change of dimensions and permutations */
/* ====================================================================== */
void elina_tcons0_array_add_dimensions_with(elina_tcons0_array_t* array,
					 elina_dimchange_t* dimchange)
{
  size_t i;
  for(i=0; i<array->size; i++){
    elina_texpr0_t* expr = array->p[i].texpr0;
    if (expr) elina_texpr0_add_dimensions_with(expr,dimchange);
  }
}
elina_tcons0_array_t elina_tcons0_array_add_dimensions(elina_tcons0_array_t* array,
						 elina_dimchange_t* dimchange)
{
  size_t i;
  elina_tcons0_array_t narray;

  narray = elina_tcons0_array_make(array->size);
  for(i=0; i<array->size; i++){
    narray.p[i] =  elina_tcons0_add_dimensions(&array->p[i],dimchange);
  }
  return narray;
}

void elina_tcons0_array_remove_dimensions_with(elina_tcons0_array_t* array,
					    elina_dimchange_t* dimchange)
{
  size_t i;
  for(i=0; i<array->size; i++){
    elina_texpr0_t* expr = array->p[i].texpr0;
    if (expr) elina_texpr0_remove_dimensions_with(expr,dimchange);
  }
}
elina_tcons0_array_t elina_tcons0_array_remove_dimensions(elina_tcons0_array_t* array,
						    elina_dimchange_t* dimchange)
{
  size_t i;
  elina_tcons0_array_t narray;

  narray = elina_tcons0_array_make(array->size);
  for(i=0; i<array->size; i++){
    narray.p[i] =  elina_tcons0_remove_dimensions(&array->p[i],dimchange);
  }
  return narray;
}

void elina_tcons0_array_permute_dimensions_with(elina_tcons0_array_t* array,
					     elina_dimperm_t* perm)
{
  size_t i;
  for(i=0; i<array->size; i++){
    elina_texpr0_t* expr = array->p[i].texpr0;
    if (expr) elina_texpr0_permute_dimensions_with(expr,perm);
  }
}
elina_tcons0_array_t elina_tcons0_array_permute_dimensions(elina_tcons0_array_t* array,
						     elina_dimperm_t* perm)
{
  size_t i;
  elina_tcons0_array_t narray;

  narray = elina_tcons0_array_make(array->size);
  for(i=0; i<array->size; i++){
    narray.p[i] =  elina_tcons0_permute_dimensions(&array->p[i],perm);
  }
  return narray;
}

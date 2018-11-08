/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
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

#include "elina_box_representation.h"

/* ********************************************************************** */
/* Internal functions */
/* ********************************************************************** */
elina_box_t* elina_box_alloc(size_t intdim, size_t realdim)
{
  elina_box_t* itv = (elina_box_t*)malloc(sizeof(elina_box_t));
  itv->inf = NULL;
  itv->sup = NULL;
  itv->intdim = intdim;
  itv->realdim = realdim;
  return itv;
}

void elina_box_init(elina_box_t* a)
{
  size_t i;
  size_t nbdims = a->intdim + a->realdim;
  assert(a->inf==NULL && a->sup==NULL);
  a->inf = (double*)malloc((nbdims+1)*sizeof(double));
  a->sup = (double*)malloc((nbdims+1)*sizeof(double));
  /* Add an unused dimension to differentiate
     empty and top values in dimension 0+0 */
}

void elina_box_set_bottom(elina_box_t* a)
{
  if (a->inf){
    free(a->inf);
    a->inf = NULL;
  }
  if(a->sup){
     free(a->sup);
     a->sup = NULL;
  }
}

void elina_box_set_top(elina_box_t* a)
{
  size_t i;
  size_t nbdims;
  
  nbdims = a->intdim + a->realdim;
  if ((a->inf==NULL)&&(a->sup==NULL)){
    elina_box_init(a);
  };
  for (i=0; i<nbdims; i++){
    a->inf[i] = INFINITY;
    a->sup[i] = INFINITY;
  }
}

void elina_box_set(elina_box_t* a, elina_box_t* b)
{
  size_t i;
  size_t nbdims;
  
  if ((b->inf==NULL)&&(b->sup==NULL))
    return;

  nbdims = b->intdim + b->realdim;
  if ((a->inf==NULL) && (a->sup==NULL)){
    elina_box_init(a);
  };
  for (i=0; i<nbdims; i++){
    a->inf[i] = b->inf[i];
    a->sup[i] = b->sup[i];
  }
}

/* ********************************************************************** */
/* 1. Memory */
/* ********************************************************************** */

/* Return a copy of an abstract value, on
   which destructive update does not affect the initial value. */
elina_box_t* elina_box_copy(elina_manager_t* man, elina_box_t* a)
{
  size_t i;
  size_t nbdims = a->intdim+a->realdim;

  elina_box_t* b = elina_box_alloc(a->intdim,a->realdim);
  if (a->inf && a->sup){
    b->inf = (double *)malloc((nbdims+1)*sizeof(double));
    b->sup = (double *)malloc((nbdims+1)*sizeof(double));
    for (i=0; i<nbdims; i++){
      b->inf[i] = a->inf[i];
      b->sup[i] = a->sup[i];
    }
    
    b->inf[nbdims] = 0.0;
    b->sup[nbdims] = 0.0;  
    /* Add an unused dimension to differentiate
       empty and top values in dimension 0+0 */ 
  }
  man->result.flag_best = true;
  man->result.flag_exact = true;
  return b;
}

/* Free all the memory used by the abstract value */
void elina_box_free(elina_manager_t* man, elina_box_t* a)
{
  if (a->inf){
    free(a->inf);
    a->inf = NULL;
  }
  if(a->sup){
    free(a->sup);
    a->sup = NULL;
  }
  free(a);
}


/* ********************************************************************** */
/* 2. Printing */
/* ********************************************************************** */

/* Print the abstract value in a pretty way, using function
   name_of_dim to name dimensions */
void elina_box_fprint(FILE* stream,
	       elina_manager_t* man,
	       elina_box_t* a,
	       char** name_of_dim)
{
  size_t i;
  size_t nbdims = a->intdim + a->realdim;

  fprintf(stream,"interval of dim (%ld,%ld):",
	  (long)a->intdim,(long)a->realdim);
  if (a->inf && a->sup){
    fprintf(stream,"\n");
    for(i=0; i<nbdims; i++){
      if (name_of_dim){
	fprintf(stream,"%8s in ", name_of_dim[i]);
      } else {
	fprintf(stream,"x%ld in ", (long)i);
      }
      fprintf(stream,"[%g,%g]",-a->inf[i],a->sup[i]);
      fprintf(stream,"\n");
    }
  }
  else {
    fprintf(stream,nbdims>0 ? " bottom\n" : "top\n");
  }
}



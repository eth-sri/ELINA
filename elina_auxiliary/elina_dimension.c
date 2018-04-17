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


/* ************************************************************************* */
/* elina_dimension.c: dimensions and related operations */
/* ************************************************************************* */


#include <string.h>
#include <limits.h>
#include <assert.h>
#include "elina_dimension.h"

/* ====================================================================== */
/* elina_dimchange_t */
/* ====================================================================== */

/* Allocating a transformation */
void elina_dimchange_init(elina_dimchange_t* dimchange, size_t intdim, size_t realdim)
{
  dimchange->dim = intdim+realdim==0 ? NULL : malloc((intdim+realdim)*sizeof(elina_dim_t));
  dimchange->intdim = intdim;
  dimchange->realdim = realdim;
}
elina_dimchange_t* elina_dimchange_alloc(size_t intdim, size_t realdim)
{
  elina_dimchange_t* res = malloc(sizeof(elina_dimchange_t));
  elina_dimchange_init(res,intdim,realdim);
  return res;
}

/* Printing a transformation */
void elina_dimchange_fprint(FILE* stream, elina_dimchange_t* dimchange)
{
  size_t i;

  fprintf(stream,"dimchange: intdim=%lu, realdim=%lu\n           ",
	  (unsigned long)dimchange->intdim,
	  (unsigned long)dimchange->realdim);
  for (i=0;i<dimchange->intdim+dimchange->realdim;i++){
    fprintf(stream,"%2lu ",(unsigned long)dimchange->dim[i]);
  }
  fprintf(stream,"\n");
}

/* Inverting in-place an adding transformation from the given dimensionality */
void elina_dimchange_add_invert(elina_dimchange_t* dimchange)
{
  size_t i;
  for (i=0;i<dimchange->intdim+dimchange->realdim; i++){
    dimchange->dim[i] += i;
  }
}

/* ====================================================================== */
/* elina_dimchange2_t */
/* ====================================================================== */

/* Clear a dimchange structure (deallocate internal arrays) */
void elina_dimchange2_clear(elina_dimchange2_t* dimchange2)
{
  if (dimchange2->add){
    elina_dimchange_free(dimchange2->add);
    dimchange2->add = NULL;
  }
  if (dimchange2->remove){
    elina_dimchange_free(dimchange2->remove);
    dimchange2->remove = NULL;
  }
}
/* Deallocate and clear a dimchange2 structure */
void elina_dimchange2_free(elina_dimchange2_t* dimchange2)
{
  elina_dimchange2_clear(dimchange2);
  free(dimchange2);
}

/* Printing */
void elina_dimchange2_fprint(FILE* stream, elina_dimchange2_t* dimchange2)
{
  fprintf(stream,"add: ");
  if (dimchange2->add){
    elina_dimchange_fprint(stream,dimchange2->add);
  }
  else {
    fprintf(stream,"NULL\n");
  }
  fprintf(stream,"remove: ");
  if (dimchange2->remove){
    elina_dimchange_fprint(stream,dimchange2->remove);
  }
  else {
    fprintf(stream,"NULL\n");
  }
}

/* ====================================================================== */
/* elina_dimperm_t */
/* ====================================================================== */

/* Allocating a permutation */
void elina_dimperm_init(elina_dimperm_t* dimperm, size_t size)
{
  dimperm->dim = size==0 ? NULL : malloc(size*sizeof(elina_dim_t));
  dimperm->size = size;
}
elina_dimperm_t* elina_dimperm_alloc(size_t size)
{
  elina_dimperm_t* dimperm = malloc(sizeof(elina_dimperm_t));
  elina_dimperm_init(dimperm,size);
  return dimperm;
}

/* Printing a permutation */
void elina_dimperm_fprint(FILE* stream, elina_dimperm_t* perm)
{
  size_t i;

  fprintf(stream,"dimperm: size=%lu\n",(unsigned long)perm->size);
  for (i=0;i<perm->size;i++){
    fprintf(stream,"%2lu -> %2lu\n",(unsigned long)i,
	    (unsigned long)perm->dim[i]);
  }
}

/* Generates the identity permutation */
void elina_dimperm_set_id(elina_dimperm_t* perm)
{
  size_t i;
  for (i=0; i<perm->size; i++) perm->dim[i] = i;
}

/* Compose 2 permutations */
void elina_dimperm_compose(elina_dimperm_t* perm, elina_dimperm_t* perm1, elina_dimperm_t* perm2)
{
  size_t i;

  assert(perm->size==perm1->size && perm->size==perm2->size);
  for (i=0; i<perm->size; i++){
    perm->dim[i] = perm2->dim[perm1->dim[i]];
  }
}
/* Invert a permutation */
void elina_dimperm_invert(elina_dimperm_t* nperm, elina_dimperm_t* perm)
{
  size_t i;

  assert(nperm->size==perm->size);
  for (i=0; i<perm->size; i++){
    nperm->dim[perm->dim[i]] = i;
  }
}


/* ====================================================================== */
/* Auxilliary Functions Definitions */
/* ====================================================================== */
void elina_dimchange_clear(elina_dimchange_t* dimchange)
{
  if (dimchange->dim) free(dimchange->dim);
  dimchange->intdim = dimchange->realdim = 0;
  dimchange->dim = NULL;
}

void elina_dimchange_free(elina_dimchange_t* dimchange)
{
  elina_dimchange_clear(dimchange);
  free(dimchange);
}


void elina_dimperm_clear(elina_dimperm_t* dimperm)
{
  if (dimperm->dim) free(dimperm->dim);
  dimperm->size = 0;
  dimperm->dim = NULL;
}


void elina_dimperm_free(elina_dimperm_t* dimperm)
{
  elina_dimperm_clear(dimperm);
  free(dimperm);
}


void elina_dimchange2_init(elina_dimchange2_t* dimchange2,
			elina_dimchange_t* add, 
			elina_dimchange_t* remove)
{ 
  dimchange2->add = add;
  dimchange2->remove = remove;
}

elina_dimchange2_t* elina_dimchange2_alloc(elina_dimchange_t* add, 
				     elina_dimchange_t* remove)
{
  elina_dimchange2_t* res = (elina_dimchange2_t*)malloc(sizeof(elina_dimchange2_t));
  elina_dimchange2_init(res,add,remove);
  return res;
}

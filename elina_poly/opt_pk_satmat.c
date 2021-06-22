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


/* ********************************************************************** */
/* opt_pk_satmat.c: operations on saturation matrices */
/* ********************************************************************** */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "opt_pk_satmat.h"

/* ********************************************************************** */
/* I. basic operations: creation, destruction, copying and printing */
/* ********************************************************************** */

/* Set all bits to zero. */
void opt_satmat_clear(opt_satmat_t* os)
{
  size_t i,j;
  for (i=0; i<os->nbrows; i++)
    for (j=0; j<os->nbcolumns; j++)
      os->p[i][j] = 0;
}

/* Standard allocation function, with initialization of the elements. */
opt_satmat_t* opt_satmat_alloc(size_t nbrows, size_t nbcols)
{
  size_t i,j;

  opt_satmat_t* os = (opt_satmat_t*)malloc(sizeof(opt_satmat_t));
  os->nbrows = os->_maxrows = nbrows;
  os->nbcolumns = nbcols;
  os->p = (opt_bitstring_t**)malloc(nbrows*sizeof(opt_bitstring_t*));
  for (i=0; i<nbrows; i++){
    os->p[i] = opt_bitstring_alloc(nbcols);
    for (j=0; j<os->nbcolumns; j++)
      os->p[i][j] = 0;
  }
  return os;
}

/* Deallocation function. */
void opt_satmat_free(opt_satmat_t* os)
{
  size_t i;
  for (i=0;i<os->_maxrows;i++){
    opt_bitstring_free(os->p[i]);
  }
  free(os->p);
  free(os);
  os=NULL;
}

/* Reallocation function, to scale up or to downsize a matrix */
void opt_satmat_resize_rows(opt_satmat_t* os, size_t nbrows)
{
  size_t i;

  if (nbrows > os->_maxrows){
    os->p = (opt_bitstring_t**)realloc(os->p, nbrows * sizeof(opt_bitstring_t*));
    for (i=os->nbrows; i<nbrows; i++){
      os->p[i] = opt_bitstring_alloc(os->nbcolumns);
    }
  }
  else if (nbrows < os->_maxrows){
    for (i=nbrows; i<os->_maxrows; i++){
      opt_bitstring_free(os->p[i]);
    }
    os->p = (opt_bitstring_t**)realloc(os->p,nbrows * sizeof(opt_bitstring_t*));
  } 
  os->nbrows = nbrows;
  os->_maxrows = nbrows;
}

/* Reallocation function, to scale up or to downsize a matrix */
void opt_satmat_resize_cols(opt_satmat_t* os, size_t nbcols)
{
  size_t i,j;

  if (nbcols!=os->nbcolumns){
    for (i=0; i<os->_maxrows; i++){
      os->p[i] = opt_bitstring_realloc(os->p[i],nbcols);
      for (j=os->nbcolumns; j<nbcols; j++){
	os->p[i][j] = 0;
      }
    }
    os->nbcolumns = nbcols;
  }
}

/* Create a copy of the matrix of size nbrows (and not
   _maxrows) and extends columns. Only ``used'' rows are copied. */
opt_satmat_t* opt_satmat_copy_resize_cols(opt_satmat_t* os, size_t nbcols)
{
  size_t i,j;
  opt_satmat_t* nos;

  assert(nbcols>=os->nbcolumns);
  nos = opt_satmat_alloc(os->nbrows,nbcols);
  for (i=0; i<os->nbrows; i++){
    for (j=0; j<os->nbcolumns; j++){
      nos->p[i][j] = os->p[i][j];
    }
    for (j=os->nbcolumns; j<nbcols; j++)
      nos->p[i][j] = 0;
  }
  return nos;
}

/* Create a copy of the matrix of size nbrows (and not
   _maxrows). Only ``used'' rows are copied. */
opt_satmat_t* opt_satmat_copy(opt_satmat_t* os)
{
  size_t i,j;
  opt_satmat_t* nos = opt_satmat_alloc(os->nbrows,os->nbcolumns);
  for (i=0; i<os->nbrows; i++){
    for (j=0; j<os->nbcolumns; j++){
      nos->p[i][j] = os->p[i][j];
    }
  }
  return nos;
}

/* Raw printing function. */
void opt_satmat_fprint(FILE* stream, opt_satmat_t* os)
{
  size_t i;

  fprintf(stream,"%lu %lu\n",
	  (unsigned long)os->nbrows,(unsigned long)os->nbcolumns);
  for (i=0; i<os->nbrows; i++){
    opt_bitstring_fprint(stream,os->p[i],os->nbcolumns);
    fprintf(stream,"\n");
  }
}

void opt_satmat_print(opt_satmat_t* os)
{
  opt_satmat_fprint(stdout,os);
}

/* ********************************************************************** */
/* II. Bit operations */
/* ********************************************************************** */

/* These function allow to read and to clear or set individual bits. i
   indicates the row and jx the column. */

opt_bitstring_t opt_satmat_get(opt_satmat_t* os, size_t i, opt_bitindex_t ojx){
  return opt_bitstring_get(os->p[i],ojx);
}
void opt_satmat_set(opt_satmat_t* os, size_t i, opt_bitindex_t ojx){
  opt_bitstring_set(os->p[i],ojx);
}
void opt_satmat_clr(opt_satmat_t* os, size_t i, opt_bitindex_t ojx){
  opt_bitstring_clr(os->p[i],ojx);
}

/* ********************************************************************** */
/* III. Matrix operations */
/* ********************************************************************** */

/* Transposition.

nbcols indicates the number of bits to be transposed (the number of columns of
the matrix is the size of the row of opt_bitstring_t, not the number of bits really
used). */

opt_satmat_t* opt_satmat_transpose(opt_satmat_t* src, size_t nbcols)
{
  opt_bitindex_t i,j;
  opt_satmat_t* dst = opt_satmat_alloc(nbcols,opt_bitindex_size(src->nbrows));

  for (i = opt_bitindex_init(0); i.index < src->nbrows; opt_bitindex_inc(&i) ) {
    for (j = opt_bitindex_init(0); j.index < nbcols; opt_bitindex_inc(&j) ){
      if (opt_satmat_get(src,i.index,j)) opt_satmat_set(dst,j.index,i);
    }
  }
  return dst;
}

/* Row exchange. */
void opt_satmat_exch_rows(opt_satmat_t* os, size_t l1, size_t l2)
{
  opt_bitstring_t* tmp=os->p[l1];
  os->p[l1]=os->p[l2];
  os->p[l2]=tmp;
}

void opt_satmat_exch_cols(opt_satmat_t* os, size_t l1, size_t l2)
{
  size_t i;
  opt_bitindex_t j1 = opt_bitindex_init(l1);
  opt_bitindex_t j2 = opt_bitindex_init(l2);
  for(i=0; i < os->nbrows; i++){
		opt_bitstring_t tmp1 = opt_satmat_get(os,i,j1);
		opt_bitstring_t tmp2 = opt_satmat_get(os,i,j2);
		if(tmp1 && !tmp2){
			opt_satmat_clr(os,i,j1);
			opt_satmat_set(os,i,j2);
		}
		else if(!tmp1 && tmp2){
			opt_satmat_set(os,i,j1);
			opt_satmat_clr(os,i,j2);
		}
  }
  
}

void opt_satmat_move_rows(opt_satmat_t* os, size_t destrow, size_t orgrow, size_t size)
{
  int i,offset;

  offset = destrow-orgrow;
  if (offset>0){
    assert(destrow+size<=os->_maxrows);
    for (i=(int)(destrow+size)-1; i>=(int)destrow; i--){
      opt_satmat_exch_rows(os,(size_t)i,(size_t)(i-offset));
    }
  } else {
    assert(orgrow+size<=os->_maxrows);
    for(i=(int)destrow; i<(int)(destrow+size); i++){
      opt_satmat_exch_rows(os,(size_t)i,(size_t)(i-offset));
    }
  }
}

size_t opt_satmat_serialize_common(void* dst, opt_satmat_t* oc, int dry_run){
  size_t idx = 0;

  if(!dry_run){
    *(size_t*)(dst + idx) = oc->nbrows;
  }
  idx += sizeof(size_t);

  if(!dry_run){
    *(unsigned short int*)(dst + idx) = oc->nbcolumns;
  }
  idx += sizeof(unsigned short int);

  if(!dry_run){
    *(size_t*)(dst + idx) = oc->_maxrows;
  }
  idx += sizeof(size_t);

  for(size_t i = 0; i < oc->_maxrows; i++){
    for(unsigned short int j = 0; j < oc->nbcolumns; j++){
      if(!dry_run){
	*(opt_bitstring_t*)(dst + idx) = oc->p[i][j];
      }
      idx += sizeof(opt_bitstring_t);
    }
  }

  return idx;
}

opt_satmat_t* opt_satmat_deserialize(void* p, size_t* size){
  opt_satmat_t* oc = (opt_satmat_t*)malloc(sizeof(opt_satmat_t));
  size_t idx = 0;

  oc->nbrows = *(size_t*)(p + idx);
  idx += sizeof(size_t);

  oc->nbcolumns = *(unsigned short int*)(p + idx);
  idx += sizeof(unsigned short int);

  oc->_maxrows = *(size_t*)(p + idx);
  idx += sizeof(size_t);

  oc->p = (opt_bitstring_t**)malloc(sizeof(opt_bitstring_t*) * oc->_maxrows);
  for(size_t i = 0; i < oc->_maxrows; i++){
    oc->p[i] = (opt_bitstring_t*)malloc(sizeof(opt_bitstring_t) * oc->nbcolumns);
    for(unsigned short int j = 0; j < oc->nbcolumns; j++){
      oc->p[i][j] = *(opt_bitstring_t*)(p + idx);
      idx += sizeof(opt_bitstring_t);
    }
  }

  *size = idx;
  return oc;
}

/*
	Copyright 2016 Software Reliability Lab, ETH Zurich

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

		http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
*/


/* ********************************************************************** */
/* opt_pk_expandfold.c: expanding and folding dimensions */
/* ********************************************************************** */


#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"
#include "opt_pk_representation.h"
#include "opt_pk_user.h"
#include "opt_pk_constructor.h"
#include "opt_pk_meetjoin.h"


/* ********************************************************************** */
/* I. Expand */
/* ********************************************************************** */

/* ---------------------------------------------------------------------- */
/* Matrix */
/* ---------------------------------------------------------------------- */

/* Expand the dimension dim of the matrix into (dimsup+1)
   dimensions, with dimsup new dimensions inserted just before
   offset. */
static
opt_matrix_t* opt_matrix_expand(opt_pk_internal_t* opk,
			bool destructive,
			opt_matrix_t* C,
			elina_dim_t dim,
			size_t offset,
			size_t dimsup)
{
  elina_dimchange_t* dimchange;
  size_t i,j,row,col,nb;
  size_t nbrows, nbcols;
  opt_numint_t** p;
  opt_matrix_t* nC;

  if (dimsup==0){
    return destructive ? C : opt_matrix_copy(C);
  }
  nbrows = C->nbrows;
  nbcols = C->nbcolumns;
  col = opk->dec + dim;
  /* Count the number of constraints to duplicate */
  nb=0;
  p = C->p;
  for (i=0; i<nbrows; i++){
    if (p[i][col]!=0)
      nb++;
  }
  /* Redimension matrix */
  dimchange = elina_dimchange_alloc(0,dimsup);
  for (i=0;i<dimsup;i++){
    dimchange->dim[i]=offset;
  }
  nC = opt_matrix_add_dimensions(opk,destructive,C,dimchange, false);
  elina_dimchange_free(dimchange);
  opt_matrix_resize_rows(nC,nbrows+nb*dimsup);
  if (nb==0)
    return nC;

  /* Duplicate constraints */
  p = nC->p;
  row = nbrows;
  for (i=0; i<nbrows; i++){
    if (p[i][col]!=0){
      for (j=offset;j < offset+dimsup; j++){
	opt_vector_copy(p[row],
		    p[i],
		    nbcols+dimsup);
	p[row][opk->dec+j] = p[row][col];
	p[row][col] = 0;
	row++;
      }
    }
  }
  nC->_sorted = false;
  return nC;
}

/* ---------------------------------------------------------------------- */
/* Polyhedra */
/* ---------------------------------------------------------------------- */

opt_pk_t* opt_pk_expand(elina_manager_t* man,
		bool destructive, opt_pk_t* oa,
		elina_dim_t dim, size_t dimsup)
{
  size_t intdimsup,realdimsup;
  size_t nintdim,nrealdim;
  opt_pk_t* op;

  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_EXPAND);
  opt_pk_internal_realloc_lazy(opk,oa->intdim+oa->realdim+dimsup);
  man->result.flag_best = man->result.flag_exact = true;   

  if (dim<oa->intdim){
    intdimsup = dimsup;
    realdimsup = 0;
  } else {
    intdimsup = 0;
    realdimsup = dimsup;
  }
  nintdim = oa->intdim + intdimsup;
  nrealdim = oa->realdim + realdimsup;

  if (dimsup==0){
    return (destructive ? oa : opt_pk_copy(man,oa));
  }

  /* Get the constraints system, and possibly minimize */
  if (opk->funopt->algorithm>0)
    opt_poly_minimize(man,oa);
  //else
    //poly_chernikova(man,pa,"of the argument");

  if (destructive){
    op = oa;
    op->intdim+=intdimsup;
    op->realdim+=realdimsup;
    op->status &= ~opt_pk_status_consgauss & ~opt_pk_status_gengauss & ~opt_pk_status_minimaleps;
  }
  else {
    op = opt_poly_alloc(nintdim,nrealdim);
  }

  if (opk->exn){
    opk->exn = ELINA_EXC_NONE;
    if (!oa->C){
      man->result.flag_best = man->result.flag_exact = false;   
      opt_poly_set_top(opk,op);
      return op;
    }
    /* We can still proceed, although it is likely 
       that the problem is only delayed
    */
  }
  /* if empty, return empty */
  if (!oa->C){
    opt_poly_set_bottom(opk,op);
    return op;
  }
  /* Prepare resulting matrix */
  if (destructive){
    op->nbeq  = 0;
    op->status &= ~opt_pk_status_consgauss & ~opt_pk_status_gengauss & ~opt_pk_status_minimaleps;
  }
  op->C = opt_matrix_expand(opk, destructive, oa->C, 
			dim, 
			(dim + dimsup < op->intdim ?
			 op->intdim-dimsup :
			 op->intdim+op->realdim-dimsup),
			dimsup);
  /* Minimize the result */
  if (opk->funopt->algorithm>0){
    opt_poly_minimize(man,op);
    if (opk->exn){
      opk->exn = ELINA_EXC_NONE;
      if (!op->C){
	man->result.flag_best = man->result.flag_exact = false;   
	opt_poly_set_top(opk,op);
	return op;
      }
    }
  }
  //assert(poly_check(pk,po));
  return op;
}

/* ********************************************************************** */
/* II. Fold */
/* ********************************************************************** */

/* ---------------------------------------------------------------------- */
/* Matrix */
/* ---------------------------------------------------------------------- */

/* Fold the last dimsup dimensions with dimension dim (not in the last dimsup
   ones) in the matrix */


/* the array tdim is assumed to be sorted */

opt_pk_t* opt_pk_fold(elina_manager_t* man,
	      bool destructive, opt_pk_t* oa,
	      elina_dim_t* tdim, size_t size)
{
  size_t intdimsup,realdimsup;
  opt_pk_t* op;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_FOLD);
  man->result.flag_best = man->result.flag_exact = true;   

  if (tdim[0]<oa->intdim){
    intdimsup = size - 1;
    realdimsup = 0;
  } else {
    intdimsup = 0;
    realdimsup = size - 1;
  }
  
  opt_poly_minimize(man,oa);

  if (destructive){
    op = oa;
    op->intdim -= intdimsup;
    op->realdim -= realdimsup;
  }
  else {
    op = opt_poly_alloc(oa->intdim-intdimsup,oa->realdim-realdimsup);
  }
  if (opk->exn){
    opk->exn = ELINA_EXC_NONE;
    if (!oa->C){
      man->result.flag_best = man->result.flag_exact = false;   
      opt_poly_set_top(opk,op);
      return op;
    }
  }
  /* if empty, return empty */
  if (!oa->C){
    man->result.flag_best = man->result.flag_exact = true;   
    opt_poly_set_bottom(opk,op);
    return op;
  }

  /* Prepare resulting matrix */
  if (destructive){
    op->nbeq = 0;
    op->status &= ~opt_pk_status_consgauss & ~opt_pk_status_gengauss & ~opt_pk_status_minimaleps;
  }
  
  //op->C = matrix_fold(pk, destructive, pa->F, 
	//	      tdim, size);
  /* Minimize the result */
  if (opk->funopt->algorithm>0){
    opt_poly_minimize(man,op);
    if (opk->exn){
      opk->exn = ELINA_EXC_NONE;
      if (!op->C){
	man->result.flag_best = man->result.flag_exact = false;   
	opt_poly_set_top(opk,op);
	return op;
      }
    }
  }

  man->result.flag_best = (intdimsup==0);
  man->result.flag_exact = (intdimsup==0);
  //assert(poly_check(pk,po));
  return op;
}

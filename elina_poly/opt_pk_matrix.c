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

/* ********************************************************************** */
/* opt_matrix.c: operations on matrices */
/* ********************************************************************** */


#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_satmat.h"
#include "opt_pk_matrix.h"
#include "opt_mf_qsort.h"

/* ********************************************************************** */
/* I. basic operations: creation, destruction, copying and printing */
/* ********************************************************************** */

/* Internal allocation function: the elements are not initialized.
   mr is the maximum number of rows, and nc the number of
   columns. By default, nbrows is initialized to mr . */
opt_matrix_t* _opt_matrix_alloc_int(size_t nbrows, unsigned short int nbcols, bool s)
{
  size_t i;

  assert(nbcols>0 || nbrows==0);

  opt_matrix_t* mat = (opt_matrix_t*)malloc(sizeof(opt_matrix_t));
  mat->nbrows = mat->_maxrows = nbrows;
  mat->nbcolumns = nbcols;
  mat->_sorted = s;
  mat->p = (opt_numint_t**)malloc(nbrows * sizeof(opt_numint_t*));
  for (i=0;i<nbrows;i++){
    mat->p[i] = _opt_vector_alloc_int(nbcols);
  }
  return mat;
}


/* Standard allocation function, with initialization of the elements. */
opt_matrix_t* opt_matrix_alloc(size_t nbrows, unsigned short int nbcols, bool s)
{
  size_t i;

  assert(nbcols>0 || nbrows==0);

  opt_matrix_t* mat = (opt_matrix_t*)malloc(sizeof(opt_matrix_t));
  mat->nbrows = mat->_maxrows = nbrows;
  mat->nbcolumns = nbcols;
  mat->_sorted = s;
  mat->p = (opt_numint_t**)malloc(nbrows * sizeof(opt_numint_t*));
  for (i=0;i<nbrows;i++){
    mat->p[i] = opt_vector_alloc(nbcols);
  }
  return mat;
}

/* Reallocation function, to scale up or to downsize a matrix */
void opt_matrix_resize_rows(opt_matrix_t* mat, size_t nbrows)
{
  size_t i;

  assert (nbrows>0);
  if (nbrows > mat->_maxrows){
    mat->p = (opt_numint_t**)realloc(mat->p, nbrows * sizeof(opt_numint_t*));
    for (i=mat->_maxrows; i<nbrows; i++){
      mat->p[i] = opt_vector_alloc(mat->nbcolumns);
    }
    mat->_sorted = false;
  }
  else if (nbrows < mat->_maxrows){
    for (i=nbrows; i<mat->_maxrows; i++){
      opt_vector_free(mat->p[i],mat->nbcolumns);
    }
    mat->p = (opt_numint_t**)realloc(mat->p,nbrows * sizeof(opt_numint_t*));
  }
  mat->_maxrows = nbrows;
  mat->nbrows = nbrows;
}

/* Ensures a minimum size */
void opt_matrix_resize_rows_lazy(opt_matrix_t* mat, size_t nbrows)
{
  if (nbrows>mat->_maxrows)
    opt_matrix_resize_rows(mat,nbrows);
  else {
    mat->_sorted = mat->_sorted && nbrows<mat->nbrows;
    mat->nbrows = nbrows;
  }
}

/* Minimization */
//void matrix_minimize(matrix_t* mat)
//{
//  matrix_resize_rows(mat,mat->nbrows);
//}

/* Deallocation function. */
void opt_matrix_free(opt_matrix_t* omat)
{
  size_t i;
  for (i=0;i<omat->_maxrows;i++){
    opt_vector_free(omat->p[i],omat->nbcolumns);
  }
  if(omat->p){
  	free(omat->p);
  }
  free(omat);
  omat = NULL;
}

/* Set all elements to zero. */
void opt_matrix_clear(opt_matrix_t* mat)
{
  size_t i,j;
  for (i=0; i<mat->nbrows; i++){
    for (j=0; j<mat->nbcolumns; j++){
      mat->p[i][j] = 0;
    }
  }
}

/* Create a copy of the matrix of size nbrows (and not
   _maxrows). Only ``used'' rows are copied. */
opt_matrix_t* opt_matrix_copy(opt_matrix_t* src)
{
  size_t i;
  opt_matrix_t* dst = _opt_matrix_alloc_int(src->nbrows,src->nbcolumns,src->_sorted);
  size_t nbcons = src->nbrows;
  size_t nbcolumns = src->nbcolumns;
  for (i=0;i < nbcons;i++){
	 opt_numint_t *si = src->p[i];
	 opt_numint_t *di = dst->p[i];
	 opt_vector_copy(di, si, nbcolumns);
  }
  return dst;
}

/* Return true iff the matrices are equal, coeff by coeff */
bool matrix_equal(opt_matrix_t* oma, opt_matrix_t* omb)
{
  int i;
  size_t j;
  bool res;

  res = oma->nbrows==omb->nbrows && oma->nbcolumns==omb->nbcolumns;
  if (!res) return res;
  for (i=(int)oma->nbrows-1;i>=0;i--){
    for (j=0; j<oma->nbcolumns; j++){
      res = (oma->p[i][j]==omb->p[i][j]);
      if (!res) return res;
    }
  }
  return res;
}

/* Raw printing function. */
void opt_matrix_fprint(FILE* stream, opt_matrix_t* mat)
{
  size_t i,j;
  fprintf(stream,"%lu %lu\n", 
	  (unsigned long)mat->nbrows, (unsigned long)mat->nbcolumns);
  for (i=0;i<mat->nbrows;i++) {
    for (j=0;j<mat->nbcolumns;j++){
	//if((j>2) && (mat->p[i][j] !=0))
		//fprintf(stdout,"%d ",j-3);
      //numint_fprint(stream,mat->p[i][j]);
      fprintf(stream,"%lld ",mat->p[i][j]);
    }
    fprintf(stream,"\n");
  }
}


void opt_matrix_print(opt_matrix_t* mat)
{
  opt_matrix_fprint(stdout,mat);
}


/* ********************************************************************** */
/* II.Operation on rows */
/* ********************************************************************** */

/* compare_rows compares rows of matrix, exch_rows exchanges
   two rows; normalize_row normalizes a row of a matrix but without
   considering the first coefficient; combine_rows combine rows
   l1 and l2 and puts the result in l3 such that
   l3[k] is zero. */

int opt_matrix_compare_rows(opt_pk_internal_t* opk,
			opt_matrix_t* oc, size_t l1, size_t l2)
{
  return opt_vector_compare(opk, 
			oc->p[l1],
			oc->p[l2],oc->nbcolumns);
}
void opt_matrix_normalize_row(opt_pk_internal_t* opk,
			  opt_matrix_t* oc, size_t l)
{
  opt_vector_normalize(opk, oc->p[l],oc->nbcolumns);
}
void opt_matrix_combine_rows(opt_pk_internal_t* opk,
			 opt_matrix_t* oc, size_t l1, size_t l2, size_t l3, size_t k, bool add)
{
  opt_vector_combine(opk, 
		 oc->p[l1],
		 oc->p[l2],
		 oc->p[l3],k,oc->nbcolumns,add);
}
void opt_matrix_exch_rows(opt_matrix_t* om, size_t l1, size_t l2)
{
  opt_numint_t* tmp =om->p[l1];
  om->p[l1]=om->p[l2];
  om->p[l2]=tmp;
}

void opt_matrix_move_rows(opt_matrix_t* oc, size_t destrow, size_t orgrow, size_t size)
{
  int offset;
  int i;

  offset = destrow-orgrow;
  if (offset>0){
    assert(destrow+size<=oc->_maxrows);
    for (i=(int)(destrow+size)-1; i>=(int)destrow; i--){
      opt_matrix_exch_rows(oc,(size_t)i,(size_t)(i-offset));
    }
  } else {
    assert(orgrow+size<=oc->_maxrows);
    for(i=(int)destrow; i<(int)(destrow+size); i++){
      opt_matrix_exch_rows(oc,(size_t)i,(size_t)(i-offset));
    }
  }
}

void opt_matrix_union(opt_matrix_t *op1, opt_matrix_t *op2){
	size_t nbcons1 = op1->nbrows;
	size_t nbcons2 = op2->nbrows;
	size_t nbcons = nbcons1 + nbcons2;
	unsigned short int nbcolumns = op1->nbcolumns;
	size_t i;
	assert(nbcolumns == op2->nbcolumns);
	opt_matrix_resize_rows(op1,nbcons);
	opt_numint_t **p1 = op1->p;
	opt_numint_t **p2 = op2->p;
	for(i=0; i < nbcons2; i++){
		size_t ind = nbcons1 + i;
		opt_vector_copy(p1[ind],p2[i],nbcolumns);
	}
}

/* ********************************************************************** */
/* */
/* ********************************************************************** */

bool opt_matrix_normalize_constraint(opt_pk_internal_t* opk,
				 opt_matrix_t* oc, 
				 size_t intdim, size_t realdim)
{
  bool change1, change2;
  size_t i;

  if ( opk->strict && realdim>0 ){
    change2=false;
    for (i=0; i<oc->nbrows; i++){
      change1 = opt_vector_normalize_constraint(opk,oc->p[i],intdim,realdim);
      change2 = change2 || change1;
    }
    if (change2){
      oc->_sorted = false;
      /* Add again \xi-\epsilon<=1 */
      size_t nbrows= oc->nbrows;
      opt_matrix_resize_rows_lazy(oc,nbrows+1);
      opt_vector_clear(oc->p[nbrows],oc->nbcolumns);
      oc->p[nbrows][0] = 1;
      oc->p[nbrows][opt_polka_cst] = 1;
      oc->p[nbrows][opt_polka_eps] = -1;
    }
    return change2;
  }
  else 
    return false;
}



/* ********************************************************************** */
/* III. Sorting and merging */
/* ********************************************************************** */

/* ====================================================================== */
/* III.1 Sorting */
/* ====================================================================== */

/* We use here the quick sort. There is here no handling of doublons */
typedef struct opt_qsort_man_t {
  opt_pk_internal_t* opk;
  unsigned short int size;
} opt_qsort_man_t;

static int opt_qsort_rows_compar(void* opt_qsort_man, const void* pq1, const void* pq2)
{
  opt_qsort_man_t* qm = (opt_qsort_man_t*)opt_qsort_man;
  opt_numint_t* q1 = *((opt_numint_t**)pq1);
  opt_numint_t* q2 = *((opt_numint_t**)pq2);
  return opt_vector_compare(qm->opk,q1,q2,qm->size);
}


void opt_matrix_sort_rows_from(opt_pk_internal_t* opk,
		      opt_matrix_t* om, size_t start, size_t num){
    opt_qsort_man_t opt_qsort_man;
    opt_qsort_man.opk = opk;
    opt_qsort_man.size = om->nbcolumns;
    opt_qsort2(om->p + start, num, sizeof(opt_numint_t*),
	   opt_qsort_rows_compar,
	   &opt_qsort_man);	
}

void opt_matrix_sort_rows(opt_pk_internal_t* opk,
		      opt_matrix_t* om)
{
  opt_qsort_man_t opt_qsort_man;
  if (!om->_sorted){
    opt_qsort_man.opk = opk;
    opt_qsort_man.size = om->nbcolumns;
    opt_qsort2(om->p, om->nbrows, sizeof(opt_numint_t*),
	   opt_qsort_rows_compar,
	   &opt_qsort_man);
    om->_sorted = true;
  }
}

/* This variant permutes also the saturation matrix together with the matrix.
   There is here no handling of doublons. */

typedef struct opt_qsort_t {
  opt_numint_t* p;
  opt_bitstring_t* satp;
} opt_qsort_t;

static int opt_qsort_rows_with_sat_compar(void* qsort_man, const void* q1, const void* q2)
{
  opt_qsort_man_t* qm = (opt_qsort_man_t*)qsort_man;
  const opt_qsort_t* qs1 = (const opt_qsort_t*)q1;
  const opt_qsort_t* qs2 = (const opt_qsort_t*)q2;
  return opt_vector_compare( qm->opk, 
			 qs1->p, 
			 qs2->p,
			 qm->size );
}

void opt_matrix_sort_rows_with_sat(opt_pk_internal_t* opk,
			       opt_matrix_t* mat, opt_satmat_t* sat)
{
  size_t i;
  opt_qsort_t* qsort_tab;
  opt_qsort_man_t qsort_man;

  if (!mat->_sorted){
    qsort_man.opk = opk;
    qsort_man.size = mat->nbcolumns;
    qsort_tab = (opt_qsort_t*)malloc(mat->nbrows * sizeof(opt_qsort_t));
    for (i=0; i<mat->nbrows; i++){
      qsort_tab[i].p = mat->p[i];
      qsort_tab[i].satp = sat->p[i];
    }
    opt_qsort2(qsort_tab,
	   mat->nbrows, sizeof(opt_qsort_t),
	   opt_qsort_rows_with_sat_compar,
	   &qsort_man);
    for (i=0; i<mat->nbrows; i++){
      mat->p[i] = qsort_tab[i].p;
      sat->p[i] = qsort_tab[i].satp;
    }
    free(qsort_tab);
    mat->_sorted = true;
  }
}

/* ====================================================================== */
/* III.2 Append */
/* ====================================================================== */

/* Appending matrices */
opt_matrix_t* opt_matrix_append(opt_matrix_t* oma, opt_matrix_t* omb)
{
  opt_matrix_t* mat;
  size_t i,l;

  assert (oma->nbcolumns == omb->nbcolumns);
  mat = _opt_matrix_alloc_int(oma->nbrows+omb->nbrows,oma->nbcolumns,false);
  for (i=0;i<oma->nbrows; i++){
    opt_numint_t * src = oma->p[i];
    opt_numint_t *dst = mat->p[i];
    opt_vector_copy(dst, src, oma->nbcolumns);
  }
  for (i=0;i<omb->nbrows; i++){
    opt_numint_t * src = omb->p[i];
    opt_numint_t *dst = mat->p[oma->nbrows+i];
    opt_vector_copy(dst, src, omb->nbcolumns);
  }
  return mat;
}

void opt_matrix_append_with(opt_matrix_t* oma, opt_matrix_t* omb)
{
  size_t i,l;
  size_t nbrows;

  assert (oma->nbcolumns == omb->nbcolumns);
  
  nbrows = oma->nbrows;
  opt_matrix_resize_rows_lazy(oma,nbrows+omb->nbrows);
  for (i=0;i<omb->nbrows; i++){
    opt_vector_copy(oma->p[nbrows+i], omb->p[i], oma->nbcolumns);
  }
  oma->_sorted = false;
}

/* Given matrices with rows p1,p2,... and q1,q2,...., fills the initial matrix
   with rows q1,q2,...,p1,p2,.... */

void opt_matrix_revappend_with(opt_matrix_t* oma, opt_matrix_t* omb)
{
  int i;
  unsigned short int l;
  size_t nbrows;

  assert(oma->nbcolumns == omb->nbcolumns);
  nbrows = oma->nbrows;
  opt_matrix_resize_rows_lazy(oma,nbrows+omb->nbrows);
  for (i=nbrows-1; i>=0; i--){
    /* exchanging rows i and i+cmat->nbrows */
    opt_numint_t* tmp = oma->p[i+omb->nbrows];
    oma->p[i+omb->nbrows] = oma->p[i];
    oma->p[i] = tmp;
  }
  for (i=0; i<(int)omb->nbrows; i++){
    for (l=0;l<omb->nbcolumns; l++){
      oma->p[i][l] = omb->p[i][l];
    }
  }
}


/* ====================================================================== */
/* III.3 Addition of sorted rows */
/* ====================================================================== */

/* Merging with sorting */

opt_matrix_t* opt_matrix_merge_sort(opt_pk_internal_t* opk,
			    opt_matrix_t* oma, opt_matrix_t* omb)
{
  size_t i,ia,ib,l;
  opt_matrix_t* mat;
  size_t nbrows;
  assert (oma->nbcolumns == omb->nbcolumns);
  if (!oma->_sorted || !omb->_sorted){
    mat = opt_matrix_append(oma,omb);
    opt_matrix_sort_rows(opk,mat);
  }
  else {
    mat = _opt_matrix_alloc_int(oma->nbrows+omb->nbrows,oma->nbcolumns,true);
    i = 0;
    ia = 0;
    ib = 0;
    while (ia < oma->nbrows && ib < omb->nbrows) {
      int res = opt_vector_compare(opk,
			       oma->p[ia],
			       omb->p[ib],
			       mat->nbcolumns);
      if (res<=0){
	opt_vector_copy(mat->p[i], oma->p[ia], mat->nbcolumns);
	ia++;
	if (res==0) ib++;
      }
      else {
	opt_vector_copy(mat->p[i], omb->p[ib], mat->nbcolumns);
	ib++;
      }
      i++;
    }
    /* does some constraint remain ? */
    if (ia < oma->nbrows) {
      do {
	opt_vector_copy(mat->p[i], oma->p[ia], mat->nbcolumns);
	ia++; i++;
      } while (ia < oma->nbrows);
    } else {
      while (ib < omb->nbrows){
	opt_vector_copy(mat->p[i], omb->p[ib], mat->nbcolumns);
	ib++; i++;
      }
    }
    nbrows = (size_t)i;
    /* initialize last rows of mat to zero */
    while (i<mat->nbrows){
      opt_vector_clear(mat->p[i], mat->nbcolumns);
      i++;
    }
    mat->nbrows = nbrows;
  }

  return mat;
}

/* This function adds to a sorted matrix the rows of another sorted matrix
   and leaves the resulting matrix sorted. Identical rows are eliminated. The
   modified matrix is supposed to be big enough to store the new rows. */

void opt_matrix_merge_sort_with(opt_pk_internal_t* opk,
			    opt_matrix_t* oma, opt_matrix_t* omb)
{
  size_t i,ia,ib,j,k,nbrows,nbrowsa;
  unsigned short int nbcols;
  opt_numint_t** pp;
  
  assert (oma->nbcolumns == omb->nbcolumns);
  assert (oma->_sorted && omb->_sorted);

  nbrowsa = oma->nbrows;
  nbcols = oma->nbcolumns;
  opt_matrix_resize_rows_lazy(oma, nbrowsa + omb->nbrows);
  
  /* one adds the coefficients of omb to oma */
  for (i=0; i<oma->nbrows; i++){
      opt_vector_copy(oma->p[nbrowsa + i], omb->p[i], nbcols);
  }
  /* now we fill pp, which will contain the unsorted rows */
  nbrows = nbrowsa + omb->nbrows;
  pp = malloc(nbrows*sizeof(opt_numint_t*));
  for (i=0; i<nbrows; i++){
    pp[i] = oma->p[i];
  }
  
  /* Now we fill oma->p from pp */
  
  ia = 0;
  ib = nbrowsa;
  i = 0;
  k = 0;
  while (ia < nbrowsa && ib < nbrows){
    int res = opt_vector_compare(opk,
			     pp[ia],
			     pp[ib],nbcols);
    if (res<=0){
      oma->p[i] = pp[ia]; ia++;
      if (res==0){
	k++;
	oma->p[nbrows-k] = pp[ib]; ib++;
      }
    }
    else {
      oma->p[i] = pp[ib]; ib++;
    }
    i++;
  }
  /* Are there still constraints ? */
  while (ia < nbrowsa){
    oma->p[i] = pp[ia];
    i++; ia++;
  }
  while (ib < nbrows){
    oma->p[i] = pp[ib];
    i++; ib++;
  }
  oma->nbrows -= k;
  oma->_sorted = true;
  free(pp);
}



/* ************************************************** */
/* Resize Matrices  ******/
/* ************************************************** */

void opt_matrix_project_var(opt_pk_internal_t *opk, opt_matrix_t * oc, size_t start, unsigned short int nbcolumns, elina_dim_t *tdim, size_t dimsup){
	//size_t dimsup = dimchange->intdim+dimchange->realdim;
	size_t end = start + dimsup;
	size_t *map = (size_t *)calloc(dimsup, sizeof(size_t));
	size_t k = 0;
	unsigned short int i;
	for(i = 0; i <= nbcolumns-opk->dec; i++){
		while((k < dimsup) && (i==tdim[k])){
			map[k] = i + k;
			k++;
		}
	}

	opt_numint_t **p = oc->p;
	k = 0;
	for(i = start; i < end; i++){
		opt_numint_t *pi = p[i];
		pi[0] = 0;
		pi[1] = 0;
		size_t ind = map[k];
		pi[ind] = 1;
		k++;
	}
	free(map);
}

void opt_matrix_resize_diffcols(opt_matrix_t* oc, int diff)
{
  if (diff != 0){
    size_t i;
    for(i=0; i<oc->_maxrows; i++){
      opt_vector_realloc(&oc->p[i],
		     oc->nbcolumns,
		     oc->nbcolumns+diff);
    }
    oc->nbcolumns += diff;
  }
}

opt_matrix_t* opt_matrix_add_dimensions(opt_pk_internal_t* opk,
				bool destructive,
				opt_matrix_t* oc,
				elina_dimchange_t* dimchange, bool project)
{
  opt_matrix_t* noc;
  size_t i,dimsup;
  size_t nbrows = oc->nbrows;
  unsigned short int nbcolumns = oc->nbcolumns;
  dimsup = dimchange->intdim+dimchange->realdim;
  if (destructive){
    noc = oc;
    if(project){
	if(nbrows+dimsup > oc->_maxrows){
		opt_matrix_resize_rows(noc,nbrows+dimsup);
	}
    }
    opt_matrix_resize_diffcols(noc,(int)dimsup);
  }
  else {
    if(project){
	noc = opt_matrix_alloc(nbrows+dimsup,nbcolumns+dimsup,oc->_sorted);
    }
    else{
    	noc = opt_matrix_alloc(nbrows,nbcolumns+dimsup,oc->_sorted);
    }
  }
  for (i=0; i<nbrows; i++){
    opt_vector_add_dimensions(opk,noc->p[i],oc->p[i],noc->nbcolumns-dimsup,dimchange);
  }
  if(project){
	opt_matrix_project_var(opk,noc,nbrows,nbcolumns, dimchange->dim,dimchange->intdim+dimchange->realdim);
  }
  return noc;
}

opt_matrix_t* opt_matrix_remove_dimensions(opt_pk_internal_t* opk,
				   bool destructive,
				   opt_matrix_t* oc,
				   elina_dimchange_t* dimchange)
{
  opt_matrix_t* noc;
  size_t i,dimsup;

  dimsup = dimchange->intdim + dimchange->realdim;
  noc = 
    destructive ? 
    oc : 
    opt_matrix_alloc(oc->nbrows, oc->nbcolumns-dimsup, false);
  for (i=0; i<oc->nbrows; i++){
    opt_vector_remove_dimensions(opk,
			     noc->p[i],
			     oc->p[i],
			     oc->nbcolumns,
			     dimchange);
    opt_vector_normalize(opk,noc->p[i],oc->nbcolumns-dimsup);
  }
  if (destructive){
    opt_matrix_resize_diffcols(noc, -(int)dimsup);
  }
  noc->_sorted = false;
  return noc;
}

opt_matrix_t* opt_matrix_permute_dimensions(opt_pk_internal_t* opk,
				    bool destructive,
				    opt_matrix_t* oc,
				    elina_dim_t* permutation)
{
  opt_matrix_t* noc;
  size_t i;

  noc = destructive ? oc : opt_matrix_alloc(oc->nbrows,oc->nbcolumns,false);
  for (i=0; i<oc->nbrows; i++){
    opt_vector_permute_dimensions(opk,noc->p[i],oc->p[i],oc->nbcolumns,permutation);
  }
  noc->_sorted = false;
  return noc;
}

/*********************************************************
	Rearranges the matrix so all the equalities are at the start
**********************************************************/
void opt_matrix_rearrange(opt_matrix_t *oc, size_t nbeq){
	size_t i=0,k=0;
	size_t nbcons = oc->nbrows;
        unsigned short int nbcolumns = oc->nbcolumns;
        opt_numint_t **p = oc->p;
	while((i<nbcons) && (k < nbeq)){
		while((i < nbcons) && p[i][0]){
			i++;
		}
		if(i< nbcons){
			if(i!=k){
				opt_matrix_exch_rows(oc,i,k);
			}
			i++;
			k++;
		}
	}
}


/*********************************************************
	Rearrange the generators so all the vertices are at the start, satF can be NULL
**********************************************************/
size_t opt_generator_rearrange(opt_matrix_t *F, opt_satmat_t * satF){
	size_t i=0,k=0;
	size_t num_vertex = 0;
	size_t nbgens = F->nbrows;
	char * vmap = (char *)calloc(nbgens, sizeof(char));
	for(i=0; i < nbgens; i++){
		if(is_vertex(F->p[i])){
			vmap[i] = 1;
			num_vertex++;
		}
	}
	i=0;
        opt_numint_t **p = F->p;
        while((i<nbgens) && (k < num_vertex)){
		
		while((i < nbgens) && !vmap[i]){
			i++;
		}
		
		if(i< nbgens){
			if(i!=k){
				
				opt_matrix_exch_rows(F,i,k);
				if(satF!=NULL){
					opt_satmat_exch_cols(satF,i,k);
				}
			}
			i++;
			k++;
		}
	}
	free(vmap);
	return num_vertex;
}


/**************************************************
Assumes inequalities are between 0 and nbeq 
**************************************************/
size_t opt_matrix_gauss_elimination(opt_pk_internal_t *opk, opt_matrix_t *oc, size_t nbeq){
  size_t nbcons = oc->nbrows;
  unsigned short int nbcolumns = oc->nbcolumns;
  size_t i, j;
  unsigned short int k = 0;
  size_t rank = 0;
  opt_numint_t **p = oc->p;

  /******************************************
          Perform  Gaussian Elimination now
  ******************************************/
  for (k = opk->dec; k < nbcolumns; k++) {
    int s;
    for (i = rank; i < nbeq; i++) {
      s = elina_int_sgn(p[i][k]);
      if (s) {
        break;
      }
    }
    if (i < nbeq) {

      if (i > rank) {
        opt_matrix_exch_rows(oc, i, rank);
      }
      // if(s > 0){
      // for(j = 1; j < nbcolumns; j++){
      //	p[rank][j] = -p[rank][j];
      //}
      //}
      p[rank][0] = 0;
      for (j = i + 1; j < nbeq; j++) {
        int sj = elina_int_sgn(p[j][k]);
        if (s * sj < 0) {
          opt_matrix_combine_rows(opk, oc, j, rank, j, k, true);
        } else if (s * sj > 0) {
          opt_matrix_combine_rows(opk, oc, j, rank, j, k, false);
        }
      }
      opk->cherni_intp[rank] = k;
      rank++;
    }
	}
	return rank;
}


/* ====================================================================== */
/* III.2 Backsubstitution */
/* ====================================================================== */

/* This function backsubstitute the coefficients according to the system of
   equations and the array opk->cherni_intp properly set by
   gauss. */


void gauss_backsubstitute(opt_pk_internal_t* opk, opt_matrix_t* oc, size_t rank)
{
  size_t i,j;
  int k;
  //printf("Before substitution %d\n",rank);
  //opt_matrix_fprint(stdout,oc);
  for (k=(int)rank-1; k>=0; k--) {
    j = opk->cherni_intp[k];
    //if(oc->p[k][j]< 0){
	// size_t l;
	 //for(l=1; l < oc->nbcolumns; l++){
	//	oc->p[k][l] = -oc->p[k][l];
	 //}
    //}
    int sk = elina_int_sgn(oc->p[k][j]);
    for (i=0; i<oc->nbrows; i++) {
      int si = elina_int_sgn(oc->p[i][j]);
      if (i != (size_t)k){
        if(sk*si < 0){
		opt_matrix_combine_rows(opk,oc,i,(size_t)k,i,j,true);
	}else if(sk*si > 0){
		opt_matrix_combine_rows(opk,oc,i,(size_t)k,i,j,false);
	}
      }
    }
  }
  //printf("After substitution\n");
  //opt_matrix_fprint(stdout,oc);
}


void opt_matrix_backsubstitute(opt_pk_internal_t *opk, opt_matrix_t *oc, size_t rank){
	size_t i, j;
	int k;
	size_t nbcons = oc->nbrows;
	opt_numint_t **p = oc->p;
	for(k=(int)rank-1; k >=0; k++){
		j = opk->cherni_intp[k];
		for(i=0; i < nbcons; i++){
			if((i!=(size_t)k) && elina_int_sgn(p[i][j])){
				opt_matrix_combine_rows(opk,oc,i,(size_t)k,i,j,true);
			}
		}
	}
}

/* *************************************************************************
Extract Functionality
************************************************************************* */

/* ---------------------------------------------------------------------- */
/* Substitution by an array of equations */
/* ---------------------------------------------------------------------- */



/* ---------------------------------------------------------------------- */
/* Assignement of an expression to a variable */
/* ---------------------------------------------------------------------- */

size_t  opt_matrix_assign_variable(opt_pk_internal_t* opk,opt_matrix_t *nmat,
				 opt_matrix_t* mat,
				 elina_dim_t dim, opt_numint_t* tab)
{
  size_t i,j,var,nbline = 0;
  bool den;
 

  var = opk->dec + dim;
  den = (tab[0] > 1);
  bool flag = false;

  size_t i1 = 0;
	
  //nmat->_sorted = false;
  for (i=0; i<mat->nbrows; i++){
    /* product for var column */
    opk->matrix_prod = opt_vector_product(opk, mat->p[i],
		   tab,mat->nbcolumns);
    /* columns != var */
    opt_numint_t abs_tab = opt_numint_abs(tab[0]);
    opt_numint_t x1 = INT64_MAX / abs_tab;
    opt_numint_t x2 = INT64_MIN / abs_tab;
    /* Side-effect */
    for (j = 0; j < mat->nbcolumns; j++) {
      if (j != var) {
        if (den)
          flag = flag ||
                 opt_int64_mult(mat->p[i][j], tab[0], x1, x2, &mat->p[i][j]);
        else
          mat->p[i][j] = mat->p[i][j];
      }
      }
	
    /* var column */
      mat->p[i][var] = opk->matrix_prod;

    opt_matrix_normalize_row(opk,nmat,i);
    if(!opt_vector_is_null(opk,mat->p[i],mat->nbcolumns)){
		if(mat->p[i][0]==0){
			nbline++;
		}
		opt_vector_copy(nmat->p[i1],mat->p[i],mat->nbcolumns);
		i1++;
    }
  }
  nmat->nbrows = i1;
  if (flag) {
    printf("exception assign variable\n");
    fflush(stdout);
    opk->exn = ELINA_EXC_OVERFLOW;
    return nbline;
  }
  return nbline;
}


/* ---------------------------------------------------------------------- */
/* Assignement by an array of equations */
/* ---------------------------------------------------------------------- */

opt_matrix_t* opt_matrix_assign_variables(opt_pk_internal_t* opk,
				  opt_matrix_t* mat,
				  elina_dim_t* tdim,
				  opt_numint_t** tvec,
				  size_t size)
{
  size_t i,j,eindex;
  opt_matrix_t* nmat = _opt_matrix_alloc_int(mat->nbrows, mat->nbcolumns,false);
  opt_numint_t den;

  /* Computing common denominator */
  den = tvec[0][0];
  for (i=1; i<size; i++){
    den = den*tvec[i][0];
  }

  if (den != 1){
    /* General case */
    opt_numint_t* vden = opt_vector_alloc(size);
    for (i=0; i<size; i++){
      vden[i]=den/tvec[i][0];
    }
    /* Column 0: copy */
    for (i=0; i<mat->nbrows; i++){
      nmat->p[i][0] = mat->p[i][0];
    }
    /* Other columns */
    eindex = 0;
    for (j=1; j<mat->nbcolumns; j++){
      if (eindex < size && opk->dec + tdim[eindex] == j){
	/* We are on an assigned column */
	for (i=0; i<mat->nbrows; i++){ /* For each row */
	  opk->matrix_prod = opt_vector_product(opk, mat->p[i], tvec[eindex],mat->nbcolumns);
	  opk->matrix_prod = opk->matrix_prod*vden[eindex];
	  /* Put the result */
	  nmat->p[i][j] = opk->matrix_prod;
	}
	eindex++;
      }
      else {
	/* We are on a normal column */
	for (i=0; i<mat->nbrows; i++){ /* For each row */
	  nmat->p[i][j]=0;
	  nmat->p[i][j]=mat->p[i][j]*den;
	}
      }
    }
    opt_vector_free(vden,size);
  }
  else {
    /* Special case: all denominators are 1 */
    /* Column 0: copy */
    for (i=0; i<mat->nbrows; i++){
      nmat->p[i][0] = mat->p[i][0];
    }
    /* Other columns */
    eindex = 0;
    for (j=1; j<mat->nbcolumns; j++){
      if (eindex < size && opk->dec + tdim[eindex] == j){
	/* We are on a assigned column */
	for (i=0; i<mat->nbrows; i++){ /* For each row */
	  opk->matrix_prod = opt_vector_product(opk, mat->p[i],
			 tvec[eindex],mat->nbcolumns);
	  nmat->p[i][j] = opk->matrix_prod;
	}
	eindex++;
      }
      else {
	/* We are on a normal column */
	for (i=0; i<mat->nbrows; i++){ /* For each row */
	  nmat->p[i][j] = mat->p[i][j];
	}
      }
    }
  }
  
  for (i=0; i<mat->nbrows; i++){
    opt_matrix_normalize_row(opk,nmat,i);
  }

  return nmat;
}



/* ---------------------------------------------------------------------- */
/* Substitution of a variable by an expression */
/* ---------------------------------------------------------------------- */

/* Hypothesis:

  - either nmat is a matrix allocated with _matrix_alloc_int,
    and his coefficients are not initialized,

  - or nmat==mat
*/
opt_matrix_t* opt_matrix_substitute_variable(opt_pk_internal_t* opk,
				     bool destructive,
				     opt_matrix_t* mat,
				     elina_dim_t dim, opt_numint_t* tab)
{
  size_t i,j,var;
  bool den;
  opt_matrix_t* nmat;
  bool flag = false;
  //matrix_fprint(stdout,mat);
  var = opk->dec + dim;
  den = tab[0] > 1;
  nmat =
    destructive ?
    mat :
    _opt_matrix_alloc_int(mat->nbrows,mat->nbcolumns,false);

  nmat->_sorted = false;

  for (i=0; i<mat->nbrows; i++) {
    if (elina_int_sgn(mat->p[i][var])) {
      /* The substitution must be done */
      if (!destructive){
	/* Functional */
	nmat->p[i][0] = mat->p[i][0];
	/* columns != var */
        opt_numint_t abs_tab, x1, x2;
        abs_tab = opt_numint_abs(tab[0]);
        x1 = INT64_MAX / abs_tab;
        x2 = INT64_MIN / abs_tab;
        for (j=1; j<mat->nbcolumns; j++) {
	  if (j!=var){
	    if (den){
	      nmat->p[i][j] = 0;
              flag = flag || opt_int64_mult(mat->p[i][j], tab[0], x1, x2,
                                            &nmat->p[i][j]);
            }
	    else {
	      nmat->p[i][j] = mat->p[i][j];
	    }
            if (tab[j]) {
              opt_numint_t abs_tab_j = opt_numint_abs(tab[j]);
              opt_numint_t x1_j = INT64_MAX / abs_tab_j;
              opt_numint_t x2_j = INT64_MIN / abs_tab_j;
              flag = flag || opt_int64_mult(mat->p[i][var], tab[j], x1_j, x2_j,
                                            &opk->matrix_prod);
              // opk->matrix_prod = mat->p[i][var]*tab[j];
              flag = flag || opt_int64_add(nmat->p[i][j], opk->matrix_prod,
                                           &nmat->p[i][j]);
              // nmat->p[i][j] = nmat->p[i][j] + opk->matrix_prod;
            }
          }
	}
	/* var column */
	nmat->p[i][var]  = 0;
        opt_numint_t abs_tab_var = opt_numint_abs(tab[var]);
        opt_numint_t x1_var = INT64_MAX / abs_tab_var;
        opt_numint_t x2_var = INT64_MIN / abs_tab_var;
        flag = flag || opt_int64_mult(mat->p[i][var], tab[var], x1_var, x2_var,
                                      &nmat->p[i][var]);
        // nmat->p[i][var] = mat->p[i][var] * tab[var];
      }
      else {
	/* Side-effect */
	/* columns != var */
        opt_numint_t abs_tab = opt_numint_abs(tab[0]);
        opt_numint_t x1 = INT64_MAX / abs_tab;
        opt_numint_t x2 = INT64_MIN / abs_tab;
        for (j=1; j<mat->nbcolumns; j++) {
	  if (j!=var){
	    if (den){
              flag = flag || opt_int64_mult(nmat->p[i][j], tab[0], x1, x2,
                                            &nmat->p[i][j]);
              // nmat->p[i][j] = nmat->p[i][j]*tab[0];
            }
            if (tab[j]) {
              opt_numint_t abs_tab_j = opt_numint_abs(tab[j]);
              opt_numint_t x1_j = INT64_MAX / abs_tab_j;
              opt_numint_t x2_j = INT64_MIN / abs_tab_j;
              flag = flag || opt_int64_mult(mat->p[i][var], tab[j], x1_j, x2_j,
                                            &opk->matrix_prod);
              // opk->matrix_prod = mat->p[i][var]*tab[j];
              flag = flag || opt_int64_add(nmat->p[i][j], opk->matrix_prod,
                                           &nmat->p[i][j]);
              // nmat->p[i][j] = nmat->p[i][j] + opk->matrix_prod;
            }
          }
	}
	/* var column */
        opt_numint_t abs_tab_var = opt_numint_abs(tab[var]);
        opt_numint_t x1_var = INT64_MAX / abs_tab_var;
        opt_numint_t x2_var = INT64_MIN / abs_tab_var;
        flag = flag || opt_int64_mult(nmat->p[i][var], tab[var], x1_var, x2_var,
                                      &nmat->p[i][var]);
        // nmat->p[i][var] = nmat->p[i][var]*tab[var];
      }
      opt_matrix_normalize_row(opk,nmat,i);
    }
    else {
      /* No substitution */
      if (!destructive){
	for (j=0; j<mat->nbcolumns; j++) {
	  nmat->p[i][j] = mat->p[i][j];
	}
      }
    }
  }
  if (flag) {
    printf("exception substiture variable\n");
    fflush(stdout);
    opk->exn = ELINA_EXC_OVERFLOW;
    return nmat;
  }
  return nmat;
}


/* ---------------------------------------------------------------------- */
/* Substitution by an array of equations */
/* ---------------------------------------------------------------------- */

opt_matrix_t* opt_matrix_substitute_variables(opt_pk_internal_t* opk,
				      opt_matrix_t* mat,
				      elina_dim_t* tdim,
				      opt_numint_t** tvec,
				      size_t size)
{
  size_t i,j,eindex;
  opt_matrix_t* nmat = opt_matrix_alloc(mat->nbrows, mat->nbcolumns,false);
  opt_numint_t den;

  /* Computing common denominator */
  den = tvec[0][0];
  for (i=1; i<size; i++){
       den = den*tvec[i][0];
  }

  if (den!=1){
    /* General case */
    opt_numint_t* vden = opt_vector_alloc(size);
    for (i=0; i<size; i++){
         vden[i] = den/tvec[i][0];
    }
    /* For each row */
    for (i=0; i<mat->nbrows; i++) {
      /* Column 0 */
        nmat->p[i][0] = mat->p[i][0];
      /* Other columns */
      /* First, copy the row and sets to zero substituted variables */
      eindex = 0;
      for (j=1; j<mat->nbcolumns; j++){
	if (eindex < size && opk->dec + tdim[eindex] == j)
	  eindex++;
	else
	  nmat->p[i][j] = mat->p[i][j]*den;
      }
      /* Second, add things coming from substitution */
      for (j=1; j<mat->nbcolumns; j++){
	for (eindex=0; eindex<size; eindex++){
	  if (opt_numint_sgn(mat->p[i][opk->dec + tdim[eindex]])) {
	      opk->matrix_prod = mat->p[i][opk->dec + tdim[eindex]]*
		       tvec[eindex][j];
	      opk->matrix_prod = opk->matrix_prod * vden[eindex];
	      nmat->p[i][j] = nmat->p[i][j] + opk->matrix_prod;
	  }
	}
      }
    }
    opt_vector_free(vden,size);
  }
  else {
    /* Special case: all denominators are 1 */
    /* For each row */
    for (i=0; i<mat->nbrows; i++) {
      /* Column 0 */
      nmat->p[i][0] = mat->p[i][0];
      /* Other columns */
      /* First, copy the row and sets to zero substituted variables */
      eindex = 0;
      for (j=1; j<mat->nbcolumns; j++){
	if (eindex < size && opk->dec + tdim[eindex] == j)
	  eindex++;
	else
	  nmat->p[i][j]= mat->p[i][j];
      }
	
      /* Second, add things coming from substitution */
      for (j=1; j<mat->nbcolumns; j++){
	for (eindex=0; eindex<size; eindex++){
	  if (opt_numint_sgn(mat->p[i][opk->dec + tdim[eindex]])) {
	      opk->matrix_prod = mat->p[i][opk->dec + tdim[eindex]]*tvec[eindex][j];
		
	      nmat->p[i][j] = nmat->p[i][j] + opk->matrix_prod;
	  }
	}
      }
    }
  }

  for (i=0; i<mat->nbrows; i++){
    opt_matrix_normalize_row(opk,nmat,i);
  }

  return nmat;
}

void opt_generator_init(opt_pk_internal_t *opk, opt_matrix_t * mat, unsigned short int comp_size, size_t start){
	//assert(nbrows>comp_size);
	unsigned short int i;
	size_t end = start + comp_size;
	unsigned short int j = opk->dec;
	for(i=start; i < end; i++){
		mat->p[i][j] = 1;
		j++;
	}
	mat->p[end][0] = 1;
	mat->p[end][1] = 1;
	mat->nbrows = end+1;
}


/************************************
	Remove common generators
***********************************/

void remove_common_gen(opt_pk_internal_t *opk, opt_matrix_t *F, size_t start) {
  size_t nbrows = F->nbrows;
  size_t nb = nbrows;
  assert(start < nbrows);
  unsigned short int nbcolumns;
  nbcolumns = F->nbcolumns;
  // size_t nb = nbcons;
  // size_t nbeq = o->nbeq;
  char *rmap = (char *)calloc(nbrows - start, sizeof(char));
  char *map = (char *)calloc(nbrows, sizeof(char));
  int i, j;
  // opt_matrix_fprint(stdout,F);
  for (i = start; i < nbrows; i++) {
    opt_numint_t *pi = F->p[i];
    unsigned short int ind = 0;
    for (j = 0; j < start; j++) {
      opt_numint_t *pj = F->p[j];
      if (!map[j] && opt_vector_equal(opk, pi, pj, nbcolumns, &ind)) {
        rmap[i - start] = 1;
        map[j] = 1;
        break;
      }
    }
    for (j = i + 1; j < nbrows; j++) {
      opt_numint_t *pj = F->p[j];
      if (!map[j] && opt_vector_equal(opk, pi, pj, nbcolumns, &ind)) {
        rmap[i - start] = 1;
        map[j] = 1;
        break;
      }
    }
  }
  j = nbrows - 1;
  for (i = start; i < nb; i++) {
    if (rmap[i - start]) {
      nbrows--;
      while ((j > i) && rmap[j - start]) {
        j--;
      }
      if (j > i) {
        opt_matrix_exch_rows(F, i, j);
      }
      j--;
    }
  }
  free(map);
  free(rmap);
  F->nbrows = nbrows;
}

void remove_common_gen_upto(opt_pk_internal_t *opk, opt_matrix_t *F,
                            size_t start, size_t end) {
  size_t nbrows = F->nbrows;
  size_t nb = end;
  assert(start < nbrows && end < nbrows);
  unsigned short int nbcolumns;
  nbcolumns = F->nbcolumns;
  // size_t nb = nbcons;
  // size_t nbeq = o->nbeq;
  char *rmap = (char *)calloc(end - start, sizeof(char));
  char *map = (char *)calloc(end, sizeof(char));
  int i, j;
  // opt_matrix_fprint(stdout,F);
  for (i = start; i < end; i++) {
    opt_numint_t *pi = F->p[i];
    unsigned short int ind = 0;
    for (j = 0; j < start; j++) {
      opt_numint_t *pj = F->p[j];
      if (!map[j] && opt_vector_equal(opk, pi, pj, nbcolumns, &ind)) {
        rmap[i - start] = 1;
        map[j] = 1;
        break;
      }
    }
    for (j = i + 1; j < end; j++) {
      opt_numint_t *pj = F->p[j];
      if (!map[j] && opt_vector_equal(opk, pi, pj, nbcolumns, &ind)) {
        printf("equal %d %d\n", i, j);
        rmap[i - start] = 1;
        map[j] = 1;
        break;
      }
    }
  }
  j = end - 1;
  for (i = start; i < nb; i++) {
    if (rmap[i - start]) {

      nbrows--;
      while ((j > i) && rmap[j - start]) {
        j--;
      }
      if (j > i) {
        opt_matrix_exch_rows(F, i, j);
      }
      j--;
    }
  }
  free(map);
  free(rmap);
  F->nbrows = nbrows;
}

/***********************
	Compute the bounds for a variable
************************/
void opt_generator_bound_dimension(opt_pk_internal_t* opk,
			    elina_interval_t * interval,
			    elina_dim_t dim,
			    opt_matrix_t* of)
{
  size_t i, index;
  int sgn;
  assert(opk->dec+dim<of->nbcolumns);
  elina_rat_t *inf = (elina_rat_t *)malloc(sizeof(elina_rat_t));
  elina_rat_set_infty(inf,1);
  elina_rat_t *sup = (elina_rat_t *)malloc(sizeof(elina_rat_t));
  elina_rat_set_infty(sup,-1);
  index = opk->dec+dim;
  for (i=0; i<of->nbrows; i++){
    if (!opk->strict || elina_int_sgn(of->p[i][opt_polka_eps])==0){
      sgn = elina_int_sgn(of->p[i][index]);
      if (!elina_int_sgn(of->p[i][0])){
	/* line: result should be zero */
	if (sgn){
	  elina_interval_set_top(interval);
	  return;
	}
      }
      else if (!elina_int_sgn(of->p[i][opt_polka_cst])){
	/* ray */
	if (sgn > 0){
	  elina_rat_set_infty(sup,+1);
	  if (elina_rat_infty(inf) && elina_rat_sgn(inf)<0){
		elina_scalar_set_elina_rat(interval->inf,inf);
  		elina_scalar_set_elina_rat(interval->sup,sup);
		return;
	  }
	    
	}
	else if (sgn < 0){
	  elina_rat_set_infty(inf,-1);
	  if (elina_rat_infty(sup) && elina_rat_sgn(sup)>0){
	    elina_scalar_set_elina_rat(interval->inf,inf);
  	    elina_scalar_set_elina_rat(interval->sup,sup);
	    return;
	  }
	}
      }
      else {
	/* point */
	
	elina_rat_set_int2(opk->poly_numrat,
			   of->p[i][index],
			   of->p[i][opt_polka_cst]);
	if (elina_rat_cmp(sup,opk->poly_numrat)<0){
	  sup->n = opk->poly_numrat->n;
	  sup->d = opk->poly_numrat->d;
	}
	//numrat_neg(opk->poly_numrat,opk->poly_numrat);
	if (elina_rat_cmp(inf,opk->poly_numrat)>0){
	  inf->n = opk->poly_numrat->n;
	  inf->d = opk->poly_numrat->d;
	}
      }	  
    }
  }
  elina_scalar_set_elina_rat(interval->inf,inf);
  elina_scalar_set_elina_rat(interval->sup,sup);
  free(inf);
  free(sup);
}

elina_interval_t ** opt_generator_to_box(opt_pk_internal_t* opk,
		     opt_matrix_t* of)
{
  unsigned short int i,dim;
  elina_interval_t** res;

  assert(of);
  assert(of->nbcolumns>=opk->dec);
  dim = of->nbcolumns - opk->dec;
  res = elina_interval_array_alloc(dim);
  for (i=0;i<dim;i++){
    opt_generator_bound_dimension(opk,res[i],i,of);
  }
  return res;
}

elina_interval_t** opt_matrix_to_box(opt_pk_internal_t* opk,
		     opt_matrix_t* F)
{
  size_t i,dim;
  elina_interval_t** res;

  assert(F);
  assert(F->nbcolumns>=opk->dec);
  dim = F->nbcolumns - opk->dec;
  res = elina_interval_array_alloc(dim);
  for (i=0;i<dim;i++){
    opt_generator_bound_dimension(opk,res[i],i,F);
  }
  return res;
}


/************************
	Compute the bounds for a linear expression
*************************/
void opt_generator_bound_elina_linexpr0(opt_pk_internal_t *opk, elina_rat_t *inf, elina_rat_t *sup,
				     elina_linexpr0_t *expr, opt_matrix_t * F){
  size_t i;
  int sgn;
  elina_rat_set_infty(inf,1);
  elina_rat_set_infty(sup,-1);
  elina_rat_t *prod1 = (elina_rat_t *)malloc(sizeof(elina_rat_t));
  elina_rat_t *prod2 = (elina_rat_t *)malloc(sizeof(elina_rat_t));
  for (i=0; i<F->nbrows; i++){
    if (!opk->strict || opt_numint_sgn(F->p[i][opt_polka_eps])==0 ){
      opt_vector_bound_elina_linexpr0(opk, prod1, prod2, expr, F->p[i], F->nbcolumns);
      if (opt_numint_sgn(F->p[i][0])==0){
	// line: result should be zero 
	if (elina_rat_sgn(prod1)!=0 || elina_rat_sgn(prod2)!=0){
	  elina_rat_set_infty(inf,-1);
	  elina_rat_set_infty(sup,1);
	  goto opt_generator_bound_elina_linexpr0_exit;
	}
      }
      else if (opt_numint_sgn(F->p[i][opt_polka_cst])==0){
	// ray 
	if (elina_rat_sgn(prod1)!=0 || elina_rat_sgn(prod2)!=0){
	  if (elina_rat_sgn(prod1)>0){
	    // [inf,sup]>0 
	    elina_rat_set_infty(sup,+1);
	    if (elina_rat_infty(inf) && elina_rat_sgn(inf)<0)
	      goto opt_generator_bound_elina_linexpr0_exit;
	  }
	  else if (elina_rat_sgn(prod2)<0){
	    // [inf,sup]<0 
	    elina_rat_set_infty(inf,-1);
	    if (elina_rat_infty(sup) && elina_rat_sgn(sup)>0)
	      goto opt_generator_bound_elina_linexpr0_exit;
	  }
	  else {
	    elina_rat_set_infty(inf,-1);
	    elina_rat_set_infty(sup,1);
	    goto opt_generator_bound_elina_linexpr0_exit;
	  }
	}
      }
      else {
        elina_rat_min(inf,inf,prod1);
	elina_rat_max(sup,sup,prod2);
      }
    }
  }
 opt_generator_bound_elina_linexpr0_exit:
  free(prod1);
  free(prod2);
  return;

}


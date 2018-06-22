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

/* ************************************************************************* */
/* elina_linexpr0.c: linear expressions */
/* ************************************************************************* */

#include "elina_linexpr0.h"

#include <stdarg.h>

/* ====================================================================== */
/* I. Auxilliary functions */
/* ====================================================================== */

size_t elina_linexpr0_size(elina_linexpr0_t* expr)
  { return expr->size; }



elina_coeff_t* elina_linexpr0_cstref(elina_linexpr0_t* expr)
  { return &expr->cst; }


void elina_linexpr0_get_cst(elina_coeff_t* coeff, elina_linexpr0_t* expr)
  { elina_coeff_set(coeff,&expr->cst); }


void elina_linexpr0_set_cst(elina_linexpr0_t* expr, elina_coeff_t* cst)
  { elina_coeff_set(&expr->cst,cst); }


void elina_linexpr0_set_cst_scalar(elina_linexpr0_t* expr, elina_scalar_t* scalar)
  { elina_coeff_set_scalar(&expr->cst, scalar); }


void elina_linexpr0_set_cst_scalar_int(elina_linexpr0_t* expr, int num)
  { elina_coeff_set_scalar_int(&expr->cst, num); }


void elina_linexpr0_set_cst_scalar_frac(elina_linexpr0_t* expr, long int num, unsigned long int den)
  { elina_coeff_set_scalar_frac(&expr->cst, num, den); }


void elina_linexpr0_set_cst_scalar_double(elina_linexpr0_t* expr, double num)
  { elina_coeff_set_scalar_double(&expr->cst, num); }


void elina_linexpr0_set_cst_interval(elina_linexpr0_t* expr, elina_interval_t* itv)
  { elina_coeff_set_interval(&expr->cst, itv); }


void elina_linexpr0_set_cst_interval_int(elina_linexpr0_t* expr, int inf, int sup)
  { elina_coeff_set_interval_int(&expr->cst, inf,sup); }


void elina_linexpr0_set_cst_interval_scalar(elina_linexpr0_t* expr, elina_scalar_t* inf, elina_scalar_t* sup)
  { elina_coeff_set_interval_scalar(&expr->cst, inf,sup); }


void elina_linexpr0_set_cst_interval_frac(elina_linexpr0_t* expr,
				  long int numinf, unsigned long int deninf,
				  long int numsup, unsigned long int densup)
  { elina_coeff_set_interval_frac(&expr->cst, numinf,deninf, numsup,densup); }


void elina_linexpr0_set_cst_interval_double(elina_linexpr0_t* expr, double inf, double sup)
  { elina_coeff_set_interval_double(&expr->cst, inf,sup); }


bool elina_linexpr0_set_coeff(elina_linexpr0_t* expr, elina_dim_t dim, elina_coeff_t* coeff)
  { elina_coeff_t* ecoeff = elina_linexpr0_coeffref(expr,dim); if (ecoeff){elina_coeff_set(ecoeff,coeff); return false;} else return true; }


bool elina_linexpr0_set_coeff_scalar(elina_linexpr0_t* expr, elina_dim_t dim, elina_scalar_t* scalar)
  { elina_coeff_t* ecoeff = elina_linexpr0_coeffref(expr,dim); if (ecoeff){ elina_coeff_set_scalar(ecoeff,scalar); return false; } else return true; }


bool elina_linexpr0_set_coeff_scalar_int(elina_linexpr0_t* expr, elina_dim_t dim, int num)
  { elina_coeff_t* ecoeff = elina_linexpr0_coeffref(expr,dim); if (ecoeff){ elina_coeff_set_scalar_int(ecoeff,num); return false; } else return true; }


bool elina_linexpr0_set_coeff_scalar_frac(elina_linexpr0_t* expr, elina_dim_t dim, long int num, unsigned long int den)
  { elina_coeff_t* ecoeff = elina_linexpr0_coeffref(expr,dim); if (ecoeff){ elina_coeff_set_scalar_frac(ecoeff,num, den); return false; } else return true; }


bool elina_linexpr0_set_coeff_scalar_double(elina_linexpr0_t* expr, elina_dim_t dim, double num)
  { elina_coeff_t* ecoeff = elina_linexpr0_coeffref(expr,dim); if (ecoeff){ elina_coeff_set_scalar_double(ecoeff,num); return false; } else return true; }


bool elina_linexpr0_set_coeff_interval(elina_linexpr0_t* expr, elina_dim_t dim, elina_interval_t* itv)
  { elina_coeff_t* ecoeff = elina_linexpr0_coeffref(expr,dim); if (ecoeff){ elina_coeff_set_interval(ecoeff,itv); return false; } else return true; }


bool elina_linexpr0_set_coeff_interval_int(elina_linexpr0_t* expr, elina_dim_t dim, int inf, int sup)
  { elina_coeff_t* ecoeff = elina_linexpr0_coeffref(expr,dim); if (ecoeff){ elina_coeff_set_interval_int(ecoeff,inf,sup); return false; } else return true; }


bool elina_linexpr0_set_coeff_interval_scalar(elina_linexpr0_t* expr, elina_dim_t dim, elina_scalar_t* inf, elina_scalar_t* sup)
  { elina_coeff_t* ecoeff = elina_linexpr0_coeffref(expr,dim); if (ecoeff){ elina_coeff_set_interval_scalar(ecoeff,inf,sup); return false; } else return true; }


bool elina_linexpr0_set_coeff_interval_frac(elina_linexpr0_t* expr, elina_dim_t dim,
				  long int numinf, unsigned long int deninf,
				  long int numsup, unsigned long int densup)
  { elina_coeff_t* ecoeff = elina_linexpr0_coeffref(expr,dim); if (ecoeff){ elina_coeff_set_interval_frac(ecoeff,numinf,deninf, numsup,densup); return false; } else return true; }


bool elina_linexpr0_set_coeff_interval_double(elina_linexpr0_t* expr, elina_dim_t dim, double inf, double sup)
  { elina_coeff_t* ecoeff = elina_linexpr0_coeffref(expr,dim); if (ecoeff){ elina_coeff_set_interval_double(ecoeff,inf,sup); return false; } else return true; }

/* ====================================================================== */
/* II. Memory management and printing */
/* ====================================================================== */

elina_linexpr0_t* elina_linexpr0_alloc(elina_linexpr_discr_t lin_discr, size_t size)
{
  elina_linexpr0_t* e;
  size_t i;

  e = (elina_linexpr0_t*)malloc(sizeof(elina_linexpr0_t));
  elina_coeff_init(&e->cst,ELINA_COEFF_SCALAR);
  e->discr = lin_discr;
  e->size = size;
  switch (lin_discr){
  case ELINA_LINEXPR_DENSE:
    e->p.coeff = size==0 ? NULL : (elina_coeff_t*)malloc(size*sizeof(elina_coeff_t));
    for (i=0; i<size; i++) elina_coeff_init(&e->p.coeff[i],ELINA_COEFF_SCALAR);
    break;
  case ELINA_LINEXPR_SPARSE:
    e->p.linterm = size==0 ? NULL : (elina_linterm_t*)malloc(size*sizeof(elina_linterm_t));
    for (i=0; i<size; i++){
      elina_coeff_init(&e->p.linterm[i].coeff,ELINA_COEFF_SCALAR);
      e->p.linterm[i].dim = ELINA_DIM_MAX;
    }
    break;
  }
  return e;
}

void elina_linexpr0_realloc(elina_linexpr0_t* e, size_t size)
{
  size_t i;
  if (e->size==size)
    return;
  switch (e->discr){
  case ELINA_LINEXPR_DENSE:
    for (i=size; i<e->size; i++){
      elina_coeff_clear(&e->p.coeff[i]);
    }
    if (size==0){
      free(e->p.coeff);
      e->p.coeff = NULL;
    }
    else {
      e->p.coeff = realloc(e->p.coeff,size*sizeof(elina_coeff_t));
      for (i=e->size; i<size; i++){
	elina_coeff_init(&e->p.coeff[i],ELINA_COEFF_SCALAR);
      }
    }
    break;
  case ELINA_LINEXPR_SPARSE:
    for (i=size; i<e->size; i++){
      elina_coeff_clear(&e->p.linterm[i].coeff);
    }
    if (size==0){
      free(e->p.linterm);
      e->p.linterm = NULL;
    }
    else {
      e->p.linterm = realloc(e->p.linterm, size*sizeof(elina_linterm_t));
      for (i=e->size; i<size; i++){
	elina_coeff_init(&e->p.linterm[i].coeff,ELINA_COEFF_SCALAR);
	e->p.linterm[i].dim = ELINA_DIM_MAX;
      }
    }
    break;
  }
  e->size = size;
}

void elina_linexpr0_minimize(elina_linexpr0_t* e)
{
  size_t i,j,nsize;

  if (e->discr==ELINA_LINEXPR_DENSE){
    for (i=0; i<e->size; i++){
      elina_coeff_reduce(&e->p.coeff[i]);
    }
  }
  else {
    nsize = 0;
    for (i=0; i<e->size; i++){
      elina_coeff_t* p = &e->p.linterm[i].coeff;
      elina_coeff_reduce(p);
      if (!elina_coeff_zero(p) && e->p.linterm[i].dim!=ELINA_DIM_MAX)
	nsize++;
    }
    if (nsize!=e->size){
      elina_linterm_t* linterm = malloc(nsize*sizeof(elina_linterm_t));
      j = 0;
      for (i=0; i<e->size; i++){
	elina_coeff_t* p = &e->p.linterm[i].coeff;
	if (!elina_coeff_zero(p) && e->p.linterm[i].dim!=ELINA_DIM_MAX){
	  linterm[j] = e->p.linterm[i];
	  j++;
	}
	else
	  elina_coeff_clear(p);
      }
      free(e->p.linterm);
      e->p.linterm = linterm;
      e->size = nsize;
    }
  }
}

elina_linexpr0_t* elina_linexpr0_copy(elina_linexpr0_t* a)
{
  elina_linexpr0_t* e;
  size_t i;

  e = (elina_linexpr0_t*)malloc(sizeof(elina_linexpr0_t));
  elina_coeff_init_set(&e->cst,&a->cst);
  e->discr = a->discr;
  e->size = a->size;
  switch (a->discr){
  case ELINA_LINEXPR_DENSE:
    e->p.coeff = a->size==0 ? NULL : (elina_coeff_t*)malloc(a->size*sizeof(elina_coeff_t));
    for (i=0; i<a->size; i++)
      elina_coeff_init_set(&e->p.coeff[i],&a->p.coeff[i]);
    break;
  case ELINA_LINEXPR_SPARSE:
    e->p.linterm = a->size==0 ? NULL : (elina_linterm_t*)malloc(a->size*sizeof(elina_linterm_t));
    for (i=0; i<a->size; i++){
      elina_coeff_init_set(&e->p.linterm[i].coeff,&a->p.linterm[i].coeff);
      e->p.linterm[i].dim = a->p.linterm[i].dim;
    }
    break;
  }
  e->size = a->size;
  return e;
}

static
void elina_linexpr0_clear(elina_linexpr0_t* e)
{
  size_t i;

  if(e!=NULL){
    elina_coeff_clear(&e->cst);
    switch (e->discr){
    case ELINA_LINEXPR_DENSE:
      if (e->p.coeff!=NULL){
	for (i=0; i<e->size; i++) elina_coeff_clear(&e->p.coeff[i]);
	free(e->p.coeff);
      }
      break;
    case ELINA_LINEXPR_SPARSE:
      if (e->p.linterm!=NULL){
	for (i=0; i<e->size; i++) elina_coeff_clear(&e->p.linterm[i].coeff);
	free(e->p.linterm);
      }
      break;
    }
  }
}

void elina_linexpr0_free(elina_linexpr0_t* e)
{
  elina_linexpr0_clear(e);
  free(e);
}

void elina_linexpr0_print(elina_linexpr0_t* a, char** name_of_dim)
{ elina_linexpr0_fprint(stdout,a,name_of_dim); }



void elina_linexpr0_fprint(FILE* stream, elina_linexpr0_t* a, char** name_of_dim)
{
  size_t i;
  elina_scalar_t* pscalar = 0;
  elina_scalar_t* scalar;
  elina_coeff_t* coeff;
  elina_dim_t dim;
  bool first;
  int sgn;

  scalar = elina_scalar_alloc();

  first = true;
  elina_linexpr0_ForeachLinterm(a,i,dim,coeff){
    if (! elina_coeff_zero(coeff)){
      switch(coeff->discr){
      case ELINA_COEFF_SCALAR:
	pscalar = coeff->val.scalar;
	sgn = elina_scalar_sgn(pscalar);
	if (sgn > 0){
	  elina_scalar_set(scalar,pscalar);
	  if (!first)
	    fprintf(stream," + ");
	} else {
	  elina_scalar_neg(scalar,pscalar);
	  fprintf(stream, first ? "-" : " - ");
	}
	if (!elina_scalar_equal_int(scalar,1))
	  elina_scalar_fprint(stream,scalar);
	break;
      case ELINA_COEFF_INTERVAL:
	if (!first)
	  fprintf(stream," + ");
	elina_interval_fprint(stream,coeff->val.interval);
	break;
      }
      if (name_of_dim)
	fprintf(stream,"%s",name_of_dim[dim]);
      else
	fprintf(stream,"x%lu",(unsigned long)dim);
      first = false;
    }
  }
  /* Constant */
  if (first || !elina_coeff_zero(&a->cst)){
    switch (a->cst.discr){
    case ELINA_COEFF_SCALAR:
      pscalar = a->cst.val.scalar;
      sgn = elina_scalar_sgn(pscalar);
      if (sgn >= 0){
	elina_scalar_set(scalar,pscalar);
	if (!first)
	  fprintf(stream," + ");
      } else {
	elina_scalar_neg(scalar,pscalar);
	fprintf(stream, first ? "-" : " - ");
      }
      elina_scalar_fprint(stream,scalar);
      break;
    case ELINA_COEFF_INTERVAL:
      if (!first)
	fprintf(stream," + ");
      elina_interval_fprint(stream, a->cst.val.interval);
      break;
    }
  }
  elina_scalar_free(scalar);
}

/* ====================================================================== */
/* III. Tests */
/* ====================================================================== */

/* Does the linear expression involve only real variables ? */
bool elina_linexpr0_is_real(elina_linexpr0_t* a, size_t intdim)
{
  elina_coeff_t* coeff;
  elina_dim_t dim;
  size_t i;
  elina_linexpr0_ForeachLinterm(a,i,dim,coeff){
    if (dim>=intdim) return true;
    if (! elina_coeff_zero(coeff)) return false;
  }
  return true;
}

/* Does the linear expression involve only integer variables ? */
bool elina_linexpr0_is_integer(elina_linexpr0_t* a, size_t intdim)
{
  size_t i;
  elina_coeff_t* coeff;
  elina_dim_t dim;

  switch (a->discr){
  case ELINA_LINEXPR_DENSE:
    for (dim=intdim; dim<a->size; dim++){
      if (! elina_coeff_zero(&a->p.coeff[dim]))
	return false;
    }
    break;
  case ELINA_LINEXPR_SPARSE:
    elina_linexpr0_ForeachLinterm(a,i,dim,coeff){
      if (dim>=intdim && ! elina_coeff_zero(coeff))
	return false;
    }
    break;
  default:
    abort();
  }
  return true;
}

/* Return the rtype of an expression */
elina_linexpr_type_t elina_linexpr0_type(elina_linexpr0_t* expr)
{
  size_t i;
  elina_dim_t dim;
  elina_coeff_t* coeff;
  elina_linexpr_type_t res = ELINA_LINEXPR_LINEAR;
  elina_linexpr0_ForeachLinterm(expr,i,dim,coeff){
    if (coeff->discr!=ELINA_COEFF_SCALAR){
      res = ELINA_LINEXPR_INTLINEAR;
      break;
    }
  }
  if (res==ELINA_LINEXPR_LINEAR && expr->cst.discr!=ELINA_COEFF_SCALAR){
    res = ELINA_LINEXPR_QUASILINEAR;
  }
  return res;
}
/* Return true iff all involved coefficients are scalars */
bool elina_linexpr0_is_linear(elina_linexpr0_t* expr)
{
  size_t i;
  elina_dim_t dim;
  elina_coeff_t* coeff;
  bool res;
  res = (expr->cst.discr==ELINA_COEFF_SCALAR);
  if (res){
    elina_linexpr0_ForeachLinterm(expr,i,dim,coeff){
      res = (coeff->discr==ELINA_COEFF_SCALAR);
      if (!res) break;
    }
  }
  return res;
}
/* Return true iff all involved coefficients but the constant are scalars */
bool elina_linexpr0_is_quasilinear(elina_linexpr0_t* expr)
{
  size_t i;
  elina_dim_t dim;
  elina_coeff_t* coeff;
  bool res;

  res = true;
  elina_linexpr0_ForeachLinterm(expr,i,dim,coeff){
    res = (coeff->discr==ELINA_COEFF_SCALAR);
    if (!res) break;
  }
  return res;
}
elina_linexpr_type_t elina_linexpr0_array_type(elina_linexpr0_t** texpr, size_t size)
{
  size_t i;
  elina_linexpr_type_t res = ELINA_LINEXPR_LINEAR;
  for (i=0; i<size; i++){
    elina_linexpr_type_t type = elina_linexpr0_type(texpr[i]);
    if (type<res) res = type;
    if (res==ELINA_LINEXPR_INTLINEAR) break;
  }
  return res;
}
bool elina_linexpr0_array_is_linear(elina_linexpr0_t** texpr, size_t size)
{
  size_t i;
  bool res = true;
  for (i=0; i<size; i++){
    res = elina_linexpr0_is_linear(texpr[i]);
    if (!res) break;
  }
  return res;
}
bool elina_linexpr0_array_is_quasilinear(elina_linexpr0_t** texpr, size_t size)
{
  size_t i;
  bool res = true;
  for (i=0; i<size; i++){
    res = elina_linexpr0_is_quasilinear(texpr[i]);
    if (!res) break;
  }
  return res;
}


/* ====================================================================== */
/* IV. Access */
/* ====================================================================== */

static size_t index_of_or_after_dim(elina_dim_t dim, elina_linterm_t* linterm, size_t size)
{
  if (size<=10){
    size_t i;
    for (i=0; i<size; i++){
      if (dim<=linterm[i].dim)
	return i;
    }
    return size;
  }
  else {
    size_t mid = size/2;
    if (dim==linterm[mid].dim)
      return mid;
    else if (dim<linterm[mid].dim)
      return index_of_or_after_dim(dim, linterm, mid);
    else
      return mid + 1 + index_of_or_after_dim(dim, linterm+mid+1, size-mid-1);
  }
}

elina_coeff_t* elina_linexpr0_coeffref(elina_linexpr0_t* expr, elina_dim_t dim)
{
  size_t index;

  switch(expr->discr){
  case ELINA_LINEXPR_DENSE:
    return (dim<expr->size) ? &expr->p.coeff[dim] : NULL;
  case ELINA_LINEXPR_SPARSE:
    if (dim==ELINA_DIM_MAX) return NULL;
    index = index_of_or_after_dim(dim,expr->p.linterm,expr->size);
    if (index>=expr->size || dim!=expr->p.linterm[index].dim){
      if (index<expr->size && expr->p.linterm[index].dim==ELINA_DIM_MAX){
	/* We have a free linterm at the right place */
	expr->p.linterm[index].dim=dim;
	return &expr->p.linterm[index].coeff;
      }
      if (expr->size==0 || expr->p.linterm[expr->size-1].dim!=ELINA_DIM_MAX){
	/* We have to insert a new linterm at the end */
	elina_linexpr0_realloc(expr, expr->size+1);
      }
      /* We insert a linterm with ELINA_DIM_MAX at the right place */
      if (index<expr->size-1){
	elina_linterm_t tmp = expr->p.linterm[expr->size-1];
	memmove(&expr->p.linterm[index+1], &expr->p.linterm[index],
		(expr->size-index-1)*sizeof(elina_linterm_t));
	expr->p.linterm[index] = tmp;
      }
      expr->p.linterm[index].dim = dim;
    }
    return &expr->p.linterm[index].coeff;
  default:
    abort();
  }
}

/* If dense representation, coefficients should be already present, otherwise
   undefined behaviour */
bool elina_linexpr0_get_coeff(elina_coeff_t* coeff, elina_linexpr0_t* expr, elina_dim_t dim)
{
  size_t index;
  switch(expr->discr){
  case ELINA_LINEXPR_DENSE:
    if (dim<expr->size){
      elina_coeff_set(coeff,&expr->p.coeff[dim]);
      return false;
    } else {
      return true;
    }
  case ELINA_LINEXPR_SPARSE:
    if (dim==ELINA_DIM_MAX)
      return true;
    else {
      index = index_of_or_after_dim(dim,expr->p.linterm,expr->size);
      if (index<expr->size && dim==expr->p.linterm[index].dim){
	elina_coeff_set(coeff, &expr->p.linterm[index].coeff);
      }
      else {
	elina_coeff_set_scalar_double(coeff,0.0);
      }
      return false;
    }
  default:
    abort();
  }
}

bool elina_linexpr0_set_list_generic(elina_coeff_t* (*get_pcoeff)(void* expr, bool cst, va_list* va),
				  void* expr, va_list* va)
{
  bool cst;
  elina_coeff_t* coeff;
  elina_coeff_t* pcoeff;
  elina_scalar_t *scalar,*scalar1,*scalar2;
  elina_interval_t* interval;
  int num,num1,num2,den,den1,den2;
  double k,k1,k2;
  MP_RAT *mpq,*mpq1,*mpq2;
  __mpfr_struct *mpfr,*mpfr1,*mpfr2;
  elina_coefftag_t tag;

  while (true){
    tag = va_arg(*va,elina_coefftag_t);
    if (tag==ELINA_END)
      break;

    cst = (tag>=ELINA_CST);
    if (cst)
      tag-=ELINA_CST;

    switch (tag){
    case ELINA_COEFF:
      coeff = va_arg(*va,elina_coeff_t*);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set(pcoeff,coeff);
      break;
    case ELINA_COEFF_S:
      scalar = va_arg(*va,elina_scalar_t*);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set_scalar(pcoeff,scalar);
      break;
    case ELINA_COEFF_S_MPQ:
      mpq = va_arg(*va,MP_RAT*);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set_scalar_mpq(pcoeff,mpq);
      break;
    case ELINA_COEFF_S_MPFR:
      mpfr = va_arg(*va,__mpfr_struct*);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set_scalar_mpfr(pcoeff,mpfr);
      break;
    case ELINA_COEFF_S_INT:
      num = va_arg(*va,int);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set_scalar_int(pcoeff,num);
      break;
    case ELINA_COEFF_S_FRAC:
      num = va_arg(*va,int);
      den = va_arg(*va,int);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set_scalar_frac(pcoeff,num,den);
      break;
    case ELINA_COEFF_S_DOUBLE:
      k = va_arg(*va,double);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set_scalar_double(pcoeff,k);
      break;
    case ELINA_COEFF_I:
      interval = va_arg(*va,elina_interval_t*);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set_interval(pcoeff,interval);
      break;
    case ELINA_COEFF_I_SCALAR:
      scalar1 = va_arg(*va,elina_scalar_t*);
      scalar2 = va_arg(*va,elina_scalar_t*);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set_interval_scalar(pcoeff,scalar1,scalar2);
      break;
    case ELINA_COEFF_I_MPQ:
      mpq1 = va_arg(*va,MP_RAT*);
      mpq2 = va_arg(*va,MP_RAT*);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set_interval_mpq(pcoeff,mpq1,mpq2);
      break;
    case ELINA_COEFF_I_MPFR:
      mpfr1 = va_arg(*va,__mpfr_struct*);
      mpfr2 = va_arg(*va,__mpfr_struct*);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set_interval_mpfr(pcoeff,mpfr1,mpfr2);
      break;
    case ELINA_COEFF_I_INT:
      num1 = va_arg(*va,int);
      num2 = va_arg(*va,int);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set_interval_int(pcoeff,num1,num2);
      break;
    case ELINA_COEFF_I_FRAC:
      num1 = va_arg(*va,int);
      den1 = va_arg(*va,int);
      num2 = va_arg(*va,int);
      den2 = va_arg(*va,int);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set_interval_frac(pcoeff,num1,den1,num2,den2);
      break;
    case ELINA_COEFF_I_DOUBLE:
      k1 = va_arg(*va,double);
      k2 = va_arg(*va,double);
      pcoeff = get_pcoeff(expr,cst,va);
      if (pcoeff==NULL) return true;
      elina_coeff_set_interval_double(pcoeff,k1,k2);
      break;
    default:
      fprintf(stderr,
	      "elina_linexpr0_set_list_generic: probably bad structure for the argument list\n");
      abort();
    }
  }
  return false;
}

elina_coeff_t* elina_linexpr0_set_list_get_pcoeff(void* expr, bool cst, va_list* va)
{
  elina_coeff_t* pcoeff;
  if (cst){
    pcoeff = elina_linexpr0_cstref(expr);
  } else {
    elina_dim_t dim = va_arg(*va,elina_dim_t);
    pcoeff = elina_linexpr0_coeffref(expr,dim);
  }
  return pcoeff;
}

bool elina_linexpr0_set_list(elina_linexpr0_t* expr, ...)
{
  va_list va;
  bool res;
  va_start(va,expr);
  res = elina_linexpr0_set_list_generic(elina_linexpr0_set_list_get_pcoeff,
				     expr,&va);
  va_end(va);
  return res;
}

/* ====================================================================== */
/* V. Change of dimensions and permutations */
/* ====================================================================== */

/* Change current environment with a super-environment */

static int elina_linterm_cmp(const void* a, const void* b)
{
  const elina_linterm_t* aa = (const elina_linterm_t*)a;
  const elina_linterm_t* bb = (const elina_linterm_t*)b;
  return
    (aa->dim > bb->dim) ? 1 :
    ( (aa->dim < bb->dim) ? -1 : 0 );
}

elina_linexpr0_t*
elina_linexpr0_add_dimensions(elina_linexpr0_t* expr,
			   elina_dimchange_t* dimchange)
{
  elina_linexpr0_t* nexpr;

  if (expr==NULL) return NULL;
  nexpr = elina_linexpr0_copy(expr);
  elina_linexpr0_add_dimensions_with(nexpr,dimchange);
  return nexpr;
}

void
elina_linexpr0_add_dimensions_with(elina_linexpr0_t* expr,
				elina_dimchange_t* dimchange)
{
  if (expr==NULL) return;
  switch(expr->discr){
  case ELINA_LINEXPR_SPARSE:
    {
      size_t i,k,dimsup;
      dimsup = dimchange->intdim+dimchange->realdim;
      k=0;
      for (i=0; i<expr->size; i++){
	elina_dim_t* pdim = &expr->p.linterm[i].dim;
	if (*pdim==ELINA_DIM_MAX)
	  break;
	while (k<dimsup && *pdim>=dimchange->dim[k]){
	  k++;
	}
	*pdim += k;
      }
    }
    break;
  case ELINA_LINEXPR_DENSE:
    {
      int i,k;
      size_t size,dimsup;

      size = expr->size;
      dimsup = dimchange->intdim+dimchange->realdim;
      elina_linexpr0_realloc(expr,
			  size+dimsup);
      k = dimsup;
      for (i=size; i>=0; i--){
	if (i<(int)size){
	  elina_coeff_set(&expr->p.coeff[i+k],&expr->p.coeff[i]);
	}
	while (k>=1 && dimchange->dim[k-1]==(elina_dim_t)i){
	  k--;
	  elina_coeff_set_scalar_double(&expr->p.coeff[i+k],0.0);
	}
      }
    }
    break;
  default:
    abort();
  }
}

elina_linexpr0_t*
elina_linexpr0_permute_dimensions(elina_linexpr0_t* expr,
			       elina_dimperm_t* perm)
{
  if (expr==NULL) return NULL;
  elina_linexpr0_t* nexpr = elina_linexpr0_copy(expr);
  switch(nexpr->discr){
  case ELINA_LINEXPR_SPARSE:
    elina_linexpr0_permute_dimensions_with(nexpr,perm);
    break;
  case ELINA_LINEXPR_DENSE:
    {
      size_t i;
      for (i=0; i<nexpr->size; i++){
	elina_coeff_set(&nexpr->p.coeff[perm->dim[i]],&expr->p.coeff[i]);
      }
    }
    break;
  default:
    abort();
  }
  return nexpr;
}
void
elina_linexpr0_permute_dimensions_with(elina_linexpr0_t* expr,
				    elina_dimperm_t* perm)
{
  if (expr==NULL) return;
  switch(expr->discr){
  case ELINA_LINEXPR_SPARSE:
    {
      size_t i;
      for (i=0; i<expr->size; i++){
	elina_dim_t dim = expr->p.linterm[i].dim;
	if (dim==ELINA_DIM_MAX) continue;
	expr->p.linterm[i].dim = perm->dim[dim];
      }
      qsort(expr->p.linterm,
	    expr->size,
	    sizeof(elina_linterm_t),
	    &elina_linterm_cmp);
    }
    break;
  case ELINA_LINEXPR_DENSE:
    {
      elina_linexpr0_t* nexpr = elina_linexpr0_permute_dimensions(expr,perm);
      elina_linexpr0_clear(expr);
      *expr = *nexpr;
      free(nexpr);
    }
    break;
  default:
    abort();
  }
}

/* ====================================================================== */
/* VI. Hashing, comparison */
/* ====================================================================== */

long elina_linexpr0_hash(elina_linexpr0_t* expr)
{
  if (expr->size==0){
    return elina_coeff_hash(&expr->cst);
  }
  else {
    elina_coeff_t* pcoeff;
    size_t i,dec;
    long res,res1;
    res = expr->size << 8;
    dec = 0;
    for (i=0; i<expr->size; i += (expr->size+7)/8){
      pcoeff = elina_linexpr0_coeffref(expr,i);
      res1 = (pcoeff==NULL) ? 0 : elina_coeff_hash(pcoeff);
      res += res1<<dec;
      dec++;
    }
    return res;
  }
}
bool elina_linexpr0_equal(elina_linexpr0_t* expr1,
		       elina_linexpr0_t* expr2)
{
  bool res;
  size_t i1,i2;
  elina_dim_t dim1,dim2;
  elina_coeff_t* coeff1;
  elina_coeff_t* coeff2;

  elina_coeff_reduce(&expr1->cst);
  elina_coeff_reduce(&expr2->cst);
  res = elina_coeff_equal(&expr1->cst,&expr2->cst);
  i1 = i2 = 0;
  while (res && (i1<expr1->size || i2<expr2->size)){
    if (i1<expr1->size){
      if (expr1->discr==ELINA_LINEXPR_DENSE){
	dim1 = i1; coeff1 = &expr1->p.coeff[i1];
      } else {
	dim1 = expr1->p.linterm[i1].dim; coeff1 = &expr1->p.linterm[i1].coeff;
      }
      elina_coeff_reduce(coeff1);
    }
    else {
      dim1 = ELINA_DIM_MAX;
      coeff1 = NULL;
    }
    if (i2<expr2->size){
      if (expr2->discr==ELINA_LINEXPR_DENSE){
	dim2 = i2; coeff2 = &expr2->p.coeff[i2];
      } else {
	dim2 = expr2->p.linterm[i2].dim; coeff2 = &expr2->p.linterm[i2].coeff;
      }
      elina_coeff_reduce(coeff2);
    }
    else {
      dim2 = ELINA_DIM_MAX;
      coeff2 = NULL;
    }
    if (dim1==dim2){
      i1++; i2++;
      res = elina_coeff_equal(coeff1,coeff2);
    }
    else if (dim1<dim2){
      i1++;
      res = elina_coeff_zero(coeff1);
    }
    else { /* if (dim2<dim1) */
      i2++;
      res = elina_coeff_zero(coeff2);
    }
  }
  return res;
}

int elina_linexpr0_compare(elina_linexpr0_t* expr1,
			elina_linexpr0_t* expr2)
{
  bool res;
  size_t i1,i2;
  elina_dim_t dim1,dim2;
  elina_coeff_t* coeff1;
  elina_coeff_t* coeff2;
  elina_coeff_t* coeffzero;

  coeffzero = elina_coeff_alloc(ELINA_COEFF_SCALAR);
  elina_coeff_set_scalar_double(coeffzero,0.0);
  res = 0;
  i1 = i2 = 0;
  while (res==0 && (i1<expr1->size || i2<expr2->size)){
    if (i1<expr1->size){
      if (expr1->discr==ELINA_LINEXPR_DENSE){
	dim1 = i1; coeff1 = &expr1->p.coeff[i1];
      } else {
	dim1 = expr1->p.linterm[i1].dim; coeff1 = &expr1->p.linterm[i1].coeff;
      }
      elina_coeff_reduce(coeff1);
    }
    else {
      dim1 = ELINA_DIM_MAX;
      coeff1 = NULL;
    }
    if (i2<expr2->size){
      if (expr2->discr==ELINA_LINEXPR_DENSE){
	dim2 = i2; coeff2 = &expr2->p.coeff[i2];
      } else {
	dim2 = expr2->p.linterm[i2].dim; coeff2 = &expr2->p.linterm[i2].coeff;
      }
      elina_coeff_reduce(coeff2);
    }
    else {
      dim2 = ELINA_DIM_MAX;
      coeff2 = NULL;
    }
    if (dim1==dim2){
      i1++; i2++;
      res = elina_coeff_cmp(coeff1,coeff2);
    }
    else if (dim1<dim2){
      i1++;
      res = elina_coeff_cmp(coeff1,coeffzero);
    }
    else { /* if (dim2<dim1) */
      i2++;
      res = elina_coeff_cmp(coeffzero,coeff2);
    }
  }
  if (res==0){
    elina_coeff_reduce(&expr1->cst);
    elina_coeff_reduce(&expr2->cst);
    res = elina_coeff_cmp(&expr1->cst,&expr2->cst);
  }
  elina_coeff_free(coeffzero);
  return res;
}

/* ====================================================================== */
/* VII. Array of expressions */
/* ====================================================================== */
elina_linexpr0_t ** elina_linexpr0_array_alloc(size_t size){
	elina_linexpr0_t ** res = (elina_linexpr0_t **)malloc(size*sizeof(elina_linexpr0_t*));
	return res;
}


void elina_linexpr0_array_free(elina_linexpr0_t** texpr, size_t size)
{
  size_t i;
  for (i=0;i<size;i++){
    if (texpr[i]){
      elina_linexpr0_free(texpr[i]);
      texpr[i] = NULL;
    }
  }
  free(texpr);
}

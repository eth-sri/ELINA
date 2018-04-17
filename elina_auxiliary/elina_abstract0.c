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
/* elina_abstract0.c: generic interface */
/* ************************************************************************* */


#include "elina_abstract0.h"

/* ********************************************************************** */
/* 0. Utility and checking functions */
/* ********************************************************************** */

/* ====================================================================== */
/* 0.1 Checking typing w.r.t. manager */
/* ====================================================================== */

/*
  These functions return true if everything is OK, otherwise they raise an
  exception in the manager and return false.
*/

/* One abstract value */

void elina_abstract0_checkman1_raise(elina_funid_t funid, elina_manager_t* man, elina_abstract0_t* a)
{
  char str[160];
  snprintf(str,159,"\
The abstract value of type %s is not of the type %s expected by the manager\
",
	   a->man->library,man->library);
  elina_manager_raise_exception(man,
			     ELINA_EXC_INVALID_ARGUMENT,
			     funid,
			     str);
}

/* Two abstract values */
bool elina_abstract0_checkman2(elina_funid_t funid,
			    elina_manager_t* man, elina_abstract0_t* a1, elina_abstract0_t* a2)
{
  bool res;
  char str[160];
  res = true;
  if (man->library != a1->man->library){
    snprintf(str,159,"\
The first abstract value of type %s is not of the type %s expected by the manager\
",
	     a1->man->library,man->library);
    res = false;
  }
  else if (man->library != a2->man->library){
    snprintf(str,159,"\
The second abstract value of type %s is not of the type %s expected by the manager\
",
	     a2->man->library,man->library);
    res = false;
  }
  if (!res){
    elina_manager_raise_exception(man,
			       ELINA_EXC_INVALID_ARGUMENT,
			       funid,
			       str);
  }
  return res;
}
/* Array of abstract values */
bool elina_abstract0_checkman_array(elina_funid_t funid,
				 elina_manager_t* man, elina_abstract0_t** tab, size_t size)
{
  size_t i;
  for (i=0;i<size;i++){
    if (man->library != tab[i]->man->library){
      char str[160];
      snprintf(str,159,"\
The %luth abstract value of the array is of type %s and not of the type %s expected by the manager\
",
	       (unsigned long)i,tab[i]->man->library,man->library);
      elina_manager_raise_exception(man,
				 ELINA_EXC_INVALID_ARGUMENT,
				 funid,
				 str);
      return false;
    }
  }
  return true;
}

/* ====================================================================== */
/* 0.2 Checking compatibility of arguments: abstract values */
/* ====================================================================== */

/* Check that the 2 abstract values have the same dimensionality */
bool elina_abstract0_check_abstract2(elina_funid_t funid, elina_manager_t* man,
				  elina_abstract0_t* a1, elina_abstract0_t* a2)
{
  elina_dimension_t dim1 = _elina_abstract0_dimension(a1);
  elina_dimension_t dim2 = _elina_abstract0_dimension(a2);
  if ( (dim1.intdim != dim2.intdim) || (dim1.realdim != dim2.realdim) ){
    char str[160];

    snprintf(str,159,"\
incompatible dimensions for the two arguments:\n\
first abstract0:  (%3lu,%3lu)\n\
second abstract0: (%3lu,%3lu)",
	     (unsigned long)dim1.intdim,
	     (unsigned long)dim1.realdim,
	     (unsigned long)dim2.intdim,
	     (unsigned long)dim2.realdim);
    elina_manager_raise_exception(man,
			       ELINA_EXC_INVALID_ARGUMENT,
			       funid,str);
    return false;
  } else {
    return true;
  }
}

/* Check that the array of abstract values have the same dimensionality.*/
bool elina_abstract0_check_abstract_array(elina_funid_t funid, elina_manager_t* man,
				       elina_abstract0_t** tab, size_t size)
{
  size_t i=0;
  elina_dimension_t dim;
  elina_dimension_t dim0;
  bool res;

  res = true;

  if (size==0){
    elina_manager_raise_exception(man,
			       ELINA_EXC_INVALID_ARGUMENT,
			       funid,"array of abstract values of size 0");
    res = false;
  }
  else {
    dim0 = _elina_abstract0_dimension(tab[0]);
    for (i=1; i<size; i++){
      dim = _elina_abstract0_dimension(tab[i]);
      if ( (dim.intdim != dim0.intdim) || (dim.realdim != dim0.realdim) ){
	res = false;
	break;
      }
    }
  }
  if (!res){
    char str[160];
    if (size==0){
      snprintf(str,159,"array of size 0");
    }
    else {
      snprintf(str,159,"\
incompatible dimensions for the array of polyhedra:\n\
abstract0 0: (%3lu,%3lu)\n\
abstract0 %lu: (%3lu,%3lu)\
",
	       (unsigned long)dim0.intdim,(unsigned long)dim0.realdim,
	       (unsigned long)i,
	       (unsigned long)dim.intdim,(unsigned long)dim.realdim
	       );
    }
    elina_manager_raise_exception(man,
			       ELINA_EXC_INVALID_ARGUMENT,
			       funid,str);
  }
  return res;
}

/* ====================================================================== */
/* 0.3 Checking compatibility of arguments: dimensions */
/* ====================================================================== */

/* Check that the dimension makes sense in the given dimensionality */
void elina_abstract0_check_dim_raise(elina_funid_t funid, elina_manager_t* man,
				  elina_dimension_t dimension, elina_dim_t dim,
				  const char* prefix)
{
  char str[160];

  snprintf(str,159,"\
%s:\n\
abstract0:  (%3lu,%3lu)\n\
dimension:  %3lu\n",
	   prefix,
	   (unsigned long)dimension.intdim,
	   (unsigned long)dimension.realdim,
	   (unsigned long)dim);
  elina_manager_raise_exception(man,
			     ELINA_EXC_INVALID_ARGUMENT,
			     funid,str);
}

/* Check that the array of dimensions make sense in the given dimensionality */
bool elina_abstract0_check_dim_array(elina_funid_t funid, elina_manager_t* man,
				  elina_dimension_t dimension, elina_dim_t* tdim, size_t size)
{
  size_t i;
  for (i=0;i<size;i++){
    elina_dim_t dim = tdim[i];
    if (dim>=dimension.intdim+dimension.realdim){
      char str[80];
      sprintf(str,"incompatible %luth dimension in the array for the abstract value",(unsigned long)i);
      elina_abstract0_check_dim_raise(funid,man,dimension,dim,str);
      return false;
    }
  }
  return true;
}

/* ====================================================================== */
/* 0.4 Checking compatibility of arguments: expressions */
/* ====================================================================== */

void elina_abstract0_check_expr_raise(elina_funid_t funid, elina_manager_t* man,
				   elina_dimension_t dimension,
				   elina_dim_t dim,
				   char* prefix)
{
  char str[160];
  snprintf(str,159,"\
%s:\n\
abstract0: (%3lu,%3lu)\n\
dimension: %3lu\
",
	   prefix,
	   (unsigned long)dimension.intdim,(unsigned long)dimension.realdim,
	   (unsigned long)dim);
  elina_manager_raise_exception(man,
			     ELINA_EXC_INVALID_ARGUMENT,
			     funid,str);
}

/* Check that the linear expression makes sense in the given dimensionality */
elina_dim_t elina_abstract0_check_linexpr_check(elina_dimension_t dimension,
					  elina_linexpr0_t* expr)
{
  int i;
  size_t nbdims;
  elina_dim_t dim;

  nbdims = dimension.intdim+dimension.realdim;
  dim = 0;
  switch (expr->discr){
  case ELINA_LINEXPR_DENSE:
    if (expr->size > nbdims){
      dim = expr->size-1;
    }
    else {
      dim = ELINA_DIM_MAX;
    }
    break;
  case ELINA_LINEXPR_SPARSE:
    /* Assumed to be sorted (this is not checked */
    for (i=expr->size-1; i>=0; i--){
      dim = expr->p.linterm[i].dim;
      if (dim!=ELINA_DIM_MAX){
	if (dim < nbdims)
	  dim = ELINA_DIM_MAX;
	;
	break;
      }
    }
    if (!dim) dim = ELINA_DIM_MAX;
    break;
  default:
    abort();
  }
  return dim;
}

bool elina_abstract0_check_linexpr(elina_funid_t funid, elina_manager_t* man,
				elina_dimension_t dimension,
				elina_linexpr0_t* expr)
{
  elina_dim_t dim = elina_abstract0_check_linexpr_check(dimension,expr);
  if (dim!=ELINA_DIM_MAX){
    elina_abstract0_check_expr_raise(funid,man,dimension,dim,
				  "incompatible dimension in the linear expression for the abstract value");
    return false;
  } else {
    return true;
  }
}

/* Check that the tree expression makes sense in the given dimensionality */
elina_dim_t elina_abstract0_check_texpr_check(elina_dimension_t dimension,
					elina_texpr0_t* expr)
{
  elina_dim_t dim;

  dim = elina_texpr0_max_dim(expr);
  if (dim <= dimension.intdim+dimension.realdim)
    return ELINA_DIM_MAX;
  else
    return dim;
}

bool elina_abstract0_check_texpr(elina_funid_t funid, elina_manager_t* man,
			      elina_dimension_t dimension,
			      elina_texpr0_t* expr)
{
  elina_dim_t dim = elina_abstract0_check_texpr_check(dimension,expr);
  if (dim!=ELINA_DIM_MAX){
    elina_abstract0_check_expr_raise(funid,man,dimension,dim,
				  "incompatible dimension in the tree expression for the abstract value");
    return false;
  } else {
    return true;
  }
}

/* ====================================================================== */
/* 0.5 Checking compatibility of arguments: array of expressions/constraints/generators */
/* ====================================================================== */

/* Check that array of linear expressions makes sense in the given dimensionality */
bool elina_abstract0_check_linexpr_array(elina_funid_t funid, elina_manager_t* man,
				      elina_dimension_t dimension,
				      elina_linexpr0_t** texpr, size_t size)
{
  size_t i;

  for (i=0;i<size; i++){
    if (texpr[i]==NULL){
      char str[80];
      sprintf(str,"null pointer in the %luth expression of the array",(unsigned long)i);
      elina_manager_raise_exception(man,
				 ELINA_EXC_INVALID_ARGUMENT,
				 funid,str);
      return false;
    }
    elina_dim_t dim = elina_abstract0_check_linexpr_check(dimension,texpr[i]);
    if (dim!=ELINA_DIM_MAX){
      char str[80];
      sprintf(str,"incompatible dimension in the %luth expression of the array",(unsigned long)i);
      elina_abstract0_check_expr_raise(funid,man,dimension,dim,str);
      return false;
    }
  }
  return true;
}
/* Check that array of linear constraint makes sense in the given dimensionality */
bool elina_abstract0_check_lincons_array(elina_funid_t funid, elina_manager_t* man,
				      elina_dimension_t dimension,
				      elina_lincons0_array_t* array)
{
  size_t i;

  for (i=0;i<array->size; i++){
    if (array->p[i].linexpr0==NULL){
      char str[80];
      sprintf(str,"null pointer in the %luth constraint of the array",(unsigned long)i);
      elina_manager_raise_exception(man,
				 ELINA_EXC_INVALID_ARGUMENT,
				 funid,str);
      return false;
    }
    elina_dim_t dim = elina_abstract0_check_linexpr_check(dimension,array->p[i].linexpr0);
    if (dim!=ELINA_DIM_MAX){
      char str[80];
      sprintf(str,"incompatible dimension in the %luth constraint of the array",(unsigned long)i);
      elina_abstract0_check_expr_raise(funid,man,dimension,dim,str);
      return false;
    }
  }
  return true;
}


/* Check that array of tree expressions makes sense in the given dimensionality */
bool elina_abstract0_check_texpr_array(elina_funid_t funid, elina_manager_t* man,
				    elina_dimension_t dimension,
				    elina_texpr0_t** texpr, size_t size)
{
  size_t i;

  for (i=0;i<size; i++){
    if (texpr[i]==NULL){
      char str[80];
      sprintf(str,"null pointer in the %luth expression of the array",(unsigned long)i);
      elina_manager_raise_exception(man,
				 ELINA_EXC_INVALID_ARGUMENT,
				 funid,str);
      return false;
    }
    elina_dim_t dim = elina_abstract0_check_texpr_check(dimension,texpr[i]);
    if (dim!=ELINA_DIM_MAX){
      char str[80];
      sprintf(str,"incompatible dimension in the %luth expression of the array",(unsigned long)i);
      elina_abstract0_check_expr_raise(funid,man,dimension,dim,str);
      return false;
    }
  }
  return true;
}
/* Check that array of tree constraint makes sense in the given dimensionality */
bool elina_abstract0_check_tcons_array(elina_funid_t funid, elina_manager_t* man,
				      elina_dimension_t dimension,
				      elina_tcons0_array_t* array)
{
  size_t i;

  for (i=0;i<array->size; i++){
    if (array->p[i].texpr0==NULL){
      char str[80];
      sprintf(str,"null pointer in the %luth constraint of the array",(unsigned long)i);
      elina_manager_raise_exception(man,
				 ELINA_EXC_INVALID_ARGUMENT,
				 funid,str);
      return false;
    }
    elina_dim_t dim = elina_abstract0_check_texpr_check(dimension,array->p[i].texpr0);
    if (dim!=ELINA_DIM_MAX){
      char str[80];
      sprintf(str,"incompatible dimension in the %luth constraint of the array",(unsigned long)i);
      elina_abstract0_check_expr_raise(funid,man,dimension,dim,str);
      return false;
    }
  }
  return true;
}

/* ====================================================================== */
/* 0.6 Checking compatibility of arguments: dimchange */
/* ====================================================================== */

bool elina_abstract0_check_elina_dimchange_add(elina_funid_t funid, elina_manager_t* man,
					 elina_dimension_t dimension, elina_dimchange_t* dimchange)
{
  size_t i,size;
  bool res = true;
  /* Check consistency between intdim and the dimensions in the array */
  for (i=0; i<dimchange->intdim; i++){
    if (dimchange->dim[i]>dimension.intdim){
      
      res = false;
      break;
    }
  }
  size = dimchange->intdim+dimchange->realdim;
  if (res && size>0){
    /* Check sortedness */
    elina_dim_t dim = 0;
    for (i=1; i<size; i++){
      if (dim>dimchange->dim[i]){
	res = false;
	
	break;
      } else {
	dim = dimchange->dim[i];
      }
    }
    res = res && dimchange->dim[size-1]<=dimension.intdim+dimension.realdim;
  }
  if (!res){
    elina_manager_raise_exception(man,
			       ELINA_EXC_INVALID_ARGUMENT,
			       funid,
			       "dimchange->dim is not sorted or is inconsistent wrt dimchange->intdim or abstract0");
  }
  return res;
}

bool elina_abstract0_check_elina_dimchange_remove(elina_funid_t funid, elina_manager_t* man,
					    elina_dimension_t dimension, elina_dimchange_t* dimchange)
{
  size_t i,size;
  bool res = true;
  /* Check consistency between intdim and the dimensions in the array */
  for (i=0; i<dimchange->intdim; i++){
    if (dimchange->dim[i]>=dimension.intdim){
      res = false;
      break;
    }
  }
  size = dimchange->intdim+dimchange->realdim;
  if (res && size>0){
    /* Check sortedness */
    elina_dim_t dim = 0;
    for (i=1; i<dimchange->intdim+dimchange->realdim; i++){
      if (dim>=dimchange->dim[i]){
	res = false;
	break;
      } else {
	dim = dimchange->dim[i];
      }
    }
    res = res && dimchange->dim[size-1]<dimension.intdim+dimension.realdim;
  }
  if (!res){
    elina_manager_raise_exception(man,
			       ELINA_EXC_INVALID_ARGUMENT,
			       funid,
			       "dimchange->dim is not sorted, contains duplicatas, or is inconsistent wrt dimchange->intdim or abstract0");
  }
  return res;
}
bool elina_abstract0_check_dimperm(elina_funid_t funid, elina_manager_t* man,
				elina_dimension_t dimension, elina_dimperm_t* perm)
{
  size_t i;
  size_t size = dimension.intdim+dimension.realdim;
  bool res = (dimension.intdim+dimension.realdim==perm->size);
  if (res){
    for (i=0; i<perm->size; i++){
      if (perm->dim[i]>=size){
	res = false;
	break;
      }
    }
  }
  if (!res){
    elina_manager_raise_exception(man,
			       ELINA_EXC_INVALID_ARGUMENT,
			       funid,
			       "permutation is not valid");
  }
  return res;
}

/* ********************************************************************** */
/* I. General management */
/* ********************************************************************** */

/* ============================================================ */
/* I.1 Memory */
/* ============================================================ */

elina_abstract0_t* elina_abstract0_copy(elina_manager_t* man, elina_abstract0_t* a)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_COPY,man,a)){
    void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_COPY];
    return elina_abstract0_cons(man,ptr(man,a->value));
  }
  else {
    elina_dimension_t dimension = _elina_abstract0_dimension(a);
    return elina_abstract0_top(man,
			    dimension.intdim,
			    dimension.realdim);
  }
}

void elina_abstract0_free(elina_manager_t* man, elina_abstract0_t* a)
{
  if (a->man==NULL){
    fprintf(stderr,"elina_abstract0_c: elina_abstract0_free: the abstract value has probably already been deallocated !\n");
    abort();
  }
  if (elina_abstract0_checkman1(ELINA_FUNID_FREE,man,a)){
    void (*ptr)(elina_manager_t*,elina_abstract0_t*) = man->funptr[ELINA_FUNID_FREE];
    ptr(man,a->value);
  }
  else {
    void (*ptr)(elina_manager_t*,elina_abstract0_t*) = a->man->funptr[ELINA_FUNID_FREE];
    ptr(a->man,a->value);
  }
  elina_manager_free(a->man);
  a->man = NULL;
  a->value = NULL;
  free(a);
}
size_t elina_abstract0_size(elina_manager_t* man, elina_abstract0_t* a)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_ASIZE,man,a)){
    size_t (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_ASIZE];
    return ptr(man,a->value);
  }
  else {
    return 0;
  }
}

/* ============================================================ */
/* I.2 Control of internal representation */
/* ============================================================ */
void elina_abstract0_minimize(elina_manager_t* man, elina_abstract0_t* a)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_MINIMIZE,man,a)){
    void (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_MINIMIZE];
    ptr(man,a->value);
  }
}
void elina_abstract0_canonicalize(elina_manager_t* man, elina_abstract0_t* a)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_CANONICALIZE,man,a)){
    void (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_CANONICALIZE];
    ptr(man,a->value);
  }
}
int elina_abstract0_hash(elina_manager_t* man, elina_abstract0_t* a)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_CANONICALIZE,man,a)){
    int (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_HASH];
    return ptr(man,a->value);
  }
  else
    return 0;
}
void elina_abstract0_approximate(elina_manager_t* man, elina_abstract0_t* a, int n)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_APPROXIMATE,man,a)){
    void (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_APPROXIMATE];
    ptr(man,a->value,n);
  }
}

/* ============================================================ */
/* I.3 Printing */
/* ============================================================ */
void elina_abstract0_fprint(FILE* stream, elina_manager_t* man,
			 elina_abstract0_t* a,
			 char** name_of_dim)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_FPRINT,man,a)){
    void (*ptr)(FILE*,elina_manager_t*,...) = man->funptr[ELINA_FUNID_FPRINT];
    ptr(stream,man,a->value,name_of_dim);
  }
  else {
    fprintf(stream,"elina_abstract0_c: elina_abstract0_fprint: INVALID_ARGUMENT\n");
  }
}

void elina_abstract0_fprintdiff(FILE* stream, elina_manager_t* man,
			     elina_abstract0_t* a, elina_abstract0_t* b,
			     char** name_of_dim)
{
  if (elina_abstract0_checkman2(ELINA_FUNID_FPRINTDIFF,man,a,b) &&
      elina_abstract0_check_abstract2(ELINA_FUNID_FPRINTDIFF,man,a,b) ){
    void (*ptr)(FILE*,elina_manager_t*,...) = man->funptr[ELINA_FUNID_FPRINTDIFF];
    ptr(stream,man,a->value,b->value,name_of_dim);
  }
  else {
    fprintf(stream,"elina_abstract0_c: elina_abstract0_fprintdiff: INVALID ARGUMENT\n");
  }
}
void elina_abstract0_fdump(FILE* stream, elina_manager_t* man, elina_abstract0_t* a)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_FDUMP,man,a)){
    void (*ptr)(FILE*,elina_manager_t*,...) = man->funptr[ELINA_FUNID_FDUMP];
    ptr(stream,man,a->value);
  }
  else {
    fprintf(stream,"elina_abstract0_c: elina_abstract0_fdump: INVALID_ARGUMENT\n");
  }
}

/* ============================================================ */
/* I.4 Serialization */
/* ============================================================ */
elina_membuf_t elina_abstract0_serialize_raw(elina_manager_t* man, elina_abstract0_t* a)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_SERIALIZE_RAW,man,a)){
    elina_membuf_t (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_SERIALIZE_RAW];
    return ptr(man,a->value);
  }
  else {
    elina_membuf_t res = { NULL, 0 };
    return res;
  }
}
elina_abstract0_t* elina_abstract0_deserialize_raw(elina_manager_t* man, void* p, size_t* size)
{
  void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_DESERIALIZE_RAW];
  return elina_abstract0_cons(man,ptr(man,p,size));
}

/* ********************************************************************** */
/* II. Constructor, accessors, tests and property extraction */
/* ********************************************************************** */
/* ============================================================ */
/* II.1 Basic constructors */
/* ============================================================ */
elina_abstract0_t* elina_abstract0_bottom(elina_manager_t* man, size_t intdim, size_t realdim)
{
  void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_BOTTOM];
  return elina_abstract0_cons(man,ptr(man,intdim,realdim));
}
elina_abstract0_t* elina_abstract0_top(elina_manager_t* man, size_t intdim, size_t realdim){
  void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_TOP];
  return elina_abstract0_cons(man,ptr(man,intdim,realdim));
}
elina_abstract0_t* elina_abstract0_of_box(elina_manager_t* man,
				    size_t intdim, size_t realdim,
				    elina_interval_t** tinterval){
  void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_OF_BOX];
  return elina_abstract0_cons(man,ptr(man,intdim,realdim,tinterval));
}

/* ============================================================ */
/* II.2 Accessors */
/* ============================================================ */
elina_dimension_t elina_abstract0_dimension(elina_manager_t* man, elina_abstract0_t* a)
{
  elina_abstract0_checkman1(ELINA_FUNID_DIMENSION,man,a);
  elina_dimension_t (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_DIMENSION];
  return ptr(a->man,a->value);
}

/* ============================================================ */
/* II.3 Tests */
/* ============================================================ */
bool elina_abstract0_is_bottom(elina_manager_t* man, elina_abstract0_t* a)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_IS_BOTTOM,man,a)){
    bool (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_IS_BOTTOM];
    return ptr(man,a->value);
  }
  else {
    man->result.flag_exact = false;
    return false;
  }
}
bool elina_abstract0_is_top(elina_manager_t* man, elina_abstract0_t* a)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_IS_TOP,man,a)){
    bool (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_IS_TOP];
    return ptr(man,a->value);
  }
  else {
    man->result.flag_exact = false;
    return false;
  }
}
bool elina_abstract0_is_leq(elina_manager_t* man, elina_abstract0_t* a1, elina_abstract0_t* a2)
{
  if (a1==a2){
    return true;
  }
  else if (elina_abstract0_checkman2(ELINA_FUNID_IS_LEQ,man,a1,a2) &&
	   elina_abstract0_check_abstract2(ELINA_FUNID_IS_EQ,man,a1,a2)){
    if (a1->value==a2->value) return true;
    bool (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_IS_LEQ];
    return ptr(man,a1->value,a2->value);
  }
  else {
    man->result.flag_exact = false;
    return false;
  }
}
bool elina_abstract0_is_eq(elina_manager_t* man, elina_abstract0_t* a1, elina_abstract0_t* a2)
{
  if (a1==a2){
    return true;
  }
  if (elina_abstract0_checkman2(ELINA_FUNID_IS_EQ,man,a1,a2) &&
      elina_abstract0_check_abstract2(ELINA_FUNID_IS_EQ,man,a1,a2)){
    if (a1->value==a2->value) return true;
    bool (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_IS_EQ];
    return ptr(man,a1->value,a2->value);
  }
  else {
    man->result.flag_exact = false;
    return false;
  }
}
bool elina_abstract0_sat_lincons(elina_manager_t* man, elina_abstract0_t* a, elina_lincons0_t* lincons)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_SAT_LINCONS,man,a) &&
      elina_abstract0_check_linexpr(ELINA_FUNID_SAT_LINCONS,man,_elina_abstract0_dimension(a),lincons->linexpr0) ){
    bool (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_SAT_LINCONS];
    return ptr(man,a->value,lincons);
  }
  else {
    man->result.flag_exact = false;
    return false;
  }
}
bool elina_abstract0_sat_tcons(elina_manager_t* man, elina_abstract0_t* a, elina_tcons0_t* tcons)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_SAT_TCONS,man,a) &&
      elina_abstract0_check_texpr(ELINA_FUNID_SAT_TCONS,man,_elina_abstract0_dimension(a),tcons->texpr0) ){
    bool (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_SAT_TCONS];
    return ptr(man,a->value,tcons);
  }
  else {
    man->result.flag_exact = false;
    return false;
  }
}
bool elina_abstract0_sat_interval(elina_manager_t* man, elina_abstract0_t* a,
			       elina_dim_t dim, elina_interval_t* interval)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_SAT_INTERVAL,man,a) &&
      elina_abstract0_check_dim(ELINA_FUNID_SAT_INTERVAL,man,_elina_abstract0_dimension(a),dim)){
    bool (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_SAT_INTERVAL];
    return ptr(man,a->value,dim,interval);
  }
  else {
    man->result.flag_exact = false;
    return false;
  }
}
bool elina_abstract0_is_dimension_unconstrained(elina_manager_t* man, elina_abstract0_t* a,
					     elina_dim_t dim)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED,man,a) &&
      elina_abstract0_check_dim(ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED,man,_elina_abstract0_dimension(a),dim)){
    bool (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED];
    return ptr(man,a->value,dim);
  }
  else {
    man->result.flag_exact = false;
    return false;
  }
}

/* ============================================================ */
/* II.4 Extraction of properties */
/* ============================================================ */
elina_interval_t* elina_abstract0_bound_linexpr(elina_manager_t* man,
					  elina_abstract0_t* a, elina_linexpr0_t* expr)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_BOUND_LINEXPR,man,a) &&
      elina_abstract0_check_linexpr(ELINA_FUNID_BOUND_LINEXPR,man,_elina_abstract0_dimension(a),expr)){
    elina_interval_t* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_BOUND_LINEXPR];
    return ptr(man,a->value,expr);
  }
  else {
    elina_interval_t* itv = elina_interval_alloc();
    elina_interval_reinit(itv,man->option.scalar_discr);
    elina_interval_set_top(itv);
    return itv;
  }
}
elina_interval_t* elina_abstract0_bound_texpr(elina_manager_t* man,
					elina_abstract0_t* a, elina_texpr0_t* expr)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_BOUND_TEXPR,man,a) &&
      elina_abstract0_check_texpr(ELINA_FUNID_BOUND_TEXPR,man,_elina_abstract0_dimension(a),expr)){
    elina_interval_t* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_BOUND_TEXPR];
    return ptr(man,a->value,expr);
  }
  else {
    elina_interval_t* itv = elina_interval_alloc();
    elina_interval_reinit(itv,man->option.scalar_discr);
    elina_interval_set_top(itv);
    return itv;
  }
}
elina_interval_t* elina_abstract0_bound_dimension(elina_manager_t* man,
					    elina_abstract0_t* a, elina_dim_t dim)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_BOUND_DIMENSION,man,a) &&
      elina_abstract0_check_dim(ELINA_FUNID_BOUND_DIMENSION,man,_elina_abstract0_dimension(a),dim)){
    elina_interval_t* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_BOUND_DIMENSION];
    return ptr(man,a->value,dim);
  }
  else {
    elina_interval_t* itv = elina_interval_alloc();
    elina_interval_reinit(itv,man->option.scalar_discr);
    elina_interval_set_top(itv);
    return itv;
  }
}
elina_lincons0_array_t elina_abstract0_to_lincons_array(elina_manager_t* man, elina_abstract0_t* a)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_TO_LINCONS_ARRAY,man,a)){
    elina_lincons0_array_t (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_TO_LINCONS_ARRAY];
    return ptr(man,a->value);
  }
  else {
    elina_lincons0_array_t res = { NULL, 0 };
    return res;
  }
}
elina_tcons0_array_t elina_abstract0_to_tcons_array(elina_manager_t* man, elina_abstract0_t* a)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_TO_TCONS_ARRAY,man,a)){
    elina_tcons0_array_t (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_TO_TCONS_ARRAY];
    return ptr(man,a->value);
  }
  else {
    elina_tcons0_array_t res = { NULL, 0 };
    return res;
  }
}
elina_interval_t** elina_abstract0_to_box(elina_manager_t* man, elina_abstract0_t* a)
{
  if (elina_abstract0_checkman1(ELINA_FUNID_TO_BOX,man,a)){
    elina_interval_t** (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_TO_BOX];
    return ptr(man,a->value);
  }
  else {
    size_t i;
    elina_dimension_t d = _elina_abstract0_dimension(a);
    elina_interval_t** titv = elina_interval_array_alloc(d.intdim+d.realdim);
    for (i=0; i<d.intdim+d.realdim; i++){
      elina_interval_reinit(titv[i],man->option.scalar_discr);
      elina_interval_set_top(titv[i]);
    }
    return titv;
  }
}


/* ********************************************************************** */
/* III. Operations */
/* ********************************************************************** */
/* ============================================================ */
/* III.1 Meet and Join */
/* ============================================================ */

elina_abstract0_t* elina_abstract0_meetjoin(elina_funid_t funid,
				      elina_manager_t* man, bool destructive, elina_abstract0_t* a1, elina_abstract0_t* a2)
{
  if (elina_abstract0_checkman2(funid,man,a1,a2) &&
      elina_abstract0_check_abstract2(funid,man,a1,a2)){
    void* (*ptr)(elina_manager_t*,...) = man->funptr[funid];
    void* value = ptr(man,destructive,a1->value,a2->value);
    return elina_abstract0_cons2(man,destructive,a1,value);
  }
  else {
    elina_dimension_t dimension = _elina_abstract0_dimension(a1);
    if (destructive) _elina_abstract0_free(a1);
    return elina_abstract0_top(man,
			    dimension.intdim,
			    dimension.realdim);
  }
}
elina_abstract0_t* elina_abstract0_meet(elina_manager_t* man, bool destructive, elina_abstract0_t* a1, elina_abstract0_t* a2){
  return  elina_abstract0_meetjoin(ELINA_FUNID_MEET,man,destructive,a1,a2);
}
elina_abstract0_t* elina_abstract0_join(elina_manager_t* man, bool destructive, elina_abstract0_t* a1, elina_abstract0_t* a2){
  return  elina_abstract0_meetjoin(ELINA_FUNID_JOIN,man,destructive,a1,a2);
}

elina_abstract0_t* elina_abstract0_meetjoin_array(elina_funid_t funid, elina_manager_t* man, elina_abstract0_t** tab, size_t size)
{
  if (elina_abstract0_checkman_array(funid,man,tab,size) &&
      elina_abstract0_check_abstract_array(funid,man,tab,size)){
    size_t i;
    elina_abstract0_t* res;
    void* (*ptr)(elina_manager_t*,...) = man->funptr[funid];
    void** ntab = malloc(size*sizeof(void*));
    for (i=0;i<size;i++) ntab[i] = tab[i]->value;
    res = elina_abstract0_cons(man,ptr(man,ntab,size));
    free(ntab);
    return res;
  }
  else {
    elina_dimension_t dimension = { 0, 0};
    if (size>0){
      dimension = _elina_abstract0_dimension(tab[0]);
    }
    return elina_abstract0_top(man,
			    dimension.intdim,
			    dimension.realdim);
  }
}
elina_abstract0_t* elina_abstract0_meet_array(elina_manager_t* man, elina_abstract0_t** tab, size_t size){
  return elina_abstract0_meetjoin_array(ELINA_FUNID_MEET_ARRAY,man,tab,size);
}
elina_abstract0_t* elina_abstract0_join_array(elina_manager_t* man, elina_abstract0_t** tab, size_t size){
  return elina_abstract0_meetjoin_array(ELINA_FUNID_JOIN_ARRAY,man,tab,size);
}
elina_abstract0_t* elina_abstract0_meet_lincons_array(elina_manager_t* man,
						bool destructive,
						elina_abstract0_t* a,
						elina_lincons0_array_t* array)
{
  elina_dimension_t dimension = _elina_abstract0_dimension(a);
  if (elina_abstract0_checkman1(ELINA_FUNID_MEET_LINCONS_ARRAY,man,a) &&
      elina_abstract0_check_lincons_array(ELINA_FUNID_MEET_LINCONS_ARRAY,man,dimension,array) ){
    void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_MEET_LINCONS_ARRAY];
    void* value = ptr(man,destructive,a->value,array);
    return elina_abstract0_cons2(man,destructive,a,value);
  }
  else {
    if (destructive) _elina_abstract0_free(a);
    return elina_abstract0_top(man,
			    dimension.intdim,
			    dimension.realdim);
  }
}
elina_abstract0_t* elina_abstract0_meet_tcons_array(elina_manager_t* man,
					      bool destructive,
					      elina_abstract0_t* a,
					      elina_tcons0_array_t* array)
{
  elina_dimension_t dimension = _elina_abstract0_dimension(a);
  if (elina_abstract0_checkman1(ELINA_FUNID_MEET_TCONS_ARRAY,man,a) &&
      elina_abstract0_check_tcons_array(ELINA_FUNID_MEET_TCONS_ARRAY,man,dimension,array) ){
    void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_MEET_TCONS_ARRAY];
    void* value = ptr(man,destructive,a->value,array);
    return elina_abstract0_cons2(man,destructive,a,value);
  }
  else {
    if (destructive) _elina_abstract0_free(a);
    return elina_abstract0_top(man,
			    dimension.intdim,
			    dimension.realdim);
  }
}


/* ============================================================ */
/* III.2 Assignment and Substitutions */
/* ============================================================ */

elina_abstract0_t* elina_abstract0_asssub_linexpr_array(elina_funid_t funid,
						  elina_manager_t* man,
						  bool destructive,
						  elina_abstract0_t* a,
						  elina_dim_t* tdim,
						  elina_linexpr0_t** texpr,
						  size_t size,
						  elina_abstract0_t* dest)
{
  if (size==0){
    if (dest){
      return elina_abstract0_meet(man,destructive,a,dest);
    }
    else if (destructive){
      return a;
    }
    else {
      return elina_abstract0_copy(man,a);
    }
  }
  else {
    elina_dimension_t dimension = _elina_abstract0_dimension(a);
    if (elina_abstract0_checkman1(funid,man,a) &&
	(dest!=NULL ? (elina_abstract0_checkman1(funid,man,dest) && elina_abstract0_check_abstract2(funid,man,a,dest)) : true) &&
	elina_abstract0_check_dim_array(funid,man,dimension,tdim,size) &&
	elina_abstract0_check_linexpr_array(funid,man,dimension,texpr,size) ){
      void* (*ptr)(elina_manager_t*,...) = man->funptr[funid];
      void* value = ptr(man,destructive,a->value,tdim,texpr,size,dest ? dest->value : NULL);
      return elina_abstract0_cons2(man,destructive,a,value);
    }
    else {
      if (destructive) _elina_abstract0_free(a);
      return elina_abstract0_top(man,
			      dimension.intdim,
			      dimension.realdim);
    }
  }
}
elina_abstract0_t* elina_abstract0_assign_linexpr_array(elina_manager_t* man,
						  bool destructive,
						  elina_abstract0_t* a,
						  elina_dim_t* tdim,
						  elina_linexpr0_t** texpr,
						  size_t size,
						  elina_abstract0_t* dest)
{
  return elina_abstract0_asssub_linexpr_array(ELINA_FUNID_ASSIGN_LINEXPR_ARRAY,
					   man,destructive,a,tdim,texpr,size,dest);
}
elina_abstract0_t* elina_abstract0_substitute_linexpr_array(elina_manager_t* man,
						      bool destructive,
						      elina_abstract0_t* a,
						      elina_dim_t* tdim,
						      elina_linexpr0_t** texpr,
						      size_t size,
						      elina_abstract0_t* dest)
{
  return elina_abstract0_asssub_linexpr_array(ELINA_FUNID_SUBSTITUTE_LINEXPR_ARRAY,
					   man,destructive,a,tdim,texpr,size,dest);
}
elina_abstract0_t* elina_abstract0_asssub_texpr_array(elina_funid_t funid,
						elina_manager_t* man,
						bool destructive,
						elina_abstract0_t* a,
						elina_dim_t* tdim,
						elina_texpr0_t** texpr,
						size_t size,
						elina_abstract0_t* dest)
{
  if (size==0){
    if (dest){
      return elina_abstract0_meet(man,destructive,a,dest);
    }
    else if (destructive){
      return a;
    }
    else {
      return elina_abstract0_copy(man,a);
    }
  }
  else {
    elina_dimension_t dimension = _elina_abstract0_dimension(a);
    if (elina_abstract0_checkman1(funid,man,a) &&
	(dest!=NULL ? (elina_abstract0_checkman1(funid,man,dest) && elina_abstract0_check_abstract2(funid,man,a,dest)) : true) &&
	elina_abstract0_check_dim_array(funid,man,dimension,tdim,size) &&
	elina_abstract0_check_texpr_array(funid,man,dimension,texpr,size) ){
      void* (*ptr)(elina_manager_t*,...) = man->funptr[funid];
      void* value = ptr(man,destructive,a->value,tdim,texpr,size,dest ? dest->value : NULL);
      return elina_abstract0_cons2(man,destructive,a,value);
    }
    else {
      if (destructive) _elina_abstract0_free(a);
      return elina_abstract0_top(man,
			      dimension.intdim,
			      dimension.realdim);
    }
  }
}
elina_abstract0_t* elina_abstract0_assign_texpr_array(elina_manager_t* man,
						bool destructive,
						elina_abstract0_t* a,
						elina_dim_t* tdim,
						elina_texpr0_t** texpr,
						size_t size,
						elina_abstract0_t* dest)
{
  return elina_abstract0_asssub_texpr_array(ELINA_FUNID_ASSIGN_TEXPR_ARRAY,
					 man,destructive,a,tdim,texpr,size,dest);
}
elina_abstract0_t* elina_abstract0_substitute_texpr_array(elina_manager_t* man,
						    bool destructive,
						    elina_abstract0_t* a,
						    elina_dim_t* tdim,
						    elina_texpr0_t** texpr,
						    size_t size,
						    elina_abstract0_t* dest)
{
  return elina_abstract0_asssub_texpr_array(ELINA_FUNID_SUBSTITUTE_TEXPR_ARRAY,
					 man,destructive,a,tdim,texpr,size,dest);
}

/* ============================================================ */
/* III.3 Projections */
/* ============================================================ */

elina_abstract0_t* elina_abstract0_forget_array(elina_manager_t* man,
					  bool destructive,
					  elina_abstract0_t* a,
					  elina_dim_t* tdim, size_t size,
					  bool project)
{
  if (size==0){
    if (destructive){
      return a;
    }
    else {
      return elina_abstract0_copy(man,a);
    }
  }
  else {
    elina_dimension_t dimension = _elina_abstract0_dimension(a);
    if (elina_abstract0_checkman1(ELINA_FUNID_FORGET_ARRAY,man,a) &&
	elina_abstract0_check_dim_array(ELINA_FUNID_FORGET_ARRAY,man,dimension,tdim,size)){
      void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_FORGET_ARRAY];
      void* value = ptr(man,destructive,a->value,tdim,size,project);
      return elina_abstract0_cons2(man,destructive,a,value);
    }
    else {
      if (destructive) _elina_abstract0_free(a);
      return elina_abstract0_top(man,
			      dimension.intdim,
			      dimension.realdim);
    }
  }
}

/* ============================================================ */
/* III.4 Change and permutation of dimensions */
/* ============================================================ */

elina_abstract0_t* elina_abstract0_add_dimensions(elina_manager_t* man,
					    bool destructive,
					    elina_abstract0_t* a,
					    elina_dimchange_t* dimchange,
					    bool project)
{
  if (dimchange->intdim+dimchange->realdim==0){
    if (destructive){
      return a;
    }
    else {
      return elina_abstract0_copy(man,a);
    }
  }
  else {
    elina_dimension_t dimension = _elina_abstract0_dimension(a);
    if (elina_abstract0_checkman1(ELINA_FUNID_ADD_DIMENSIONS,man,a) &&
	elina_abstract0_check_elina_dimchange_add(ELINA_FUNID_ADD_DIMENSIONS,man,dimension,dimchange)){
      void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_ADD_DIMENSIONS];
      void* value = ptr(man,destructive,a->value,dimchange,project);
      return elina_abstract0_cons2(man,destructive,a,value);
    }
    else {
      if (destructive) _elina_abstract0_free(a);
      return elina_abstract0_top(man,
			      dimension.intdim+dimchange->intdim,
			      dimension.realdim+dimchange->realdim);
    }
  }
}
elina_abstract0_t* elina_abstract0_remove_dimensions(elina_manager_t* man,
					       bool destructive,
					       elina_abstract0_t* a,
					       elina_dimchange_t* dimchange)
{
  if (dimchange->intdim+dimchange->realdim==0){
    if (destructive){
      return a;
    }
    else {
      return elina_abstract0_copy(man,a);
    }
  }
  else {
    elina_dimension_t dimension = _elina_abstract0_dimension(a);
    if (elina_abstract0_checkman1(ELINA_FUNID_REMOVE_DIMENSIONS,man,a) &&
	elina_abstract0_check_elina_dimchange_remove(ELINA_FUNID_REMOVE_DIMENSIONS,man,dimension,dimchange)){
      void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_REMOVE_DIMENSIONS];
      void* value = ptr(man,destructive,a->value,dimchange);
      return elina_abstract0_cons2(man,destructive,a,value);
    }
    else {
      if (destructive) _elina_abstract0_free(a);
      return elina_abstract0_top(man,
			      dimension.intdim>=dimchange->intdim ? dimension.intdim-dimchange->intdim : 0,
			      dimension.realdim>=dimchange->realdim ? dimension.realdim-dimchange->realdim : 0);
    }
  }
}
elina_abstract0_t* elina_abstract0_permute_dimensions(elina_manager_t* man,
						bool destructive,
						elina_abstract0_t* a,
						elina_dimperm_t* perm)
{
  elina_dimension_t dimension = _elina_abstract0_dimension(a);
  if (elina_abstract0_checkman1(ELINA_FUNID_PERMUTE_DIMENSIONS,man,a) &&
      elina_abstract0_check_dimperm(ELINA_FUNID_PERMUTE_DIMENSIONS,man,dimension,perm)){
    void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_PERMUTE_DIMENSIONS];
    void* value = ptr(man,destructive,a->value,perm);
    return elina_abstract0_cons2(man,destructive,a,value);
  }
  else {
    if (destructive) _elina_abstract0_free(a);
    return elina_abstract0_top(man,
			    dimension.intdim,
			    dimension.realdim);
  }
}

/* ============================================================ */
/* III.5 Expansion and folding of dimensions */
/* ============================================================ */
elina_abstract0_t* elina_abstract0_expand(elina_manager_t* man,
				    bool destructive,
				    elina_abstract0_t* a,
				    elina_dim_t dim,
				    size_t n)
{
 
  if (n==0){
    if (destructive){
	return a;
    }
    else {
      return elina_abstract0_copy(man,a);
    }
  }
  else {
    elina_dimension_t dimension = _elina_abstract0_dimension(a);
    if (elina_abstract0_checkman1(ELINA_FUNID_EXPAND,man,a) &&
	elina_abstract0_check_dim(ELINA_FUNID_EXPAND,man,dimension,dim)){
      void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_EXPAND];
      void* value = ptr(man,destructive,a->value,dim,n);
      return elina_abstract0_cons2(man,destructive,a,value);
    }
    else {
      if (destructive) _elina_abstract0_free(a);
      return elina_abstract0_top(man,
			      dimension.intdim + (dim<dimension.intdim ? n : 0),
			      dimension.realdim + (dim<dimension.intdim ? 0 : n));
    }
  }
}
elina_abstract0_t* elina_abstract0_fold(elina_manager_t* man,
				  bool destructive,
				  elina_abstract0_t* a,
				  elina_dim_t* tdim,
				  size_t size)
{
  elina_dimension_t dimension = _elina_abstract0_dimension(a);
  if (elina_abstract0_checkman1(ELINA_FUNID_FOLD,man,a) &&
      elina_abstract0_check_dim_array(ELINA_FUNID_FOLD,man,dimension,tdim,size)){
    if (size==0){
	elina_manager_raise_exception(man,
				   ELINA_EXC_INVALID_ARGUMENT,
				   ELINA_FUNID_FOLD,
				   "The array of dimension is empty");
	goto _elina_abstract0_fold_exc;
    }
    /* Check also that the array is sorted and contains only integer or real
       dimensions */
    size_t i;
    for (i=1;i<size; i++){
      if (tdim[i-1]>=tdim[i]){
	elina_manager_raise_exception(man,
				   ELINA_EXC_INVALID_ARGUMENT,
				   ELINA_FUNID_FOLD,
				   "The array of dimension is not sorted");
	goto _elina_abstract0_fold_exc;
      }
    }
    if (tdim[0]<dimension.intdim && tdim[size-1]>=dimension.intdim){
      elina_manager_raise_exception(man,
				 ELINA_EXC_INVALID_ARGUMENT,
				 ELINA_FUNID_FOLD,
				 "Mixed integer and real dimensions in the array");
      goto _elina_abstract0_fold_exc;
    }
    /* OK now */
    if (size==1){
      if (destructive){
	return a;
      }
      else {
	void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_COPY];
	return elina_abstract0_cons(man,ptr(man,a->value));
      }
    }
    else {
      void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_FOLD];
      void* value = ptr(man,destructive,a->value,tdim,size);
      return elina_abstract0_cons2(man,destructive,a,value);
    }
  }
  else {
  _elina_abstract0_fold_exc:
    if (destructive) _elina_abstract0_free(a);
    return elina_abstract0_top(man,
			    dimension.intdim - ( (size>0 && tdim[0]<dimension.intdim) ? (size-1) : 0),
			    dimension.realdim - ( (size>0 && tdim[0]<dimension.intdim) ? 0 : (size-1)));
  }
}
/* ============================================================ */
/* III.6 Widening */
/* ============================================================ */
elina_abstract0_t* elina_abstract0_widening(elina_manager_t* man,
				      elina_abstract0_t* a1, elina_abstract0_t* a2)
{
  if (elina_abstract0_checkman2(ELINA_FUNID_WIDENING,man,a1,a2) &&
      elina_abstract0_check_abstract2(ELINA_FUNID_WIDENING,man,a1,a2)){
    void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_WIDENING];
    void* value = ptr(man,a1->value,a2->value);
    return elina_abstract0_cons(man,value);
  }
  else {
    elina_dimension_t dimension = _elina_abstract0_dimension(a1);
    return elina_abstract0_top(man,
			    dimension.intdim,
			    dimension.realdim);
  }
}

/* ============================================================ */
/* III.7 Closure operation */
/* ============================================================ */
elina_abstract0_t* elina_abstract0_closure(elina_manager_t* man, bool destructive, elina_abstract0_t* a)
{
  elina_dimension_t dimension = _elina_abstract0_dimension(a);
  if (elina_abstract0_checkman1(ELINA_FUNID_CLOSURE,man,a)){
    void* (*ptr)(elina_manager_t*,...) = man->funptr[ELINA_FUNID_CLOSURE];
    void* value = ptr(man,destructive,a->value);
    return elina_abstract0_cons2(man,destructive,a,value);
  }
  else {
    if (destructive) _elina_abstract0_free(a);
    return elina_abstract0_top(man,
			    dimension.intdim,
			    dimension.realdim);
  }
}

/* ********************************************************************** */
/* IV. Functions offered by the ELINA interface */
/* ********************************************************************** */

/* These functions do not correspond to functions in the underlying library. */


/*
   These two functions implement of_lincons/tcons_array constructors
   using top and meet_lincons/tcons_array operations.
*/
elina_abstract0_t* elina_abstract0_of_lincons_array(elina_manager_t* man,
					      size_t intdim, size_t realdim,
					      elina_lincons0_array_t* array)
{
  elina_abstract0_t* res = elina_abstract0_top(man,intdim,realdim);
  res = elina_abstract0_meet_lincons_array(man,true,res,array);
  return res;
}
elina_abstract0_t* elina_abstract0_of_tcons_array(elina_manager_t* man,
				  size_t intdim, size_t realdim,
				  elina_tcons0_array_t* array)
{
  elina_abstract0_t* res = elina_abstract0_top(man,intdim,realdim);
  res = elina_abstract0_meet_tcons_array(man,true,res,array);
  return res;
}

/*
  These four functions implement assignment and substitution of a single
  dimension by a single expression.
*/
elina_abstract0_t* elina_abstract0_assign_linexpr(elina_manager_t* man,
					    bool destructive,
					    elina_abstract0_t* a,
					    elina_dim_t dim, elina_linexpr0_t* expr,
					    elina_abstract0_t* dest)
{
  return elina_abstract0_asssub_linexpr_array(ELINA_FUNID_ASSIGN_LINEXPR_ARRAY,
					   man,destructive,a,&dim,&expr,1,dest);
}
elina_abstract0_t* elina_abstract0_substitute_linexpr(elina_manager_t* man,
						bool destructive,
						elina_abstract0_t* a,
						elina_dim_t dim, elina_linexpr0_t* expr,
						elina_abstract0_t* dest)
{
  return elina_abstract0_asssub_linexpr_array(ELINA_FUNID_SUBSTITUTE_LINEXPR_ARRAY,
					   man,destructive,a,&dim,&expr,1,dest);
}
elina_abstract0_t* elina_abstract0_assign_texpr(elina_manager_t* man,
					  bool destructive,
					  elina_abstract0_t* a,
					  elina_dim_t dim, elina_texpr0_t* expr,
					  elina_abstract0_t* dest)
{
  return elina_abstract0_asssub_texpr_array(ELINA_FUNID_ASSIGN_TEXPR_ARRAY,
					 man,destructive,a,&dim,&expr,1,dest);
}
elina_abstract0_t* elina_abstract0_substitute_texpr(elina_manager_t* man,
					      bool destructive,
					      elina_abstract0_t* a,
					      elina_dim_t dim, elina_texpr0_t* expr,
					      elina_abstract0_t* dest)
{
  return elina_abstract0_asssub_texpr_array(ELINA_FUNID_SUBSTITUTE_TEXPR_ARRAY,
					 man,destructive,a,&dim,&expr,1,dest);
}

/* Applying an elina_dimchange2_t transformation (dimension adding followed by
   dimension removal/projection).  If project is true, the newly added
   dimension are projected on their 0-hyperplane. */
elina_abstract0_t* elina_abstract0_apply_dimchange2(elina_manager_t* man,
					      bool destructive,
					      elina_abstract0_t* a, elina_dimchange2_t* dimchange2,
					      bool project)
{
  elina_abstract0_t* res;
  if (dimchange2->add){
    res = elina_abstract0_add_dimensions(man,destructive,a,dimchange2->add,project);
    if (dimchange2->remove){
      res = elina_abstract0_remove_dimensions(man,true,res,dimchange2->remove);
    }
  }
  else if (dimchange2->remove){
    res = elina_abstract0_remove_dimensions(man,destructive,a,dimchange2->remove);
  }
  else {
    res = destructive ? a : elina_abstract0_copy(man,a);
  }
  return res;
}

/* This function implements widening with threshold, relying on the
   widening, sat_lincons and meet_lincons_array operations.
*/
elina_abstract0_t* elina_abstract0_widening_threshold(elina_manager_t* man,
						elina_abstract0_t* a1,
						elina_abstract0_t* a2,
						elina_lincons0_array_t* array)
{
  void* (*ptr)(elina_manager_t*,...);
  bool (*ptr2)(elina_manager_t*,...);
  void* value;
  size_t i,j,size;
  elina_lincons0_t tmp;

  elina_dimension_t dimension = _elina_abstract0_dimension(a1);
  if (elina_abstract0_checkman2(ELINA_FUNID_WIDENING,man,a1,a2) &&
      elina_abstract0_check_abstract2(ELINA_FUNID_WIDENING,man,a1,a2) &&
      elina_abstract0_check_lincons_array(ELINA_FUNID_WIDENING,man,dimension,array) ){
    ptr = man->funptr[ELINA_FUNID_WIDENING];
    value = ptr(man,a1->value,a2->value);

    ptr2 = man->funptr[ELINA_FUNID_SAT_LINCONS];
    size = array->size;
    i = j = 0;
    while (i<size-j){
      if (ptr2(man,a2->value,&array->p[i])){
	i++;
      }
      else {
	j++;
	tmp = array->p[i];
	array->p[i] = array->p[array->size-j];
	array->p[array->size-j] = tmp;
      }
    }
    if (i>0){
      array->size = i;
      ptr = man->funptr[ELINA_FUNID_MEET_LINCONS_ARRAY];
      value = ptr(man,true,value,array);
      array->size = size;
    }
    return elina_abstract0_cons(man,value);
  }
  else {
    return elina_abstract0_top(man,
			    dimension.intdim,
			    dimension.realdim);
  }
}

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
/* opt_pk_user.c: conversions with interface datatypes */
/* ********************************************************************** */


#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"
#include "opt_pk_internal.h"
#include "opt_pk_user.h"

static void elina_coeff_set_scalar_numint(elina_coeff_t* coeff, opt_numint_t *num){
  elina_coeff_reinit(coeff,ELINA_COEFF_SCALAR,ELINA_SCALAR_MPQ);
  mpq_set_numint(coeff->val.scalar->val.mpq,num);
}

/* ********************************************************************** */
/* From ITV to PK */
/* ********************************************************************** */

/* Fills the vector with the constraint:
   dim <= num if sgn>0,
   dim = num if sgn=0
   dim >= -num if sgn<0 (an inferior bound is implicitly negated, as in itv).

   Returns false if equality of an integer dimension with a non-integer numbers
*/

bool opt_vector_set_dim_bound(opt_pk_internal_t* opk,
			  opt_numint_t* vec,
			  elina_dim_t dim,
			  numrat_t numrat,
			  int mode,
			  size_t intdim, size_t realdim,
			  bool integer)
{
  numrat_t bound;
  size_t size;

  size = opk->dec+intdim+realdim;

  numrat_init(bound);
  if (integer && dim<intdim){
    if (mode>0){
      numint_fdiv_q(numrat_numref(bound),
		    numrat_numref(numrat),numrat_denref(numrat));
      numint_set_int(numrat_denref(bound),1);
    }
    else if (mode<0){
      numint_cdiv_q(numrat_numref(bound),
		    numrat_numref(numrat),numrat_denref(numrat));
      numint_set_int(numrat_denref(bound),1);
    }
    else {
      if (numint_cmp_int(numrat_denref(numrat),1)!=0){
	numrat_clear(bound);
	return false;
      }
      numrat_set(bound,numrat);
    }
  } else {
    numrat_set(bound,numrat);
  }
  /* Write the constraint num + den*x >= 0 */
  opt_vector_clear(vec,size);
  vec[0] =  (mode ? 1 : 0);
  numint_set(vec + opt_polka_cst, numrat_numref(bound));
  numint_set(vec+opk->dec+dim, numrat_denref(bound));
  numrat_clear(bound);
  /* put the right sign now */
  if (mode>=0){
    vec[opk->dec+dim] = -vec[opk->dec+dim];
  }
  return true;
}


/* Fills the vector with the constraint:
   expr <= bound if sgn>0,
   expr = bound if sgn=0
   expr >= -bound if sgn<0
   (an inferior bound is implicitly negated, as in itv).

   Returns false if equality with a non-integer numbers (like 2x=1).
*/

/*bool vector_set_linexpr_bound(pk_internal_t* pk,
			      numint_t* vec,
			      numint_t* vec2,
			      numrat_t numrat,
			      int mode,
			      size_t intdim, size_t realdim,
			      bool integer)
{
  size_t i;
  size_t size;
  numint_t cst;

  assert (numint_sgn(numrat_denref(numrat)));

  size = pk->dec+intdim+realdim;

  if (vec!=vec2){
    vector_copy(vec,vec2,size);
  }

  if (numint_cmp_int(numrat_denref(numrat),1) != 0){
    for (i=1; i<size; i++){
      numint_mul(vec[i],vec[i],numrat_denref(numrat));
    }
  }
  numint_init(cst);
  numint_mul(cst,numrat_numref(numrat),vec[0]);
  if (mode>=0){
    numint_sub(vec[polka_cst],vec[polka_cst],cst);
    for (i=pk->dec; i<size; i++){
      numint_neg(vec[i],vec[i]);
    }
  }
  else {
    numint_add(vec[polka_cst],vec[polka_cst],cst);
  }
  numint_set_int(vec[0], mode ? 1 : 0);
  vector_normalize(pk,vec,size);
  if (integer) vector_normalize_constraint_int(pk,vec,intdim,realdim);

  return true;
}*/

/* Fills the vector with the quasi-linear expression (itv_linexpr) */
void opt_vector_set_itv_linexpr(opt_pk_internal_t* opk,
			    opt_numint_t* ov,
			    itv_linexpr_t* expr,
			    size_t dim,
			    int mode)
{
  size_t i,d;
  bool* peq;
  itv_ptr pitv;
  /* compute lcm of denominators, in vec[0] */
  if (mode>=0){
    assert(!bound_infty(expr->cst->sup));
    if (bound_sgn(expr->cst->sup)){
      numint_set(ov,
		 numrat_denref(bound_numref(expr->cst->sup)));
    }
    else{
      ov[0] = 1;
    }
    
  } else {
    assert(!bound_infty(expr->cst->inf));
    if (bound_sgn(expr->cst->inf))
      ov[0] =  numrat_denref(bound_numref(expr->cst->inf));
    else
      ov[0] = 1;
  }
  
  itv_linexpr_ForeachLinterm(expr,i,d,pitv,peq){
    assert(*peq);
    numint_lcm(ov,ov,numrat_denref(bound_numref(pitv->sup)));
  }
  
  /* Fill the vector */
  if (opk->strict) ov[opt_polka_eps] = 0;
  /* constant coefficient */
  if (mode>=0){
    numint_divexact(ov + opt_polka_cst,
		    ov,
		    numrat_denref(bound_numref(expr->cst->sup)));
    numint_mul(ov + opt_polka_cst,
	       ov + opt_polka_cst,
	       numrat_numref(bound_numref(expr->cst->sup)));
  } else {
    numint_divexact(ov + opt_polka_cst,
		    ov,
		    numrat_denref(bound_numref(expr->cst->inf)));
    numint_mul(ov + opt_polka_cst,
	       ov + opt_polka_cst,
	       numrat_numref(bound_numref(expr->cst->inf)));
    numint_neg(ov + opt_polka_cst,ov + opt_polka_cst);
  }
  
  /* Other coefficients */
  for (i=opk->dec;i<opk->dec+dim; i++){
       ov[i] = 0;
  }
  itv_linexpr_ForeachLinterm(expr,i,d,pitv,peq){
     size_t index = opk->dec + d;
     numint_divexact(ov + index,ov,numrat_denref(bound_numref(pitv->sup)));
     numint_mul(ov+index,ov+index,numrat_numref(bound_numref(pitv->sup)));
  }
  //opt_vector_print(ov,opk->dec+dim);
  return;
}

/* Fills the vector(s) with the fully linear constraint cons */
void opt_vector_set_itv_lincons(opt_pk_internal_t* opk,
			    opt_numint_t* ov,
			    itv_lincons_t* cons,
			    size_t intdim, size_t realdim,
			    bool integer)
{
  size_t i,nb;
  assert(cons->constyp == ELINA_CONS_EQ ||
	 cons->constyp == ELINA_CONS_SUPEQ ||
	 cons->constyp == ELINA_CONS_SUP);
  assert(itv_linexpr_is_scalar(&cons->linexpr));

  opt_vector_set_itv_linexpr(opk, ov, &cons->linexpr, intdim+realdim,1);
  opt_vector_normalize(opk,ov,opk->dec+intdim+realdim);
  if (cons->constyp == ELINA_CONS_EQ){
    ov[0] = 0;
  }
  else {
    ov[0] = 1;
  }
  if (cons->constyp == ELINA_CONS_SUP){
    if (opk->strict){
      ov[opt_polka_eps] = -1;
    }
    else if (integer && opt_vector_is_integer(opk, ov, intdim, realdim)){
      ov[opt_polka_cst]= ov[opt_polka_cst] - 1;
    }
  }
  if (integer)
    opt_vector_normalize_constraint_int(opk,ov,intdim,realdim);
  
}

/* Fills the vector(s) with the fully linear constraint cons for testing
   satisfiability.

   Returns false if unsatisfiable
 */
bool opt_vector_set_itv_lincons_sat(opt_pk_internal_t* opk,
				opt_numint_t* ov,
				itv_lincons_t* cons,
				size_t intdim, size_t realdim,
				bool integer)
{
  bool sat;
  size_t i;

  if (cons->constyp == ELINA_CONS_EQ && cons->linexpr.equality != true){
    return false;
  }

  assert(cons->constyp == ELINA_CONS_EQ ||
	 cons->constyp == ELINA_CONS_SUPEQ ||
	 cons->constyp == ELINA_CONS_SUP);

  if (!bound_infty(cons->linexpr.cst->inf)){
    opt_vector_set_itv_linexpr(opk, ov, &cons->linexpr, intdim+realdim,-1);
    opt_vector_normalize(opk,ov,opk->dec+intdim+realdim);
    if (cons->constyp == ELINA_CONS_EQ && cons->linexpr.equality){
      ov[0] = 0;
    }
    else {
      ov[0] = 1;
    }
    if (cons->constyp == ELINA_CONS_SUP){
      if (opk->strict){
	ov[opt_polka_eps] = -1;
      }
      else if (integer && opt_vector_is_integer(opk, ov, intdim, realdim)){
	ov[opt_polka_cst] =  ov[opt_polka_cst] - 1;
      }
    }
    if (integer)
      opt_vector_normalize_constraint_int(opk,ov,intdim,realdim);
    return true;
  }
  else {
    return false;
  }
 }

/* ********************************************************************** */
/* From APRON to PK */
/* ********************************************************************** */


/*
static
bool matrix_append_elina_lincons0_array(pk_internal_t* pk,
				     matrix_t* mat,
				     size_t** tabindex, size_t* size,
				     elina_lincons0_array_t* array,
				     size_t intdim, size_t realdim,
				     bool integer)
{
  bool exact,res;
  size_t nbrows,i,j,k;
  size_t* tab;

  nbrows = mat->nbrows;
  matrix_resize_rows_lazy(mat,nbrows+2*array->size);

  res = true;
  tab = NULL;
  j = nbrows;
  k = 0;
  for (i=0; i<array->size; i++){
    if (elina_lincons0_is_unsat(&array->p[i])){
      if (tab){
	free(tab); tab=NULL;
	k=0;
      }
      vector_clear(mat->p[0],mat->nbcolumns);
      numint_set_int(mat->p[0][0],1);
      numint_set_int(mat->p[0][polka_cst],-1);
      j = 1;
      res = true;
      break;
    }
    switch (array->p[i].constyp){
    case ELINA_CONS_EQ:
    case ELINA_CONS_SUPEQ:
    case ELINA_CONS_SUP:
      if (elina_linexpr0_is_quasilinear(array->p[i].linexpr0)){
	exact = itv_lincons_set_elina_lincons0(pk->itv,
					    &pk->poly_itv_lincons,
					    &array->p[i]);
	res = res && exact;
	j += vector_set_itv_lincons(pk,
				    &mat->p[j], &pk->poly_itv_lincons,
				    intdim,realdim,integer);
      }
      else {
	if (tab==NULL){
	  tab = (size_t*)malloc(array->size*sizeof(size_t));
	}
	tab[k] = i;
	k++;
      }
      break;
    default:
      res = false;
      break;
    }
  }
  mat->nbrows = j;
  if (tab){
    tab = (size_t*)realloc(tab,k*sizeof(size_t));
  }
  *tabindex = tab;
  *size = k;
  return res;
}

bool matrix_set_elina_lincons0_array(pk_internal_t* pk,
				  matrix_t** mat,
				  size_t** tabindex, size_t* size,
				  elina_lincons0_array_t* array,
				  size_t intdim, size_t realdim,
				  bool integer)
{
  *mat = matrix_alloc(2*array->size,pk->dec+intdim+realdim,false);
  (*mat)->nbrows = 0;
  return matrix_append_elina_lincons0_array(pk,
					 *mat,
					 tabindex,size,
					 array,
					 intdim,realdim,integer);
}
static
bool matrix_append_elina_intlincons0_array(pk_internal_t* pk,
					matrix_t* mat,
					itv_t* titv,
					elina_lincons0_array_t* array,
					size_t* tab, size_t size,
					size_t intdim, size_t realdim,
					bool integer)
{
  bool exact;
  size_t nbrows,i,j;

  nbrows = mat->nbrows;
  matrix_resize_rows_lazy(mat,nbrows+2*array->size);
  exact = true;
  j = nbrows;
  for (i=0; i<size; i++){
    size_t index = tab[i];

    exact = itv_lincons_set_elina_lincons0(pk->itv,
					&pk->poly_itv_lincons,
					&array->p[index])
      && exact;
    exact = itv_quasilinearize_lincons(pk->itv,&pk->poly_itv_lincons,titv,true)
      && exact;
    j += vector_set_itv_lincons(pk,&mat->p[j],&pk->poly_itv_lincons,
				intdim,realdim,integer);
  }
  mat->nbrows = j;
  return exact;
}

bool matrix_set_elina_intlincons0_array(pk_internal_t* pk,
				     matrix_t** mat,
				     itv_t* titv,
				     elina_lincons0_array_t* array,
				     size_t* tab, size_t size,
				     size_t intdim, size_t realdim,
				     bool integer)
{
  *mat = matrix_alloc(2*array->size,pk->dec+intdim+realdim,false);
  (*mat)->nbrows = 0;
  return matrix_append_elina_intlincons0_array(pk,
					    *mat,
					    titv,
					    array,tab,size,
					    intdim,realdim,integer);
}
*/
static
bool opt_matrix_append_itv_lincons_array(opt_pk_internal_t* opk,
				     opt_matrix_t* mat,
				     itv_lincons_array_t* array,
				     size_t intdim, size_t realdim,
				     bool integer)
{
  bool exact,res;
  size_t nbrows,i,j;
  size_t* tab;

  nbrows = mat->nbrows;
  opt_matrix_resize_rows_lazy(mat,nbrows+array->size);
  
  res = true;
  j = nbrows;
  for (i=0; i<array->size; i++){
    assert(itv_linexpr_is_scalar(&array->p[i].linexpr));
    switch (array->p[i].constyp){
    case ELINA_CONS_EQ:
    case ELINA_CONS_SUPEQ:
    case ELINA_CONS_SUP:
      opt_vector_set_itv_lincons(opk,
			     mat->p[j], &array->p[i],
			     intdim,realdim,integer);
      j++;
      break;
    default:
      res = false;
      break;
    }
  }
  mat->nbrows = j;
  return res;
}

bool opt_matrix_set_itv_lincons_array(opt_pk_internal_t* opk,
				  opt_matrix_t** mat,
				  itv_lincons_array_t* array,
				  size_t intdim, size_t realdim,
				  bool integer)
{


  *mat = opt_matrix_alloc(array->size,opk->dec+intdim+realdim,false);
  (*mat)->nbrows = 0;
  return opt_matrix_append_itv_lincons_array(opk,
					 *mat,array,
					 intdim,realdim,integer);
}

/* ********************************************************************** */
/* From PK to APRON */
/* ********************************************************************** */

elina_lincons0_t opt_lincons0_of_vector(opt_pk_internal_t* opk,
				 opt_numint_t* ov, unsigned short int * ca,
				 unsigned short int vsize, unsigned short int size)
{
  elina_lincons0_t lincons;
  elina_linexpr0_t* linexpr;
  size_t i;
  linexpr = elina_linexpr0_alloc(ELINA_LINEXPR_DENSE, size - opk->dec);
  elina_coeff_set_scalar_numint(&linexpr->cst, ov + opt_polka_cst);
  for(i=opk->dec; i < size; i++){
    elina_dim_t dim = i - opk->dec;
    elina_coeff_set_scalar_int(&linexpr->p.coeff[dim], 0);
  }
  for (i=0; i<vsize; i++){
    elina_dim_t dim = ca[i] - opk->dec;
    elina_coeff_set_scalar_numint(&linexpr->p.coeff[dim], ov + i+ opk->dec);
  }
 
  if (ov[0]){
    if (opk->strict && numint_sgn(ov + opt_polka_eps)<0)
      lincons.constyp = ELINA_CONS_SUP;
    else
      lincons.constyp = ELINA_CONS_SUPEQ;
  }
  else {
    lincons.constyp = ELINA_CONS_EQ;
  }
  lincons.linexpr0 = linexpr;
  lincons.scalar = NULL;
  return lincons;
}



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
/* opt_pk_user.c: conversions with interface datatypes */
/* ********************************************************************** */


#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"
#include "opt_pk_internal.h"
#include "opt_pk_user.h"

static void elina_coeff_set_scalar_numint(elina_coeff_t* coeff, opt_numint_t num){
  elina_coeff_reinit(coeff,ELINA_COEFF_SCALAR,ELINA_SCALAR_MPQ);
  mpq_set_si(coeff->val.scalar->val.mpq,num,1);
}



bool opt_vector_set_dim_bound(opt_pk_internal_t* opk,
			  opt_numint_t* vec,
			  elina_dim_t dim,
			  elina_scalar_t *scalar,
			  int mode,
			  size_t intdim, size_t realdim,
			  bool integer)
{
  elina_rat_t * rat = elina_scalar_set_rat(scalar);
  size_t size;
  if(mode < 0){
	rat->n= -rat->n;
  }
  size = opk->dec+intdim+realdim;

  if (integer && dim<intdim){
    if (mode>0){
      rat->n = elina_int_fdiv_q(rat->n,rat->d);
      rat->d = 1;
    }
    else if (mode<0){
      rat->n = elina_int_cdiv_q(rat->n,rat->d);
      rat->d = 1;
    }
    else {
      if (rat->d!=1){
	free(rat);
	return false;
      }
    }
  } 
  /* Write the constraint num + den*x >= 0 */
  opt_vector_clear(vec,size);
  vec[0] =  (mode ? 1 : 0);
  vec[opt_polka_cst] = rat->n;
  vec[opk->dec+dim] = rat->d;
  free(rat);
  /* put the right sign now */
  if (mode>=0){
    vec[opk->dec+dim]= -vec[opk->dec+dim];
  }
  /* put the right sign now */
  return true;
}


/* Fills the vector with the quasi-linear expression (elina_linexpr0) */
void opt_vector_set_elina_linexpr0(opt_pk_internal_t* opk,
			    opt_numint_t* ov,
			    elina_linexpr0_t* expr,
			    size_t dim,
			    int mode)
{
  size_t i;
  bool* peq;
  elina_dim_t d;
  elina_coeff_t* coeff;
  elina_scalar_t * scalar=NULL;
  elina_coeff_t * cst = &(expr->cst);
  /* compute lcm of denominators, in vec[0] */
  switch(cst->discr){
	case ELINA_COEFF_SCALAR:
	    scalar = cst->val.scalar;
	    assert(!elina_scalar_infty(scalar));
	    if (elina_scalar_sgn(scalar)){
		elina_rat_set_elina_scalar(opk->poly_numrat, scalar);
		ov[0] = opk->poly_numrat->d;
	    }
	    else{
	      ov[0] = 1;
	    }
	    break;  
	case ELINA_COEFF_INTERVAL:	
	   if (mode>=0){  
	       scalar = cst->val.interval->sup;
	       assert(!elina_scalar_infty(scalar));
	       if (elina_scalar_sgn(scalar)){
		   elina_rat_set_elina_scalar(opk->poly_numrat,scalar);
		   ov[0] = opk->poly_numrat->d;
	       }
	       else{
	          ov[0] = 1;
	       }
	    } 
	  else {
	       scalar = cst->val.interval->inf;
	       assert(!elina_scalar_infty(scalar));
	       if (elina_scalar_sgn(scalar)){
	           elina_rat_set_elina_scalar(opk->poly_numrat,scalar);
	           ov[0] =  opk->poly_numrat->d;
	       }
	       else{
	           ov[0] = 1;
	       }
	  }
	 break;
  }
  
   elina_linexpr0_ForeachLinterm(expr,i,d,coeff){
     assert(coeff.discr == ELINA_COEFF_SCALAR);
     scalar = coeff->val.scalar;
     elina_rat_set_elina_scalar(opk->poly_numrat, scalar);
     ov[0] = elina_int_lcm(ov[0], opk->poly_numrat->d);
  }
  
  /* Fill the vector */
  if (opk->strict) ov[opt_polka_eps] = 0;
  /* constant coefficient */
  switch(cst->discr){
	  case ELINA_COEFF_SCALAR:
		scalar = cst->val.scalar;
	  break;
	  case ELINA_COEFF_INTERVAL:
		  if (mode>=0){
		      scalar = cst->val.interval->sup;
		  } 
		  else {
		      scalar = cst->val.interval->inf;
		  }
	  break;
  }
  elina_rat_set_elina_scalar(opk->poly_numrat,scalar);
  ov[opt_polka_cst] = ov[0] /opk->poly_numrat->d;
  ov[opt_polka_cst] = ov[opt_polka_cst] * opk->poly_numrat->n;
  /* Other coefficients */
  for (i=opk->dec;i<opk->dec+dim; i++){
       ov[i] = 0;
  }
  elina_linexpr0_ForeachLinterm(expr,i,d,coeff){
     size_t index = opk->dec + d;
     scalar = coeff->val.scalar;    
     elina_rat_set_elina_scalar(opk->poly_numrat,scalar);
     ov[index] = ov[0]/opk->poly_numrat->d;
     ov[index] = ov[index]*opk->poly_numrat->n;
  }
  //opt_vector_print(ov,opk->dec+dim);
  return;
}

/* Fills the vector(s) with the fully linear constraint cons */
void opt_vector_set_elina_lincons0(opt_pk_internal_t* opk,
			    opt_numint_t* ov,
			    elina_lincons0_t* cons,
			    size_t intdim, size_t realdim,
			    bool integer)
{
  size_t i,nb;
  assert(cons->constyp == ELINA_CONS_EQ ||
	 cons->constyp == ELINA_CONS_SUPEQ ||
	 cons->constyp == ELINA_CONS_SUP);
  assert(elina_linexpr0_is_linear(&cons->linexpr));

  opt_vector_set_elina_linexpr0(opk, ov, cons->linexpr0, intdim+realdim,1);
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
bool opt_vector_set_elina_lincons0_sat(opt_pk_internal_t* opk,
				opt_numint_t* ov,
				elina_lincons0_t* cons,
				size_t intdim, size_t realdim,
				bool integer)
{
  bool sat;
  size_t i;
  elina_coeff_t *cst = &(cons->linexpr0->cst);
  if (cons->constyp == ELINA_CONS_EQ && cst->discr!=ELINA_COEFF_SCALAR){
    return false;
  }

  assert(cons->constyp == ELINA_CONS_EQ ||
	 cons->constyp == ELINA_CONS_SUPEQ ||
	 cons->constyp == ELINA_CONS_SUP);
  elina_scalar_t *scalar;
  
  if(cst->discr==ELINA_COEFF_SCALAR){
	scalar = cst->val.scalar;
  }
  else{
	scalar = cst->val.interval->inf;
  }
  if (!elina_scalar_infty(scalar)){
    opt_vector_set_elina_linexpr0(opk, ov, cons->linexpr0, intdim+realdim,-1);
    opt_vector_normalize(opk,ov,opk->dec+intdim+realdim);
    if (cons->constyp == ELINA_CONS_EQ && cst->discr==ELINA_COEFF_SCALAR){
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


static
bool opt_matrix_append_elina_lincons0_array(opt_pk_internal_t* opk,
				     opt_matrix_t* mat,
				     elina_lincons0_array_t* array,
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
    assert(elina_linexpr0_is_linear(&array->p[i].linexpr0));
    switch (array->p[i].constyp){
    case ELINA_CONS_EQ:
    case ELINA_CONS_SUPEQ:
    case ELINA_CONS_SUP:
      opt_vector_set_elina_lincons0(opk,
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

bool opt_matrix_set_elina_lincons0_array(opt_pk_internal_t* opk,
				  opt_matrix_t** mat,
				  elina_lincons0_array_t* array,
				  size_t intdim, size_t realdim,
				  bool integer)
{


  *mat = opt_matrix_alloc(array->size,opk->dec+intdim+realdim,false);
  (*mat)->nbrows = 0;
  return opt_matrix_append_elina_lincons0_array(opk,
					 *mat,array,
					 intdim,realdim,integer);
}

elina_lincons0_t opt_lincons0_of_vector(opt_pk_internal_t* opk,
				 opt_numint_t* ov, unsigned short int * ca,
				 unsigned short int vsize, unsigned short int size)
{
  elina_lincons0_t lincons;
  elina_linexpr0_t* linexpr;
  size_t i;
  linexpr = elina_linexpr0_alloc(ELINA_LINEXPR_DENSE, size - opk->dec);
  elina_coeff_set_scalar_numint(&linexpr->cst, ov[opt_polka_cst]);
  for(i=opk->dec; i < size; i++){
    elina_dim_t dim = i - opk->dec;
    elina_coeff_set_scalar_int(&linexpr->p.coeff[dim], 0);
  }
  for (i=0; i<vsize; i++){
    elina_dim_t dim = ca[i] - opk->dec;
    elina_coeff_set_scalar_numint(&linexpr->p.coeff[dim], ov[opk->dec+i]);
  }
 
  if (ov[0]){
    if (opk->strict && elina_int_sgn(ov[opt_polka_eps])<0)
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




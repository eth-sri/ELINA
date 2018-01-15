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

/* *******************************************************************************    */
/* elina_linearize_texpr.c: functions for (quasi)linearisation of tree expressions   */
/* *******************************************************************************   */


#include "elina_linearize_texpr.h"


/***************************************
	Allocation and deallocation of auxiliary data structures
****************************************/

static
bool elina_intlinearize_alloc(elina_manager_t* man, elina_abstract0_t* abs,
			elina_interval_t*** penv, elina_dimension_t* pdim, elina_scalar_discr_t discr)
{
  size_t size,i;
  
  assert(!elina_abstract0_is_bottom(man,abs));
  *penv = elina_abstract0_to_box(man,abs);
  elina_interval_t **box = *penv;
  *pdim = elina_abstract0_dimension(man,abs);
  size = pdim->intdim+pdim->realdim;
  for(i=0; i < size; i++){
	elina_interval_convert(box[i],discr);
  }
  bool exact = man->result.flag_exact;
  return exact;
}

static inline
void elina_intlinearize_free(elina_interval_t* env, elina_dimension_t dim)
{
  if (env){
	size_t i, size = dim.intdim + dim.realdim;
	for(i=0; i <size; i++){
	     elina_interval_free(env+ i);
	}
  }
}


/* ********************************************************************** */
/* Evaluation of tree expressions */
/* ********************************************************************** */



/* Evaluates node into interval res, assuming operator arguments are arg1 and
   (for binary operator) arg2
*/

/* ********************************************************************** */
/* General Rounding */
/* ********************************************************************** */
void elina_interval_round(elina_interval_t *res, elina_interval_t* arg,
		      elina_texpr_rtype_t t, elina_texpr_rdir_t d, elina_scalar_discr_t discr)
{
  switch (t) {

  case ELINA_RTYPE_REAL:
    if (&res!=&arg) elina_interval_set(res,arg);
    break;

  case ELINA_RTYPE_INT:
    switch (d) {
    case ELINA_RDIR_ZERO:
      elina_interval_trunc(res,arg,discr);
      break;
    case ELINA_RDIR_UP:
      elina_interval_ceil(res,arg,discr);
      break;
    case ELINA_RDIR_DOWN:
      elina_interval_floor(res,arg,discr);
      break;
    case ELINA_RDIR_RND:
    case ELINA_RDIR_NEAREST: /* 'to nearest' case could be improved */
      elina_interval_to_int(res,arg,discr);
      break;
    default:
      assert(0);
    }
    break;

  case ELINA_RTYPE_SINGLE:
    /* directed rounding cases (+oo, -oo, 0) could be improved */
    elina_interval_to_float(res,arg,discr);
    break;

  case ELINA_RTYPE_QUAD:     /* 'round to quad' could be improved */
  case ELINA_RTYPE_EXTENDED: /* 'round to extended' could be improved */
  case ELINA_RTYPE_DOUBLE:
    /* directed rounding cases (+oo, -oo, 0) could be improved */
    elina_interval_to_double(res,arg,discr);
    break;

  default:
    assert(0);
  }
}


static void
elina_interval_eval_elina_texpr0_node(elina_texpr0_node_t* n,
			elina_interval_t* res, elina_interval_t *op1, elina_interval_t * op2,
			elina_scalar_discr_t discr)
{
  switch (n->op) {
  case ELINA_TEXPR_NEG:
    elina_interval_neg(res, op1);
    return; /* no rounding */
  case ELINA_TEXPR_CAST:
    elina_interval_set(res, op1);
    break;
  case ELINA_TEXPR_SQRT:
    elina_interval_sqrt(res, op1,discr);
    break;
  case ELINA_TEXPR_ADD:
    elina_interval_add(res, op1, op2,discr);
    break;
  case ELINA_TEXPR_SUB:
    elina_interval_sub(res, op1, op2, discr);
    break;
  case ELINA_TEXPR_MUL:
    elina_interval_mul(res, op1, op2,discr);
    break;
  case ELINA_TEXPR_DIV:
    elina_interval_div(res, op1, op2,discr);
    break;
  case ELINA_TEXPR_POW:
    elina_interval_pow(res, op1, op2, discr);
    break;
  case ELINA_TEXPR_MOD:
    elina_interval_mod(res, op1, op2, n->type==ELINA_RTYPE_INT,discr);
    return; /* no rounding */
  default:
    assert(0);
  }
  elina_interval_round(res,res,n->type,n->dir,discr);
}

/* evaluates expr into intervalres,
   assuming env maps dimensions to interval values */
void elina_interval_eval_elina_texpr0(elina_interval_t *res,
				elina_texpr0_t* expr,elina_scalar_discr_t discr,
				elina_interval_t** env)
{
  if (!expr) {
    elina_interval_set_bottom(res);
    return;
  }

  switch(expr->discr){
  case ELINA_TEXPR_CST:
    elina_interval_set_elina_coeff(res,&expr->val.cst);
    break;
  case ELINA_TEXPR_DIM:
    elina_interval_set(res,env[expr->val.dim]);
    break;
  case ELINA_TEXPR_NODE:
    if (expr->val.node->exprB) {
      /* binary */
      elina_interval_t * x = elina_interval_alloc();
      elina_interval_eval_elina_texpr0(x,expr->val.node->exprA,discr,env);
      elina_interval_eval_elina_texpr0(res,expr->val.node->exprB,discr,env);
      if (elina_interval_is_bottom(x) || elina_interval_is_bottom(res)){
	elina_interval_set_bottom(res);
      }
      else {
	elina_interval_eval_elina_texpr0_node(expr->val.node,res,x,res,discr);
      }
      elina_interval_free(x);
    }
    else {
      /* unary */
      elina_interval_eval_elina_texpr0(res,expr->val.node->exprA,discr,env);
      if (!elina_interval_is_bottom(res)){
	elina_interval_eval_elina_texpr0_node(expr->val.node,res,res,res,discr);
      }
    }
    break;
  default:
    assert(false);
  }
}

elina_interval_t* elina_eval_texpr0(elina_manager_t* man,
			      elina_abstract0_t* abs,
			      elina_texpr0_t* expr,
			      elina_scalar_discr_t discr,
			      bool* pexact)
{
  elina_interval_t** env = NULL;
  elina_dimension_t dim = {0,0};
  elina_interval_t* r = elina_interval_alloc();
  size_t i, size;
  if (pexact) *pexact = false;
  if (elina_texpr0_is_interval_cst(expr)){
    elina_interval_eval_elina_texpr0(r,expr,discr,NULL);
  }
  else {
    //elina_intlinearize_alloc(man,abs,env,&dim,discr);
    assert(!elina_abstract0_is_bottom(man,abs));
    env = elina_abstract0_to_box(man,abs);
  
    dim = elina_abstract0_dimension(man,abs);
    size = dim.intdim+dim.realdim;
    for(i=0; i < size; i++){
	elina_interval_convert(env[i],discr);
    }
    elina_interval_eval_elina_texpr0(r,expr,discr,env);
  }
  if(env)elina_interval_array_free(env,dim.intdim+dim.realdim);
  return r;
}

/* ********************************************************************** */
/*  Interval Linearization of linear tree expressions */
/* ********************************************************************** */

/* Linearize a tree expression that is (syntaxically) interval linear with
   exact arithmetic.

   Compared to elina_intlinearize_texpr0() function below, this functions does
   not require a bounding box for dimensions.

   If the precondition is violated, returns NULL.
*/

static bool elina_texpr0_node_intlinearize_linear(elina_texpr0_node_t* n, elina_linexpr0_t** dres, elina_scalar_discr_t discr)
{
  elina_linexpr0_t *res = *dres;
  elina_interval_t *i1;
  elina_linexpr0_t *l1;
  bool exc = n->type!=ELINA_RTYPE_REAL && n->op!=ELINA_TEXPR_NEG;
  if (exc) return exc;
  switch (n->op) {
  case ELINA_TEXPR_NEG:
    exc = elina_intlinearize_elina_texpr0_intlinear(dres,n->exprA,discr);
    elina_linexpr0_neg(res);
    break;

  case ELINA_TEXPR_CAST:
    exc = elina_intlinearize_elina_texpr0_intlinear(dres,n->exprA,discr);
    break;

  case ELINA_TEXPR_MOD:
  case ELINA_TEXPR_SQRT:
  case ELINA_TEXPR_POW:
    exc = true;
    break;

  case ELINA_TEXPR_ADD:
  case ELINA_TEXPR_SUB:
    //elina_linexpr0_init(&l1,0);
    l1 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
    /* intlinearize arguments */
    exc = elina_intlinearize_elina_texpr0_intlinear(&l1,n->exprA,discr);
    exc = elina_intlinearize_elina_texpr0_intlinear(dres,n->exprB,discr) || exc;
    /* add/sub linear form & interval */
    if (n->op==ELINA_TEXPR_ADD){
      elina_linexpr0_add(dres,&l1,dres,discr);
    }
    else{
      elina_linexpr0_sub(dres,&l1,dres,discr);
    }
    elina_linexpr0_free(l1);
    break;

  case ELINA_TEXPR_MUL:
  case ELINA_TEXPR_DIV:
    i1 = elina_interval_alloc();
    if (elina_texpr0_is_interval_cst(n->exprB)){
      exc = elina_intlinearize_elina_texpr0_intlinear(dres,n->exprA,discr);
      elina_interval_eval_elina_texpr0(i1,n->exprB,discr,NULL);
    }
    else if (n->op == ELINA_TEXPR_MUL && elina_texpr0_is_interval_cst(n->exprA)){
      exc = elina_intlinearize_elina_texpr0_intlinear(dres,n->exprB,discr);
      elina_interval_eval_elina_texpr0(i1,n->exprA,discr,NULL);
    }
    else {
      exc = true;
      break;
    }
    if (n->op==ELINA_TEXPR_DIV){
      elina_interval_t *i2 = elina_interval_alloc();
      elina_interval_set_to_int(i2,1,1,discr);
      elina_interval_div(i1,i2,i1,discr);
      elina_interval_free(i2);
    }
	
    elina_linexpr0_scale(res,i1,discr);
    elina_interval_free(i1);
    break;

  default:
    assert(0);
  }
  return exc;
}



bool elina_intlinearize_elina_texpr0_intlinear(elina_linexpr0_t** dres, elina_texpr0_t* expr, elina_scalar_discr_t discr)
{
  elina_linexpr0_t *res = *dres;
  bool exc = false;
  assert(expr);
  elina_coeff_t * cst;
  elina_interval_t *interval;
  elina_scalar_t * scalar,*inf,*sup;
  elina_linterm_t * term;
  switch(expr->discr){
  case ELINA_TEXPR_CST:
    elina_linexpr0_reinit(res,0);
    elina_coeff_set(&res->cst,&expr->val.cst);
    cst = &res->cst;
    break;
  case ELINA_TEXPR_DIM:
    elina_linexpr0_reinit(res,1);
    elina_coeff_reinit(&res->cst,ELINA_COEFF_SCALAR,discr);
    elina_scalar_set_to_int(res->cst.val.scalar,0,discr);
    term = &res->p.linterm[0];
    term->dim = expr->val.dim;
    elina_coeff_reinit(&term->coeff,ELINA_COEFF_SCALAR,discr);
    elina_scalar_set_to_int(term->coeff.val.scalar,1,discr);
    break;
  case ELINA_TEXPR_NODE:
    if (elina_texpr0_is_interval_cst(expr)){
      elina_linexpr0_reinit(res,0);
      cst = &res->cst;
      if(cst->discr==ELINA_COEFF_SCALAR){
	scalar = elina_scalar_alloc();
	elina_scalar_set(scalar,cst->val.scalar);
	elina_coeff_reinit(cst,ELINA_COEFF_INTERVAL,scalar->discr);
	elina_coeff_set_interval_scalar(cst,scalar,scalar);
	elina_scalar_free(scalar);
      }
       interval = cst->val.interval;
       elina_interval_eval_elina_texpr0(interval,expr,discr,NULL);
       inf = interval->inf;
       sup = interval->sup;	
	if(!elina_scalar_infty(inf) && elina_scalar_equal(inf,sup)){
		elina_coeff_reduce(cst);	
	}
    }
    else {
      exc = elina_texpr0_node_intlinearize_linear(expr->val.node,dres,discr);
    }
    break;
  default:
    exc = true;
    assert(false);
  }
  return exc;
}


elina_linexpr0_t* elina_intlinearize_texpr0_intlinear(elina_manager_t* man,
						elina_texpr0_t* expr,
						elina_scalar_discr_t discr)
{
  bool exc;
  elina_linexpr0_t* res;
  res = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
  exc = elina_intlinearize_elina_texpr0_intlinear(&res,expr,discr);
  if (exc){
    elina_linexpr0_free(res);
    return NULL;
  }
  return res;
}

/* ====================================================================== */
/*  Linearization of tree expressions */
/* ====================================================================== */

/* transform in-place
	[A0,B0] + sum Xi [Ai,Bi]
   into
       ([A0,B0] + [A0,B0][-ulp,ulp]) + [-mf,mf] +
       sum Xi ([Ai,Bi] + [Ai,Bi][-ulp,ulp])

   i.e., add a relative error of magnitude ulp as an interval linear form
*/
static void elina_make_float_const(int frac_bits, int exp_bits, int exp_bias,
			     elina_float_const* cst, elina_scalar_discr_t discr)
{
  elina_scalar_t *b,*c;
  b = elina_scalar_alloc();
  c = elina_scalar_alloc();
  cst->ulp = elina_interval_alloc();
  cst->min = elina_interval_alloc();
  cst->min_normal = elina_interval_alloc();
  cst->max = elina_interval_alloc();
  cst->max_exact = elina_interval_alloc();

  elina_scalar_set_to_int(b,1,discr);
  elina_scalar_mul_2exp(b,b,-frac_bits,discr);
  elina_scalar_neg(cst->ulp->inf,b);
  elina_scalar_set(cst->ulp->sup,b);

  elina_scalar_set_to_int(b,1,discr);
  elina_scalar_mul_2exp(b,b,1-exp_bias-frac_bits,discr);
  elina_scalar_neg(cst->min->inf,b);
  elina_scalar_set(cst->min->sup,b);

  elina_scalar_set_to_int(b,1,discr);
  elina_scalar_mul_2exp(b,b,1-exp_bias,discr);
  elina_scalar_neg(cst->min_normal->inf,b);
  elina_scalar_set(cst->min_normal->sup,b);

  elina_scalar_set_to_int(b,2,discr);
  elina_scalar_set_to_int(c,1,discr);
  elina_scalar_mul_2exp(c,c,-frac_bits,discr);
  elina_scalar_sub(b,b,c,discr);
  elina_scalar_mul_2exp(b,b,(1<<exp_bits)-2-exp_bias,discr);
  elina_scalar_neg(cst->max->inf,b);
  elina_scalar_set(cst->max->sup,b);

  elina_scalar_set_to_int(b,1,discr);
  elina_scalar_mul_2exp(b,b,frac_bits,discr);
  elina_scalar_neg(cst->max_exact->inf,b);
  elina_scalar_set(cst->max_exact->sup,b);
  elina_scalar_free(b); 
  elina_scalar_free(c);
}

static void elina_float_const_clear(elina_float_const* cst)
{
  elina_interval_free(cst->ulp); 
  elina_interval_free(cst->min); 
  elina_interval_free(cst->min_normal); 
  elina_interval_free(cst->max); 
  elina_interval_free(cst->max_exact);
}

static void elina_linexpr0_round_float_lin(elina_linexpr0_t* l ,elina_float_const* f, elina_scalar_discr_t discr)
{
  size_t i;
  elina_dim_t dim;
  elina_coeff_t *coeff;
  elina_coeff_t *cst = &l->cst;
  elina_scalar_t *scalar;
  if(cst->discr==ELINA_COEFF_SCALAR){
	scalar = elina_scalar_alloc();
	elina_scalar_set(scalar,cst->val.scalar);
	elina_coeff_reinit(cst,ELINA_COEFF_INTERVAL,discr);
	elina_coeff_set_interval_scalar(cst,scalar,scalar);
	elina_scalar_free(scalar);
  }
  scalar = elina_scalar_alloc();
  elina_interval_t *interval = cst->val.interval;
  elina_interval_magnitude(scalar,interval);
  elina_scalar_mul(scalar,scalar,f->ulp->sup,discr);
  elina_scalar_add(scalar,scalar,f->min->sup,discr);
  elina_interval_enlarge_bound(interval,interval,scalar,discr);
  elina_linexpr0_ForeachLinterm(l,i,dim,coeff) {
    if(coeff->discr==ELINA_COEFF_SCALAR){
	elina_scalar_set(scalar,coeff->val.scalar);
	elina_coeff_reinit(coeff,ELINA_COEFF_INTERVAL,discr);
	elina_coeff_set_interval_scalar(coeff,scalar,scalar);
    }
    interval = coeff->val.interval;
    elina_interval_magnitude(scalar,interval);
    elina_scalar_mul(scalar,scalar,f->ulp->sup,discr);
    elina_interval_enlarge_bound(interval,interval,scalar,discr);
  }
  elina_scalar_free(scalar);
}


/* adds an absolute error to l corresponding to a conversion to int
   assumes that i overapproximates the values of l before conversion
 */
static void
elina_texpr0_to_int(elina_linexpr0_t* l /* in/out */, elina_interval_t *i /* in */,
		 elina_texpr_rdir_t d, elina_scalar_discr_t discr)
{
  elina_coeff_t *cst = &l->cst;
  elina_scalar_t *scalar;
  if(cst->discr==ELINA_COEFF_SCALAR){
	scalar = elina_scalar_alloc();
	elina_scalar_set(scalar,cst->val.scalar);
	elina_coeff_reinit(cst,ELINA_COEFF_INTERVAL,discr);
	elina_coeff_set_interval_scalar(cst,scalar,scalar);
	elina_scalar_free(scalar);
  }
  elina_interval_t *interval = cst->val.interval;
  elina_interval_t * half;
  switch (d) {
  case ELINA_RDIR_UP:
    /* add [0,1] */
    elina_scalar_add_uint(interval->sup,interval->sup,1,discr);
    break;
  case ELINA_RDIR_DOWN:
    /* add [-1,0] */
    elina_scalar_sub_uint(interval->inf,interval->inf,1,discr);
    break;
  case ELINA_RDIR_RND:
    /* add [-1,1] */
    elina_scalar_add_uint(interval->sup,interval->sup,1,discr);
    elina_scalar_sub_uint(interval->inf,interval->inf,1,discr);
    break;
  case ELINA_RDIR_ZERO:
    /* UP or DOWN or RND, depending on sign of i */
    if (elina_scalar_sgn(i->inf)<0) elina_scalar_add_uint(interval->sup,interval->sup,1,discr);
    if (elina_scalar_sgn(i->sup)>0) elina_scalar_sub_uint(interval->inf,interval->inf,1,discr);
    break;
  case ELINA_RDIR_NEAREST:
    /* [-0.5,0.5] */
    half = elina_interval_alloc();
    elina_interval_set_to_int(half,-1,1,discr);
    elina_interval_mul_2exp(half,half,-1,discr);
    elina_interval_add(interval,interval,half,discr);
    elina_interval_free(half);
    break;
  default:
    assert(0);
  }
}

/* adds rounding error to both l and i to go from type org to type dst */
elina_texpr_rtype_t elina_texpr0_round(elina_linexpr0_t* l /* in/out */, elina_interval_t *i /* in/out */,
		elina_texpr_rtype_t org, elina_texpr_rtype_t dst, elina_texpr_rdir_t d, elina_scalar_discr_t discr)
{
  
  elina_float_const* float_cst;
  if (dst==ELINA_RTYPE_REAL) return org;
  switch (dst) {
  case ELINA_RTYPE_INT:
    if (org==ELINA_RTYPE_INT) return org;
    elina_texpr0_to_int(l,i,d,discr);
    break;
  case ELINA_RTYPE_SINGLE:
    if (org==ELINA_RTYPE_SINGLE) return org;
    float_cst = (elina_float_const*)malloc(sizeof(elina_float_const));
    elina_make_float_const(23,8,127,float_cst,discr);
    elina_linexpr0_round_float_lin(l,float_cst,discr);
    elina_float_const_clear(float_cst);
    break;
  case ELINA_RTYPE_DOUBLE:
    if (org==ELINA_RTYPE_SINGLE || org==ELINA_RTYPE_DOUBLE) return org;
    float_cst = (elina_float_const*)malloc(sizeof(elina_float_const));
    elina_make_float_const(52,11,1023,float_cst,discr);
    elina_linexpr0_round_float_lin(l,float_cst,discr);
    elina_float_const_clear(float_cst);
    break;
  case ELINA_RTYPE_EXTENDED:
    if (org==ELINA_RTYPE_SINGLE && org==ELINA_RTYPE_DOUBLE && org==ELINA_RTYPE_EXTENDED)
      return org;
    float_cst = (elina_float_const*)malloc(sizeof(elina_float_const));
    elina_make_float_const(63,15,16383,float_cst,discr);
    
    elina_linexpr0_round_float_lin(l,float_cst,discr);
    elina_float_const_clear(float_cst);
    break;
  case ELINA_RTYPE_QUAD:
    if (org==ELINA_RTYPE_SINGLE || org==ELINA_RTYPE_DOUBLE ||
	org==ELINA_RTYPE_EXTENDED || org==ELINA_RTYPE_QUAD) return org;
    float_cst = (elina_float_const*)malloc(sizeof(elina_float_const));
    elina_make_float_const(112,15,16383,float_cst,discr);
    elina_linexpr0_round_float_lin(l,float_cst,discr);
    elina_float_const_clear(float_cst);
    break;
  default:
    assert(0);
  }
  
  elina_interval_round(i,i,dst,d,discr);
  return dst;
}


/* reduce l and i:
   - intersect i with the interval evaluation of l
   - if l is constant, replace it with i
   - check for emptiness
 */
void elina_texpr0_reduce(elina_interval_t** env,
		 elina_linexpr0_t* l /* in/out */, elina_interval_t* i /* in/out */, elina_scalar_discr_t discr)
{
  elina_interval_t *tmp = elina_interval_alloc();
  elina_interval_eval_elina_linexpr0(tmp,l,env,discr);
  elina_scalar_max(i->inf,i->inf,tmp->inf);
  elina_scalar_min(i->sup,i->sup,tmp->sup);
  elina_coeff_t *cst = &l->cst;
  bool flag = false;
  if(cst->discr==ELINA_COEFF_INTERVAL && elina_interval_is_bottom(cst->val.interval)){
	flag = true;
  }
  if (elina_interval_is_bottom(i) || flag) {
    elina_interval_set_bottom(i);
    if(cst->discr==ELINA_COEFF_SCALAR){
	elina_coeff_reinit(cst,ELINA_COEFF_INTERVAL,discr);
    }
    elina_interval_set_bottom(cst->val.interval);
    if (l->size>0) elina_linexpr0_reinit(l,0);
  }
  else if (l->size==0){
    if(!elina_scalar_infty(i->inf) && elina_scalar_equal(i->inf,i->sup)){
	elina_coeff_set_scalar(cst,i->inf);
    }
    else{
	elina_coeff_set_interval(cst,i);
    }
  }
  elina_interval_free(tmp);
}


/* multiplication heuristic: choose which interval to keep (0=a, 1=b) */
int elina_texpr0_cmp_range(elina_linexpr0_t* la, elina_interval_t *ia,
		    elina_linexpr0_t* lb, elina_interval_t *ib, elina_scalar_discr_t discr)
{
  int sgn_a,sgn_b;
  /* if one linear form is an interval keep it */
  if (la->size==0) return 0;
  if (lb->size==0) return 1;
  /* if only one interval has constant sign, keep it */
  sgn_a = elina_scalar_sgn(ia->inf) >= 0 || elina_scalar_sgn(ia->sup)<=0;
  sgn_b = elina_scalar_sgn(ib->inf)>=0 || elina_scalar_sgn(ib->sup)<=0;
  if (sgn_a!=sgn_b) return sgn_a ? 0 : 1;
  /* otherwise, keep the interval with the smallest relative range */
  elina_scalar_t *tmp1 = elina_scalar_alloc();
  elina_scalar_t *tmp2 = elina_scalar_alloc();
  elina_interval_range_rel(tmp1,ia,discr);
  elina_interval_range_rel(tmp2,ib,discr);
  int res = 0;
  if (elina_scalar_cmp(tmp1,tmp2)<0){
	 res = 0;
  }
  else{
	 res = 1;
  }
  elina_scalar_free(tmp1);
  elina_scalar_free(tmp2);
  return res;
}




/* ********************************************************************** */
/*  Interval Linearization of tree expressions */
/* ********************************************************************** */

static elina_texpr_rtype_t elina_texpr0_node_intlinearize(elina_texpr0_node_t* n,
			    elina_interval_t** env, size_t intdim,
			    elina_linexpr0_t** dlres /* out */, elina_interval_t *ires /* out */,
			    elina_scalar_discr_t discr)
{
  elina_linexpr0_t * lres = *dlres;
  elina_interval_t *i1,*i2;
  elina_linexpr0_t* l1;
  elina_texpr_rtype_t t1,t2;
  switch (n->op) {
  case ELINA_TEXPR_NEG:
    /* negate linear form & interval, no rounding */
    t1 = elina_interval_intlinearize_texpr0_rec(n->exprA,env,intdim,dlres,ires,discr);
    elina_linexpr0_neg(lres);
    elina_interval_neg(ires,ires);
    return t1;

  case ELINA_TEXPR_CAST:
    /* round linear form & interval */
    t1 = elina_interval_intlinearize_texpr0_rec(n->exprA,env,intdim,dlres,ires,discr);
    elina_texpr0_round(lres,ires,t1,n->type,n->dir,discr);
    elina_texpr0_reduce(env,lres,ires,discr);
    break;

  case ELINA_TEXPR_SQRT:
    /* intlinearize argument, lres is not used */
	
    elina_interval_intlinearize_texpr0_rec(n->exprA,env,intdim,dlres,ires,discr);
    /* interval square root */
    elina_interval_sqrt(ires,ires,discr);
    elina_interval_round(ires,ires,n->type,n->dir,discr);
    elina_linexpr0_reinit(lres,0);
    if(!elina_scalar_infty(ires->inf) && elina_scalar_equal(ires->inf,ires->sup)){
	elina_coeff_set_scalar(&lres->cst,ires->inf);
    }
    else{
	 elina_coeff_set_interval(&lres->cst,ires);
    }
    break;

  case ELINA_TEXPR_ADD:
  case ELINA_TEXPR_SUB:
    i1 = elina_interval_alloc();
    l1 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
    //elina_linexpr0_init(&l1,0);
    /* intlinearize arguments */
    t1 = elina_interval_intlinearize_texpr0_rec(n->exprA,env,intdim,&l1,i1,discr);
    t2 = elina_interval_intlinearize_texpr0_rec(n->exprB,env,intdim,dlres,ires,discr);
	
    if (elina_interval_is_bottom(i1) || elina_interval_is_bottom(ires)){
      elina_interval_set_bottom(ires);
      elina_linexpr0_reinit(lres,0);
      elina_coeff_set_interval(&lres->cst,ires);
    }
    else {
      /* add/sub linear form & interval */
      if (n->op==ELINA_TEXPR_ADD) {
	
	elina_linexpr0_add(dlres,&l1,dlres,discr);
	elina_interval_add(ires,i1,ires,discr);
	
      }
      else {
	elina_linexpr0_sub(dlres,&l1,dlres,discr);
	elina_interval_sub(ires,i1,ires,discr);
      }
	
      /* round */
      elina_texpr0_round(lres,ires,
		      (t1==ELINA_RTYPE_INT && t2==ELINA_RTYPE_INT) ?
		      ELINA_RTYPE_INT : ELINA_RTYPE_REAL,
		      n->type,n->dir,discr);
	
      /* reduce */
      elina_texpr0_reduce(env,lres,ires,discr);
	
    }
    elina_interval_free(i1);
    elina_linexpr0_free(l1);
    
    break;

  case ELINA_TEXPR_DIV:
	
    i1 = elina_interval_alloc();
    //elina_linexpr0_init(&l1,0);
    l1 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
    /* intlinearize arguments, l1 is not used */
    elina_interval_intlinearize_texpr0_rec(n->exprA,env,intdim,dlres,ires,discr);
    elina_interval_intlinearize_texpr0_rec(n->exprB,env,intdim,&l1,i1,discr);
    if (elina_interval_is_bottom(i1) || elina_interval_is_bottom(ires)){
      elina_interval_set_bottom(ires);
      elina_linexpr0_reinit(lres,0);
      elina_coeff_set_interval(&lres->cst,ires);
    }
    else {
      /* divide linear form & interval */
	
	
      elina_linexpr0_div(lres,i1,discr);
      elina_interval_div(ires,ires,i1,discr);
      /* round */
      elina_texpr0_round(lres,ires,ELINA_RTYPE_REAL,n->type,n->dir,discr);
      /* reduce */
	
      elina_texpr0_reduce(env,lres,ires,discr);
	
    }
    elina_interval_free(i1);
    elina_linexpr0_free(l1);
    break;

  case ELINA_TEXPR_MUL:
    i1 = elina_interval_alloc();
    //elina_linexpr0_init(&l1,0);
    l1 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
    /* intlinearize arguments */
    t1 = elina_interval_intlinearize_texpr0_rec(n->exprA,env,intdim,&l1,i1,discr);
    t2 = elina_interval_intlinearize_texpr0_rec(n->exprB,env,intdim,dlres,ires,discr);
    if (elina_interval_is_bottom(i1) || elina_interval_is_bottom(ires)){
      elina_interval_set_bottom(ires);
      elina_linexpr0_reinit(lres,0);
      elina_coeff_set_interval(&lres->cst,ires);
      elina_linexpr0_free(l1);
    }
    else {
      /* multiply one linear form with the other interval */
	
      if (elina_texpr0_cmp_range(l1,i1,lres,ires,discr))  {
	/* res = ires * l1 */
	elina_linexpr0_clear(lres);
	elina_linexpr0_scale(l1,ires,discr);
	**dlres = *l1;
	lres = *dlres;
      }
      else {
	/* res = i1 * lres */
	elina_linexpr0_free(l1);
	elina_linexpr0_scale(lres,i1,discr);
      }
	
      elina_interval_mul(ires,i1,ires,discr);
      /* round */
	
      elina_texpr0_round(lres,ires,
		      (t1==ELINA_RTYPE_INT && t2==ELINA_RTYPE_INT) ?
		      ELINA_RTYPE_INT : ELINA_RTYPE_REAL,
		      n->type,n->dir,discr);
      /* reduce */
      elina_texpr0_reduce(env,lres,ires,discr);
    }
    elina_interval_free(i1);
	
    break;

  case ELINA_TEXPR_MOD:
   i1 = elina_interval_alloc();
   //elina_linexpr0_init(&l1,0);
   l1 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
    /* intlinearize arguments, lres & l1 are not used */
    elina_interval_intlinearize_texpr0_rec(n->exprA,env,intdim,dlres,ires,discr);
    elina_interval_intlinearize_texpr0_rec(n->exprB,env,intdim,&l1,i1,discr);
    if (elina_interval_is_bottom(i1) || elina_interval_is_bottom(ires)){
      elina_interval_set_bottom(ires);
      elina_linexpr0_reinit(lres,0);
      elina_coeff_set_interval(&lres->cst,ires);
    }
    else {
      /* interval modulo, no rounding */
      elina_interval_mod(ires,ires,i1,n->type==ELINA_RTYPE_INT,discr);
      elina_linexpr0_reinit(lres,0);
      if(!elina_scalar_infty(ires->inf) && elina_scalar_equal(ires->inf,ires->sup)){
	 elina_coeff_set_scalar(&lres->cst,ires->inf);
      }
      else{
	 elina_coeff_set_interval(&lres->cst,ires);
      }
    }
    elina_interval_free(i1);
    elina_linexpr0_free(l1);
    break;

  case ELINA_TEXPR_POW:
    i1 = elina_interval_alloc();
    //elina_linexpr0_init(&l1,0);
    l1 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
    /* intlinearize arguments, lres & l1 are not used */
    elina_interval_intlinearize_texpr0_rec(n->exprA,env,intdim,dlres,ires,discr);
    elina_interval_intlinearize_texpr0_rec(n->exprB,env,intdim,&l1,i1,discr);
    if (elina_interval_is_bottom(i1) || elina_interval_is_bottom(ires)){
      elina_interval_set_bottom(ires);
      elina_linexpr0_reinit(lres,0);
      elina_coeff_set_interval(&lres->cst,ires);
    }
    else {
      elina_interval_pow(ires,ires,i1,discr);
      elina_linexpr0_reinit(lres,0);
      if(!elina_scalar_infty(ires->inf) && elina_scalar_equal(ires->inf,ires->sup)){
	 elina_coeff_set_scalar(&lres->cst,ires->inf);
      }
      else{
	 elina_coeff_set_interval(&lres->cst,ires);
      }
    }
    elina_interval_free(i1);
    elina_linexpr0_free(l1);
    break;

  default:
    assert(0);
  }
  return n->type;
}

elina_texpr_rtype_t elina_interval_intlinearize_texpr0_rec(elina_texpr0_t* expr,
			    elina_interval_t** env, size_t intdim,
			    elina_linexpr0_t** dlres /* out */, elina_interval_t *ires /* out */
			    ,elina_scalar_discr_t discr)
{
  elina_linexpr0_t *lres = *dlres;
  elina_texpr_rtype_t t;
  assert(expr);
  elina_coeff_t * cst, *coeff;
  elina_scalar_t *scalar;
  elina_linterm_t * term;
  switch(expr->discr){
  case ELINA_TEXPR_CST:
    cst = &expr->val.cst;
    if(cst->discr==ELINA_COEFF_SCALAR){
	scalar = cst->val.scalar;
	elina_scalar_set(ires->inf,scalar);
	elina_scalar_set(ires->sup,scalar);
    }
    else{
	elina_interval_set(ires,cst->val.interval);
    }
    elina_linexpr0_reinit(lres,0);
    elina_coeff_set(&lres->cst,&expr->val.cst);
    t = elina_interval_is_int(ires,discr) ? ELINA_RTYPE_INT : ELINA_RTYPE_REAL;
    break;
  case ELINA_TEXPR_DIM:
    elina_interval_set(ires,env[expr->val.dim]);
    elina_linexpr0_reinit(lres,1);
    cst = &lres->cst;
    elina_coeff_reinit(cst,ELINA_COEFF_SCALAR,discr);
    elina_scalar_set_to_int(cst->val.scalar,0,discr);
    term = &lres->p.linterm[0];
    term->dim = expr->val.dim;
    coeff = &term->coeff;
    elina_coeff_reinit(coeff,ELINA_COEFF_SCALAR,discr);
    elina_scalar_set_to_int(coeff->val.scalar,1,discr);
    t = (expr->val.dim<intdim) ? ELINA_RTYPE_INT : ELINA_RTYPE_REAL;
    break;
  case ELINA_TEXPR_NODE:
	
    t = elina_texpr0_node_intlinearize(expr->val.node,env,intdim,dlres,ires,discr);
    break;
  default:
    t = 0;
    assert(false);
  }
  return t;
}


bool elina_interval_intlinearize_elina_texpr0(elina_linexpr0_t** dres, elina_texpr0_t* expr, elina_interval_t** env, size_t intdim, elina_scalar_discr_t discr)
{
  elina_linexpr0_t * res = *dres;
  bool exc;
  elina_interval_t *i = elina_interval_alloc();
  elina_interval_intlinearize_texpr0_rec(expr,env,intdim,dres,i,discr);
  elina_coeff_t *cst = &res->cst;
  
  if(cst->discr==ELINA_COEFF_SCALAR){
	elina_scalar_t *scalar = elina_scalar_alloc();
	elina_scalar_set(scalar,cst->val.scalar);
	elina_coeff_reinit(cst,ELINA_COEFF_INTERVAL,discr);
	elina_coeff_set_interval_scalar(cst,scalar,scalar);
	elina_scalar_free(scalar);
  }
  if (!elina_interval_is_bottom(i) && !elina_interval_is_bottom(cst->val.interval)) {
    if (res->size==0){
      elina_scalar_t *inf, *sup;
      inf = cst->val.interval->inf;
      sup = cst->val.interval->sup;	
      elina_scalar_max(inf,inf,i->inf);
      elina_scalar_min(sup,sup,i->sup);
      if(!elina_scalar_infty(inf) && elina_scalar_equal(inf,sup)){
	elina_coeff_reduce(cst);
      }
    }
    exc = false;
  }
  else {
    exc = true;
  }
  elina_interval_free(i);
  return exc;
}


elina_linexpr0_t* elina_intlinearize_texpr0(elina_manager_t* man,
				      elina_abstract0_t* abs,
				      elina_texpr0_t* expr,
				      bool* pexact,
				      elina_scalar_discr_t discr,
				      bool quasilinearize)
{
  elina_interval_t** env = NULL;
  elina_dimension_t dim = {0,0};
  elina_linexpr0_t* res = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
  size_t i, size;
  if (pexact) *pexact = false;
  if (elina_texpr0_is_interval_linear(expr)){
    elina_intlinearize_elina_texpr0_intlinear(&res, expr, discr);
	
  }
  else {
    //elina_intlinearize_alloc(man,abs,env,&dim,discr);
    assert(!elina_abstract0_is_bottom(man,abs));
    env = elina_abstract0_to_box(man,abs);
  
    dim = elina_abstract0_dimension(man,abs);
    size = dim.intdim+dim.realdim;
    for(i=0; i < size; i++){
	 elina_interval_convert(env[i],discr);
    }
	
    elina_interval_intlinearize_elina_texpr0(&res, expr, env,  dim.intdim,discr);
	
  }
  if (quasilinearize && !elina_linexpr0_is_quasilinear(res)){
   if (!env){
	assert(!elina_abstract0_is_bottom(man,abs));
        env = elina_abstract0_to_box(man,abs);
  
        dim = elina_abstract0_dimension(man,abs);
        size = dim.intdim+dim.realdim;
        for(i=0; i < size; i++){
	    elina_interval_convert(env[i],discr);
        } 
	//elina_intlinearize_alloc(man,abs,env,&dim,discr);
    }
    quasilinearize_elina_linexpr0(res,env,false,discr);
  }
  if(env)elina_interval_array_free(env,dim.intdim+dim.realdim);
  return res;
}

elina_linexpr0_t** elina_intlinearize_texpr0_array(elina_manager_t* man,
					     elina_abstract0_t* abs,
					     elina_texpr0_t** texpr0, size_t size,
					     bool* pexact,
					     elina_scalar_discr_t discr,
					     bool quasilinearize)
{
  elina_dimension_t dim = {0,0};
  elina_interval_t** env = NULL;
  elina_linexpr0_t** res;
  size_t i, abs_size;

  if (pexact) *pexact = false;
  
  res = malloc(size*sizeof(elina_linexpr0_t*));
  for (i=0; i<size; i++){
    res[i] = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
    if (elina_texpr0_is_interval_linear(texpr0[i])){
      elina_intlinearize_elina_texpr0_intlinear(&res[i], texpr0[i],discr);
    }
    else {
      if (!env){
	    assert(!elina_abstract0_is_bottom(man,abs));
	    env = elina_abstract0_to_box(man,abs);
	  
	    dim = elina_abstract0_dimension(man,abs);
	    abs_size = dim.intdim+dim.realdim;
	    for(i=0; i < abs_size; i++){
		 elina_interval_convert(env[i],discr);
	    }
	 //elina_intlinearize_alloc(man,abs,env,&dim,discr);
      }
      elina_interval_intlinearize_elina_texpr0(&res[i], texpr0[i], env, dim.intdim,discr);
    }
    if (quasilinearize && !elina_linexpr0_is_quasilinear(res[i])){
      if (!env){ 
	    assert(!elina_abstract0_is_bottom(man,abs));
	    env = elina_abstract0_to_box(man,abs);
	  
	    dim = elina_abstract0_dimension(man,abs);
	    abs_size = dim.intdim+dim.realdim;
	    for(i=0; i < abs_size; i++){
		 elina_interval_convert(env[i],discr);
	    }
		//elina_intlinearize_alloc(man,abs,env,&dim,discr);
      }
      quasilinearize_elina_linexpr0(res[i],env,false,discr);
    }
  }
  if(env)elina_interval_array_free(env,dim.intdim+dim.realdim);
  return res;
}

/********************************************/
/**** Linearize constraints 		*****/
/*********************************************/

bool elina_intlinearize_elina_tcons0_intlinear(elina_lincons0_t* res, elina_tcons0_t* cons, elina_scalar_discr_t discr)
{
  bool exc = elina_intlinearize_elina_texpr0_intlinear(&res->linexpr0,cons->texpr0,discr);
  if (exc){
    elina_lincons0_set_bool(res,false, discr);
  }
  else {
    res->constyp = cons->constyp;
    if (cons->scalar){
      elina_scalar_set(res->scalar,cons->scalar);
    }
    else {
      elina_scalar_set_to_int(res->scalar,0,discr);
    }
  }
  return exc;
}

bool elina_intlinearize_elina_tcons0(elina_lincons0_t* res,
				   elina_tcons0_t* cons,
				   elina_interval_t** env, size_t intdim, elina_scalar_discr_t discr)
{
  bool exc;
  elina_interval_t *i, *bound;
  i = elina_interval_alloc();
  bound = elina_interval_alloc();
  elina_interval_reinit(bound,discr);
  elina_interval_reinit(i,discr);
  elina_interval_intlinearize_texpr0_rec(cons->texpr0,env,intdim,&res->linexpr0,i,discr);
  /* checks that the contraint is satisfiable */
  switch (cons->constyp){
  case ELINA_CONS_EQ:
    elina_scalar_max(i->inf,i->inf,bound->inf);
    elina_scalar_min(i->sup,i->sup,bound->sup);
    break;
  case ELINA_CONS_SUPEQ: 
  case ELINA_CONS_SUP:
    elina_interval_set_top(bound);
    elina_scalar_set_to_int(bound->inf,0,discr);
    elina_scalar_set_infty(bound->sup,1);
    elina_scalar_max(i->inf,i->inf,bound->inf);
    elina_scalar_min(i->sup,i->sup,bound->sup);
    break;
  default:
    break;
  }
  elina_coeff_t * cst = &res->linexpr0->cst;
  bool was_scalar = false;
  if(cst->discr==ELINA_COEFF_SCALAR){
	elina_scalar_t *tmp = elina_scalar_alloc();
	elina_scalar_set(tmp,cst->val.scalar);
	elina_coeff_reinit(cst,ELINA_COEFF_INTERVAL,discr);
	elina_coeff_set_interval_scalar(cst,tmp,tmp);
	was_scalar = true;
	elina_scalar_free(tmp);
  } 
  elina_interval_t *interval = cst->val.interval;
  if (!elina_interval_is_bottom(i) && !elina_interval_is_bottom(interval)) {
    if (res->linexpr0->size==0){
	elina_scalar_max(interval->inf,interval->inf,i->inf);
	elina_scalar_min(interval->sup,interval->sup,i->sup);
	if(!elina_scalar_infty(interval->inf) && elina_scalar_equal(interval->inf,interval->sup)){
		elina_coeff_reduce(cst);
	}
    }
    if(was_scalar){
	elina_coeff_reduce(cst);
    }
    res->constyp = cons->constyp;
    if (cons->scalar){
      elina_scalar_set(res->scalar,cons->scalar);
    }
    else {
      elina_scalar_set_to_int(res->scalar,0,discr);
    }
    exc = false;
  }
  else {
	
    elina_lincons0_set_bool(res,false, discr);
    exc = true;
  }
  elina_interval_free(i);
  elina_interval_free(bound);
   
  return exc;
}


elina_lincons0_t elina_intlinearize_tcons0(elina_manager_t* man,
				     elina_abstract0_t* abs,
				     elina_tcons0_t* cons,
				     bool* pexact,
				     elina_scalar_discr_t discr,
				     bool quasilinearize, bool meet)
{
  elina_dimension_t dim = {0,0};
  elina_interval_t** env = NULL;
  elina_lincons0_t res;
  size_t i, size;
  elina_linexpr0_t *linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
  elina_scalar_t *scalar = elina_scalar_alloc();
  res.constyp = cons->constyp;
  res.linexpr0 = linexpr0;
  res.scalar = scalar;
  if (pexact) *pexact = false;
  if (elina_tcons0_is_interval_linear(cons)){
    elina_intlinearize_elina_tcons0_intlinear(&res, cons,discr);
  }
  else {
	assert(!elina_abstract0_is_bottom(man,abs));
  	env = elina_abstract0_to_box(man,abs);
  
  	dim = elina_abstract0_dimension(man,abs);
  	size = dim.intdim+dim.realdim;
  	for(i=0; i < size; i++){
		elina_interval_convert(env[i],discr);
  	}
	
    //elina_intlinearize_alloc(man,abs,env,&dim,discr);
    elina_intlinearize_elina_tcons0(&res, cons, env,
			       dim.intdim,discr);
  }
  if (quasilinearize && !elina_linexpr0_is_quasilinear(res.linexpr0)){
  if (!env){
	assert(!elina_abstract0_is_bottom(man,abs));
  	env = elina_abstract0_to_box(man,abs);
  
  	dim = elina_abstract0_dimension(man,abs);
  	size = dim.intdim+dim.realdim;
  	for(i=0; i < size; i++){
		elina_interval_convert(env[i],discr);
  	}
	 //elina_intlinearize_alloc(man,abs,env,&dim,discr);
    }
    quasilinearize_elina_lincons0(&res,env,meet,discr);
  }
  if(env)elina_interval_array_free(env,dim.intdim+dim.realdim);
  return res;
}


/* ********************************************************************** */
/*  Boxization of interval linear expressions */
/* ********************************************************************** */

static bool elina_boxize_lincons0(elina_interval_t** res,
			       bool* tchange,
			       elina_lincons0_t* cons,
			       elina_interval_t** env,
			       size_t intdim,
			       bool intervalonly, elina_scalar_discr_t discr)
{
  size_t i;
  elina_linexpr0_t* expr;
  bool change,globalchange;
  bool exc;
  assert(cons->constyp == ELINA_CONS_EQ ||
	 cons->constyp == ELINA_CONS_SUPEQ ||
	 cons->constyp == ELINA_CONS_SUP);

  expr = cons->linexpr0;
  globalchange = false;

  /* Iterates on coefficients */
  elina_interval_t * interval = elina_interval_alloc();
  elina_scalar_t * scalar2 = elina_scalar_alloc();
  for (i=0; i<expr->size; i++){
    elina_linterm_t *term = &expr->p.linterm[i];
    elina_dim_t dim = term->dim;
    elina_coeff_t * coeff = &term->coeff;
    elina_coeff_t * tmp = elina_coeff_alloc(coeff->discr);
    /* 1. We decompose the expression e = ax+e' */
    elina_coeff_swap(tmp,coeff);
    elina_scalar_t *inf, *sup;
    if(tmp->discr==ELINA_COEFF_SCALAR){
	inf = tmp->val.scalar;
	sup = inf;
	//elina_scalar_neg(sup,sup);
    }
    else{
	inf = tmp->val.interval->inf;
	//elina_scalar_neg(inf,inf);
	sup = tmp->val.interval->sup;
	//elina_scalar_neg(sup,sup);
    }
    /* 2. evaluate e' */
    elina_interval_eval_elina_linexpr0(interval,expr,env,discr);
    /* 3. Perform deduction from ax+e' = [-m,M]x + [-e,E] >= 0
	  we can deduce [-m,M]x + E >= 0
	  (for equality, [-m,M]x - e <= 0)
    */
    bool equality = tmp->discr==ELINA_COEFF_SCALAR;
    change = false;
    if (!elina_interval_is_top(interval)){
      if (equality && !intervalonly){
	/* [-m,M]=[a,a] */
	
	int sgn = elina_scalar_sgn(sup);
	if (sgn!=0){
	  /* From ax+E >= 0, we deduce
	     (1) If a>0, then x>=-E/a
	     (2) If a<0, then x<=-E/a
	     From ax-e <= 0, we deduce
	     (3) If a>0, then x<=e/a
	     (4) If a<0, then x>=e/a
	  */
	  if (sgn>0 || cons->constyp == ELINA_CONS_EQ){
	    /* 1,4: inf bound
	       If we have a>0 (1), we compute E/a
	       If we have a<0 (4), we compute -e/a
	    */
		
	    if (sgn>0){
		elina_scalar_neg(sup,sup);
	      elina_scalar_div(scalar2,
			interval->sup,
			sup,discr);
		elina_scalar_neg(sup,sup);
	    }
	    else {
	      elina_scalar_neg(inf,inf);
	      elina_scalar_div(scalar2,
			interval->inf,
			inf,discr);
		elina_scalar_neg(inf,inf);
	    }
	    if (dim<intdim && !elina_scalar_infty(scalar2)){
	      if (cons->constyp==ELINA_CONS_SUP &&
		  elina_scalar_is_integer(scalar2)){
		elina_scalar_add_uint(scalar2,scalar2,1,discr);
	      }
	      else {
		elina_scalar_ceil(scalar2,
			    scalar2,discr);
	      }
	    }
	    /* We update the interval */
	    if (elina_scalar_cmp(scalar2, res[dim]->inf)>0){
	      change = true;
	      if (tchange) tchange[2*dim] = true;
	      elina_scalar_set(res[dim]->inf, scalar2);
	    }
	  }
	  if (sgn<0 || cons->constyp == ELINA_CONS_EQ){
	    /* 2,3: sup bound
	       If we have a<0 (2), we compute -E/a
	       If we have a>0 (3), we compute e/a
	    */
	    if (sgn<0){
	      elina_scalar_neg(inf,inf);
	      elina_scalar_div(scalar2,
			interval->sup,
			inf,discr);
		elina_scalar_neg(inf,inf);
	    }
	    else {
		elina_scalar_neg(sup,sup);
	      elina_scalar_div(scalar2,
			interval->inf,
			sup,discr);
		elina_scalar_neg(sup,sup);
	    }
	    if (dim<intdim && !elina_scalar_infty(scalar2)){
	      if (cons->constyp==ELINA_CONS_SUP &&
		  elina_scalar_is_integer(scalar2)){
		elina_scalar_sub_uint(scalar2,scalar2,1,discr);
	      }
	      else {
		elina_scalar_floor(scalar2,
			    scalar2,discr);
	      }
	    }
	    /* We update the interval */
	    if (elina_scalar_cmp(scalar2, res[dim]->sup)<0){
	      change = true;
	      if (tchange) tchange[2*dim+1] = true;
	      elina_scalar_set(res[dim]->sup, scalar2);
	    }
	  }
	}
      }
      else if (!equality){
	/* We have a real interval [-m,M] */
	/* Case [-m,M]x+E >= 0
	  (1) If -m>0, we rewrite [-m,M]x >= -E, and we have -1/m > 1/M > 0
	  (1A) If E<=0 <=> -E>=0 then x >= -E/M
	  (1B) If E>0  <=> -E<0  then x >= E/m
	  (2) If M<0, we rewrite [-M,m]x <= E, and we have 0 < 1/m < -1/M
	  (2A) If E<=0           then x <= E/m
	  (2B) If E>0            then x <= -E/M
	  Case [-m,M]x-e <= 0
	  (3) If -m>0, we rewrite [-m,M]x <= e, and we have 0 < 1/M < -1/m
	  (3A) If e<=0           then x <= e/M
	  (3B) If e>0            then x <= -e/m
	  (4) If M<0, we rewrite [-M,m]x >= -e, and we have -1/M > 1/m > 0
	  (4A) If e<=0 <=> -e>=0 then x >= -e/m
	  (4B) If e>0  <=> -e<0  then x >= e/M
	*/
	int sgnitv =
	  elina_scalar_sgn(inf)>0 ?
	  1 :
	  ( elina_scalar_sgn(sup)<0 ?
	    -1 :
	    0 );
	if (sgnitv != 0){
	  int sgne = elina_scalar_sgn(interval->inf);
	  int sgnE = elina_scalar_sgn(interval->sup);
	  if (sgnitv>0 || (cons->constyp==ELINA_CONS_EQ && sgnitv<0)){
	    /* 1,4: inf bound */
	    if (sgnitv>0){ /* 1 */
	      if (sgnE<=0){ /* 1A */
		/* We compute E/M */
		elina_scalar_neg(sup,sup);
		elina_scalar_div(scalar2,
			  interval->sup,
			  sup,discr);
		elina_scalar_neg(sup,sup);
		
	      } else { /* 1B */
		/* We compute -E/m */
		elina_scalar_neg(inf,inf);
		elina_scalar_div(scalar2,
			  interval->sup,
			  inf,discr);
		elina_scalar_neg(inf,inf);
		//bound_neg(scalar2,
		//	  scalar2);
	      }
	    }
	    else { /* 4 */
	      if (sgne>=0){ /* 4A */
		/* We compute e/m */
		elina_scalar_neg(inf,inf);
		elina_scalar_div(scalar2,
			  interval->inf,
			  inf,discr);
		elina_scalar_neg(inf,inf);
	      } else { /* 4B */
		/* We compute -e/M */
		elina_scalar_neg(sup,sup);
		elina_scalar_div(scalar2,
			  interval->inf,
			  sup,discr);
		elina_scalar_neg(sup,sup);
		//bound_neg(scalar2,
		//	  scalar2);
	      }
	    }
	    if (dim<intdim && !elina_scalar_infty(scalar2)){
	      if (cons->constyp==ELINA_CONS_SUP &&
		  elina_scalar_is_integer(scalar2)){
		elina_scalar_add_uint(scalar2,scalar2,1,discr);
	      }
	      else {
		elina_scalar_ceil(scalar2,scalar2,discr);
	      }
	    }
	    /* We update the interval */
	    if (elina_scalar_cmp(scalar2, res[dim]->inf)>0){
	      change = true;
	      if (tchange) tchange[2*dim] = true;
	      elina_scalar_set(res[dim]->inf, scalar2);
	    }
	  }
	  if (sgnitv<0 || (cons->constyp==ELINA_CONS_EQ && sgnitv>0)){
	    /* 2,3: sup bound */
	    if (sgnitv<0){ /* 2 */
	      if (sgnE>=0){ /* 2B */
		/* We compute -E/M */
		elina_scalar_neg(sup,sup);
		elina_scalar_div(scalar2,
			  interval->sup,
			  sup,discr);
		elina_scalar_neg(sup,sup);
		//bound_neg(scalar2,
		//	  scalar2);
	      } else { /* 2A */
		/* We compute E/m */
		elina_scalar_neg(inf,inf);
		elina_scalar_div(scalar2,
			  interval->sup,
			  inf,discr);
		elina_scalar_neg(inf,inf);
	      }
	    }
	    else { /* 3 */
	      if (sgne<=0){ /* 3B */
		/* We compute -e/m */
		elina_scalar_neg(inf,inf);
		elina_scalar_div(scalar2,
			  interval->inf,
			  inf,discr);
		elina_scalar_neg(inf,inf);
		//bound_neg(scalar2,
		//	  scalar2);
	      }
	      else { /* 3A */
		/* We compute e/M */
		elina_scalar_neg(sup,sup);
		elina_scalar_div(scalar2,
			  interval->inf,
			  sup,discr);
		elina_scalar_neg(sup,sup);
	      }
	    }
	    if (dim<intdim && !elina_scalar_infty(scalar2)){
	      if (cons->constyp==ELINA_CONS_SUP &&
		  elina_scalar_is_integer(scalar2)){
		elina_scalar_sub_uint(scalar2,scalar2,1,discr);
	      }
	      else {
		elina_scalar_floor(scalar2,
			    scalar2,discr);
	      }
	    }
	    /* We update the interval */
	    if (elina_scalar_cmp(scalar2, res[dim]->sup)<0){
	      change = true;
	      if (tchange) tchange[2*dim+1] = true;
	      elina_scalar_set(res[dim]->sup, scalar2);
	    }
	  }
	}
      }
    }
    elina_coeff_swap(tmp,coeff);
    elina_coeff_free(tmp);
    if (change){
      globalchange = true;
      exc = elina_interval_canonicalize(res[dim],dim<intdim,discr);
      if (exc){
	elina_interval_set_bottom(res[0]);
	return true;
      }
    }
  }
  elina_interval_free(interval);
  elina_scalar_free(scalar2);
  if (expr->size==0 &&
      eval_elina_cstlincons0(cons)==0){
    elina_interval_set_bottom(res[0]);
    globalchange = true;
  }
  return globalchange;
}

/* This function deduces interval constraints from a set of interval linear
   constraints.

   Return true if some (better if res==env) bounds have been inferred.

   - The inferred bounds are stored in res (which may be equal to env)
   - If tchange!=NULL *and initialized to false*,
     tchange[2dim] (resp. 2dim+1) set to true indicates
     that the inf (resp. sup) bound of dimension dim has been improved.
   - env is the current bounds for variables
   - kmax specifies the maximum number of iterations, when res==env
   - if intervalonly is true, deduces bounds from a constraint only when the
     coefficient associated to the current dimension is an interval.
*/
bool elina_boxize_lincons0_array(elina_interval_t** res,bool* tchange,
				      elina_lincons0_array_t* array,
				      elina_interval_t** env, size_t intdim,
				      size_t kmax,
				      bool intervalonly, elina_scalar_discr_t discr)
{
  size_t i,k;
  bool change,globalchange;

  if (kmax<1) kmax=1;
  if (res!=env) kmax=1;

  globalchange = false;
  /* we possibly perform kmax passes */
  for (k=0;k<(size_t)kmax;k++){
    change = false;
    for (i=0; i<array->size; i++){
      if (array->p[i].constyp==ELINA_CONS_EQ ||
	  array->p[i].constyp==ELINA_CONS_SUPEQ ||
	  array->p[i].constyp==ELINA_CONS_SUP){
	
	change =
	  elina_boxize_lincons0(res,tchange,&array->p[i],env,intdim,intervalonly,discr)
	  ||
	  change
	  ;
	globalchange = globalchange || change;
	if (elina_interval_is_bottom(res[0])){
	  return true;
	}
      }
    }
    if (!change) break;
  }
  return globalchange;
}


/********************************************/
/* Linearize array of constraints 	    */
/*********************************************/




bool elina_intlinearize_elina_tcons0_array(elina_lincons0_array_t* res,
					 elina_tcons0_array_t* array,
					 elina_interval_t** env, size_t intdim,  elina_scalar_discr_t discr)
{
  bool exc;
  elina_interval_t *itv,*bound;
  size_t i,index;
 
  itv = elina_interval_alloc();
  elina_interval_reinit(itv,discr);
  bound  = elina_interval_alloc();
  elina_interval_reinit(bound,discr);
  elina_lincons0_array_reinit(res,array->size);
  exc = false;
  for (i=0; i<array->size;i++){
    elina_interval_intlinearize_texpr0_rec(array->p[i].texpr0,env,intdim,&res->p[i].linexpr0,itv,discr);
    res->p[i].constyp = array->p[i].constyp;
	
    if (array->p[i].scalar){
      elina_scalar_set(res->p[i].scalar,array->p[i].scalar);
    }
    else {
      elina_scalar_set_to_int(res->p[i].scalar,0,discr);
    }

    /* checks that the contraint is satisfiable */
    switch (array->p[i].constyp){
    case ELINA_CONS_EQ:
      elina_scalar_max(itv->inf,itv->inf,bound->inf);
      elina_scalar_min(itv->sup,itv->sup,bound->sup);
      break;
    case ELINA_CONS_SUPEQ:
    case ELINA_CONS_SUP:
      elina_interval_set_top(bound);
      elina_scalar_set_to_int(bound->inf,0,discr);
      elina_scalar_set_infty(bound->sup,1);
      elina_scalar_max(itv->inf,itv->inf,bound->inf);
      elina_scalar_min(itv->sup,itv->sup,bound->sup);
      break;
    default:
      break;
    }
     elina_coeff_t *cst = &res->p[i].linexpr0->cst;
    if (elina_interval_is_bottom(itv) ||
	(cst->discr==ELINA_COEFF_INTERVAL && elina_interval_is_bottom(cst->val.interval)) ||
	(res->p[i].linexpr0->size==0 ?
	 eval_elina_cstlincons0(&res->p[i])==0 :
	 false)){
      exc = true;
      elina_lincons0_array_reinit(res,1);
      elina_lincons0_set_bool(&res->p[0],false, discr);
      break;
    }
  }
  elina_interval_free(itv);
  elina_interval_free(bound);
  return exc;
}


elina_lincons0_array_t elina_intlinearize_tcons0_array(elina_manager_t* man,
						 elina_abstract0_t* abs,
						 elina_tcons0_array_t* array,
						 bool* pexact,
						 elina_scalar_discr_t discr,
						 elina_linexpr_type_t linearize, bool meet,
						 bool boxize, size_t kmax, bool intervalonly)
{
	
  elina_dimension_t dim = {0,0};
  elina_interval_t** env = NULL;
  elina_lincons0_array_t res;
  size_t i,size;
  bool change = false;
  bool* tchange = NULL;

  if (pexact) *pexact = false;
  assert(!elina_abstract0_is_bottom(man,abs));
  
  env = elina_abstract0_to_box(man,abs);

  dim = elina_abstract0_dimension(man,abs);
  size = dim.intdim+dim.realdim;
  for(i=0; i < size; i++){
	elina_interval_convert(env[i],discr);
  }
  //bool exact = man->result.flag_exact;
  //elina_intlinearize_alloc(man,abs,&env,&dim,discr);
  res.p = (elina_lincons0_t*)malloc(array->size*sizeof(elina_lincons0_t));
  for(i=0; i < array->size; i++){
	elina_linexpr0_t *linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
	//elina_linexpr0_init(&linexpr0,0);
	elina_scalar_t * scalar = elina_scalar_alloc();
	res.p[i].linexpr0 = linexpr0;
	res.p[i].scalar = scalar;
  }
  res.size = array->size;
  
  elina_intlinearize_elina_tcons0_array(&res,
				   array,env, dim.intdim,discr);
 
  if (res.size==1 && 
      res.p[0].linexpr0->size==0)
    goto elina_intlinearize_tcons0_array_exit;
  if (meet && boxize && 
      (intervalonly ? !elina_lincons0_array_is_quasilinear(&res) : true)){
    tchange = malloc((dim.intdim+dim.realdim)*2);
    for (i=0;i<(dim.intdim+dim.realdim)*2;i++) tchange[i]=false;
    change = elina_boxize_lincons0_array(env,tchange,&res,env,dim.intdim,kmax,intervalonly,discr);
    //change = false;
  }
  switch(linearize){
  case ELINA_LINEXPR_INTLINEAR:
    break;
  case ELINA_LINEXPR_QUASILINEAR:
    quasilinearize_elina_lincons0_array(&res,env,meet,discr);
    break;
  case ELINA_LINEXPR_LINEAR:
    quasilinearize_elina_lincons0_array(&res,env,meet,discr);
    linearize_elina_lincons0_array(&res,meet, discr);
    break;
  }
  if (res.size==1 && 
      res.p[0].linexpr0->size==0 &&
      eval_elina_cstlincons0(&res.p[0])==0){
    goto elina_intlinearize_tcons0_array_exit;
   }
  if (change){
    if (elina_interval_is_bottom(env[0])){
      elina_lincons0_array_reinit(&res,1);
      elina_lincons0_set_bool(&res.p[0],false, discr);
      goto elina_intlinearize_tcons0_array_exit;
    }
    size = res.size;
    for (i=0;i<dim.intdim+dim.realdim; i++){
      if (tchange[2*i] || tchange[2*i+1]){
	/* we add a new row and prepare it */
	if (size>=res.size) 
	  elina_lincons0_array_reinit(&res,1+(5*res.size)/4);
	elina_linexpr0_reinit(res.p[size].linexpr0,1);
	res.p[size].linexpr0->p.linterm[0].dim = i;
      }	
      if ((tchange[2*i] || tchange[2*i+1]) && !elina_scalar_infty(env[i]->inf) && elina_scalar_equal(env[i]->inf,env[i]->sup)){
	/* We have a point */
	res.p[size].constyp = ELINA_CONS_EQ;
	elina_coeff_t * coeff = &res.p[size].linexpr0->p.linterm[0].coeff;
	elina_coeff_reinit(coeff,ELINA_COEFF_SCALAR,discr);
	elina_scalar_set_to_int(coeff->val.scalar,-1,discr);
	elina_coeff_set_scalar(&res.p[size].linexpr0->cst, env[i]->sup);
	size++;
      }
      else {
	if (tchange[2*i]){
	  /* inf bound */
	  assert(!elina_scalar_infty(env[i]->inf));
	  res.p[size].constyp = ELINA_CONS_SUPEQ;
	  elina_coeff_t * coeff = &res.p[size].linexpr0->p.linterm[0].coeff;
	  elina_coeff_reinit(coeff,ELINA_COEFF_SCALAR,discr);
	  elina_scalar_set_to_int(coeff->val.scalar,1,discr);
	  elina_scalar_neg(env[i]->inf,env[i]->inf);
	  elina_coeff_set_scalar(&res.p[size].linexpr0->cst,env[i]->inf);
	  elina_scalar_neg(env[i]->inf,env[i]->inf);
	  size ++;
	}
	if (tchange[2*i+1]){
	  /* sup bound */
	  if (tchange[2*i]){
	    /* we add a new row and prepare it */
	    if (size>=res.size) 
	      elina_lincons0_array_reinit(&res,1+(5*res.size)/4);
	    elina_linexpr0_reinit(res.p[size].linexpr0,1);
	    res.p[size].linexpr0->p.linterm[0].dim = i;
	  }	    
	  assert(!elina_scalar_infty(env[i]->sup));
	  res.p[size].constyp = ELINA_CONS_SUPEQ;
	  elina_coeff_t * coeff = &res.p[size].linexpr0->p.linterm[0].coeff;
          elina_coeff_reinit(coeff,ELINA_COEFF_SCALAR,discr);
	  elina_scalar_set_to_int(coeff->val.scalar,-1,discr);
	  elina_coeff_set_scalar(&res.p[size].linexpr0->cst, env[i]->sup);
	  size++;
	}
      }
    }
    elina_lincons0_array_reinit(&res,size);
  }
  
 elina_intlinearize_tcons0_array_exit:
  if (tchange) free(tchange);
  elina_interval_array_free(env,dim.intdim+dim.realdim);
  
  return res; 
}

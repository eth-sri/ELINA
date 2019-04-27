/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
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
/* elina_texpr0.c: tree expressions */
/* ************************************************************************* */


#include "elina_texpr0.h"
//#include "elina_linearize.h"

#include <stdarg.h>

/* ====================================================================== */
/* I. Constructors and Destructors */
/* ====================================================================== */

elina_texpr0_t* elina_texpr0_cst(elina_coeff_t* coeff)
{
  elina_texpr0_t* res = malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init_set(&res->val.cst,coeff);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_scalar(elina_scalar_t* scalar)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_SCALAR);
  elina_coeff_set_scalar(&res->val.cst, scalar);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_scalar_mpq(mpq_t mpq)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_SCALAR);
  elina_coeff_set_scalar_mpq(&res->val.cst, mpq);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_scalar_mpfr(mpfr_t mpfr)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_SCALAR);
  elina_coeff_set_scalar_mpfr(&res->val.cst, mpfr);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_scalar_int(long int num)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_SCALAR);
  elina_coeff_set_scalar_int(&res->val.cst, num);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_scalar_frac(long int num, unsigned long int den)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_SCALAR);
  elina_coeff_set_scalar_frac(&res->val.cst, num, den);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_scalar_double(double num)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_SCALAR);
  elina_coeff_set_scalar_double(&res->val.cst, num);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_interval(elina_interval_t* itv)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_INTERVAL);
  elina_coeff_set_interval(&res->val.cst, itv);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_interval_scalar(elina_scalar_t* inf, elina_scalar_t* sup)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_INTERVAL);
  elina_coeff_set_interval_scalar(&res->val.cst, inf, sup);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_interval_mpq(mpq_t inf, mpq_t sup)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_INTERVAL);
  elina_coeff_set_interval_mpq(&res->val.cst, inf, sup);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_interval_mpfr(mpfr_t inf, mpfr_t sup)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_INTERVAL);
  elina_coeff_set_interval_mpfr(&res->val.cst, inf, sup);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_interval_int(long int inf, long int sup)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_INTERVAL);
  elina_coeff_set_interval_int(&res->val.cst, inf, sup);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_interval_frac(long int numinf, unsigned long int deninf,
					 long int numsup, unsigned long int densup)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_INTERVAL);
  elina_coeff_set_interval_frac(&res->val.cst, numinf, deninf, numsup, densup);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_interval_double(double inf, double sup)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_INTERVAL);
  elina_coeff_set_interval_double(&res->val.cst, inf, sup);
  return res;
}
elina_texpr0_t* elina_texpr0_cst_interval_top(void)
{
  elina_texpr0_t* res = (elina_texpr0_t*) malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_CST;
  elina_coeff_init(&res->val.cst, ELINA_COEFF_INTERVAL);
  elina_interval_set_top(res->val.cst.val.interval);
  return res;
}
elina_texpr0_t* elina_texpr0_dim(elina_dim_t dim)
{
  elina_texpr0_t* res = malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_DIM;
  res->val.dim = dim;
  return res;
}
elina_texpr0_t* elina_texpr0_node(elina_texpr_op_t op, elina_texpr_rtype_t type, elina_texpr_rdir_t dir, elina_texpr0_t* opA, elina_texpr0_t* opB)
{
  elina_texpr0_node_t* node = malloc(sizeof(elina_texpr0_node_t));
  node->op = op;
  node->type = type;
  node->dir = dir;
  node->exprA = opA;
  node->exprB = opB;
  elina_texpr0_t* res = malloc(sizeof(elina_texpr0_t));
  res->discr = ELINA_TEXPR_NODE;
  res->val.node = node;
  return res;
}
elina_texpr0_t* elina_texpr0_unop(elina_texpr_op_t op, elina_texpr0_t* opA, elina_texpr_rtype_t type, elina_texpr_rdir_t dir)
{
  if (!elina_texpr_is_unop(op)){
    fprintf(stderr,"elina_texpr0.c: elina_texpr0_unop: unary operator expected\n");
    abort();
  }
  return elina_texpr0_node(op,type,dir,opA,NULL);
}
elina_texpr0_t* elina_texpr0_binop(elina_texpr_op_t op, elina_texpr0_t* opA, elina_texpr0_t* opB, elina_texpr_rtype_t type, elina_texpr_rdir_t dir)
{
  if (!elina_texpr_is_binop(op)){
    fprintf(stderr,"elina_texpr0.c: elina_texpr0_binop: binary operator expected\n");
    abort();
  }
  return elina_texpr0_node(op,type,dir,opA,opB);
}
elina_texpr0_t* elina_texpr0_node_copy(elina_texpr0_node_t* node)
{
  elina_texpr0_t* res = malloc(sizeof(elina_texpr0_t));
  elina_texpr0_node_t* n = malloc(sizeof(elina_texpr0_node_t));
  res->discr = ELINA_TEXPR_NODE;
  res->val.node = n;
  n->op = node->op;
  n->type = node->type;
  n->dir  = node->dir;
  n->exprA = elina_texpr0_copy(node->exprA);
  n->exprB = elina_texpr0_copy(node->exprB);
  return res;
}
elina_texpr0_t* elina_texpr0_copy(elina_texpr0_t* expr)
{
  if (!expr) return NULL;
  switch (expr->discr){
  case ELINA_TEXPR_CST:
    return elina_texpr0_cst(&expr->val.cst);
  case ELINA_TEXPR_DIM:
    return elina_texpr0_dim(expr->val.dim);
  case ELINA_TEXPR_NODE:
    return elina_texpr0_node_copy(expr->val.node);
  default:
    assert(false);
    return NULL;
  }
}
void elina_texpr0_node_free(elina_texpr0_node_t* node)
{
  elina_texpr0_free(node->exprA);
  elina_texpr0_free(node->exprB);
  free(node);
}
void elina_texpr0_clear(elina_texpr0_t* expr)
{
  switch(expr->discr){
  case ELINA_TEXPR_CST:
    elina_coeff_clear(&expr->val.cst);
    break;
  case ELINA_TEXPR_DIM:
    break;
  case ELINA_TEXPR_NODE:
    elina_texpr0_node_free(expr->val.node);
    break;
  default:
    assert(false);
  }
}
void elina_texpr0_free(elina_texpr0_t* expr)
{
  if (!expr) return;
  elina_texpr0_clear(expr);
  free(expr);
}
elina_texpr0_t* elina_texpr0_from_linexpr0(elina_linexpr0_t* e)
{
  elina_texpr0_t* res = elina_texpr0_cst(&e->cst);
  size_t i;
  elina_dim_t d;
  elina_coeff_t* c;
  elina_linexpr0_ForeachLinterm(e, i, d, c) {
    res = elina_texpr0_binop(ELINA_TEXPR_ADD,
			  res,
			  elina_texpr0_binop(ELINA_TEXPR_MUL,
					  elina_texpr0_cst(c), elina_texpr0_dim(d),
					  ELINA_RTYPE_REAL, ELINA_RDIR_RND),
			  ELINA_RTYPE_REAL, ELINA_RDIR_RND);
  }
  return res;
}

/* ====================================================================== */
/* II. Printing */
/* ====================================================================== */

static const char* elina_texpr_op_name[] =
  { "+", "-", "*", "/", "%", "^", /* binary */
    "-", "cast", "sqrt",          /* unary */
  };

static const int elina_texpr_op_precedence[] =
  { 1, 1, 2, 2, 2, 3, /* binary */
    4, 5, 5           /* unary */
  };

static const char* elina_texpr_rtype_name[] =
  { "", "i", "f", "d", "l", "q", };

static const char* elina_texpr_rdir_name[] =
  { "n", "0", "+oo", "-oo", "?", "", };

/* node induces some rounding (to float or integer) */
static inline bool elina_texpr0_node_exact(elina_texpr0_node_t* a)
{
  if (a->op==ELINA_TEXPR_NEG || a->op==ELINA_TEXPR_MOD ||
      a->type==ELINA_RTYPE_REAL) return true;
  return false;
}

static inline int elina_texpr0_precedence(elina_texpr0_t* a)
{
  if (!a || a->discr!=ELINA_TEXPR_NODE) return elina_texpr_op_precedence[ELINA_TEXPR_NEG];
  return elina_texpr_op_precedence[a->val.node->op];
}

static void elina_texpr0_node_fprint(FILE* stream, elina_texpr0_node_t* a,
				  char** name_of_dim)
{
  int prec = elina_texpr_op_precedence[a->op];

  /* left argument (if binary) */
  if (a->exprB) {
    int prec2 = elina_texpr0_precedence(a->exprA);
    if (prec2<prec) fprintf(stream, "(");
    elina_texpr0_fprint(stream, a->exprA, name_of_dim);
    if (prec2<prec) fprintf(stream, ")");
  }

  /* operator & rounding mode */
  if (a->exprB) fprintf(stream, " ");
  fprintf(stream, "%s", elina_texpr_op_name[a->op]);
  if (!elina_texpr0_node_exact(a))
    fprintf(stream, "_%s,%s",
	    elina_texpr_rtype_name[a->type], elina_texpr_rdir_name[a->dir]);

  /* right argument */
  {
    elina_texpr0_t* arg = a->exprB ? a->exprB : a->exprA;
    int prec2 = elina_texpr0_precedence(arg);
    if (a->exprB) fprintf(stream, " ");
    if (prec2<=prec) fprintf(stream, "(");
    elina_texpr0_fprint(stream,arg,name_of_dim);
    if (prec2<=prec) fprintf(stream, ")");
  }
}

void elina_texpr0_fprint(FILE* stream, elina_texpr0_t* a, char** name_of_dim)
{
  if (!a) return;
  switch (a->discr) {
  case ELINA_TEXPR_CST:
    elina_coeff_fprint(stream, &a->val.cst);
    break;
  case ELINA_TEXPR_DIM:
    if (name_of_dim) fprintf(stream, "%s", name_of_dim[a->val.dim]);
    else             fprintf(stream, "x%lu", (unsigned long)a->val.dim);
    break;
  case ELINA_TEXPR_NODE:
    elina_texpr0_node_fprint(stream, a->val.node, name_of_dim);
    break;
  default:
    assert(false);
  }
}

void elina_texpr0_print(elina_texpr0_t* a, char** name_of_dim)
{ elina_texpr0_fprint(stdout, a, name_of_dim); }


/* ====================================================================== */
/* III. Tests, size */
/* ====================================================================== */

size_t elina_texpr0_depth(elina_texpr0_t* a)
{
  int l,r;
  if (!a) return 0;
  switch(a->discr) {
  case ELINA_TEXPR_CST:
  case ELINA_TEXPR_DIM:
    return 0;
  case ELINA_TEXPR_NODE:
    l = elina_texpr0_depth(a->val.node->exprA);
    r = elina_texpr0_depth(a->val.node->exprB);
    return 1 + (l>r ? l : r);
  default:
    assert(0);
    return 0;
  }
}

size_t elina_texpr0_size(elina_texpr0_t* a)
{
  if (!a) return 0;
  switch(a->discr) {
  case ELINA_TEXPR_CST:
  case ELINA_TEXPR_DIM:
    return 0;
  case ELINA_TEXPR_NODE:
    return 1 + elina_texpr0_size(a->val.node->exprA) + elina_texpr0_size(a->val.node->exprB);
  default:
    assert(0);
    return 0;
  }
}

/* maximum between all dimensions and max */
static elina_dim_t elina_texpr0_max_dim_internal(elina_texpr0_t* a, elina_dim_t max)
{
  if (!a) return max;
  switch(a->discr) {
  case ELINA_TEXPR_CST:
    return max;
  case ELINA_TEXPR_DIM:
    return (a->val.dim+1) > max ? (a->val.dim+1) : max;
  case ELINA_TEXPR_NODE:
    return elina_texpr0_max_dim_internal(a->val.node->exprB,
				      elina_texpr0_max_dim_internal(a->val.node->exprA,max));
  default:
    assert(0);
    return max;
  }
}

elina_dim_t elina_texpr0_max_dim(elina_texpr0_t* a)
{
  return elina_texpr0_max_dim_internal(a, 0);
}

bool elina_texpr0_has_dim(elina_texpr0_t* a, elina_dim_t d)
{
  if (!a) return false;
  switch(a->discr) {
  case ELINA_TEXPR_CST:
    return false;
  case ELINA_TEXPR_DIM:
    return a->val.dim == d;
  case ELINA_TEXPR_NODE:
    return
      elina_texpr0_has_dim(a->val.node->exprA, d) ||
      elina_texpr0_has_dim(a->val.node->exprB, d);
  default:
    assert(0);
    return false;
  }
}

/* fill in v, v should be pre-allocated with size max_dim */
static void elina_texpr0_dimlist_internal(elina_texpr0_t* a, char* v)
{
  if (!a) return;
  switch(a->discr) {
  case ELINA_TEXPR_CST:
    break;
  case ELINA_TEXPR_DIM:
    v[a->val.dim] = 1;
    break;
  case ELINA_TEXPR_NODE:
    elina_texpr0_dimlist_internal(a->val.node->exprA, v);
    elina_texpr0_dimlist_internal(a->val.node->exprB, v);
    break;
  default:
    assert(0);
  }
}

elina_dim_t* elina_texpr0_dimlist(elina_texpr0_t* a)
{
  elina_dim_t max,i,nb;
  elina_dim_t* d;
  char* v;

  /* compute occurence vector */
  max = elina_texpr0_max_dim(a);
  if (max==0){
    /* constant expression */
    d = malloc(sizeof(elina_dim_t));
    d[0] = ELINA_DIM_MAX;
  }
  else {
    /* get number of distinct variables */
    v = malloc(max);
    memset(v, 0, max);
    elina_texpr0_dimlist_internal(a, v);
    for (i=0, nb=0; i<max; i++)
      if (v[i]) nb++;

    /* create & fill list */
    d = malloc(sizeof(elina_dim_t) * (nb+1));
    for (i=0, nb=0; i<max; i++)
      if (v[i]) { d[nb] = i; nb++; }
    d[nb] = ELINA_DIM_MAX; /* terminator */

    /* clean-up */
    free (v);
  }
  return d;
}

bool elina_texpr0_is_interval_cst(elina_texpr0_t* a)
{
  if (!a) return true;
  switch(a->discr) {
  case ELINA_TEXPR_CST:
    return true;
  case ELINA_TEXPR_DIM:
    return false;
  case ELINA_TEXPR_NODE:
    return
      elina_texpr0_is_interval_cst(a->val.node->exprA) &&
      elina_texpr0_is_interval_cst(a->val.node->exprB);
  default:
    assert(0);
    return false;
  }
}

bool elina_texpr0_is_scalar(elina_texpr0_t* a)
{
  if (!a) return true;
  switch(a->discr) {
  case ELINA_TEXPR_CST:
    return
      (a->val.cst.discr==ELINA_COEFF_SCALAR) ||
      (elina_scalar_equal(a->val.cst.val.interval->inf,
		       a->val.cst.val.interval->sup));
  case ELINA_TEXPR_DIM:
    return true;
  case ELINA_TEXPR_NODE:
    return
      elina_texpr0_is_scalar(a->val.node->exprA) &&
      elina_texpr0_is_scalar(a->val.node->exprB);
  default:
    assert(0);
    return false;
  }
}

bool elina_texpr0_is_interval_linear(elina_texpr0_t* a)
{
  if (!a) return true;
  switch(a->discr) {
  case ELINA_TEXPR_CST:
  case ELINA_TEXPR_DIM:
    return true;
  case ELINA_TEXPR_NODE:
    switch (a->val.node->op) {
    case ELINA_TEXPR_NEG:
      return elina_texpr0_is_interval_linear(a->val.node->exprA);
    case ELINA_TEXPR_CAST:
      return
	elina_texpr0_node_exact(a->val.node) &&
	elina_texpr0_is_interval_linear(a->val.node->exprA);
    case ELINA_TEXPR_ADD:
    case ELINA_TEXPR_SUB:
      return
	elina_texpr0_node_exact(a->val.node) &&
	elina_texpr0_is_interval_linear(a->val.node->exprA) &&
	elina_texpr0_is_interval_linear(a->val.node->exprB);
    case ELINA_TEXPR_MUL:
      return
	elina_texpr0_node_exact(a->val.node) &&
	( (elina_texpr0_is_interval_linear(a->val.node->exprA) &&
	   elina_texpr0_is_interval_cst(a->val.node->exprB)) ||
	  (elina_texpr0_is_interval_linear(a->val.node->exprB) &&
	   elina_texpr0_is_interval_cst(a->val.node->exprA)) );
    case ELINA_TEXPR_DIV:
      return
	elina_texpr0_node_exact(a->val.node) &&
	elina_texpr0_is_interval_linear(a->val.node->exprA) &&
	elina_texpr0_is_interval_cst(a->val.node->exprB);
    default:
      return false;
    }
  default:
    assert(0);
    return 0;
  }
}

bool elina_texpr0_is_interval_polynomial(elina_texpr0_t* a)
{
  if (!a) return true;
  switch(a->discr) {
  case ELINA_TEXPR_CST:
  case ELINA_TEXPR_DIM:
    return true;
  case ELINA_TEXPR_NODE:
    switch (a->val.node->op) {
    case ELINA_TEXPR_NEG:
      return elina_texpr0_is_interval_polynomial(a->val.node->exprA);
    case ELINA_TEXPR_CAST:
      return
	elina_texpr0_node_exact(a->val.node) &&
	elina_texpr0_is_interval_polynomial(a->val.node->exprA);
    case ELINA_TEXPR_ADD:
    case ELINA_TEXPR_SUB:
    case ELINA_TEXPR_MUL:
     return
	elina_texpr0_node_exact(a->val.node) &&
	elina_texpr0_is_interval_polynomial(a->val.node->exprA) &&
	elina_texpr0_is_interval_polynomial(a->val.node->exprB);
    case ELINA_TEXPR_DIV:
      return
	elina_texpr0_node_exact(a->val.node) &&
	elina_texpr0_is_interval_polynomial(a->val.node->exprA) &&
	elina_texpr0_is_interval_cst(a->val.node->exprB);
    case ELINA_TEXPR_POW:
      return
	elina_texpr0_node_exact(a->val.node) &&
        elina_texpr0_is_interval_polynomial(a->val.node->exprA) &&
        elina_texpr0_is_interval_cst(a->val.node->exprB); /* check for positivity? */
    default:
      return false;
    }
  default:
    assert(0);
    return 0;
  }
}

bool elina_texpr0_is_interval_polyfrac(elina_texpr0_t* a)
{
  if (!a) return true;
  switch(a->discr) {
  case ELINA_TEXPR_CST:
  case ELINA_TEXPR_DIM:
    return true;
  case ELINA_TEXPR_NODE:
    switch (a->val.node->op) {
    case ELINA_TEXPR_NEG:
      return elina_texpr0_is_interval_polyfrac(a->val.node->exprA);
    case ELINA_TEXPR_CAST:
      return
	elina_texpr0_node_exact(a->val.node) &&
	elina_texpr0_is_interval_polyfrac(a->val.node->exprA);
    case ELINA_TEXPR_ADD:
    case ELINA_TEXPR_SUB:
    case ELINA_TEXPR_MUL:
    case ELINA_TEXPR_DIV:
      return
	elina_texpr0_node_exact(a->val.node) &&
	elina_texpr0_is_interval_polyfrac(a->val.node->exprA) &&
	elina_texpr0_is_interval_polyfrac(a->val.node->exprB);
    case ELINA_TEXPR_POW:
      return
	elina_texpr0_node_exact(a->val.node) &&
        elina_texpr0_is_interval_polyfrac(a->val.node->exprA) &&
        elina_texpr0_is_interval_cst(a->val.node->exprB);
    default:
      return false;
    }
  default:
    assert(0);
    return 0;
  }
}

static bool elina_texpr0_array_is_template(elina_texpr0_t** texpr, size_t size, bool (*is_template)(elina_texpr0_t* texpr))
{
  size_t i;
  bool res = true;
  for (i=0; i<size; i++){
    res = is_template(texpr[i]);
    if (!res) break;
  }
  return res;
}
bool elina_texpr0_array_is_interval_linear(elina_texpr0_t** texpr, size_t size)
{
  return elina_texpr0_array_is_template(texpr,size,elina_texpr0_is_interval_linear);
}
bool elina_texpr0_array_is_interval_polynomial(elina_texpr0_t** texpr, size_t size)
{
return elina_texpr0_array_is_template(texpr,size,elina_texpr0_is_interval_polynomial);
}
bool elina_texpr0_array_is_interval_polyfrac(elina_texpr0_t** texpr, size_t size)
{
  return elina_texpr0_array_is_template(texpr,size,elina_texpr0_is_interval_polyfrac);
}
bool elina_texpr0_array_is_scalar(elina_texpr0_t** texpr, size_t size)
{
  return elina_texpr0_array_is_template(texpr,size,elina_texpr0_is_scalar);
}

/* ====================================================================== */
/* IV. Operations */
/* ====================================================================== */

elina_texpr0_t* elina_texpr0_substitute(elina_texpr0_t* a, elina_dim_t dim, elina_texpr0_t *dst)
{
  elina_texpr0_t* res = elina_texpr0_copy(a);
  elina_texpr0_substitute_with(res, dim, dst);
  return res;
}

void elina_texpr0_substitute_with(elina_texpr0_t* a, elina_dim_t dim, elina_texpr0_t *dst)
{
  if (!a) return;
  switch(a->discr) {
  case ELINA_TEXPR_CST:
    return;
  case ELINA_TEXPR_DIM:
    if (a->val.dim!=dim) return;
    dst = elina_texpr0_copy(dst);
    *a = *dst;
    free(dst);
    return;
  case ELINA_TEXPR_NODE:
    elina_texpr0_substitute_with(a->val.node->exprA, dim, dst);
    elina_texpr0_substitute_with(a->val.node->exprB, dim, dst);
    break;
  default:
    assert(0);
  }
}

/* ====================================================================== */
/* V. Change of dimensions and permutations */
/* ====================================================================== */

elina_texpr0_t* elina_texpr0_add_dimensions(elina_texpr0_t* expr,
				      elina_dimchange_t* dimchange)
{
  elina_texpr0_t* res = elina_texpr0_copy(expr);
  elina_texpr0_add_dimensions_with(res,dimchange);
  return res;
}
void elina_texpr0_add_dimensions_with(elina_texpr0_t* expr,
				   elina_dimchange_t* dimchange)
{
  if (!expr) return;
  switch(expr->discr){
  case ELINA_TEXPR_CST:
    return;
  case ELINA_TEXPR_DIM:
    {
      size_t dimsup = dimchange->intdim+dimchange->realdim;
      size_t k = 0;
      while (k<dimsup && expr->val.dim>=dimchange->dim[k]){
	k++;
      }
      expr->val.dim += k;
    }
    return;
  case ELINA_TEXPR_NODE:
    elina_texpr0_add_dimensions_with(expr->val.node->exprA,dimchange);
    elina_texpr0_add_dimensions_with(expr->val.node->exprB,dimchange);
    return;
  default:
    assert(false);
  }
}
elina_texpr0_t* elina_texpr0_remove_dimensions(elina_texpr0_t* expr,
					 elina_dimchange_t* dimchange)
{
  elina_texpr0_t* res = elina_texpr0_copy(expr);
  elina_texpr0_remove_dimensions_with(res,dimchange);
  return res;
}
void elina_texpr0_remove_dimensions_with(elina_texpr0_t* expr,
				      elina_dimchange_t* dimchange)
{
  if (!expr) return;
  switch(expr->discr){
  case ELINA_TEXPR_CST:
    return;
  case ELINA_TEXPR_DIM:
    {
      size_t dimrem = dimchange->intdim+dimchange->realdim;
      size_t i;
      for (i=0;i<dimrem && expr->val.dim>dimchange->dim[i];i++);
      if (i<dimrem && expr->val.dim==dimchange->dim[i]) {
	/* replace variable with top */
	expr->discr = ELINA_TEXPR_CST;
	elina_coeff_init(&expr->val.cst, ELINA_COEFF_INTERVAL);
	elina_interval_set_top(expr->val.cst.val.interval);
      }
      else expr->val.dim -= i;
    }
    return;
  case ELINA_TEXPR_NODE:
    elina_texpr0_remove_dimensions_with(expr->val.node->exprA,dimchange);
    elina_texpr0_remove_dimensions_with(expr->val.node->exprB,dimchange);
    return;
  default:
    assert(false);
  }
}
elina_texpr0_t* elina_texpr0_permute_dimensions(elina_texpr0_t* expr,
					  elina_dimperm_t* dimperm)
{
  elina_texpr0_t* res = elina_texpr0_copy(expr);
  elina_texpr0_permute_dimensions_with(res,dimperm);
  return res;
}
void elina_texpr0_permute_dimensions_with(elina_texpr0_t* expr,
				       elina_dimperm_t* perm)
{
  if (!expr) return;
  switch(expr->discr){
  case ELINA_TEXPR_CST:
    return;
  case ELINA_TEXPR_DIM:
    expr->val.dim = perm->dim[expr->val.dim];
    return;
  case ELINA_TEXPR_NODE:
    elina_texpr0_permute_dimensions_with(expr->val.node->exprA,perm);
    elina_texpr0_permute_dimensions_with(expr->val.node->exprB,perm);
    return;
  default:
    assert(false);
  }
}

/* ====================================================================== */
/* VI. Hashing, comparisons */
/* ====================================================================== */

long elina_texpr0_hash(elina_texpr0_t* a)
{
  if (!a) return 0;
  switch(a->discr) {
  case ELINA_TEXPR_CST:
    return elina_coeff_hash(&a->val.cst);
  case ELINA_TEXPR_DIM:
    return a->val.dim;
  case ELINA_TEXPR_NODE:
    return
      a->val.node->op * 17 +
      a->val.node->type * 23 +
      a->val.node->dir * 4801 +
      elina_texpr0_hash(a->val.node->exprA) * 17053 +
      elina_texpr0_hash(a->val.node->exprB);
  default:
    assert(0);
    return 0;
  }
}

bool elina_texpr0_equal(elina_texpr0_t* a1, elina_texpr0_t* a2)
{
  if (!a1 && !a2) return true;
  if (!a1 || !a2) return false;
  if (a1->discr!=a2->discr) return false;
  switch(a1->discr) {
  case ELINA_TEXPR_CST:
    return elina_coeff_equal(&a1->val.cst, &a2->val.cst);
  case ELINA_TEXPR_DIM:
    return a1->val.dim==a2->val.dim;
  case ELINA_TEXPR_NODE:
    return
      (a1->val.node->op==a2->val.node->op) &&
      (a1->val.node->type==a2->val.node->type) &&
      (a1->val.node->dir==a2->val.node->dir) &&
      elina_texpr0_equal(a1->val.node->exprA, a2->val.node->exprA) &&
      elina_texpr0_equal(a1->val.node->exprB, a2->val.node->exprB);
  default:
    assert(0);
    return false;
  }
}

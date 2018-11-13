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

#ifndef __OPT_ZONES_INTERNAL_H__
#define __OPT_ZONES_INTERNAL_H__

#ifdef __cplusplus
extern "C" {
#endif

#define min(a,b) a < b ? a : b

#define max(a,b) a > b ? a : b

#if defined(VECTOR)

 #include <immintrin.h>
 #include "vector_intrin.h"

#else
 #define v_length 1

#endif

#if defined(THRESHOLD)
#define zone_sparse_threshold THRESHOLD

#else
#define zone_sparse_threshold 0.5

#endif

#define zone_num_incomplete  1

#define zones_flag_incomplete						\
  man->result.flag_exact = man->result.flag_best = false

#define zone_flag_algo zones_flag_incomplete

#define zone_flag_conv zones_flag_incomplete

#include <stdio.h>

#include <assert.h>
#include "opt_zones.h"
#include "comp_list.h"
#include "rdtsc.h"

typedef struct opt_zones_internal_t{
  /* Name of function */
  elina_funid_t funid;

  /* local parameters for the function */
  elina_funopt_t* funopt;

  /* growing temporary buffer */
  double* tmp;
  void* tmp2;
  int tmp_size;

  /* raised when a conversion from/to a user type resulted in an
     overapproximation
  */
  bool conv;

  /* pointer to elina_manager*/
  elina_manager_t* man;
}opt_zones_internal_t;

typedef struct opt_zones_mat_t{
	double *mat;
	array_comp_list_t *acl;
	int nni;
	bool is_top;
	bool is_dense;
	bool ti;
	int ind;
}opt_zones_mat_t;

typedef struct opt_zones_t{
	opt_zones_mat_t *m;
	opt_zones_mat_t *closed;
	unsigned short int intdim;
	unsigned short int dim;
}opt_zones_t;


typedef enum{
    OPT_EMPTY,    /* OPT_EMPTY domain */
    OPT_ZERO,     /* 0 */
    OPT_UNARY,    /* UNARY unit expression */
    OPT_BINARY,   /* BINARY unit expression that can be expressed in zones */
    OPT_OTHER,   /* not expressible in zones */
  } zone_expr_type;


typedef struct zone_expr{
  zone_expr_type type; 
  /* index and coefficient for OPT_UNARY / OPT_BINARY unit expressions */
  unsigned short int i,j;
  short int coef_i,coef_j; /* -1 or 1 */

  /* expression has integer value */
  bool is_int;

}zone_expr;

static const int opt_zones_pre_widening = 99;

static inline bool is_integer(double num){
	return (num != INFINITY) && (ceil(num) == num);
}

static inline void init_array_zones(double *arr, int size){
	#if defined(VECTOR)
		v_double_type OPT_ZERO = v_set1_double(0);
		for(int i = 0; i < size/v_length; i++){
			v_store_double(arr + i*v_length,OPT_ZERO);
		}
	#else
		for(int i=0;i<(size/v_length)*v_length;i++){
			arr[i] = 0;
		}
	#endif
	for(int i = (size/v_length)*v_length; i < size; i++){
		arr[i] = 0;
	}
	
}


/* called by each function to setup and get manager-local data */
static inline opt_zones_internal_t*
opt_zones_init_from_manager(elina_manager_t* man, elina_funid_t id, int size)
{
 
  opt_zones_internal_t* pr = (opt_zones_internal_t*) man->internal;
  
  pr->funid = id;
  pr->funopt = man->option.funopt+id;
  man->result.flag_exact = man->result.flag_best = true;
  pr->conv = false;
  if (pr->tmp_size<size) {
    pr->tmp = (double*)realloc(pr->tmp,sizeof(double)*size);
    assert(pr->tmp);
    pr->tmp_size = size;
    init_array_zones(pr->tmp,pr->tmp_size);
    pr->tmp2 = realloc(pr->tmp2,sizeof(long)*size);
    assert(pr->tmp2);
  }
  return pr;
}

/*****
Conversion to user types
****/

/* upper bound => scalar, with optional division by 2
   pr->conv is set if the conversion is not exact
*/
static inline void zones_scalar_of_upper_bound(elina_scalar_t* r, double d)
{
  elina_scalar_reinit(r,ELINA_SCALAR_DOUBLE);
  if (!isfinite(d)) elina_scalar_set_infty(r,1);
  else {
      r->val.dbl = d;
  }
}

/* opposite of lower bound => scalar, with optional division by 2
   pr->conv is set if the conversion is not exact
*/
static inline void zones_scalar_of_lower_bound(elina_scalar_t* r, double d)
{
  elina_scalar_reinit(r,ELINA_SCALAR_DOUBLE);
  if (d== INFINITY) elina_scalar_set_infty(r,-1);
  else {
      r->val.dbl = d;
      r->val.dbl = -r->val.dbl;
     //We don't handle MPQ and MPFR
  }
}

static inline void zones_interval_of_bounds(elina_interval_t* i, double minf, double sup)
{
  zones_scalar_of_upper_bound(i->sup, sup);
  zones_scalar_of_lower_bound(i->inf,minf);
}

static inline bool double_set_elina_scalar(double *r, elina_scalar_t *t){
	switch(t->discr){
		case ELINA_SCALAR_DOUBLE:
			*r = t->val.dbl;
			return true;
		default:
			abort();
	}
}

static inline void opt_bound_of_scalar(opt_zones_internal_t* pr,
				   double *r, elina_scalar_t* t,
				   bool neg)
{
  
  if (neg) elina_scalar_neg(t,t);
  if (!elina_double_set_scalar(r,t,GMP_RNDU)) pr->conv = true;
  if (neg) elina_scalar_neg(t,t);
}

static inline bool opt_bounds_of_coeff(opt_zones_internal_t* pr,
				   double *minf, double *sup,
				   elina_coeff_t c)
{
  switch (c.discr) {
  case ELINA_COEFF_SCALAR:
    opt_bound_of_scalar(pr,minf,c.val.scalar,true);
    opt_bound_of_scalar(pr,sup,c.val.scalar,false);
    return false;
  case ELINA_COEFF_INTERVAL:
    opt_bound_of_scalar(pr,minf,c.val.interval->inf,true);
    opt_bound_of_scalar(pr,sup,c.val.interval->sup,false);
    return elina_scalar_cmp(c.val.interval->inf,c.val.interval->sup)>0;
  default: abort();
	/*******
		TODO: handle arg_assert
		arg_assert(0,return false;);
	********/
  }
}

static inline elina_lincons0_t zones_lincons_of_bound(opt_zones_internal_t* pr,
					     int i, int j,
					     double  d)
{
  elina_linexpr0_t* e;
  if (i==j) {
    /* OPT_ZEROary constraint */
    e = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 0);
    zones_scalar_of_upper_bound(e->cst.val.scalar,d);
  }
  else if (i==0) {
    /* upper bound */
    e = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 1);
    e->p.linterm[0].dim = j-1;
    elina_scalar_set_int(e->p.linterm[0].coeff.val.scalar, -1);
    zones_scalar_of_upper_bound(e->cst.val.scalar,d);
  }
  else if(j==0){
	/* lower bound */
       e = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 1);
       e->p.linterm[0].dim = i-1;
       elina_scalar_set_int(e->p.linterm[0].coeff.val.scalar, 1);
       zones_scalar_of_upper_bound(e->cst.val.scalar,d);
  }
  else {
    /* OPT_BINARY constraint */
    e = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 2);
    e->p.linterm[0].dim = j-1;
    e->p.linterm[1].dim = i-1;
    elina_scalar_set_int(e->p.linterm[0].coeff.val.scalar, -1);
    elina_scalar_set_int(e->p.linterm[1].coeff.val.scalar, 1);
    zones_scalar_of_upper_bound(e->cst.val.scalar,d);
  }
  return elina_lincons0_make(ELINA_CONS_SUPEQ,e,NULL);
}


/**************
	Basic operators
*************/
opt_zones_t * opt_zones_alloc_internal(opt_zones_internal_t *pr, unsigned short int dim, unsigned short int intdim);
void opt_zones_fprint(FILE* stream, elina_manager_t* man, opt_zones_t * o,char** name_of_dim);
int opt_zones_size(elina_manager_t* man, opt_zones_t* o);
elina_dimension_t opt_zones_dimension(elina_manager_t* man, opt_zones_t* o);
opt_zones_t* opt_zones_top(elina_manager_t* man, unsigned short int intdim, unsigned short int realdim);
opt_zones_t* opt_zones_bottom(elina_manager_t* man, unsigned short int intdim, unsigned short int realdim);
void opt_zones_free(elina_manager_t* man, opt_zones_t* o);
opt_zones_t * opt_zones_copy_internal(opt_zones_internal_t *pr, opt_zones_t *o);
opt_zones_t* opt_zones_copy(elina_manager_t* man, opt_zones_t* o);
void opt_zones_minimize(elina_manager_t* man, opt_zones_t* o);
void opt_zones_canonicalize(elina_manager_t* man, opt_zones_t* o);
int opt_zones_hash(elina_manager_t* man, opt_zones_t* o);
void opt_zones_approximate(elina_manager_t* man, opt_zones_t* o, int algorithm);
void opt_zones_sparse_weak_closure(opt_zones_internal_t *pr,  opt_zones_t *o);
void opt_zones_cache_closure(opt_zones_internal_t *pr, opt_zones_t *o);
opt_zones_t* opt_zones_closure(elina_manager_t *man, bool destructive, opt_zones_t *o);
opt_zones_t* opt_zones_set_mat(opt_zones_internal_t* pr, opt_zones_t* o, 
			       opt_zones_mat_t* m, opt_zones_mat_t* closed, 
			       bool destructive);


/**************
	Predicate operators
***************/
bool opt_zones_is_bottom(elina_manager_t* man, opt_zones_t* o);
bool opt_zones_is_top(elina_manager_t* man, opt_zones_t* o);
bool opt_zones_is_leq(elina_manager_t* man, opt_zones_t* o1, opt_zones_t* o2);
bool opt_zones_is_eq(elina_manager_t* man, opt_zones_t* o1, opt_zones_t* o2);
elina_tcons0_array_t opt_zones_to_tcons_array(elina_manager_t* man, opt_zones_t* o);
elina_interval_t** opt_zones_to_box(elina_manager_t* man, opt_zones_t* o);
elina_lincons0_array_t opt_zones_to_lincons_array(elina_manager_t* man, opt_zones_t* o);
elina_interval_t* opt_zones_bound_dimension(elina_manager_t* man, opt_zones_t* o, elina_dim_t dim);
bool opt_zones_is_dimension_unconstrained(elina_manager_t* man, opt_zones_t* o, elina_dim_t dim);

/*****************
	Nary operators
******************/
opt_zones_t* opt_zones_meet(elina_manager_t* man, bool destructive, opt_zones_t* o1, opt_zones_t* o2);
opt_zones_t* opt_zones_join(elina_manager_t* man, bool destructive, opt_zones_t* o1, opt_zones_t* o2);
opt_zones_t* opt_zones_widening(elina_manager_t* man, opt_zones_t* o1, opt_zones_t* o2);

/********************
	Resize Operators
*********************/
opt_zones_t* opt_zones_forget_array(elina_manager_t* man, bool destructive, opt_zones_t* o, elina_dim_t* tdim, int size, bool project);
opt_zones_t* opt_zones_add_dimensions(elina_manager_t* man, bool destructive, opt_zones_t* o, elina_dimchange_t* dimchange, bool project);
opt_zones_t* opt_zones_remove_dimensions(elina_manager_t* man, bool destructive, opt_zones_t* o, elina_dimchange_t* dimchange);
opt_zones_t* opt_zones_permute_dimensions(elina_manager_t* man, bool destructive, opt_zones_t* o, elina_dimperm_t* permutation);
opt_zones_t* opt_zones_expand(elina_manager_t* man, bool destructive, opt_zones_t* o, elina_dim_t dim, size_t n);


/**********************
	Transfer functions
***********************/
opt_zones_t* opt_zones_assign_texpr_array(elina_manager_t* man, bool destructive, opt_zones_t* o, 
					  elina_dim_t* tdim, elina_texpr0_t** texpr, size_t size, opt_zones_t* dest);
opt_zones_t* opt_zones_meet_lincons_array(elina_manager_t* man, bool destructive, opt_zones_t* o, elina_lincons0_array_t* array);
opt_zones_t* opt_zones_meet_tcons_array(elina_manager_t* man, bool destructive, opt_zones_t* o, elina_tcons0_array_t* array);

#ifdef __cplusplus
}
#endif


#endif

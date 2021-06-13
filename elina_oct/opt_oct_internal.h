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



#ifndef __OPT_OCT_INTERNAL_H
#define __OPT_OCT_INTERNAL_H

#ifdef __cplusplus
extern "C" {
#endif


#if defined(NUM_DOUBLE)
#define is_int_flag 0

#else
#define is_int_flag 1

#endif


#if defined(THRESHOLD)
#define sparse_threshold THRESHOLD

#else
#define sparse_threshold 0.5

#endif
//#if defined(SPARSE)
//#define sparse_flag 1

//#else
//#define sparse_flag 0

//#endif


#if defined(VECTOR)

 #include <immintrin.h>
 #include "vector_intrin.h"

#else
 #define v_length 1

#endif

#if defined(TIMING)
  #include "rdtsc.h"
  		
#endif

#define num_incomplete  1

#define flag_incomplete						\
  man->result.flag_exact = man->result.flag_best = false

#define flag_algo flag_incomplete

#define flag_conv flag_incomplete

//#define INFINITY 1.0/0.0

#include <stdio.h>

#include <assert.h>
//#include <limits.h>
//#include <math.h>
#include "opt_oct.h"
#include "comp_list.h"

typedef struct opt_oct_internal_t{
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
}opt_oct_internal_t;

typedef struct opt_oct_mat_t{
	double *mat;
	array_comp_list_t *acl;
	int nni;
	bool is_top;
	bool is_dense;
	bool ti;
}opt_oct_mat_t;

typedef struct opt_oct_t{
	opt_oct_mat_t *m;
	opt_oct_mat_t *closed;
	int intdim;
	int dim;
}opt_oct_t;


typedef enum{
    OPT_EMPTY,    /* OPT_EMPTY domain */
    OPT_ZERO,     /* 0 */
    OPT_UNARY,    /* OPT_UNARY unit expression */
    OPT_BINARY,   /* OPT_BINARY unit expression */
    OPT_OTHER,
  } uexpr_type;


typedef struct opt_uexpr{
  uexpr_type type; 
  /* index and coefficient for OPT_UNARY / OPT_BINARY unit expressions */
  size_t i,j;
  int coef_i,coef_j; /* -1 or 1 */

  /* expression has integer value */
  int is_int;

}opt_uexpr;

static const int opt_oct_pre_widening = 99;


static inline void init_array(double *arr, int size){
	#if defined(VECTOR)
		v_double_type OPT_ZERO = v_set1_double(0);
		for(int i = 0; i < size/v_length; i++){
			v_store_double(arr + i*v_length,OPT_ZERO);
		}
	//}
	#else
		for(int i=0;i<(size/v_length)*v_length;i++){
			arr[i] = 0;
		}
	//}
	#endif
	for(int i = (size/v_length)*v_length; i < size; i++){
		arr[i] = 0;
	}
	
}

/* called by each function to setup and get manager-local data */
static inline opt_oct_internal_t*
opt_oct_init_from_manager(elina_manager_t* man, elina_funid_t id, int size)
{
  opt_oct_internal_t* pr = (opt_oct_internal_t*) man->internal;
  pr->funid = id;
  pr->funopt = man->option.funopt+id;
  man->result.flag_exact = man->result.flag_best = true;
  pr->conv = false;
  if (pr->tmp_size<size) {
    pr->tmp = (double*)realloc(pr->tmp,sizeof(double)*size);
    assert(pr->tmp);
    pr->tmp_size = size;
    init_array(pr->tmp,pr->tmp_size);
    pr->tmp2 = realloc(pr->tmp2,sizeof(long)*size);
    assert(pr->tmp2);
  }
  return pr;
}

static inline bool is_integer(double num){
	return (num != INFINITY) && (ceil(num) == num);
}


/****
	Size function
**/
static inline size_t opt_matsize(size_t dim)
{
  return 2 * dim * (dim+1);
}
/* position of (i,j) element, assuming j/2 <= i/2 */
static inline int opt_matpos(int i, int j)
{
  return j + ((i+1)*(i+1))/2;
}

/* position of (i,j) element, no assumption */
static inline int opt_matpos2(int i, int j)
{
  if (j>i) return opt_matpos(j^1,i^1);
  else return opt_matpos(i,j);
}

/*****
Conversion to user types
****/

/* upper bound => scalar, with optional division by 2
   pr->conv is set if the conversion is not exact
*/
static inline void opt_scalar_of_upper_bound(opt_oct_internal_t* pr,
					 elina_scalar_t* r,
					 double d,
					 bool div2)
{
  elina_scalar_reinit(r,ELINA_SCALAR_DOUBLE);
  if (!isfinite(d)) elina_scalar_set_infty(r,1);
  else {
      r->val.dbl = d;
      if(div2) pr->conv = 1;
      if (div2) r->val.dbl /= 2;
  }
}

/* opposite of lower bound => scalar, with optional division by 2
   pr->conv is set if the conversion is not exact
*/
static inline void opt_scalar_of_lower_bound(opt_oct_internal_t* pr,
					 elina_scalar_t* r,
					 double d,
					 bool div2)
{
  elina_scalar_reinit(r,ELINA_SCALAR_DOUBLE);
  if (d== INFINITY) elina_scalar_set_infty(r,-1);
  else {
      r->val.dbl = d;
      if(div2) pr->conv = 1;
      if (div2) r->val.dbl /= 2;
      r->val.dbl = -r->val.dbl;    
  }
}

static inline bool double_set_mpq_tmp(double *a, mpq_t b, mpfr_t mpfr)
{
  int res = mpfr_set_q(mpfr,b,GMP_RNDU);
#if defined(NUMFLT_DOUBLE)
  *a = mpfr_get_d(mpfr,GMP_RNDU);/* Normally, exact conversion here (unless overfloww) */
#else
  *a = mpfr_get_ld(mpfr,GMP_RNDU);/* Normally, exact conversion here (unless overfloww) */
#endif
  return (res==0);
}
static inline bool double_set_mpq(double *a, mpq_t b)
{
  mpfr_t mpfr;
  mpfr_init2(mpfr,DBL_MANT_DIG);
  bool res = double_set_mpq_tmp(a,b,mpfr);
  mpfr_clear(mpfr);
  return res;
}

static inline bool double_set_elina_scalar(double *r, elina_scalar_t *t){
	switch(t->discr){
		case ELINA_SCALAR_MPQ:
    			if (mpz_sgn(mpq_denref(t->val.mpq))==0){
    				if(mpz_sgn(mpq_numref(t->val.mpq))>0){
    					*r = INFINITY;
    				}
    				else{
    					*r = -INFINITY;
    				}
      				return true;
    			}
    			else {
      				return double_set_mpq(r,t->val.mpq);
    			}
    			break;
		case ELINA_SCALAR_DOUBLE:
			*r = t->val.dbl;
			return true;

		default:
			abort();
	}
}

static inline void opt_bound_of_scalar(opt_oct_internal_t* pr,
				   double *r, elina_scalar_t* t,
				   bool neg, bool mul2)
{
  if (neg) elina_scalar_neg(t,t);
  if (!elina_double_set_scalar(r,t,GMP_RNDU)) pr->conv = true;
  if (mul2) {
    *r = (*r)*2;
    pr->conv = true;
  }
  if (neg) elina_scalar_neg(t,t);
 
}

static inline bool opt_bounds_of_coeff(opt_oct_internal_t* pr,
				   double *minf, double *sup,
				   elina_coeff_t c,
				   bool mul2)
{
  
  switch (c.discr) {
  case ELINA_COEFF_SCALAR:
    opt_bound_of_scalar(pr,minf,c.val.scalar,true,mul2);
    opt_bound_of_scalar(pr,sup,c.val.scalar,false,mul2);
    return false;
  case ELINA_COEFF_INTERVAL:
    opt_bound_of_scalar(pr,minf,c.val.interval->inf,true,mul2);
    opt_bound_of_scalar(pr,sup,c.val.interval->sup,false,mul2);
    return elina_scalar_cmp(c.val.interval->inf,c.val.interval->sup)>0;
  default: abort();
	/*******
		TODO: handle arg_assert
		arg_assert(0,return false;);
	********/
  }
  
}


static inline elina_lincons0_t opt_lincons_of_bound(opt_oct_internal_t* pr,
					     int i, int j,
					     double  d)
{
  elina_linexpr0_t* e;
  if (i==j) {
    /* OPT_ZEROary constraint */
    e = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 0);
    opt_scalar_of_upper_bound(pr,e->cst.val.scalar,d,true);
  }
  else if (i==(j^1)) {
    /* OPT_UNARY constraint */
    e = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 1);
    e->p.linterm[0].dim = i/2;
    elina_scalar_set_int(e->p.linterm[0].coeff.val.scalar,(i&1) ? -1 : 1);
    opt_scalar_of_upper_bound(pr,e->cst.val.scalar,d,true);
  }
  else {
    /* OPT_BINARY constraint */
    e = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 2);
    e->p.linterm[0].dim = j/2;
    e->p.linterm[1].dim = i/2;
    elina_scalar_set_int(e->p.linterm[0].coeff.val.scalar,(j&1) ?  1 : -1);
    elina_scalar_set_int(e->p.linterm[1].coeff.val.scalar,(i&1) ? -1 :  1);
    opt_scalar_of_upper_bound(pr,e->cst.val.scalar,d,false);
  }
  return elina_lincons0_make(ELINA_CONS_SUPEQ,e,NULL);
}


/* makes an interval from [-minf,sup], with sound approximations
   pr->conv is set if the conversion is not exact
   note: may output an OPT_EMPTY interval
*/
static inline void opt_interval_of_bounds(opt_oct_internal_t* pr,
				      elina_interval_t* i,
				      double minf, double sup,
				      bool div2)
{
  opt_scalar_of_upper_bound(pr,i->sup, sup,div2);
  opt_scalar_of_lower_bound(pr,i->inf,minf,div2);
}



static inline void opt_oct_of_scalar(opt_oct_internal_t* pr,
				   double * r, elina_scalar_t* t,
				   bool neg, bool mul2)
{

  if (neg) elina_scalar_neg(t,t);
  double_set_elina_scalar(r,t);
  pr->conv = true;
  
  if (mul2) {
    *r = *r*2;
    pr->conv = true;
  }
  if (neg) elina_scalar_neg(t,t);
}


/* both bounds of an interval, the lower bound is negated
   pr->conv is set if the conversion is not exact
   returns true if the interval is empty
*/
static inline bool opt_oct_of_interval(opt_oct_internal_t* pr,
				      opt_oct_mat_t *oo, unsigned short int i,
				      elina_interval_t* itv,
				      bool mul2)
{
  double * m = oo->mat;
  
  
  if(!elina_interval_is_top(itv)){
  	
	comp_list_t * cl = create_comp_list();
	insert_comp(cl,i);
	insert_comp_list(oo->acl,cl);
	unsigned short int i1 = 2*i;
	unsigned short int i2 = 2*i+1;
	int ind1 = i2 + ((i1+1)*(i1+1))/2;
	int ind2 = i1 + ((i2+1)*(i2+1))/2;
	int ind3 = i2 + ((i2+1)*(i2+1))/2;
	int ind4 = i1 + ((i1+1)*(i1+1))/2;
	
	opt_oct_of_scalar(pr,&m[ind1],itv->inf,true,mul2);
  	opt_oct_of_scalar(pr,&m[ind2],itv->sup,false,mul2);
  	m[ind3] = 0;
  	m[ind4] = 0;
	oo->nni+=2;
	
	return elina_scalar_cmp(itv->inf,itv->sup)>0;
  }
  
  return false;
}

/**************
	Basic operations
*************/

opt_oct_t * opt_oct_alloc_internal(opt_oct_internal_t *pr, int dim, int intdim);
int opt_oct_size(elina_manager_t* man, opt_oct_t* o);
opt_oct_t * opt_oct_alloc_top(opt_oct_internal_t *pr, int dim, int intdim);
void opt_oct_free_internal(opt_oct_internal_t *pr, opt_oct_t *o);
opt_oct_t * opt_oct_copy_internal(opt_oct_internal_t *pr, opt_oct_t *o);
opt_oct_t* opt_oct_set_mat(opt_oct_internal_t* pr, opt_oct_t* o, opt_oct_mat_t* m, opt_oct_mat_t* closed, bool destructive);
opt_oct_t* opt_oct_copy(elina_manager_t* man, opt_oct_t* o);
void opt_oct_free(elina_manager_t* man, opt_oct_t* a);
opt_oct_t* opt_oct_bottom(elina_manager_t* man, int intdim, int realdim);
opt_oct_t* opt_oct_top(elina_manager_t* man, int intdim, int realdim);
opt_oct_t* opt_oct_of_box(elina_manager_t* man, size_t intdim, size_t realdim, elina_interval_t ** t);
elina_dimension_t opt_oct_dimension(elina_manager_t* man, opt_oct_t* o);
void opt_oct_cache_closure(opt_oct_internal_t *pr, opt_oct_t *o);
void opt_oct_close(opt_oct_internal_t *pr, opt_oct_t *o);
opt_oct_t* opt_oct_closure(elina_manager_t *man, bool destructive, opt_oct_t *o);
void opt_oct_internal_free(opt_oct_internal_t *pr);
opt_oct_t* opt_oct_of_abstract0(elina_abstract0_t* a);
elina_abstract0_t* abstract0_of_opt_oct(elina_manager_t* man, opt_oct_t* oct);
void opt_oct_minimize(elina_manager_t* man, opt_oct_t* o);
void opt_oct_canonicalize(elina_manager_t* man, opt_oct_t* o);
int opt_oct_hash(elina_manager_t* man, opt_oct_t* o);
void opt_oct_approximate(elina_manager_t* man, opt_oct_t* o, int algorithm);


/**************
	nary operations
*************/

opt_oct_t* opt_oct_meet(elina_manager_t* man, bool destructive, opt_oct_t* o1, opt_oct_t* o2);
opt_oct_t* opt_oct_join(elina_manager_t* man, bool destructive, opt_oct_t* o1, opt_oct_t* o2);
opt_oct_t* opt_oct_widening(elina_manager_t* man, opt_oct_t* o1, opt_oct_t* o2);
opt_oct_t* opt_oct_widening_thresholds(elina_manager_t* man, opt_oct_t* o1, opt_oct_t* o2, elina_scalar_t** array, size_t nb);
opt_oct_t* opt_oct_narrowing(elina_manager_t* man, opt_oct_t* o1, opt_oct_t* o2);
opt_oct_t* opt_oct_add_epsilon(elina_manager_t* man, opt_oct_t* o, elina_scalar_t* epsilon);
opt_oct_t* opt_oct_add_epsilon_bin(elina_manager_t* man, opt_oct_t* o1, opt_oct_t* o2, elina_scalar_t* epsilon);
opt_oct_t* opt_oct_join_array(elina_manager_t* man, opt_oct_t** tab, size_t size);
opt_oct_t* opt_oct_meet_array(elina_manager_t* man, opt_oct_t** tab, size_t size);


/**************
	predicate operations
*************/ 

bool opt_oct_is_bottom(elina_manager_t* man, opt_oct_t* o);
bool opt_oct_is_top(elina_manager_t* man, opt_oct_t* o);
bool opt_oct_is_leq(elina_manager_t* man, opt_oct_t* o1, opt_oct_t* o2);
bool opt_oct_is_eq(elina_manager_t* man, opt_oct_t* o1, opt_oct_t* o2);
elina_tcons0_array_t opt_oct_to_tcons_array(elina_manager_t* man, opt_oct_t* o);
elina_interval_t** opt_oct_to_box(elina_manager_t* man, opt_oct_t* o);
elina_interval_t* opt_oct_bound_texpr(elina_manager_t* man,opt_oct_t* o, elina_texpr0_t* expr);
elina_interval_t* opt_oct_bound_linexpr(elina_manager_t* man,opt_oct_t* o, elina_linexpr0_t* expr);
elina_interval_t* opt_oct_bound_dimension(elina_manager_t* man,opt_oct_t* o, elina_dim_t dim);
elina_lincons0_array_t opt_oct_to_lincons_array(elina_manager_t* man, opt_oct_t* o);
bool opt_oct_sat_interval(elina_manager_t* man, opt_oct_t* o, elina_dim_t dim, elina_interval_t* i);
bool opt_oct_is_dimension_unconstrained(elina_manager_t* man, opt_oct_t* o, elina_dim_t dim);
bool opt_oct_sat_lincons_timing(elina_manager_t* man, opt_oct_t* o, elina_lincons0_t* lincons);
bool opt_oct_sat_tcons(elina_manager_t* man, opt_oct_t* o, elina_tcons0_t* cons);


/**************
	resize operations
*************/ 

opt_oct_t* opt_oct_forget_array(elina_manager_t* man,bool destructive, opt_oct_t* o, elina_dim_t* tdim, int size, bool project);
opt_oct_t* opt_oct_add_dimensions(elina_manager_t* man, bool destructive, opt_oct_t* o, elina_dimchange_t* dimchange, bool project);
opt_oct_t* opt_oct_remove_dimensions(elina_manager_t* man, bool destructive, opt_oct_t* o, elina_dimchange_t* dimchange);
opt_oct_t* opt_oct_permute_dimensions(elina_manager_t* man, bool destructive, opt_oct_t* o, elina_dimperm_t* permutation);

/***************
	Transfer Operations
***************/

opt_oct_t* opt_oct_meet_lincons_array(elina_manager_t* man, bool destructive, opt_oct_t* o, elina_lincons0_array_t* array);
opt_oct_t* opt_oct_meet_tcons_array(elina_manager_t* man, bool destructive, opt_oct_t* o, elina_tcons0_array_t* array);
opt_oct_t* opt_oct_assign_linexpr_array(elina_manager_t* man, bool destructive, opt_oct_t* o, elina_dim_t* tdim, elina_linexpr0_t** texpr, size_t size, opt_oct_t* dest);
opt_oct_t* opt_oct_substitute_linexpr_array(elina_manager_t* man, bool destructive, opt_oct_t* o, elina_dim_t* tdim, elina_linexpr0_t** texpr, size_t size, opt_oct_t* dest);
opt_oct_t* opt_oct_assign_texpr_array(elina_manager_t* man, bool destructive, opt_oct_t* o, elina_dim_t* tdim, elina_texpr0_t** texpr, int size, opt_oct_t* dest);

#ifdef __cplusplus
}
#endif

#endif

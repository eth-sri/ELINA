/*
	Copyright 2015 Software Reliability Lab, ETH Zurich

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

#define num_incomplete  1

#define flag_incomplete						\
  man->result.flag_exact = man->result.flag_best = false

#define flag_algo flag_incomplete

#define flag_conv flag_incomplete

#define INFINITY 1.0/0.0

#include <stdio.h>

#include <assert.h>
//#include <limits.h>
//#include <math.h>
#include "opt_oct.h"
#include "comp_list.h"
#include "num.h"

typedef struct opt_oct_internal_t{
  /* Name of function */
  ap_funid_t funid;

  /* local parameters for the function */
  ap_funopt_t* funopt;

  /* growing temporary buffer */
  double* tmp;
  void* tmp2;
  int tmp_size;

  /* raised when a conversion from/to a user type resulted in an
     overapproximation
  */
  bool conv;

  /* pointer to ap_manager*/
  ap_manager_t* man;
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
  int i,j;
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
opt_oct_init_from_manager(ap_manager_t* man, ap_funid_t id, int size)
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
					 ap_scalar_t* r,
					 double d,
					 bool div2)
{
  ap_scalar_reinit(r,NUM_AP_SCALAR);
  if (!isfinite(d)) ap_scalar_set_infty(r,1);
  else {
    switch (NUM_AP_SCALAR) {
    case AP_SCALAR_DOUBLE:
	//The line below can cause errors
      //if (!double_set_num(&r->val.dbl,d) || div2) pr->conv = 1;
      r->val.dbl = d;
      if(div2) pr->conv = 1;
      if (div2) r->val.dbl /= 2;
      break;
   //We don't handle MPQ and MPFR
   
    default:
      abort();
    }
  }
}

/* opposite of lower bound => scalar, with optional division by 2
   pr->conv is set if the conversion is not exact
*/
static inline void opt_scalar_of_lower_bound(opt_oct_internal_t* pr,
					 ap_scalar_t* r,
					 double d,
					 bool div2)
{
  ap_scalar_reinit(r,NUM_AP_SCALAR);
  if (d== INFINITY) ap_scalar_set_infty(r,-1);
  else {
    switch (NUM_AP_SCALAR) {
    case AP_SCALAR_DOUBLE:
      r->val.dbl = d;
      if(div2) pr->conv = 1;
      if (div2) r->val.dbl /= 2;
      r->val.dbl = -r->val.dbl;
      break;
     //We don't handle MPQ and MPFR
    
    default:
      abort();
    }
  }
}

static inline bool double_set_ap_scalar(double *r, ap_scalar_t *t){
	switch(t->discr){
		case AP_SCALAR_DOUBLE:
			*r = t->val.dbl;
			return true;
		default:
			abort();
	}
}

static inline void opt_bound_of_scalar(opt_oct_internal_t* pr,
				   double *r, ap_scalar_t* t,
				   bool neg, bool mul2)
{
  if (neg) ap_scalar_neg(t,t);
  if (!double_set_ap_scalar(r,t)) pr->conv = true;
  if (mul2) {
    *r = (*r)*2;
    pr->conv = true;
  }
  if (neg) ap_scalar_neg(t,t);
}

static inline bool opt_bounds_of_coeff(opt_oct_internal_t* pr,
				   double *minf, double *sup,
				   ap_coeff_t c,
				   bool mul2)
{
  switch (c.discr) {
  case AP_COEFF_SCALAR:
    opt_bound_of_scalar(pr,minf,c.val.scalar,true,mul2);
    opt_bound_of_scalar(pr,sup,c.val.scalar,false,mul2);
    return false;
  case AP_COEFF_INTERVAL:
    opt_bound_of_scalar(pr,minf,c.val.interval->inf,true,mul2);
    opt_bound_of_scalar(pr,sup,c.val.interval->sup,false,mul2);
    return ap_scalar_cmp(c.val.interval->inf,c.val.interval->sup)>0;
  default: abort();
	/*******
		TODO: handle arg_assert
		arg_assert(0,return false;);
	********/
  }
}


static inline ap_lincons0_t opt_lincons_of_bound(opt_oct_internal_t* pr,
					     int i, int j,
					     double  d)
{
  ap_linexpr0_t* e;
  if (i==j) {
    /* OPT_ZEROary constraint */
    e = ap_linexpr0_alloc(AP_LINEXPR_SPARSE, 0);
    opt_scalar_of_upper_bound(pr,e->cst.val.scalar,d,true);
  }
  else if (i==(j^1)) {
    /* OPT_UNARY constraint */
    e = ap_linexpr0_alloc(AP_LINEXPR_SPARSE, 1);
    e->p.linterm[0].dim = i/2;
    ap_scalar_set_int(e->p.linterm[0].coeff.val.scalar,(i&1) ? -1 : 1);
    opt_scalar_of_upper_bound(pr,e->cst.val.scalar,d,true);
  }
  else {
    /* OPT_BINARY constraint */
    e = ap_linexpr0_alloc(AP_LINEXPR_SPARSE, 2);
    e->p.linterm[0].dim = j/2;
    e->p.linterm[1].dim = i/2;
    ap_scalar_set_int(e->p.linterm[0].coeff.val.scalar,(j&1) ?  1 : -1);
    ap_scalar_set_int(e->p.linterm[1].coeff.val.scalar,(i&1) ? -1 :  1);
    opt_scalar_of_upper_bound(pr,e->cst.val.scalar,d,false);
  }
  return ap_lincons0_make(AP_CONS_SUPEQ,e,NULL);
}


/* makes an interval from [-minf,sup], with sound approximations
   pr->conv is set if the conversion is not exact
   note: may output an OPT_EMPTY interval
*/
static inline void opt_interval_of_bounds(opt_oct_internal_t* pr,
				      ap_interval_t* i,
				      double minf, double sup,
				      bool div2)
{
  opt_scalar_of_upper_bound(pr,i->sup, sup,div2);
  opt_scalar_of_lower_bound(pr,i->inf,minf,div2);
}

/**************
	Basic operations
*************/

opt_oct_t * opt_oct_alloc_internal(opt_oct_internal_t *pr, int dim, int intdim);
int opt_oct_size(ap_manager_t* man, opt_oct_t* o);
opt_oct_t * opt_oct_alloc_top(opt_oct_internal_t *pr, int dim, int intdim);
void opt_oct_free_internal(opt_oct_internal_t *pr, opt_oct_t *o);
opt_oct_t * opt_oct_copy_internal(opt_oct_internal_t *pr, opt_oct_t *o);
opt_oct_t* opt_oct_set_mat(opt_oct_internal_t* pr, opt_oct_t* o, opt_oct_mat_t* m, opt_oct_mat_t* closed, bool destructive);
opt_oct_t* opt_oct_copy(ap_manager_t* man, opt_oct_t* o);
void opt_oct_free(ap_manager_t* man, opt_oct_t* a);
opt_oct_t* opt_oct_bottom(ap_manager_t* man, int intdim, int realdim);
opt_oct_t* opt_oct_top(ap_manager_t* man, int intdim, int realdim);
ap_dimension_t opt_oct_dimension(ap_manager_t* man, opt_oct_t* o);
void opt_oct_cache_closure(opt_oct_internal_t *pr, opt_oct_t *o);
void opt_oct_close(opt_oct_internal_t *pr, opt_oct_t *o);
opt_oct_t* opt_oct_closure(ap_manager_t *man, bool destructive, opt_oct_t *o);
void opt_oct_internal_free(opt_oct_internal_t *pr);
opt_oct_t* opt_oct_of_abstract0(ap_abstract0_t* a);
ap_abstract0_t* abstract0_of_opt_oct(ap_manager_t* man, opt_oct_t* oct);



/**************
	nary operations
*************/

opt_oct_t* opt_oct_meet(ap_manager_t* man, bool destructive, opt_oct_t* o1, opt_oct_t* o2);
opt_oct_t* opt_oct_join(ap_manager_t* man, bool destructive, opt_oct_t* o1, opt_oct_t* o2);
opt_oct_t* opt_oct_widening(ap_manager_t* man, opt_oct_t* o1, opt_oct_t* o2);
opt_oct_t* opt_oct_add_epsilon(ap_manager_t* man, opt_oct_t* o, ap_scalar_t* epsilon);
opt_oct_t* opt_oct_add_epsilon_bin(ap_manager_t* man, opt_oct_t* o1, opt_oct_t* o2, ap_scalar_t* epsilon);


/**************
	predicate operations
*************/ 

bool opt_oct_is_bottom(ap_manager_t* man, opt_oct_t* o);
bool opt_oct_is_top(ap_manager_t* man, opt_oct_t* o);
bool opt_oct_is_leq(ap_manager_t* man, opt_oct_t* o1, opt_oct_t* o2);
bool opt_oct_is_eq(ap_manager_t* man, opt_oct_t* o1, opt_oct_t* o2);
ap_interval_t** opt_oct_to_box(ap_manager_t* man, opt_oct_t* o);
ap_interval_t* opt_oct_bound_dimension(ap_manager_t* man,opt_oct_t* o, ap_dim_t dim);
ap_lincons0_array_t opt_oct_to_lincons_array(ap_manager_t* man, opt_oct_t* o);
bool opt_oct_sat_lincons_timing(ap_manager_t* man, opt_oct_t* o, ap_lincons0_t* lincons);
bool opt_oct_sat_tcons(ap_manager_t* man, opt_oct_t* o, ap_tcons0_t* cons);


/**************
	resize operations
*************/ 

opt_oct_t* opt_oct_forget_array(ap_manager_t* man,bool destructive, opt_oct_t* o, ap_dim_t* tdim, int size, bool project);
opt_oct_t* opt_oct_add_dimensions(ap_manager_t* man, bool destructive, opt_oct_t* o, ap_dimchange_t* dimchange, bool project);
opt_oct_t* opt_oct_remove_dimensions(ap_manager_t* man, bool destructive, opt_oct_t* o, ap_dimchange_t* dimchange);
opt_oct_t* opt_oct_permute_dimensions(ap_manager_t* man, bool destructive, opt_oct_t* o, ap_dimperm_t* permutation);

/***************
	Transfer Operations
***************/

opt_oct_t* opt_oct_meet_lincons_array(ap_manager_t* man, bool destructive, opt_oct_t* o, ap_lincons0_array_t* array);
opt_oct_t* opt_oct_meet_tcons_array(ap_manager_t* man, bool destructive, opt_oct_t* o, ap_tcons0_array_t* array);
opt_oct_t* opt_oct_assign_linexpr_array(ap_manager_t* man, bool destructive, opt_oct_t* o, ap_dim_t* tdim, ap_linexpr0_t** texpr, size_t size, opt_oct_t* dest);
opt_oct_t* opt_oct_assign_texpr_array(ap_manager_t* man, bool destructive, opt_oct_t* o, ap_dim_t* tdim, ap_texpr0_t** texpr, int size, opt_oct_t* dest);

#ifdef __cplusplus
}
#endif

#endif

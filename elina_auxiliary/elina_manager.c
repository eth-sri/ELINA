/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
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
/* elina_manager.c: global manager passed to all functions */
/* ************************************************************************* */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "elina_manager.h"

const char* elina_name_of_funid[ELINA_FUNID_SIZE2] = {
  "unknown",
  "copy",
  "free",
  "size",
  "minimize",
  "canonicalize",
  "hash",
  "approximate",
  "print",
  "printdiff",
  "dump",
  "serialize_raw",
  "deserialize_raw",
  "bottom",
  "top",
  "of_box",
  "dimension",
  "is_bottom",
  "is_top",
  "is_leq",
  "is_eq",
  "is_dimension_unconstrained",
  "sat_interval",
  "sat_lincons",
  "sat_tcons",
  "bound_dimension",
  "bound_linexpr",
  "bound_texpr",
  "to_box",
  "to_lincons_array",
  "to_tcons_array",
  "to_generator_array",
  "meet",
  "meet_array",
  "meet_lincons_array",
  "meet_tcons_array",
  "join",
  "join_array",
  "add_ray_array",
  "assign_linexpr_array",
  "substitute_linexpr_array",
  "assign_texpr_array",
  "substitute_texpr_array",
  "add_dimensions",
  "remove_dimensions",
  "permute_dimensions",
  "forget_array",
  "expand",
  "fold",
  "widening",
  "closure",
  "unknown",
  "change_environment",
  "rename"
};

const char* elina_name_of_exception[ELINA_EXC_SIZE] = {
  "NONE",
  "TIMEOUT",
  "OUT_OF_SPACE",
  "OVERFLOW",
  "INVALID_ARGUMENT",
  "NOT_IMPLEMENTED"
};

/* ********************************************************************** */
/* I. Constructor and destructor for options */
/* ********************************************************************** */

void elina_funopt_init(elina_funopt_t* opt)
{
  opt->algorithm = 0;
  opt->timeout = 0;
  opt->max_object_size = 0;
  opt->flag_exact_wanted = false;
  opt->flag_best_wanted = false;
}

void elina_option_init(elina_option_t* opt)
{
  elina_funid_t funid;
  elina_exc_t exn;

  for (funid=0; funid<ELINA_FUNID_SIZE; funid++){
    elina_funopt_init(&opt->funopt[funid]);
  }
  for (exn=0; exn<ELINA_EXC_SIZE; exn++){
    opt->abort_if_exception[exn] = true;
  }
  opt->scalar_discr = ELINA_SCALAR_DOUBLE;
}

/* ********************************************************************** */
/* II. Constructor and destructor for result */
/* ********************************************************************** */

elina_exclog_t* elina_exc_cons(elina_exc_t exn,
			 elina_funid_t funid,
			 const char* msg,
			 elina_exclog_t* tail)
{
  elina_exclog_t* head = (elina_exclog_t*)malloc(sizeof(elina_exclog_t));
  head->exn = exn;
  head->funid = funid;
  head->msg = strdup(msg ? msg : "");
  head->tail = tail;
  return head;
}

void elina_exclog_free(elina_exclog_t* head)
{
  elina_exclog_t* p = head;
  while (p!=NULL) {
    elina_exclog_t* tail = p->tail;
    free(p->msg);
    free(p);
    p = tail;
  }
}

void elina_result_add_exception(elina_result_t* result, elina_exc_t exn, elina_funid_t funid, const char* msg)
{
  result->exclog = elina_exc_cons(exn,funid,msg,result->exclog);
  result->exn = exn;
}

void elina_result_init(elina_result_t* result)
{
  result->exclog = NULL;
  result->exn = ELINA_EXC_NONE;
  result->flag_exact = false;
  result->flag_best = false;
}
void elina_result_clear(elina_result_t* result)
{
  elina_exclog_free(result->exclog);
  elina_result_init(result);
}

/* ********************************************************************** */
/* III. Constructor and destructor for manager */
/* ********************************************************************** */

/* Constructor and destructor for manager */
elina_manager_t* elina_manager_alloc(const char* library, const char* version,
			       void* internal,
			       void (*internal_free)(void*))
{
  elina_manager_t* man;

  assert(sizeof(bool)==1);

  man = (elina_manager_t*)malloc(sizeof(elina_manager_t));
  man->library = library;
  man->version = version;
  man->internal = internal;
  man->internal_free = internal_free;
  man->count = 1;
  elina_option_init(&man->option);
  elina_result_init(&man->result);
  return man;
}
void elina_manager_free(elina_manager_t* man)
{
  assert(man->count>=1);
  if (man->count>1){
    man->count--;
  }
  else {
    if (man->internal != NULL){
      man->internal_free(man->internal);
      man->internal = NULL;
    }
    elina_result_clear(&man->result);
    man->count = 0;
    free(man);
  }
}


/* ********************************************************************** */
/* IV. Other User functions */
/* ********************************************************************** */

const char* elina_manager_get_library(elina_manager_t* man)
{ return man->library; }
const char* elina_manager_get_version(elina_manager_t* man)
{ return man->version; }
elina_funopt_t elina_manager_get_funopt(elina_manager_t* man, elina_funid_t funid)
{
  if (funid<ELINA_FUNID_SIZE)
    return man->option.funopt[funid];
  else {
    fprintf(stderr,"elina_manager.c: elina_manager_get_funopt: funid should be less than ELINA_FUNID_SIZE\n");
    abort();
  }
}
bool elina_manager_get_abort_if_exception(elina_manager_t* man, elina_exc_t exn)
{ return man->option.abort_if_exception[exn]; }
bool elina_manager_get_flag_exact(elina_manager_t* man)
{ return man->result.flag_exact; }
bool elina_manager_get_flag_best(elina_manager_t* man)
{ return man->result.flag_best; }


void elina_manager_set_funopt(elina_manager_t* man, elina_funid_t funid, elina_funopt_t* funopt)
{ if (funid<ELINA_FUNID_SIZE) man->option.funopt[funid] = *funopt; }

void elina_manager_set_abort_if_exception(elina_manager_t* man, elina_exc_t exn, bool flag)
{ man->option.abort_if_exception[exn] = flag; }
void elina_manager_clear_exclog(elina_manager_t* man)
{
  elina_exclog_free(man->result.exclog);
  man->result.exclog = NULL;
}

/* ********************************************************************** */
/* V. Other Implementor functions */
/* ********************************************************************** */

void elina_manager_raise_exception(elina_manager_t* man,
				elina_exc_t exn,
				elina_funid_t funid,
				const char* msg)
{
  bool pabort;

  if (exn!=ELINA_EXC_NONE){
    pabort = man->option.abort_if_exception[exn];
    if (pabort){
      fprintf(stderr,"ELINA: Abort because of following exception:\nexception %s in function %s:\n%s\n",
	      elina_name_of_exception[exn], elina_name_of_funid[funid],
	      msg);
      abort();
    }
    else {
      elina_result_add_exception(&man->result,exn,funid,msg);
      man->result.flag_exact = man->result.flag_best = false;
    }
  }
  return;
}


/* ********************************************************************** */
/* V. FPU init */
/* ********************************************************************** */

/* simple run-time test that fpu behaves correctly */
static bool test_fpu(void)
{
  int i;
  long double d = 1., dd;
  /* find the minimal long double, as the fixpoint of x -> x/2 with rounding
     towards +oo;
     the max iteration value should be enough for 128-bit floating-point */
  for (i=0;i<5000000;i++) {
    dd = d;
    d /= 2;
    if (d==dd || d==0.) break;
  }
  /* fails if flush to 0 */
  if (d!=dd) { fprintf(stderr,"test_fpu failed test #1 after %i iterations\n",i); return false; }
  /* fails if long double rounding is not towards +oo */
  if (d*0.25!=dd) { fprintf(stderr,"test_fpu failed test #2\n"); return false; }
  /* fails if double rounding is not towards +oo */
  if ((double)d<dd) { fprintf(stderr,"test_fpu failed test #3\n"); return false; }
  /* fails if float rounding is not towards +oo */
  if ((float)d<dd) { fprintf(stderr,"test_fpu failed test #4\n"); return false; }
  return true;
}

#if defined(__ppc__)
bool elina_fpu_init(void)
{
  __asm volatile ("mtfsfi 7,2");
  return test_fpu();
}

#elif defined(__linux) || defined (__APPLE__)
#include <fenv.h>
bool elina_fpu_init(void)
{
  if (!fesetround(FE_UPWARD)) return test_fpu();
  fprintf(stderr,"could not set fpu rounding mode: fesetround failed\n");
  return false;
}

#elif defined(__FreeBSD__) || defined(sun)
#include <ieeefp.h>
bool elina_fpu_init(void)
{
  fpsetround(FP_RP);
  return test_fpu();
}

#else
bool elina_fpu_init(void)
{
  fprintf(stderr,"could not set fpu rounding mode: platform not supported\n");
  return false;
}

#endif

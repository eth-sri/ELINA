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
/* elina_manager.h: global manager passed to all functions */
/* ************************************************************************* */


#ifndef _ELINA_MANAGER_H_
#define _ELINA_MANAGER_H_

#include <stdlib.h>
#include <stdio.h>

#include "elina_coeff.h"

#ifdef __cplusplus
extern "C" {
#endif


/* ********************************************************************** */
/* I. Types */
/* ********************************************************************** */

/* ====================================================================== */
/* I.O General usage */
/* ====================================================================== */

/* For serialization */
typedef struct elina_membuf_t {
  void* ptr;
  size_t size;
} elina_membuf_t;

/* ====================================================================== */
/* I.1 Identifying functions */
/* ====================================================================== */

typedef enum elina_funid_t {
  ELINA_FUNID_UNKNOWN,
  ELINA_FUNID_COPY,
  ELINA_FUNID_FREE,
  ELINA_FUNID_ASIZE, /* For avoiding name conflict with ELINA_FUNID_SIZE */
  ELINA_FUNID_MINIMIZE,
  ELINA_FUNID_CANONICALIZE,
  ELINA_FUNID_HASH,
  ELINA_FUNID_APPROXIMATE,
  ELINA_FUNID_FPRINT,
  ELINA_FUNID_FPRINTDIFF,
  ELINA_FUNID_FDUMP,
  ELINA_FUNID_SERIALIZE_RAW,
  ELINA_FUNID_DESERIALIZE_RAW,
  ELINA_FUNID_BOTTOM,
  ELINA_FUNID_TOP,
  ELINA_FUNID_OF_BOX,
  ELINA_FUNID_DIMENSION,
  ELINA_FUNID_IS_BOTTOM,
  ELINA_FUNID_IS_TOP,
  ELINA_FUNID_IS_LEQ,
  ELINA_FUNID_IS_EQ,
  ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED,
  ELINA_FUNID_SAT_INTERVAL,
  ELINA_FUNID_SAT_LINCONS,
  ELINA_FUNID_SAT_TCONS,
  ELINA_FUNID_BOUND_DIMENSION,
  ELINA_FUNID_BOUND_LINEXPR,
  ELINA_FUNID_BOUND_TEXPR,
  ELINA_FUNID_TO_BOX,
  ELINA_FUNID_TO_LINCONS_ARRAY,
  ELINA_FUNID_TO_TCONS_ARRAY,
  ELINA_FUNID_TO_GENERATOR_ARRAY,
  ELINA_FUNID_MEET,
  ELINA_FUNID_MEET_ARRAY,
  ELINA_FUNID_MEET_LINCONS_ARRAY,
  ELINA_FUNID_MEET_TCONS_ARRAY,
  ELINA_FUNID_JOIN,
  ELINA_FUNID_JOIN_ARRAY,
  ELINA_FUNID_ADD_RAY_ARRAY,
  ELINA_FUNID_ASSIGN_LINEXPR_ARRAY,
  ELINA_FUNID_SUBSTITUTE_LINEXPR_ARRAY,
  ELINA_FUNID_ASSIGN_TEXPR_ARRAY,
  ELINA_FUNID_SUBSTITUTE_TEXPR_ARRAY,
  ELINA_FUNID_ADD_DIMENSIONS,
  ELINA_FUNID_REMOVE_DIMENSIONS,
  ELINA_FUNID_PERMUTE_DIMENSIONS,
  ELINA_FUNID_FORGET_ARRAY,
  ELINA_FUNID_EXPAND,
  ELINA_FUNID_FOLD,
  ELINA_FUNID_WIDENING,
  ELINA_FUNID_CLOSURE,
  ELINA_FUNID_SIZE,
  ELINA_FUNID_CHANGE_ENVIRONMENT,
  ELINA_FUNID_RENAME_ARRAY,
  ELINA_FUNID_SIZE2
} elina_funid_t;

extern const char* elina_name_of_funid[ELINA_FUNID_SIZE2];
/* give the name of a function identifier */


/* ====================================================================== */
/* I.2 Exceptions */
/* ====================================================================== */

/* Exceptions (public type) */
typedef enum elina_exc_t {
  ELINA_EXC_NONE,             /* no exception detected */
  ELINA_EXC_TIMEOUT,          /* timeout detected */
  ELINA_EXC_OUT_OF_SPACE,     /* out of space detected */
  ELINA_EXC_OVERFLOW,         /* magnitude overflow detected */
  ELINA_EXC_INVALID_ARGUMENT, /* invalid arguments */
  ELINA_EXC_NOT_IMPLEMENTED,  /* not implemented */
  ELINA_EXC_SIZE
} elina_exc_t;

extern const char* elina_name_of_exception[ELINA_EXC_SIZE];

/* Exception log */
typedef struct elina_exclog_t {
  elina_exc_t exn;
  elina_funid_t funid;
  char* msg;                   /* dynamically allocated */
  struct elina_exclog_t* tail;
} elina_exclog_t;

/* Exceptions and other indications (out) (opaque type) */
typedef struct elina_result_t {
  elina_exclog_t* exclog; /* history of exceptions */
  elina_exc_t exn;        /* exception for the last called function */
  bool flag_exact;  /* result is mathematically exact or don't know */
  bool flag_best;   /* result is best correct approximation or don't know */
} elina_result_t;


/* ====================================================================== */
/* I.2 Options */
/* ====================================================================== */

/* Option associated to each function (public type) */
typedef struct elina_funopt_t {
  int algorithm;
  /* Algorithm selection:
     - 0 is default algorithm;
     - MAX_INT is most accurate available;
     - MIN_INT is most efficient available;
     - otherwise, no accuracy or speed meaning
  */
  size_t timeout; /* unit !? */
  /* Above the given computation time, the function may abort with the
     exception flag flag_time_out on.
  */
  size_t max_object_size; /* in abstract object size unit. */
  /* If during the computation, the size of some object reach this limit, the
     function may abort with the exception flag flag_out_of_space on.
  */
  bool flag_exact_wanted;
  /* return information about exactitude if possible
  */
  bool flag_best_wanted;
  /* return information about best correct approximation if possible
  */
} elina_funopt_t;

/* Options (in) (opaque type) */
typedef struct elina_option_t {
  elina_funopt_t funopt[ELINA_FUNID_SIZE];
  bool abort_if_exception[ELINA_EXC_SIZE];
  elina_scalar_discr_t scalar_discr; /* Preferred type for scalars */
} elina_option_t;

/* ====================================================================== */
/* I.3 Manager */
/* ====================================================================== */

/* Manager (opaque type) */
typedef struct elina_manager_t {
  const char* library;                 /* name of the effective library */
  const char* version;                 /* version of the effective library */
  void* internal;                /* library dependent,
				    should be different for each thread
				    (working space) */
  void* funptr[ELINA_FUNID_SIZE];     /* Array of function pointers,
				   initialized by the effective library */
  elina_option_t option;            /* Options (in) */
  elina_result_t result;            /* Exceptions and other indications (out) */
  void (*internal_free)(void*);  /* deallocation function for internal */
  size_t count;                  /* reference counter */
} elina_manager_t;

/* ********************************************************************** */
/* II. User Functions */
/* ********************************************************************** */

void elina_manager_clear_exclog(elina_manager_t* man);
  /* erase the current log of exception */
void elina_manager_free(elina_manager_t* man);
  /* dereference the counter,
     and possibly free internal field if it is not yet put to NULL */

/* Reading fields */
const char* elina_manager_get_library(elina_manager_t* man);
const char* elina_manager_get_version(elina_manager_t* man);

elina_funopt_t elina_manager_get_funopt(elina_manager_t* man, elina_funid_t funid);
bool elina_manager_get_abort_if_exception(elina_manager_t* man, elina_exc_t exn);

elina_exc_t elina_manager_get_exception(elina_manager_t* man);
  /* Get the last exception raised */
elina_exclog_t* elina_manager_get_exclog(elina_manager_t* man);
  /* Get the full log of exception */
bool elina_manager_get_flag_exact(elina_manager_t* man);
bool elina_manager_get_flag_best(elina_manager_t* man);

/* Settings fields */
void elina_funopt_init(elina_funopt_t* fopt);
void elina_manager_set_funopt(elina_manager_t* man, elina_funid_t funid, elina_funopt_t* funopt);
void elina_manager_set_abort_if_exception(elina_manager_t* man, elina_exc_t exn, bool flag);

bool elina_fpu_init(void);
/* tries to set the FPU rounding-mode towards +oo, returns true if successful */


/* ********************************************************************** */
/* III. Implementor Functions */
/* ********************************************************************** */

elina_manager_t* elina_manager_alloc(const char* library, const char* version,
			       void* internal,
			       void (*internal_free)(void*));
static inline
elina_manager_t* elina_manager_copy(elina_manager_t* man);
  /* Increment the reference counter and return its argument */
void elina_manager_raise_exception(elina_manager_t* man,
				elina_exc_t exn, elina_funid_t funid, const char* msg);
  /* raise an exception and put fields
     man->result.flag_exact and man->result.flag_best to
     false
  */
elina_exclog_t* elina_exc_cons(elina_exc_t exn,
			 elina_funid_t funid, const char* msg,
			 elina_exclog_t* tail);
void elina_exclog_free(elina_exclog_t* head);

/* ********************************************************************** */
/* IV. Definition of previously declared inline functions */
/* ********************************************************************** */

static inline
elina_manager_t* elina_manager_copy(elina_manager_t* man)
{ man->count++; return man; }
#ifdef __cplusplus
}
#endif

#endif

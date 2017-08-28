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

/* ************************************************************************* */
/* elina_tcons0.h: tree expressions constraints and arrays */
/* ************************************************************************* */



#ifndef _ELINA_TCONS0_H_
#define _ELINA_TCONS0_H_

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "elina_coeff.h"
#include "elina_texpr0.h"
#include "elina_lincons0.h"

#ifdef __cplusplus
extern "C" {
#endif


/* Datatypes */
/* ====================================================================== */

/* Represents the constraint "expr constyp 0" */
typedef struct elina_tcons0_t {
  elina_texpr0_t* texpr0;  /* expression */
  elina_constyp_t constyp; /* type of constraint */
  elina_scalar_t* scalar;  /* maybe NULL.  
			   For EQMOD constraint, indicates the
			   modulo */
} elina_tcons0_t;

/* Array of constraints */
typedef struct elina_tcons0_array_t {
  elina_tcons0_t* p;
  size_t size;
} elina_tcons0_array_t;

/* ********************************************************************** */
/* I. elina_tcons0_t */
/* ********************************************************************** */

/* ====================================================================== */
/* I.1 Memory management and printing */
/* ====================================================================== */

static inline
elina_tcons0_t elina_tcons0_make(elina_constyp_t constyp, 
			   elina_texpr0_t* texpr,
			   elina_scalar_t* scalar);
  /* Create a constraint of given type with the given expression.
     The expression and the coefficient are not duplicated, just pointed to */

elina_tcons0_t elina_tcons0_make_unsat(void);
  /* Create the constraint -1>=0 */

static inline
elina_tcons0_t elina_tcons0_from_lincons0(elina_lincons0_t* cons);
  /* From linear constraint to comb-like tree expression constraint  */

static inline
elina_tcons0_t elina_tcons0_copy(elina_tcons0_t* cons);
  /* Duplication */

static inline
void elina_tcons0_clear(elina_tcons0_t* cons);
  /* Free the linear expression of the constraint and set pointer to NULL */

void elina_tcons0_fprint(FILE* stream,
		      elina_tcons0_t* cons, char** name_of_dim);
  /* Printing a linear constraint */

/* ====================================================================== */
/* I.2 Tests */
/* ====================================================================== */

static inline
bool elina_tcons0_is_interval_cst(elina_tcons0_t* a);
  /* no-variable, only constant leaves */

static inline
bool elina_tcons0_is_interval_linear(elina_tcons0_t* a);
  /* linear with possibly interval coefficients, no rounding */

static inline
bool elina_tcons0_is_interval_polynomial(elina_tcons0_t* a);
  /* polynomial with possibly interval coefficients, no rounding  */

static inline
bool elina_tcons0_is_interval_polyfrac(elina_tcons0_t* a);
  /* polynomial fraction with possibly interval coefficients, no rounding */

static inline
bool elina_tcons0_is_scalar(elina_tcons0_t* a);
  /* all coefficients are scalar (non-interval) */

/* ====================================================================== */
/* I.3 Change of dimensions and permutations */
/* ====================================================================== */

static inline
void elina_tcons0_add_dimensions_with(elina_tcons0_t* cons,
				   elina_dimchange_t* dimchange);
static inline
elina_tcons0_t elina_tcons0_add_dimensions(elina_tcons0_t* cons,
				     elina_dimchange_t* dimchange);

static inline
void elina_tcons0_remove_dimensions_with(elina_tcons0_t* cons,
				      elina_dimchange_t* dimchange);
static inline
elina_tcons0_t elina_tcons0_remove_dimensions(elina_tcons0_t* cons,
					elina_dimchange_t* dimchange);

static inline
void elina_tcons0_permute_dimensions_with(elina_tcons0_t* cons,
				       elina_dimperm_t* perm);
static inline
elina_tcons0_t elina_tcons0_permute_dimensions(elina_tcons0_t* cons,
					 elina_dimperm_t* perm);

/* ********************************************************************** */
/* II. Array of linear constraints */
/* ********************************************************************** */

elina_tcons0_array_t elina_tcons0_array_make(size_t size);
  /* Allocate an array of size constraints.
     The constraints are initialized with NULL pointers, */

void elina_tcons0_array_resize(elina_tcons0_array_t* array, size_t size);
  /* Resize an array of size constraints.
     New constraints are initialized with NULL pointers,
     Removed constraints with non-NULL pointers are deallocated */

void elina_tcons0_array_clear(elina_tcons0_array_t* array);
  /* Clear the constraints of the array, and then the array itself */

void elina_tcons0_array_fprint(FILE* stream,
			      elina_tcons0_array_t* elina_tcons0_array,
			      char** name_of_dim);
  /* Printing */

bool elina_tcons0_array_is_interval_linear(elina_tcons0_array_t* array);
  /* Are all the expressions involved interval linear ? */

/* ====================================================================== */
/* II.1 Change of dimensions and permutations */
/* ====================================================================== */

void elina_tcons0_array_add_dimensions_with(elina_tcons0_array_t* array,
					 elina_dimchange_t* dimchange);
elina_tcons0_array_t elina_tcons0_array_add_dimensions(elina_tcons0_array_t* array,
						 elina_dimchange_t* dimchange);

void elina_tcons0_array_remove_dimensions_with(elina_tcons0_array_t* array,
					    elina_dimchange_t* dimchange);
elina_tcons0_array_t elina_tcons0_array_remove_dimensions(elina_tcons0_array_t* array,
						    elina_dimchange_t* dimchange);

void elina_tcons0_array_permute_dimensions_with(elina_tcons0_array_t* array,
					     elina_dimperm_t* perm);
elina_tcons0_array_t elina_tcons0_array_permute_dimensions(elina_tcons0_array_t* array,
						     elina_dimperm_t* perm);

/* ********************************************************************** */
/* III. Inline functions definitions */
/* ********************************************************************** */

static inline elina_tcons0_t elina_tcons0_make(elina_constyp_t constyp, elina_texpr0_t* texpr, elina_scalar_t* scalar)
{
  elina_tcons0_t cons;
  cons.constyp = constyp;
  cons.texpr0 = texpr;
  cons.scalar = scalar;
  return cons;
}
static inline elina_tcons0_t elina_tcons0_from_lincons0(elina_lincons0_t* cons)
{
  elina_tcons0_t res;
  res.texpr0 = elina_texpr0_from_linexpr0(cons->linexpr0);
  res.constyp = cons->constyp;
  res.scalar = cons->scalar ? elina_scalar_alloc_set(cons->scalar) : NULL;
  return res;
}
static inline elina_tcons0_t elina_tcons0_copy(elina_tcons0_t* cons)
{
  return elina_tcons0_make(cons->constyp, 
			elina_texpr0_copy(cons->texpr0),
			cons->scalar ? elina_scalar_alloc_set(cons->scalar) : NULL);
}
static inline void elina_tcons0_clear(elina_tcons0_t* tcons)
{
  if (tcons->texpr0){
    elina_texpr0_free(tcons->texpr0);
  }
  tcons->texpr0 = NULL;
  if (tcons->scalar){
    elina_scalar_free(tcons->scalar);
  }
  tcons->scalar = NULL;
}
  
static inline
bool elina_tcons0_is_interval_cst(elina_tcons0_t* a)
{ return elina_texpr0_is_interval_cst(a->texpr0); }
static inline
bool elina_tcons0_is_interval_linear(elina_tcons0_t* a)
{ return elina_texpr0_is_interval_linear(a->texpr0); }
static inline
bool elina_tcons0_is_interval_polynomial(elina_tcons0_t* a)
{ return elina_texpr0_is_interval_polynomial(a->texpr0); }
static inline
bool elina_tcons0_is_interval_polyfrac(elina_tcons0_t* a)
{ return elina_texpr0_is_interval_polyfrac(a->texpr0); }
static inline
bool elina_tcons0_is_scalar(elina_tcons0_t* a)
{ return elina_texpr0_is_scalar(a->texpr0); }

static inline
void elina_tcons0_add_dimensions_with(elina_tcons0_t* cons,
				   elina_dimchange_t* dimchange)
{ elina_texpr0_add_dimensions_with(cons->texpr0,dimchange); }
static inline
elina_tcons0_t elina_tcons0_add_dimensions(elina_tcons0_t* cons,
				     elina_dimchange_t* dimchange)
{
  return elina_tcons0_make(cons->constyp,
			elina_texpr0_add_dimensions(cons->texpr0,dimchange),
			cons->scalar ? elina_scalar_alloc_set(cons->scalar) : NULL);
}
static inline
void elina_tcons0_remove_dimensions_with(elina_tcons0_t* cons,
				      elina_dimchange_t* dimchange)
{ elina_texpr0_remove_dimensions_with(cons->texpr0,dimchange); }
static inline
elina_tcons0_t elina_tcons0_remove_dimensions(elina_tcons0_t* cons,
					elina_dimchange_t* dimchange)
{
  return elina_tcons0_make(cons->constyp,
			elina_texpr0_remove_dimensions(cons->texpr0,dimchange),
			cons->scalar ? elina_scalar_alloc_set(cons->scalar) : NULL);
}
static inline
void elina_tcons0_permute_dimensions_with(elina_tcons0_t* cons,
				       elina_dimperm_t* perm)
{ elina_texpr0_permute_dimensions_with(cons->texpr0,perm); }
static inline
elina_tcons0_t elina_tcons0_permute_dimensions(elina_tcons0_t* cons,
					 elina_dimperm_t* perm)
{
  return elina_tcons0_make(cons->constyp,
			elina_texpr0_permute_dimensions(cons->texpr0,perm),
			cons->scalar ? elina_scalar_alloc_set(cons->scalar) : NULL);
}
  
#ifdef __cplusplus
}
#endif

#endif

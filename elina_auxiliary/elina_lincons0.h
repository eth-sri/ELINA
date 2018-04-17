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

/* ************************************************************************* */
/* elina_lincons0.h: linear constraints and arrays */
/* ************************************************************************* */


#ifndef _ELINA_LINCONS0_H_
#define _ELINA_LINCONS0_H_

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "elina_coeff.h"
#include "elina_linexpr0.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ====================================================================== */
/* Datatypes */
/* ====================================================================== */

/* Datatype for type of constraints */
typedef enum elina_constyp_t {
  ELINA_CONS_EQ,    /* equality constraint */
  ELINA_CONS_SUPEQ, /* >= constraint */
  ELINA_CONS_SUP,   /* > constraint */
  ELINA_CONS_EQMOD, /* congruence equality constraint */
  ELINA_CONS_DISEQ  /* disequality constraint */
} elina_constyp_t;

/* Represents the constraint "expr constyp 0" */
typedef struct elina_lincons0_t {
  elina_linexpr0_t* linexpr0;  /* expression */
  elina_constyp_t constyp;     /* type of constraint */
  elina_scalar_t* scalar;      /* maybe NULL.  

			       For EQMOD constraint, indicates the
			       modulo */
} elina_lincons0_t;

/* Array of constraints */
typedef struct elina_lincons0_array_t {
  elina_lincons0_t* p;
  size_t size;
} elina_lincons0_array_t;

/* ********************************************************************** */
/* I. elina_lincons0_t */
/* ********************************************************************** */

/* ====================================================================== */
/* I.1 Memory management and printing */
/* ====================================================================== */


elina_lincons0_t elina_lincons0_make(elina_constyp_t constyp, 
			       elina_linexpr0_t* linexpr,
			       elina_scalar_t* scalar);
  /* Create a constraint of given type with the given expression.
     The expression and the coefficient are not duplicated, just pointed to */

elina_lincons0_t elina_lincons0_make_unsat(void);
  /* Create the constraint -1>=0 */


elina_lincons0_t elina_lincons0_copy(elina_lincons0_t* cons);
  /* Duplication */


void elina_lincons0_clear(elina_lincons0_t* cons);
  /* Free the linear expression of the constraint and set pointer to NULL */

void elina_lincons0_print(elina_lincons0_t* cons, char** name_of_dim);
void elina_lincons0_fprint(FILE* stream,
			elina_lincons0_t* cons, char** name_of_dim);
  /* Printing a linear constraint */

/* ====================================================================== */
/* I.2 Tests */
/* ====================================================================== */

bool elina_lincons0_is_unsat(elina_lincons0_t* cons);
  /* True if the constraint is b>=0 or [a,b]>=0 with b negative */

bool elina_lincons0_is_sat(elina_lincons0_t* cons);
  /* True if the constraint is trivially satisfiable, e.g. [a,b]>=0 with a positive */

/* ====================================================================== */
/* I.3 Change of dimensions and permutations */
/* ====================================================================== */


void elina_lincons0_add_dimensions_with(elina_lincons0_t* cons,
				     elina_dimchange_t* dimchange);

elina_lincons0_t elina_lincons0_add_dimensions(elina_lincons0_t* cons,
					 elina_dimchange_t* dimchange);


void elina_lincons0_permute_dimensions_with(elina_lincons0_t* cons,
					 elina_dimperm_t* perm);

elina_lincons0_t elina_lincons0_permute_dimensions(elina_lincons0_t* cons,
					     elina_dimperm_t* perm);

/* ********************************************************************** */
/* II. Array of linear constraints */
/* ********************************************************************** */

elina_lincons0_array_t elina_lincons0_array_make(size_t size);
  /* Allocate an array of size constraints.
     The constraints are initialized with NULL pointers, */

void elina_lincons0_array_resize(elina_lincons0_array_t* array, size_t size);
  /* Resize an array of size constraints.
     New constraints are initialized with NULL pointers,
     Removed constraints with non-NULL pointers are deallocated */

void elina_lincons0_array_clear(elina_lincons0_array_t* array);
  /* Clear the constraints of the array, and then the array itself */

void elina_lincons0_array_print(elina_lincons0_array_t* elina_lincons0_array,
			     char** name_of_dim);
void elina_lincons0_array_fprint(FILE* stream,
			      elina_lincons0_array_t* elina_lincons0_array,
			      char** name_of_dim);
  /* Printing */

elina_linexpr_type_t elina_lincons0_array_type(elina_lincons0_array_t* array);
bool elina_lincons0_array_is_linear(elina_lincons0_array_t* array);
bool elina_lincons0_array_is_quasilinear(elina_lincons0_array_t* array);
  /* Are all the expressions involved linear (resp. quasilinear) */

/* ====================================================================== */
/* II.1 Change of dimensions and permutations */
/* ====================================================================== */
void elina_lincons0_array_add_dimensions_with(elina_lincons0_array_t* array,
					   elina_dimchange_t* dimchange);
elina_lincons0_array_t elina_lincons0_array_add_dimensions(elina_lincons0_array_t* array,
						     elina_dimchange_t* dimchange);

void elina_lincons0_array_permute_dimensions_with(elina_lincons0_array_t* array,
					       elina_dimperm_t* perm);
elina_lincons0_array_t elina_lincons0_array_permute_dimensions(elina_lincons0_array_t* array,
							 elina_dimperm_t* perm);

  
#ifdef __cplusplus
}
#endif

#endif

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
/* elina_linexpr0.h: linear expressions */
/* ************************************************************************* */


#ifndef _ELINA_LINEXPR0_H_
#define _ELINA_LINEXPR0_H_

#include <limits.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "elina_coeff.h"
#include "elina_dimension.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ====================================================================== */
/* Datatypes */
/* ====================================================================== */

/* Discriminant for dense or sparse representation */
typedef enum elina_linexpr_discr_t {
  ELINA_LINEXPR_DENSE,
  ELINA_LINEXPR_SPARSE
} elina_linexpr_discr_t;

/* A term, for use in sparse representation */
/* Meant to be an abstract datatype ! */
typedef struct elina_linterm_t {
  elina_dim_t dim;
  elina_coeff_t coeff;
} elina_linterm_t;

/* A linear expression. */
/* Meant to be an abstract datatype ! */
typedef struct elina_linexpr0_t {
  elina_coeff_t cst;             /* constant */
  elina_linexpr_discr_t discr;   /* discriminant for array */
  size_t size;             /* size of the array */
  union {
    elina_coeff_t* coeff;     /* array of coefficients */
    elina_linterm_t* linterm; /* array of linear terms */
  } p;
} elina_linexpr0_t;
/* Important invariant:
   If sparse representation,

   - linear terms are sorted in increasing order wrt their dimension.

   - ELINA_DIM_MAX dimensions are meaningless: they serve as free linterm when a new dimension
     is needed (this avoids to permanently reallocating the array.
     They should be ignored.

*/

/* Comment: we do not inline the array in the structure, because this allows to
   redimension (with realloc) the array in a transparent way for the user. */

/* - An interval linear expression is the more general form.
   - A quasilinear expression is such that the only non-scalar
     coefficient is the constant coefficient.

   - A linear expression contains no non-scalar coefficients
*/
typedef enum elina_linexpr_type_t {
  ELINA_LINEXPR_INTLINEAR,
  ELINA_LINEXPR_QUASILINEAR,
  ELINA_LINEXPR_LINEAR
} elina_linexpr_type_t;

/* ====================================================================== */
/* I. Memory management and printing */
/* ====================================================================== */

elina_linexpr0_t* elina_linexpr0_alloc(elina_linexpr_discr_t lin_discr, size_t size);
  /* Allocates a linear expressions with coefficients by default of type SCALAR
     and DOUBLE. If sparse representation, corresponding new dimensions are
     initialized with ELINA_DIM_MAX. */

void elina_linexpr0_realloc(elina_linexpr0_t* e, size_t size);
  /* Change the dimensions of the array in linexpr0.
     If new coefficients are added, their type is of type SCALAR and DOUBLE.
     If sparse representation, corresponding new dimensions are initialized
     with ELINA_DIM_MAX. */

void elina_linexpr0_minimize(elina_linexpr0_t* e);
  /* Reduce the coefficients (transform intervals into scalars when possible).
     In case of sparse representation, also remove zero coefficients */

void elina_linexpr0_free(elina_linexpr0_t* linexpr);
  /* Free the linear expression */

elina_linexpr0_t* elina_linexpr0_copy(elina_linexpr0_t* a);
  /* Duplication */

void elina_linexpr0_fprint(FILE* stream, elina_linexpr0_t* a, char** name_of_dim);
void elina_linexpr0_print(elina_linexpr0_t* a, char** name_of_dim);
  /* Printing a linear expression */

/* ====================================================================== */
/* II. Tests */
/* ====================================================================== */

bool elina_linexpr0_is_integer(elina_linexpr0_t* a, size_t intdim);
  /* Does the expression depends only on integer variables ? assuming
     that the first intdim dimensions are integer */
bool elina_linexpr0_is_real(elina_linexpr0_t* a, size_t intdim);
  /* Does the expression depends only on real variables ? assuming
     that the first intdim dimensions are integer */

  /* Expression classification */

elina_linexpr_type_t elina_linexpr0_type(elina_linexpr0_t* a);
  /* Return the type of the linear expression */
bool elina_linexpr0_is_linear(elina_linexpr0_t* a);
  /* Return true iff all involved coefficients are scalars */
bool elina_linexpr0_is_quasilinear(elina_linexpr0_t* a);
  /* Return true iff all involved coefficients but the constant are scalars */

elina_linexpr_type_t elina_linexpr0_array_type(elina_linexpr0_t** texpr, size_t size);
bool elina_linexpr0_array_is_linear(elina_linexpr0_t** texpr, size_t size);
bool elina_linexpr0_array_is_quasilinear(elina_linexpr0_t** texpr, size_t size);
  /* Idem for arrays */

/* ====================================================================== */
/* III. Access */
/* ====================================================================== */


size_t elina_linexpr0_size(elina_linexpr0_t* expr);
  /* Get the size of the linear expression */


elina_coeff_t* elina_linexpr0_cstref(elina_linexpr0_t* expr);
  /* Get a reference to the constant. Do not free it. */

elina_coeff_t* elina_linexpr0_coeffref(elina_linexpr0_t* expr, elina_dim_t dim);
  /* Get a reference to the coefficient associated to the dimension.
     Do not free it.
     In case of sparse representation,
     possibly induce the addition of a new linear term.
     Return NULL if:
     In case of dense representation, dim>=expr->size.
     In case of sparse representation, dim==ELINA_DIM_MAX.
 */


void elina_linexpr0_get_cst(elina_coeff_t* coeff, elina_linexpr0_t* expr);
  /* Get the constant and assign it to coeff */

bool elina_linexpr0_get_coeff(elina_coeff_t* coeff, elina_linexpr0_t* expr, elina_dim_t dim);
  /* Get coefficient of dimension dim in the expression and assign it to coeff
     Return true in case elina_linexpr0_coeffref returns NULL */

/* Set the constant of the linear expression */
 void elina_linexpr0_set_cst(elina_linexpr0_t* expr, elina_coeff_t* cst);
 void elina_linexpr0_set_cst_scalar(elina_linexpr0_t* expr, elina_scalar_t* scalar);
 void elina_linexpr0_set_cst_scalar_int(elina_linexpr0_t* expr, int num);
 void elina_linexpr0_set_cst_scalar_frac(elina_linexpr0_t* expr, long int num, unsigned long int den);
 void elina_linexpr0_set_cst_scalar_double(elina_linexpr0_t* expr, double num);
 void elina_linexpr0_set_cst_interval(elina_linexpr0_t* expr, elina_interval_t* itv);
 void elina_linexpr0_set_cst_interval_scalar(elina_linexpr0_t* expr, elina_scalar_t* inf, elina_scalar_t* sup);
 void elina_linexpr0_set_cst_interval_int(elina_linexpr0_t* expr, int inf, int sup);
 void elina_linexpr0_set_cst_interval_frac(elina_linexpr0_t* expr,
							 long int numinf, unsigned long int deninf,
							 long int numsup, unsigned long int densup);
 void elina_linexpr0_set_cst_interval_double(elina_linexpr0_t* expr, double inf, double sup);

/* Set the coefficient of dimension dim in the expression.
   Return true in case elina_linexpr0_coeffref returns NULL */
 bool elina_linexpr0_set_coeff(elina_linexpr0_t* expr, elina_dim_t dim, elina_coeff_t* coeff);
 bool elina_linexpr0_set_coeff_scalar(elina_linexpr0_t* expr, elina_dim_t dim, elina_scalar_t* scalar);
 bool elina_linexpr0_set_coeff_scalar_int(elina_linexpr0_t* expr, elina_dim_t dim, int num);
 bool elina_linexpr0_set_coeff_scalar_frac(elina_linexpr0_t* expr, elina_dim_t dim, long int num, unsigned long int den);
 bool elina_linexpr0_set_coeff_scalar_double(elina_linexpr0_t* expr, elina_dim_t dim, double num);
 bool elina_linexpr0_set_coeff_interval(elina_linexpr0_t* expr, elina_dim_t dim, elina_interval_t* itv);
 bool elina_linexpr0_set_coeff_interval_scalar(elina_linexpr0_t* expr, elina_dim_t dim, elina_scalar_t* inf, elina_scalar_t* sup);
 bool elina_linexpr0_set_coeff_interval_int(elina_linexpr0_t* expr, elina_dim_t dim, int inf, int sup);
 bool elina_linexpr0_set_coeff_interval_frac(elina_linexpr0_t* expr, elina_dim_t dim,
							  long int numinf, unsigned long int deninf,
							  long int numsup, unsigned long int densup);
 bool elina_linexpr0_set_coeff_interval_double(elina_linexpr0_t* expr, elina_dim_t dim, double inf, double sup);

/*
bool elina_linexpr0_set_format_generic(elina_coeff_t* (*get_pcoeff)(char*,va_list*,void*,bool*),
				 void* expr, char* fmt, va_list* ap);

bool elina_linexpr0_set_format(elina_linexpr0_t* expr, char* fmt, ...);
*/

typedef enum elina_coefftag_t {
  ELINA_COEFF,          /* waiting for a coeff_t* object and a dimension */
  ELINA_COEFF_S,        /* waiting for a scalar_t* object and a dimension */
  ELINA_COEFF_S_MPQ,    /* waiting for a mpq_t object and a dimension */
  ELINA_COEFF_S_MPFR,   /* waiting for a mpfr_t object and a dimension */
  ELINA_COEFF_S_INT,    /* waiting for a int object and a dimension */
  ELINA_COEFF_S_FRAC,   /* waiting for 2 int objects and a dimension */
  ELINA_COEFF_S_DOUBLE, /* waiting for a double object and a dimension */
  ELINA_COEFF_I,        /* waiting for a interval_t* object and a dimension */
  ELINA_COEFF_I_SCALAR, /* waiting for 2 scalar_t* objects and a dimension */
  ELINA_COEFF_I_MPQ,    /* waiting for 2 mpq_t objects and a dimension */
  ELINA_COEFF_I_MPFR,   /* waiting for 2 mpfr_t objects and a dimension */
  ELINA_COEFF_I_INT,    /* waiting for 2 int objects and a dimension */
  ELINA_COEFF_I_FRAC,   /* waiting for 4 int objects and a dimension */
  ELINA_COEFF_I_DOUBLE, /* waiting for 2 double objects and a dimension */
  ELINA_CST,            /* waiting for a coeff_t* object */
  ELINA_CST_S,          /* waiting for a scalar_t* object */
  ELINA_CST_S_MPQ,      /* waiting for a mpq_t object */
  ELINA_CST_S_MPFR,     /* waiting for a mpfr_t object */
  ELINA_CST_S_INT,      /* waiting for a int object */
  ELINA_CST_S_FRAC,     /* waiting for 2 int objects */
  ELINA_CST_S_DOUBLE,   /* waiting for a double object */
  ELINA_CST_I,          /* waiting for a interval_t* object */
  ELINA_CST_I_SCALAR,   /* waiting for 2 scalar_t* objects */
  ELINA_CST_I_MPQ,      /* waiting for 2 mpq_t objects */
  ELINA_CST_I_MPFR,     /* waiting for 2 mpfr_t objects */
  ELINA_CST_I_INT,      /* waiting for 2 int objects */
  ELINA_CST_I_FRAC,     /* waiting for 4 int objects */
  ELINA_CST_I_DOUBLE,   /* waiting for 2 double objects */
  ELINA_END
} elina_coefftag_t;

bool elina_linexpr0_set_list_generic(elina_coeff_t* (*get_pcoeff)(void* expr, bool cst, va_list* va),
				  void* expr, va_list* va);

bool elina_linexpr0_set_list(elina_linexpr0_t* expr, ...);

/* Iterator (Macro): use:
   elina_linexpr0_ForeachLinterm(elina_linexpr0_t* e, size_t i, elina_dim_t d, elina_coeff_t* coeff){
     ..
   }
   where
   - e is the inspected expression,
   - i is the internal iterator (of type size_t or int)
   - dim is the dimension of one linear term
   - coeff is a pointer to the corresponding coefficient

   ELINA_DIM_MAX dimensions are filtered out.

*/
#define elina_linexpr0_ForeachLinterm(_p_e_, _p_i_, _p_dim_, _p_elina_coeff) \
  for ((_p_i_)=0; \
       (_p_i_)<(_p_e_)->size ? \
	 ((_p_e_)->discr==ELINA_LINEXPR_DENSE ? \
	  ((_p_dim_) = (_p_i_), \
	   (_p_elina_coeff) = &(_p_e_)->p.coeff[_p_i_], \
	   true) :				\
	  ((_p_dim_) = (_p_e_)->p.linterm[_p_i_].dim, \
	   (_p_elina_coeff) = &(_p_e_)->p.linterm[_p_i_].coeff, \
	   (_p_dim_)!=ELINA_DIM_MAX)) :			   \
	 false; \
       (_p_i_)++)

/* ====================================================================== */
/* IV. Change of dimensions and permutations */
/* ====================================================================== */

/* These two functions add dimensions to the expressions, following the
   semantics of dimchange (see the type definition of dimchange).  */
void elina_linexpr0_add_dimensions_with(elina_linexpr0_t* expr,
				  elina_dimchange_t* dimchange);
elina_linexpr0_t* elina_linexpr0_add_dimensions(elina_linexpr0_t* expr,
					  elina_dimchange_t* dimchange);

/* These two functions apply the given permutation to the dimensions. If dense
   representation, the size of the permutation should be expr->size. If sparse
   representation, the dimensions present in the expression should just be less
   than the size of the permutation. */
void elina_linexpr0_permute_dimensions_with(elina_linexpr0_t* expr,
					 elina_dimperm_t* perm);
elina_linexpr0_t* elina_linexpr0_permute_dimensions(elina_linexpr0_t* expr,
					      elina_dimperm_t* perm);

/* ====================================================================== */
/* V. Hashing, comparison */
/* ====================================================================== */

/* Induces reduction of the coefficients */

long elina_linexpr0_hash(elina_linexpr0_t* expr);
bool elina_linexpr0_equal(elina_linexpr0_t* expr1,
		    elina_linexpr0_t* expr2);

/* Lexicographic ordering, terminating by constant coefficients */
int elina_linexpr0_compare(elina_linexpr0_t* expr1,
		     elina_linexpr0_t* expr2);

/* ====================================================================== */
/* Vb. Array of expressions */
/* ====================================================================== */

/* Free the array of expressions of size size */
elina_linexpr0_t ** elina_linexpr0_array_alloc(size_t size);

void elina_linexpr0_array_free(elina_linexpr0_t** texpr, size_t size);

#ifdef __cplusplus
}
#endif

#endif

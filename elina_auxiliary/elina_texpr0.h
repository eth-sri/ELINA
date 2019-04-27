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
/* elina_texpr0.h: tree expressions */
/* ************************************************************************* */


#include <stdarg.h>

#ifndef _ELINA_TEXPR0_H_
#define _ELINA_TEXPR0_H_

#include "elina_dimension.h"
#include "elina_coeff.h"
#include "elina_linexpr0.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ====================================================================== */
/* Datatypes */
/* ====================================================================== */


/*
  IMPORTANT NOTE
  --------------
  correct use of floating-point ELINA_RTYPE_xxx currently supposes that the
  FPU rounds towards +oo
*/



/* Operators */
typedef enum elina_texpr_op_t {

  /* Binary operators */
  ELINA_TEXPR_ADD, ELINA_TEXPR_SUB, ELINA_TEXPR_MUL, ELINA_TEXPR_DIV,
  ELINA_TEXPR_MOD,  /* either integer or real, no rounding */
  ELINA_TEXPR_POW,

  /* Unary operators */
  ELINA_TEXPR_NEG   /* no rounding */,
  ELINA_TEXPR_CAST, ELINA_TEXPR_SQRT,

} elina_texpr_op_t;

/* Numerical type (destination of rounding) */
typedef enum elina_texpr_rtype_t {
  ELINA_RTYPE_REAL,     /* real (no rounding) */
  ELINA_RTYPE_INT,      /* integer */
  ELINA_RTYPE_SINGLE,   /* IEEE 754 32-bit single precision, e.g.: C's float */
  ELINA_RTYPE_DOUBLE,   /* IEEE 754 64-bit double precision, e.g.: C's double */
  ELINA_RTYPE_EXTENDED, /* non-standard 80-bit double extended, e.g.: Intel's long double */
  ELINA_RTYPE_QUAD,     /* non-standard 128-bit quadruple precision, e.g.: Motorola's long double */
  ELINA_RTYPE_SIZE      /* Not to be used ! */
} elina_texpr_rtype_t;

/* Rounding direction */
typedef enum elina_texpr_rdir_t {
  ELINA_RDIR_NEAREST = GMP_RNDN, /* Nearest */
  ELINA_RDIR_ZERO    = GMP_RNDZ, /* Zero (truncation for integers) */
  ELINA_RDIR_UP      = GMP_RNDU, /* + Infinity */
  ELINA_RDIR_DOWN    = GMP_RNDD, /* - Infinity */
  ELINA_RDIR_RND,    /* All possible mode, non deterministically */
  ELINA_RDIR_SIZE    /* Not to be used ! */
} elina_texpr_rdir_t;

/* Internal (operator) nodes */
typedef struct elina_texpr0_node_t {
  elina_texpr_op_t    op;
  elina_texpr_rtype_t type;
  elina_texpr_rdir_t  dir;
  struct elina_texpr0_t* exprA; /* First operand */
  struct elina_texpr0_t* exprB; /* Second operand (for binary operations) or NULL */
} elina_texpr0_node_t;

/* Node types */
typedef enum elina_texpr_discr_t {
  ELINA_TEXPR_CST, ELINA_TEXPR_DIM, ELINA_TEXPR_NODE
} elina_texpr_discr_t;

typedef struct elina_texpr0_t {
  elina_texpr_discr_t discr;
  union {
    elina_coeff_t cst;
    elina_dim_t dim;
    elina_texpr0_node_t* node;
  } val;
} elina_texpr0_t;

/* ====================================================================== */
/* I. Constructors and Destructors */
/* ====================================================================== */

elina_texpr0_t* elina_texpr0_cst                 (elina_coeff_t* coeff);
elina_texpr0_t* elina_texpr0_cst_scalar          (elina_scalar_t* scalar);
elina_texpr0_t* elina_texpr0_cst_scalar_mpq      (mpq_t mpq);
elina_texpr0_t* elina_texpr0_cst_scalar_mpfr     (mpfr_t mpfr);
elina_texpr0_t* elina_texpr0_cst_scalar_int      (long int num);
elina_texpr0_t* elina_texpr0_cst_scalar_frac     (long int num, unsigned long int den);
elina_texpr0_t* elina_texpr0_cst_scalar_double   (double num);
elina_texpr0_t* elina_texpr0_cst_interval        (elina_interval_t* itv);
elina_texpr0_t* elina_texpr0_cst_interval_scalar (elina_scalar_t* inf, elina_scalar_t* sup);
elina_texpr0_t* elina_texpr0_cst_interval_mpq    (mpq_t inf, mpq_t sup);
elina_texpr0_t* elina_texpr0_cst_interval_mpfr   (mpfr_t inf, mpfr_t sup);
elina_texpr0_t* elina_texpr0_cst_interval_int    (long int inf, long int sup);
elina_texpr0_t* elina_texpr0_cst_interval_frac   (long int numinf, unsigned long int deninf,
					    long int numsup, unsigned long int densup);
elina_texpr0_t* elina_texpr0_cst_interval_double (double inf, double sup);
elina_texpr0_t* elina_texpr0_cst_interval_top    (void);
  /* Create a constant leaf expression */

elina_texpr0_t* elina_texpr0_dim(elina_dim_t dim);
  /* Create a dimension (variable) leaf expression */

elina_texpr0_t* elina_texpr0_unop(elina_texpr_op_t op,
			    elina_texpr0_t* opA,
			    elina_texpr_rtype_t type, elina_texpr_rdir_t dir);
  /* Create unary operator node */

elina_texpr0_t* elina_texpr0_binop(elina_texpr_op_t op,
			     elina_texpr0_t* opA, elina_texpr0_t* opB,
			     elina_texpr_rtype_t type, elina_texpr_rdir_t dir);
  /* Create binary operator node */

elina_texpr0_t* elina_texpr0_copy(elina_texpr0_t* expr);
  /* Recursive (deep) copy */

void elina_texpr0_free(elina_texpr0_t* expr);
  /* Recursive (deep) free */

elina_texpr0_t* elina_texpr0_from_linexpr0(elina_linexpr0_t* e);
  /* From linear expression to comb-like expression tree */


/* ====================================================================== */
/* II. Printing */
/* ====================================================================== */

void elina_texpr0_fprint(FILE* stream, elina_texpr0_t* a, char** name_of_dim);
void elina_texpr0_print(elina_texpr0_t* a, char** name_of_dim);
  /* Prints the expression, name_of_dim can be NULL */


/* ====================================================================== */
/* III. Tests, size */
/* ====================================================================== */

static inline bool elina_texpr_is_unop(elina_texpr_op_t op){
  return (op>=ELINA_TEXPR_NEG && op<=ELINA_TEXPR_SQRT);
}
static inline bool elina_texpr_is_binop(elina_texpr_op_t op){
  return (op<=ELINA_TEXPR_POW);
}
  /* Operator classification */

size_t elina_texpr0_depth(elina_texpr0_t* a);
  /* Returns the depth, in operator nodes */

size_t elina_texpr0_size(elina_texpr0_t* a);
  /* Returns the number of operator nodes */

elina_dim_t elina_texpr0_max_dim(elina_texpr0_t* a);
  /* Returns the maximum elina_dim_t PLUS ONE of all dimensions in expression, and
     0 if no dimension at all. 

     For instance, it returns 3 on the expression x2. */

bool elina_texpr0_has_dim(elina_texpr0_t* a, elina_dim_t d);
   /* Returns true if dimension d appears in the expression */

elina_dim_t* elina_texpr0_dimlist(elina_texpr0_t* a);
  /* Returns an ordered, ELINA_DIM_MAX-terminated array of occurring dimensions;
     caller should free() the array after use
   */


  /* Expression classification */

bool elina_texpr0_is_interval_cst(elina_texpr0_t* a);
  /* no-variable, only constant leaves */

bool elina_texpr0_is_interval_linear(elina_texpr0_t* a);
  /* linear with possibly interval coefficients, no rounding */
bool elina_texpr0_is_interval_polynomial(elina_texpr0_t* a);
  /* polynomial with possibly interval coefficients, no rounding  */
bool elina_texpr0_is_interval_polyfrac(elina_texpr0_t* a);
  /* polynomial fraction with possibly interval coefficients, no rounding */
bool elina_texpr0_is_scalar(elina_texpr0_t* a);
  /* all coefficients are scalar (non-interval) */

bool elina_texpr0_array_is_interval_linear(elina_texpr0_t** texpr, size_t size);
bool elina_texpr0_array_is_interval_polynomial(elina_texpr0_t** texpr, size_t size);
bool elina_texpr0_array_is_interval_polyfrac(elina_texpr0_t** texpr, size_t size);
bool elina_texpr0_array_is_scalar(elina_texpr0_t** texpr, size_t size);
  /* idem for arrays */

/* ====================================================================== */
/* IV. Operations */
/* ====================================================================== */


elina_texpr0_t* elina_texpr0_substitute(elina_texpr0_t* a, elina_dim_t dim, elina_texpr0_t *dst);
void elina_texpr0_substitute_with   (elina_texpr0_t* a, elina_dim_t dim, elina_texpr0_t *dst);
  /* Substitute every occurrence of dimension dim with a copy of dst  */

/* ====================================================================== */
/* V. Change of dimensions and permutations */
/* ====================================================================== */

elina_texpr0_t* elina_texpr0_add_dimensions(elina_texpr0_t* expr,
				      elina_dimchange_t* dimchange);
elina_texpr0_t* elina_texpr0_remove_dimensions(elina_texpr0_t* expr,
					 elina_dimchange_t* dimchange);
elina_texpr0_t* elina_texpr0_permute_dimensions(elina_texpr0_t* expr,
					  elina_dimperm_t* dimperm);
void elina_texpr0_add_dimensions_with(elina_texpr0_t* expr,
				   elina_dimchange_t* dimchange);
void elina_texpr0_remove_dimensions_with(elina_texpr0_t* expr,
					elina_dimchange_t* dimchange);
void elina_texpr0_permute_dimensions_with(elina_texpr0_t* expr,
				       elina_dimperm_t* perm);

/* Note: removed variables are replaced with a top interval */

/* ====================================================================== */
/* VI. Hashing, comparisons */
/* ====================================================================== */

long elina_texpr0_hash(elina_texpr0_t* a);
  /* Recursive hashing */

bool elina_texpr0_equal(elina_texpr0_t* a1, elina_texpr0_t* a2);
  /* Structural (recursive) equality */



/* used internally */
elina_texpr0_t* elina_texpr0_node(elina_texpr_op_t op,
			    elina_texpr_rtype_t type, elina_texpr_rdir_t dir,
			    elina_texpr0_t* opA, elina_texpr0_t* opB);
void elina_texpr0_node_free(elina_texpr0_node_t* node);
void elina_texpr0_clear(elina_texpr0_t* node);


#ifdef __cplusplus
}
#endif

#endif

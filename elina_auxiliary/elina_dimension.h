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

/* ********************************************************************** */
/* elina_dimension.h: dimensions and related operations */
/* ********************************************************************** */


#ifndef _ELINA_DIMENSION_H_
#define _ELINA_DIMENSION_H_

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#ifdef __cplusplus
extern "C" {
#endif


/* ====================================================================== */
/* Datatypes */
/* ====================================================================== */

/* Datatype for dimensions */
typedef unsigned int elina_dim_t;
#define ELINA_DIM_MAX UINT_MAX
/* Used for sparse representations (mean: to be ignored) and also as
   a result when an error occurred */

/* Datatype for specifying the dimensionality of an abstract value */
typedef struct elina_dimension_t {
  size_t intdim;
  size_t realdim;
} elina_dimension_t;

/* Datatype for specifying change of dimension (addition or removal) */
typedef struct elina_dimchange_t {
  elina_dim_t* dim;  /* Assumed to be an array of size intdim+realdim */
  size_t intdim ; /* Number of integer dimensions to add/remove */
  size_t realdim; /* Number of real dimensions to add/remove */
} elina_dimchange_t;

/* The semantics is the following:

- Addition of dimensions:

  dimchange.dim[i]=k means: add one dimension at dimension k and shift the
  already existing dimensions greater than or equal to k one step on the
  right (or increment them).

  if k is equal to the size of the vector, then it means: add a dimension at
  the end.

  Repetition are allowed, and means that one inserts more than one dimensions.

  Example:
  linexpr0_add_dimensions([i0 i1 r0 r1], { [0 1 2 2 4],3,1 }) returns
  [0 i0 0 i1 0 0 r0 r1 0], considered as a vector with 6 integer dimensions
  and 3 real dimensions.

- Removal of dimensions

  dimchange.dim[i]=k means: remove the dimension k and shift the dimensions
  greater than k one step on the left (or decrement them).

  Repetitions are meaningless (and are not correct specification)

  Example:
  linexpr0_remove_dimensions([i0 i1 i2 r0 r1 r2], { [0 2 4],2,1 }) returns
  [i1 r0 r2], considered as a vector with 1 integer dimensions
  and 2 real dimensions.

*/

/* Datatype for specifying double changes of dimensions (combination of
   addition and then removal). Used by level 1 function
   change_environment. */
typedef struct elina_dimchange2_t {
  elina_dimchange_t* add;    /* If not NULL, specifies the adding new dimensions */
  elina_dimchange_t* remove; /* If not NULL, specifies the removal of dimensions */
} elina_dimchange2_t;

/* Datatype for permutations */
typedef struct elina_dimperm_t {
  elina_dim_t* dim;    /* Array assumed to be of size size */
  size_t size;
} elina_dimperm_t;
/* Such an object represent the permutation
   i -> dimperm.p[i] for 0<=i<dimperm.size */

/* ====================================================================== */
/* Functions */
/* ====================================================================== */

/* ---------------------------------------------------------------------- */
/* elina_dimchange_t */
/* ---------------------------------------------------------------------- */

void elina_dimchange_init(elina_dimchange_t* dimchange, size_t intdim, size_t realdim);
  /* Initialize a dimchange structure (allocate internal array) */
elina_dimchange_t* elina_dimchange_alloc(size_t intdim, size_t realdim);
  /* Allocate and initialize a dimchange structure */

void elina_dimchange_clear(elina_dimchange_t* dimchange);
  /* Clear a dimchange structure (deallocate internal array) */
void elina_dimchange_free(elina_dimchange_t* dimchange);
  /* Deallocate and clear a dimchange structure */

void elina_dimchange_fprint(FILE* stream, elina_dimchange_t* dimchange);
  /* Printing */
void elina_dimchange_add_invert(elina_dimchange_t* dimchange);
  /* Assuming that dimchange is a transformation for add_dimensions,
     invert it to obtain the inverse transformation using remove_dimensions */

/* ---------------------------------------------------------------------- */
/* elina_dimchange2_t */
/* ---------------------------------------------------------------------- */

void elina_dimchange2_init(elina_dimchange2_t* dimchange2,
				      elina_dimchange_t* add, 
				      elina_dimchange_t* remove);
  /* Initialize a dimchange2 structure by filling its fields with
     arguments */
elina_dimchange2_t* elina_dimchange2_alloc(elina_dimchange_t* add, 
						   elina_dimchange_t* remove);
  /* Allocate and initialize a dimchange2 structure */

void elina_dimchange2_clear(elina_dimchange2_t* dimchange2);
  /* Clear a dimchange structure (deallocate its fields) */
void elina_dimchange2_free(elina_dimchange2_t* dimchange2);
  /* Deallocate and clear a dimchange2 structure */

void elina_dimchange2_fprint(FILE* stream, elina_dimchange2_t* dimchange2);
  /* Printing */

/* ---------------------------------------------------------------------- */
/* elina_dimperm_t */
/* ---------------------------------------------------------------------- */

void elina_dimperm_init(elina_dimperm_t* dimperm, size_t size);
  /* Initialize a dimperm structure (allocate internal array) */
elina_dimperm_t* elina_dimperm_alloc(size_t size);
  /* Allocate and initialize a dimperm structure */

void elina_dimperm_clear(elina_dimperm_t* dimperm);
  /* Clear a dimperm structure (deallocate internal arrau) */
void elina_dimperm_free(elina_dimperm_t* dimperm);
  /* Deallocate and clear a dimchange structure */

void elina_dimperm_fprint(FILE* stream, elina_dimperm_t* perm);
  /* Print a permutation under the form:
     dimperm: size=...
     0 -> perm->dim[0]
     1 -> perm->dim[1]
     ...
 */

void elina_dimperm_set_id(elina_dimperm_t* perm);
  /* Generate the identity permutation */

void elina_dimperm_compose(elina_dimperm_t* perm,
			elina_dimperm_t* perm1, elina_dimperm_t* perm2);
  /* Compose the 2 permutations perm1 and perm2 (in this order)
     and store the result the already allocated perm.
     The sizes of permutations are supposed to be equal.
     At exit, we have perm.dim[i] = perm2.dim[perm1.dim[i]]
  */
void elina_dimperm_invert(elina_dimperm_t* nperm, elina_dimperm_t* perm);
  /* Invert the permutation perm and store it in the already allocated nperm.
     The sizes of permutations are supposed to be equal.
  */


#ifdef __cplusplus
}
#endif

#endif

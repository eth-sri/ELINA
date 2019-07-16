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

/* ********************************************************************** */
/* opt_pk_constructor.h: constructors and accessors */
/* ********************************************************************** */



#ifndef _OPT_PK_CONSTRUCTOR_H_
#define _OPT_PK_CONSTRUCTOR_H_

#include "opt_pk_config.h"
#include "opt_pk_vector.h"

#include "opt_pk_matrix.h"
#include "opt_pk.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Fill the first (opk->dec-1) rows of the matrix with the constraints of the
   universe polyhedron */
void opt_matrix_fill_constraint_top(opt_pk_internal_t* opk, opt_matrix_t* oc, size_t start);


void opt_poly_set_bottom(opt_pk_internal_t* opk, opt_pk_array_t* op);
void opt_poly_set_top(opt_pk_internal_t* opk, opt_pk_array_t* op);

#ifdef __cplusplus
}
#endif

#endif

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


/* ********************************************************************** */
/* opt_pk_satmat.h: operations on saturation matrices */
/* ********************************************************************** */

#ifndef __OPT_PK_SATMAT_H__
#define __OPT_PK_SATMAT_H__

#include <stdlib.h>
#include "opt_pk_bit.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct opt_satmat_t {
  /* public part */
  opt_bitstring_t** p;
  size_t nbrows;
  size_t nbcolumns;
  /* private part */
  size_t  _maxrows;   /* number of rows allocated */
} opt_satmat_t;

opt_satmat_t* opt_satmat_alloc(size_t nbrows, size_t nbcols);
void opt_satmat_resize_rows(opt_satmat_t* os, size_t nbrows);
void opt_satmat_resize_cols(opt_satmat_t* os, size_t nbcols);
opt_satmat_t* opt_satmat_copy_resize_cols(opt_satmat_t* os, size_t nbcols);
void opt_satmat_free(opt_satmat_t* os);
void opt_satmat_clear(opt_satmat_t* os);
opt_satmat_t* opt_satmat_copy(opt_satmat_t* os);
void opt_satmat_print(opt_satmat_t* os);
void opt_satmat_fprint(FILE* stream, opt_satmat_t* os);

opt_bitstring_t opt_satmat_get(opt_satmat_t* os, size_t i, opt_bitindex_t ojx);
void opt_satmat_set(opt_satmat_t* os, size_t i, opt_bitindex_t ojx);
void opt_satmat_clr(opt_satmat_t* os, size_t i, opt_bitindex_t ojx);

opt_satmat_t* opt_satmat_transpose(opt_satmat_t* os, size_t nbcols);

void opt_satmat_exch_rows(opt_satmat_t* os, size_t l1, size_t l2);
void opt_satmat_exch_cols(opt_satmat_t* os, size_t l1, size_t l2);
void opt_satmat_move_rows(opt_satmat_t* os, size_t destrow, size_t orgrow, size_t size);

#ifdef __cplusplus
}
#endif

#endif

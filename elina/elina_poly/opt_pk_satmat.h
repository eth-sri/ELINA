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

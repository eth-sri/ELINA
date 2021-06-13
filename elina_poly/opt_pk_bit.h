/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
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

#ifndef __OPT_PK_BIT_H_
#define __OPT_PK_BIT_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned int opt_bitstring_t;

typedef struct opt_bitindex_t {
  size_t index;
  size_t word;
  opt_bitstring_t bit;
} opt_bitindex_t;

#define opt_bitstring_size (sizeof(opt_bitstring_t)*8)
#define opt_bitstring_msb (1U<<(opt_bitstring_size-1))

/* Operations on \verb-bitindex_t- */
void opt_bitindex_print(opt_bitindex_t* bi);
void opt_bitindex_fprint(FILE* stream, opt_bitindex_t* bi);
opt_bitindex_t opt_bitindex_init(size_t col);
void opt_bitindex_inc(opt_bitindex_t*);
void opt_bitindex_dec(opt_bitindex_t*);
size_t opt_bitindex_size(size_t n);

/* Operations on \verb-bitstring_t- */
opt_bitstring_t* opt_bitstring_alloc(size_t n);
opt_bitstring_t* opt_bitstring_realloc(opt_bitstring_t* b, size_t n);
void opt_bitstring_free(opt_bitstring_t* b);
void opt_bitstring_clear(opt_bitstring_t* b, size_t size);
void opt_bitstring_copy(opt_bitstring_t* b2, opt_bitstring_t* b1, size_t size);
int opt_bitstring_cmp(opt_bitstring_t* r1, opt_bitstring_t* r2, size_t size);

void opt_bitstring_print(opt_bitstring_t* b, size_t size);
void opt_bitstring_fprint(FILE* stream, opt_bitstring_t* b, size_t size);

int opt_bitstring_get(opt_bitstring_t* b, opt_bitindex_t ix);
void opt_bitstring_set(opt_bitstring_t* b, opt_bitindex_t ix);
void opt_bitstring_clr(opt_bitstring_t* b, opt_bitindex_t ix);

void opt_bitstring_move(opt_bitstring_t *dst,opt_bitstring_t *src, size_t size, size_t start);

#ifdef __cplusplus
}
#endif

#endif

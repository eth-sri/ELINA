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

/* ********************************************************************** */
/* opt_pk_bit.c: operations on bitstrings */
/* ********************************************************************** */


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "opt_pk_bit.h"

/* ********************************************************************** */
/* I. Bitindices */
/* ********************************************************************** */

void opt_bitindex_fprint(FILE* stream, opt_bitindex_t* obi)
{
  int k;
  opt_bitstring_t m = obi->bit;
  assert (m!=0);
  k=(-1);
  do {
    k++;
    m >>= 1;
  } while (m!=0);
  fprintf(stream,"index=%lu, word=%lu, bit=%d\n",
	  (unsigned long)obi->index,(unsigned long)obi->word,k);
}

void opt_bitindex_print(opt_bitindex_t* obi)
{
  opt_bitindex_fprint(stdout,obi);
}


/* opt_bitindex_init()- takes as parameter a flat index of a bit and
   returns the corresponding structured index.  opt_bitindex_inc() and
   opt_bitindex_dec()- allow to increment and decrement an index.
   opt_bitindex_size(n)- returns the size of an array of opt_bitstring_t-
   containing n- bits. */

opt_bitindex_t opt_bitindex_init(const size_t col)
{
  opt_bitindex_t res;
  res.index = col;
  res.word = col / opt_bitstring_size;
  res.bit = opt_bitstring_msb >> (col % opt_bitstring_size);
  return res;
}

void opt_bitindex_inc(opt_bitindex_t* const obi){
  obi->index++;
  obi->bit >>= 1;
  if (obi->bit==0){
    obi->bit = opt_bitstring_msb;
    obi->word++;
  }
}


void opt_bitindex_dec(opt_bitindex_t* const obi){
  obi->index--;
  if (obi->bit != opt_bitstring_msb){
     obi->bit <<= 1;
  }
  else {
    obi->bit = 1;
    obi->word--;
  }
}

size_t opt_bitindex_size(const size_t n){
  size_t size = n / opt_bitstring_size;
  if (n % opt_bitstring_size) size++;
  return size;
}

/* ********************************************************************** */
/* II. Bitstrings */
/* ********************************************************************** */

/*
  opt_bitstring_alloc- allocates a new bitstring and
  opt_bitstring_free()- frees the bitstring.

  opt_bitstring_clear- sets to 0- the bits, opt_bitstring_cmp-
  compares two bitfields; be careful, it takes also in account unused bits of
  the last word. Last, opt_bitstring_print()- writes the bits of a
  bitstring. 
*/
  
opt_bitstring_t* opt_bitstring_alloc(size_t n){
  return (opt_bitstring_t*)malloc(n*sizeof(opt_bitstring_t));
}

opt_bitstring_t* opt_bitstring_realloc(opt_bitstring_t* ob, size_t n){
  return (opt_bitstring_t*)realloc(ob, n*sizeof(opt_bitstring_t));
}

void opt_bitstring_free(opt_bitstring_t* ob){
  free(ob);
}

void opt_bitstring_clear(opt_bitstring_t* ob, size_t size){
  size_t i;
  for (i=0; i<size; i++) ob[i]=0;
}

void opt_bitstring_copy(opt_bitstring_t* ob2, opt_bitstring_t* ob1, size_t size){
  size_t i;
  for (i=0; i<size; i++) ob2[i]= ob1[i];
}

void opt_bitstring_fprint(FILE* stream, opt_bitstring_t* ob, size_t size)
{
  size_t j,k;
  opt_bitstring_t m;

  for (j=0; j<size; j++){
    m = opt_bitstring_msb; k = 1;
    while (m!=0) {
      if (ob[j] & m) fprintf(stream,"1"); else fprintf(stream,"0");
      if (k % 8 == 0) fprintf(stream," ");
      else if (k % 4 == 0) fprintf(stream,",");
      m >>= 1; k++;
    }
  }
}

void opt_bitstring_print(opt_bitstring_t* ob, size_t size)
{
  opt_bitstring_fprint(stdout,ob,size);
}

int opt_bitstring_cmp(opt_bitstring_t* const ob1, opt_bitstring_t* const ob2, size_t size){
  size_t i;
  int res=0;
  for (i=0; i<size; i++){
    if (ob1[i] < ob2[i]){ res=-1; break; }
    else if (ob1[i] > ob2[i]){ res=1; break; }
  }
  return res;
}

/* These functions allow to read, set or clear individual bits of a bitstring, 
   referenced by a bitindex. */

int opt_bitstring_get(opt_bitstring_t* const ob, opt_bitindex_t oix) { 
  return ob[oix.word] & oix.bit; 
}

void opt_bitstring_set(opt_bitstring_t* const ob, opt_bitindex_t oix){
  ob[oix.word] |= oix.bit; 
}

void opt_bitstring_clr(opt_bitstring_t* const ob, opt_bitindex_t oix){ 
  ob[oix.word] &= ~oix.bit; 
}

void opt_bitstring_move(opt_bitstring_t *dst,opt_bitstring_t *src, size_t size, size_t start){
	opt_bitindex_t ind = opt_bitindex_init(0);
	opt_bitindex_t counter = opt_bitindex_init(start);
	while(ind.index < size){
		if(opt_bitstring_get(src,ind)){
			opt_bitstring_set(dst,counter);
		}
		opt_bitindex_inc(&ind);
		opt_bitindex_inc(&counter);
	}
}



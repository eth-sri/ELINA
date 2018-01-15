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


#ifndef _ZONOTOPE_REPRESENTATION_H_
#define _ZONOTOPE_REPRESENTATION_H_

#ifdef __cplusplus
extern "C" {
#endif



zonotope_t* zonotope_alloc(elina_manager_t* man, size_t intdim, size_t realdim);

/* Return a copy of an abstract value, on
 * which destructive update does not affect the initial value. */
zonotope_t* zonotope_copy(elina_manager_t* man, zonotope_t* a);

/* free all the memory used by abstract value */
void zonotope_free(elina_manager_t* man, zonotope_t* a);

size_t zonotope_size(elina_manager_t* man, zonotope_t* a);

/* ********************************************************************** */
/* 2. Control of internal representation */
/* ********************************************************************** */
void zonotope_minimize(elina_manager_t* man, zonotope_t* a);

void zonotope_canonicalize(elina_manager_t* man, zonotope_t* a);

int zonotope_hash(elina_manager_t* man, zonotope_t* a);

void zonotope_approximate(elina_manager_t* man, zonotope_t* a, int algorithm);

/* ********************************************************************** */
/* 3. Printing */
/* ********************************************************************** */
void zonotope_fprint(FILE* stream,
		elina_manager_t* man,
		zonotope_t* a,
		char** name_of_dim);

void zonotope_fprintdiff(FILE* stream,
		elina_manager_t* man,
		zonotope_t* a1, zonotope_t* a2,
		char** name_of_dim);

//void zonotope_fdump(FILE* stream, elina_manager_t* man, zonotope_t* a);


#ifdef __cplusplus
}
#endif

#endif

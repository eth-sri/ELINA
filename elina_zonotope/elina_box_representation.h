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

#ifndef _ELINA_BOX_REPRESENTATION_H_
#define _ELINA_BOX_REPRESENTATION_H_

#include "elina_box_internal.h"
#include "elina_box_representation.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Internal functions */
elina_box_t* elina_box_alloc(size_t intdim, size_t realdim);
void elina_box_init(elina_box_t* a);
void elina_box_set_bottom(elina_box_t* a);
void elina_box_set_top(elina_box_t* a);
void elina_box_set(elina_box_t* a, elina_box_t* b);

/* 1. Memory */
elina_box_t* elina_box_copy(elina_manager_t* man, elina_box_t* a);
void elina_box_free(elina_manager_t* man, elina_box_t* a);


/* 3. Printing */
void elina_box_fprint(FILE* stream,
	       elina_manager_t* man,elina_box_t* a,char** name_of_dim);

#ifdef __cplusplus
}
#endif

#endif

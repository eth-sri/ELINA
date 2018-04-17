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

#ifndef __OPT_PK_WIDENING_H
#define __OPT_PK_WIDENING_H

#include "opt_pk_config.h"
#include "opt_pk_internal.h"
#include "opt_pk.h"

#ifdef __cplusplus
extern "C" {
#endif

bool is_vectors_equal_comp_list(opt_pk_internal_t *opk, opt_numint_t * v1, 
				opt_numint_t * v2, unsigned short int * ca1, 
				unsigned short int * ca2, unsigned short int comp_size1, 
				unsigned short int comp_size2);

void vector_copy_comp_list(opt_pk_internal_t *opk, opt_numint_t * dst, opt_numint_t * src, 
			   unsigned short int * ind_map, unsigned short int comp_size);

#ifdef __cplusplus
}
#endif

#endif

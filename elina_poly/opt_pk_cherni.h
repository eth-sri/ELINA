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

#ifndef __OPT_PK_CHERNI_H_
#define __OPT_PK_CHERNI_H_

#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_satmat.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"

#ifdef __cplusplus
extern "C"{
#endif

/* ********************************************************************** */
/* Conversion algorithm */
/* ********************************************************************** */

size_t opt_cherni_conversion(opt_pk_internal_t* opk,
			 opt_matrix_t* con, size_t start,
			 opt_matrix_t* ray, opt_satmat_t* osc, size_t nbline);


int  opt_cherni_simplify(opt_pk_internal_t* opk,
		     opt_matrix_t* con, opt_matrix_t* ray,
		     opt_satmat_t* satf, const size_t nbline);

void opt_cherni_minimize(opt_pk_internal_t* opk,
		     bool con_to_ray,
		     opt_pk_t* op);


void opt_cherni_add_and_minimize(opt_pk_internal_t* opk, 
			     bool con_to_ray,
			     opt_pk_t* op,
			     size_t start);


#ifdef __cplusplus
}
#endif

#endif



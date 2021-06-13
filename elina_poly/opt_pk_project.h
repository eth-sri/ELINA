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


#ifndef __OPT_PK_PROJECT_H__
#define __OPT_PK_PROJECT_H__

#include "opt_pk_config.h"
#include "opt_pk_internal.h"
#include "opt_pk_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

void opt_poly_projectforget_array(bool project,
			      elina_manager_t* man,	
			      opt_pk_t* op, opt_pk_t* oa, 
			      elina_dim_t* tdim, size_t size, bool destructive);

opt_matrix_t * extreme_projection(opt_pk_internal_t *opk, 
				  opt_matrix_t *oc, 
				  elina_dim_t * tdim, size_t size, size_t proj);


#ifdef __cplusplus
}
#endif

#endif

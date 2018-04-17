/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
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


#ifndef __OPT_ZONES_INCR_CLOSURE_H__
#define __OPT_ZONES_INCR_CLOSURE_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "opt_zones.h"
#include "opt_mat.h"
#include "opt_zones_internal.h"
#include "opt_zones_closure.h"

bool incr_closure_dense_scalar(opt_zones_mat_t *oz, unsigned short int dim, unsigned short int v);

bool incr_closure_dense(opt_zones_mat_t *oz, unsigned short int dim, unsigned short int v);

bool incr_closure_comp_sparse(opt_zones_mat_t *oz, unsigned short int dim, unsigned short int v);

#ifdef __cplusplus
}
#endif

#endif

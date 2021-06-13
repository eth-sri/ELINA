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


#ifndef _ZONOML_FUN_H_
#define _ZONOML_FUN_H_

#include "zonoml.h"
#include "zonoml_internal.h"
#include "zonoml_reduced_product.h"


#ifdef __cplusplus
extern "C" {
#endif


zonotope_aff_t * zonotope_aff_from_dense_weights_bias(zonotope_internal_t* pr, double * weights, double bias, size_t offset, size_t size, zonotope_t *z);

zonotope_aff_t* zonotope_aff_mul_weight(zonotope_internal_t* pr, zonotope_aff_t* src, double lambda);

zonotope_aff_t * zonotope_aff_from_sparse_weights_bias(zonotope_internal_t* pr, double * weights, double bias, elina_dim_t *dim, size_t size, zonotope_t *z);

#ifdef __cplusplus
}
#endif

#endif

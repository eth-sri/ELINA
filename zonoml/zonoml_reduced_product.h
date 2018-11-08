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


#ifndef _ZONOML_REDUCED_PRODUCT_H_
#define _ZONOML_REDUCED_PRODUCT_H_

#include "zonoml.h"
#include "zonoml_internal.h"
#include "zonoml_fun.h"


#ifdef __cplusplus
extern "C" {
#endif

elina_abstract0_t *relu_zono(elina_manager_t* man, bool destructive, elina_abstract0_t * abs, elina_dim_t x);

elina_abstract0_t *relu_zono_refined(elina_manager_t* man, bool destructive, elina_abstract0_t * abs,  elina_dim_t x, double new_inf, double new_sup);

elina_abstract0_t * sigmoid_zono(elina_manager_t *man, bool destructive,
				 elina_abstract0_t *abs, elina_dim_t x);

elina_abstract0_t * tanh_zono(elina_manager_t *man, bool destructive,
				 elina_abstract0_t *abs, elina_dim_t x);

elina_abstract0_t* maxpool_zono(elina_manager_t *man, bool destructive, elina_abstract0_t *abs, 
			   size_t *pool_size, size_t *input_size, size_t src_offset, size_t* strides, 
			   size_t dimensionality, size_t dst_offset, bool is_valid_padding);
//void reduced_product_zono_to_oct(elina_manager_t* man, zonoml_t *zo, elina_dim_t *tdim,elina_linexpr0_t ** lexpr_arr, size_t size, elina_abstract0_t* dest);

//void reduced_product_oct_to_zono(elina_manager_t *man, zonoml_t *zo, elina_dim_t y, elina_dim_t x);

#ifdef __cplusplus
}
#endif

#endif

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

#ifndef _ELINA_BOX_RESIZE_H_
#define _ELINA_BOX_RESIZE_H_

#include "elina_box.h"

#ifdef __cplusplus
extern "C" {
#endif

elina_box_t* elina_box_forget_array(elina_manager_t* man,
                                    bool destructive,
                                    elina_box_t* a,
                                    elina_dim_t* tdim,
                                    size_t size,
                                    bool project);
    
elina_box_t* elina_box_add_dimensions(elina_manager_t* man,
			  bool destructive, elina_box_t* a,
			  elina_dimchange_t* dimchange,
			  bool project);

elina_box_t* elina_box_remove_dimensions(elina_manager_t* man,
			     bool destructive, elina_box_t* a,
			     elina_dimchange_t* dimchange);

elina_box_t* elina_box_permute_dimensions(elina_manager_t* man,
                                  bool destructive,
                                  elina_box_t* a,
                                  elina_dimperm_t* permutation);

#ifdef __cplusplus
}
#endif

#endif

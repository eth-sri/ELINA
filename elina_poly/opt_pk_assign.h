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


/* ********************************************************************** */
/* opt_pk_assign.h: Assignements and Substitutions */
/* ********************************************************************** */



#ifndef _OPT_PK_ASSIGN_H_
#define _OPT_PK_ASSIGN_H_

#include "opt_pk_config.h"
#include "opt_pk.h"

#ifdef __cplusplus
extern "C" {
#endif

opt_pk_array_t* opt_poly_asssub_linexpr_det(bool assign, elina_manager_t* man,
			      bool destructive,
			      opt_pk_array_t* oa,
			      elina_dim_t dim, elina_linexpr0_t* linexpr);

opt_pk_array_t* opt_poly_asssub_linexpr_array_det(bool assign, elina_manager_t* man,
				    bool destructive,
				    opt_pk_array_t* pa,
				    elina_dim_t* tdim, elina_linexpr0_t** texpr, 
				    size_t size);

#ifdef __cplusplus
}
#endif

#endif

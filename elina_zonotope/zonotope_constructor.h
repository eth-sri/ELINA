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

#ifndef _ZONOTOPE_CONSTRUCTOR_H_
#define _ZONOTOPE_CONSTRUCTOR_H_



#ifdef __cplusplus
extern "C" {
#endif

/****************/
/* Constructors */
/****************/
/* 1.Basic constructors */
zonotope_t* zonotope_bottom(elina_manager_t* man, size_t intdim, size_t realdim);
zonotope_t* zonotope_top(elina_manager_t* man, size_t intdim, size_t realdim);
zonotope_t* zonotope_of_box(elina_manager_t* man, size_t intdim, size_t realdim, elina_interval_t** tinterval);

/* 2.Accessors */
elina_dimension_t zonotope_dimension(elina_manager_t* man, zonotope_t* a);

/* 3.Tests */
bool zonotope_is_bottom(elina_manager_t* man, zonotope_t* z);
bool zonotope_is_top(elina_manager_t* man, zonotope_t* z);

//tbool_t zonotope_is_leq(elina_manager_t* man, zonotope_t* a, zonotope_t* b);
bool zonotope_is_eq(elina_manager_t* man, zonotope_t* z1, zonotope_t* z2);
//tbool_t zonotope_is_dimension_unconstrained(elina_manager_t* man, zonotope_t* a, elina_dim_t dim);
//tbool_t zonotope_sat_lincons(elina_manager_t* man, zonotope_t* a, elina_lincons0_t* lincons);
//tbool_t zonotope_sat_interval(elina_manager_t* man, zonotope_t* a, elina_interval_t* interval);
//tbool_t zonotope_sat_tcons(elina_manager_t* man, zonotope_t* a, elina_tcons0_t* tcons);

/* 4.Extraction of properties */
elina_interval_t* zonotope_bound_texpr(elina_manager_t* man, zonotope_t* a, elina_texpr0_t* expr);
elina_interval_t* zonotope_bound_dimension(elina_manager_t* man, zonotope_t* a, elina_dim_t dim);
elina_interval_t* zonotope_bound_linexpr(elina_manager_t* man, zonotope_t* a, elina_linexpr0_t* expr);

elina_interval_t** zonotope_to_box(elina_manager_t* man, zonotope_t* z);
elina_tcons0_array_t zonotope_to_tcons_array(elina_manager_t* man, zonotope_t* a);
elina_lincons0_array_t zonotope_to_lincons_array(elina_manager_t* man, zonotope_t* a);
//elina_generator0_array_t zonotope_to_generator_array(elina_manager_t* man, zonotope_t* a);

#ifdef __cplusplus
}
#endif

#endif

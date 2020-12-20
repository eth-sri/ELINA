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



#ifndef __COMPUTE_BOUNDS_H_INCLUDED__
#define __COMPUTE_BOUNDS_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#include "backsubstitute.h"

expr_t * replace_input_poly_cons_in_lexpr(fppoly_internal_t *pr, expr_t * expr, fppoly_t * fp);

expr_t * replace_input_poly_cons_in_uexpr(fppoly_internal_t *pr, expr_t * expr, fppoly_t * fp);

double compute_lb_from_expr(fppoly_internal_t *pr, expr_t * expr, fppoly_t * fp, int layerno);

double compute_ub_from_expr(fppoly_internal_t *pr, expr_t * expr, fppoly_t * fp, int layerno);

double get_lb_using_previous_layers(elina_manager_t *man, fppoly_t *fp, expr_t *expr, size_t layerno);

double get_ub_using_previous_layers(elina_manager_t *man, fppoly_t *fp, expr_t *expr, size_t layerno);

elina_linexpr0_t *get_output_lexpr_defined_over_previous_layers(elina_manager_t *man, elina_abstract0_t *element, size_t neuron_no, int prev_layer);

elina_linexpr0_t *get_output_uexpr_defined_over_previous_layers(elina_manager_t *man, elina_abstract0_t *element, size_t neuron_no, int prev_layer);

#ifdef __cplusplus
 }
#endif

#endif

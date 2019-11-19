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



#ifndef __S_CURVE_APPROX_H_INCLUDED__
#define __S_CURVE_APPROX_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#include "backsubstitute.h"

expr_t * uexpr_replace_sigmoid_bounds(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons);

expr_t * lexpr_replace_sigmoid_bounds(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons);

expr_t * lexpr_replace_tanh_bounds(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons);

expr_t * uexpr_replace_tanh_bounds(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons);

void compute_slope_and_intercept_s_curve_lexpr(fppoly_internal_t * pr, double * slope_inf, double *slope_sup, 
						double *intercept_inf, double *intercept_sup, double inf_coeff, 
						double sup_coeff, double lb, double ub, bool is_sigmoid, bool *boxify);

void compute_slope_and_intercept_s_curve_uexpr(fppoly_internal_t * pr, double * slope_inf, double *slope_sup, 
						double *intercept_inf, double *intercept_sup, double inf_coeff, 
						double sup_coeff, double lb, double ub, bool is_sigmoid, bool *boxify);

double apply_sigmoid_lexpr(fppoly_internal_t *pr, expr_t **lexpr_p, neuron_t * neuron);

double apply_sigmoid_uexpr(fppoly_internal_t *pr, expr_t **lexpr_p, neuron_t * neuron);

double apply_tanh_lexpr(fppoly_internal_t *pr, expr_t **lexpr_p, neuron_t * neuron);

double apply_tanh_uexpr(fppoly_internal_t *pr, expr_t **lexpr_p, neuron_t * neuron);

#ifdef __cplusplus
 }
#endif

#endif

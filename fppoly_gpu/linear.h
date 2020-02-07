/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
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
 *  @file linear.h
 *  @author Christoph Müller
 *  @brief Provides all linear and affine transformation functions.
 */

#ifndef __LINEAR_H_INCLUDED__
#define __LINEAR_H_INCLUDED__

#include "interval.h"

__device__ void add_coeff(float_type &lhs_inf_coeff, float_type &lhs_sup_coeff,
                          float_type &maxRes, float_type &maxMul,
                          const float_type rhs_inf_coeff,
                          const float_type rhs_sup_coeff) {
  maxRes = max(abs(lhs_inf_coeff), abs(lhs_sup_coeff));
  maxMul = max(abs(rhs_inf_coeff), abs(rhs_sup_coeff));

  lhs_inf_coeff =
      add_rd(lhs_inf_coeff,
             add_rd(rhs_inf_coeff, -mul_ru(add_ru(maxRes, maxMul), ulp)));
  lhs_sup_coeff =
      add_ru(lhs_sup_coeff,
             add_ru(rhs_sup_coeff, mul_ru(add_ru(maxRes, maxMul), ulp)));
}

__device__ void add_cst(float_type &lhs_inf_cst, float_type &lhs_sup_cst,
                        float_type &maxRes, float_type &maxMul,
                        const float_type rhs_inf_cst,
                        const float_type rhs_sup_cst) {
  maxRes = max(abs(lhs_inf_cst), abs(lhs_sup_cst));
  maxMul = max(abs(rhs_inf_cst), abs(rhs_sup_cst));

  lhs_inf_cst =
      add_rd(lhs_inf_cst,
             add_rd(rhs_inf_cst, -add_rd(mul_ru(add_ru(maxRes, maxMul), ulp),
                                         -min_denormal)));
  lhs_sup_cst =
      add_ru(lhs_sup_cst,
             add_ru(rhs_sup_cst,
                    add_ru(mul_ru(add_ru(maxRes, maxMul), ulp), min_denormal)));
}

__device__ void affine_trans_coeff(float_type &inf_coeff, float_type &sup_coeff,
                                   float_type &tmp1, float_type &tmp2,
                                   float_type &maxRes, float_type &maxMul,
                                   const float_type prev_inf_coeff,
                                   const float_type prev_sup_coeff,
                                   const float_type aux_coeff) {
  interval_mul_expr_coeff_const_expr(tmp1, tmp2, prev_inf_coeff, prev_sup_coeff,
                                     aux_coeff);

  maxRes = max(abs(inf_coeff), abs(sup_coeff));
  maxMul = max(abs(tmp1), abs(tmp2));

  inf_coeff =
      add_rd(inf_coeff, add_rd(tmp1, -mul_ru(add_ru(maxRes, maxMul), ulp)));
  sup_coeff =
      add_ru(sup_coeff, add_ru(tmp2, mul_ru(add_ru(maxRes, maxMul), ulp)));
}

__device__ void affine_trans_cst(float_type &inf_cst, float_type &sup_cst,
                                 float_type &tmp1, float_type &tmp2,
                                 float_type &maxRes, float_type &maxMul,
                                 const float_type inf_coeff,
                                 const float_type sup_coeff,
                                 const float_type aux_cst) {
  interval_mul_cst_coeff_const_expr(tmp1, tmp2, inf_coeff, sup_coeff, aux_cst);

  maxRes = max(abs(inf_cst), abs(sup_cst));
  maxMul = max(abs(tmp1), abs(tmp2));

  inf_cst =
      add_rd(inf_cst, add_rd(tmp1, -add_rd(mul_ru(add_ru(maxRes, maxMul), ulp),
                                           -min_denormal)));
  sup_cst = add_ru(
      sup_cst,
      add_ru(tmp2, add_ru(mul_ru(add_ru(maxRes, maxMul), ulp), min_denormal)));
}

__device__ void lb_component(float_type &res_inf, float_type &tmp1,
                             float_type &tmp2, const float_type inf_coeff,
                             const float_type sup_coeff,
                             const float_type input_inf,
                             const float_type input_sup) {
  interval_mul(tmp1, tmp2, inf_coeff, sup_coeff, input_inf, input_sup);

  res_inf = add_rd(res_inf, tmp1);
}

__device__ void ub_component(float_type &res_sup, float_type &tmp1,
                             float_type &tmp2, const float_type inf_coeff,
                             const float_type sup_coeff,
                             const float_type input_inf,
                             const float_type input_sup) {
  interval_mul(tmp1, tmp2, inf_coeff, sup_coeff, input_inf, input_sup);

  res_sup = add_ru(res_sup, tmp2);
}

__device__ void lcst_input_poly_neutral(
    float_type &inf_cst, float_type &sup_cst, float_type &tmp1,
    float_type &tmp2, float_type &maxRes, float_type &maxMul,
    const float_type prev_inf_coeff, const float_type prev_sup_coeff,
    const float_type input_inf, const float_type input_sup) {
  interval_mul(tmp1, tmp2, prev_inf_coeff, prev_sup_coeff, input_inf,
               input_sup);

  maxRes = max(abs(inf_cst), abs(sup_cst));
  maxMul = max(abs(-tmp1), abs(-tmp1));

  inf_cst =
      add_rd(inf_cst, -add_rd(tmp1, -add_rd(mul_ru(add_ru(maxRes, maxMul), ulp),
                                            -min_denormal)));
  sup_cst = add_ru(
      sup_cst,
      -add_ru(tmp1, add_ru(mul_ru(add_ru(maxRes, maxMul), ulp), min_denormal)));
}

__device__ void ucst_input_poly_neutral(
    float_type &inf_cst, float_type &sup_cst, float_type &tmp1,
    float_type &tmp2, float_type &maxRes, float_type &maxMul,
    const float_type prev_inf_coeff, const float_type prev_sup_coeff,
    const float_type input_inf, const float_type input_sup) {
  interval_mul(tmp1, tmp2, prev_inf_coeff, prev_sup_coeff, input_inf,
               input_sup);

  maxRes = max(abs(inf_cst), abs(sup_cst));
  maxMul = max(abs(tmp2), abs(tmp2));

  inf_cst =
      add_rd(inf_cst, add_rd(tmp2, -add_rd(mul_ru(add_ru(maxRes, maxMul), ulp),
                                           -min_denormal)));
  sup_cst = add_ru(
      inf_cst,
      add_ru(tmp2, add_ru(mul_ru(add_ru(maxRes, maxMul), ulp), min_denormal)));
}

#endif //__LINEAR_H_INCLUDED__

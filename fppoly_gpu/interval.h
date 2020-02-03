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
 *  @file interval.h
 *  @author Christoph Müller
 *  @brief Provides basic interval-interval or interval-scalar operations.
 */

#pragma once

__device__ void interval_mul(float_type &a_inf, float_type &a_sup,
                             const float_type b_inf, const float_type b_sup,
                             const float_type c_inf, const float_type c_sup) {
  float_type inf_inf = b_inf * c_inf;
  float_type inf_sup = b_inf * c_sup;
  float_type sup_inf = b_sup * c_inf;
  float_type sup_sup = b_sup * c_sup;

  if (inf_inf < inf_sup) {
    if (sup_inf < sup_sup) {
      if (inf_inf < sup_inf) {
        a_inf = mul_rd(b_inf, c_inf);
      } else {
        a_inf = mul_rd(b_sup, c_inf);
      }

      if (inf_sup < sup_sup) {
        a_sup = mul_ru(b_sup, c_sup);
      } else {
        a_sup = mul_ru(b_inf, c_sup);
      }
    } else {
      if (inf_inf < sup_sup) {
        a_inf = mul_rd(b_inf, c_inf);
      } else {
        a_inf = mul_rd(b_sup, c_sup);
      }

      if (inf_sup < sup_inf) {
        a_sup = mul_ru(b_sup, c_inf);
      } else {
        a_sup = mul_ru(b_inf, c_sup);
      }
    }
  } else {
    if (sup_inf < sup_sup) {
      if (inf_sup < sup_inf) {
        a_inf = mul_rd(b_inf, c_sup);
      } else {
        a_inf = mul_rd(b_sup, c_inf);
      }

      if (inf_inf < sup_sup) {
        a_sup = mul_ru(b_sup, c_sup);
      } else {
        a_sup = mul_ru(b_inf, c_inf);
      }
    } else {
      if (inf_sup < sup_sup) {
        a_inf = mul_rd(b_inf, c_sup);
      } else {
        a_inf = mul_rd(b_sup, c_sup);
      }

      if (inf_inf < sup_inf) {
        a_sup = mul_ru(b_sup, c_inf);
      } else {
        a_sup = mul_ru(b_inf, c_inf);
      }
    }
  }
}

__device__ void interval_mul_symmetric_c(float_type &a_inf, float_type &a_sup,
                                         const float_type b_inf,
                                         const float_type b_sup,
                                         const float_type c) {
  if (-b_inf < b_sup) {
    a_inf = mul_rd(b_sup, -c);
    a_sup = mul_ru(b_sup, c);
  } else {
    a_inf = mul_rd(b_inf, c);
    a_sup = mul_ru(b_inf, -c);
  }
}

__device__ void interval_mul_const_c(float_type &a_inf, float_type &a_sup,
                                     const float_type b_inf,
                                     const float_type b_sup,
                                     const float_type c) {
  if (c >= 0) {
    a_inf = mul_rd(b_inf, c);
    a_sup = mul_ru(b_sup, c);
  } else {
    a_inf = mul_rd(b_sup, c);
    a_sup = mul_ru(b_inf, c);
  }
}

__device__ void
interval_mul_expr_coeff(float_type &res_inf, float_type &res_sup,
                        const float_type inf, const float_type sup,
                        const float_type inf_expr, const float_type sup_expr) {
  interval_mul(res_inf, res_sup, inf, sup, inf_expr, sup_expr);

  const float_type maxA = max(abs(inf_expr), abs(sup_expr));
  float_type tmp1, tmp2;

  interval_mul_symmetric_c(tmp1, tmp2, inf, sup, maxA * ulp);

  res_inf = add_rd(res_inf, tmp1);
  res_sup = add_ru(res_sup, tmp2);
}

__device__ void interval_mul_cst_coeff(float_type &res_inf, float_type &res_sup,
                                       const float_type inf,
                                       const float_type sup,
                                       const float_type inf_expr,
                                       const float_type sup_expr) {
  interval_mul_expr_coeff(res_inf, res_sup, inf, sup, inf_expr, sup_expr);

  res_inf = add_rd(res_inf, -min_denormal);
  res_sup = add_ru(res_sup, min_denormal);
}

__device__ void interval_mul_expr_coeff_const_expr(float_type &res_inf,
                                                   float_type &res_sup,
                                                   const float_type inf,
                                                   const float_type sup,
                                                   const float_type expr) {
  interval_mul_const_c(res_inf, res_sup, inf, sup, expr);

  float_type tmp1, tmp2;

  interval_mul_symmetric_c(tmp1, tmp2, inf, sup, abs(expr) * ulp);

  res_inf = add_rd(res_inf, tmp1);
  res_sup = add_ru(res_sup, tmp2);
}

__device__ void interval_mul_cst_coeff_const_expr(float_type &res_inf,
                                                  float_type &res_sup,
                                                  const float_type inf,
                                                  const float_type sup,
                                                  const float_type expr) {
  interval_mul_expr_coeff_const_expr(res_inf, res_sup, inf, sup, expr);

  res_inf = add_rd(res_inf, -min_denormal);
  res_sup = add_ru(res_sup, min_denormal);
}

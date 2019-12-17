#pragma once

#include "interval.h"

__device__
void add_coeff(float_type& lhs_inf_coeff, float_type& lhs_sup_coeff, float_type& maxRes, float_type& maxMul, const float_type rhs_inf_coeff, const float_type rhs_sup_coeff)
{
    maxRes = max(abs(lhs_inf_coeff), abs(lhs_sup_coeff));
    maxMul = max(abs(rhs_inf_coeff), abs(rhs_sup_coeff));

    lhs_inf_coeff = lhs_inf_coeff + rhs_inf_coeff - (maxRes + maxMul)*ulp;
    lhs_sup_coeff = lhs_sup_coeff + rhs_sup_coeff + (maxRes + maxMul)*ulp;
}


__device__
void add_cst(float_type& lhs_inf_cst, float_type& lhs_sup_cst, float_type& maxRes, float_type& maxMul, const float_type rhs_inf_cst, const float_type rhs_sup_cst)
{
    maxRes = max(abs(lhs_inf_cst), abs(lhs_sup_cst));
    maxMul = max(abs(rhs_inf_cst), abs(rhs_sup_cst));

    lhs_inf_cst = lhs_inf_cst + rhs_inf_cst - (maxRes + maxMul)*ulp - min_denormal;
    lhs_sup_cst = lhs_sup_cst + rhs_sup_cst + (maxRes + maxMul)*ulp + min_denormal;
}


__device__
void affine_trans_coeff(float_type& inf_coeff, float_type& sup_coeff, float_type& tmp1, float_type& tmp2, float_type& maxRes, float_type& maxMul, const float_type prev_inf_coeff, const float_type prev_sup_coeff, const float_type aux_coeff)
{
    interval_mul_expr_coeff_const_expr(tmp1, tmp2, prev_inf_coeff, prev_sup_coeff, aux_coeff);

    maxRes = max(abs(inf_coeff), abs(sup_coeff));
    maxMul = max(abs(tmp1), abs(tmp2));

    inf_coeff = inf_coeff + tmp1 - (maxRes + maxMul)*ulp;
    sup_coeff = sup_coeff + tmp2 + (maxRes + maxMul)*ulp;
}


__device__
void affine_trans_cst(float_type& inf_cst, float_type& sup_cst, float_type& tmp1, float_type& tmp2, float_type& maxRes, float_type& maxMul, const float_type inf_coeff, const float_type sup_coeff, const float_type aux_cst)
{
    interval_mul_cst_coeff_const_expr(tmp1, tmp2, inf_coeff, sup_coeff, aux_cst);

    maxRes = max(abs(inf_cst), abs(sup_cst));
    maxMul = max(abs(tmp1), abs(tmp2));

    inf_cst = inf_cst + tmp1 - (maxRes + maxMul)*ulp - min_denormal;
    sup_cst = sup_cst + tmp2 + (maxRes + maxMul)*ulp + min_denormal;
}


__device__
void lb_component(float_type& res_inf, float_type& tmp1, float_type& tmp2, const float_type inf_coeff, const float_type sup_coeff, const float_type input_inf, const float_type input_sup)
{
    interval_mul(tmp1, tmp2, inf_coeff, sup_coeff, input_inf, input_sup);

    res_inf = res_inf + tmp1;
}


__device__
void ub_component(float_type& res_sup, float_type& tmp1, float_type& tmp2, const float_type inf_coeff, const float_type sup_coeff, const float_type input_inf, const float_type input_sup)
{
    interval_mul(tmp1, tmp2, inf_coeff, sup_coeff, input_inf, input_sup);

    res_sup = res_sup + tmp2;
}


__device__
void lcst_input_poly_neutral(float_type& inf_cst, float_type& sup_cst, float_type& tmp1, float_type& tmp2, float_type& maxRes, float_type& maxMul, const float_type prev_inf_coeff, const float_type prev_sup_coeff, const float_type input_inf, const float_type input_sup)
{
    interval_mul(tmp1, tmp2, prev_inf_coeff, prev_sup_coeff, input_inf, input_sup);

    maxRes = max(abs(inf_cst), abs(sup_cst));
    maxMul = max(abs(-tmp1), abs(-tmp1));

    inf_cst = inf_cst - tmp1 - (maxRes + maxMul)*ulp - min_denormal;
    sup_cst = sup_cst - tmp1 + (maxRes + maxMul)*ulp + min_denormal;
}


__device__
void ucst_input_poly_neutral(float_type& inf_cst, float_type& sup_cst, float_type& tmp1, float_type& tmp2, float_type& maxRes, float_type& maxMul, const float_type prev_inf_coeff, const float_type prev_sup_coeff, const float_type input_inf, const float_type input_sup)
{
    interval_mul(tmp1, tmp2, prev_inf_coeff, prev_sup_coeff, input_inf, input_sup);

    maxRes = max(abs(inf_cst), abs(sup_cst));
    maxMul = max(abs(tmp2), abs(tmp2));

    inf_cst = inf_cst + tmp2 - (maxRes + maxMul)*ulp - min_denormal;
    sup_cst = inf_cst + tmp2 + (maxRes + maxMul)*ulp + min_denormal;
}

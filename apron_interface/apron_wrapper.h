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

#ifndef _APRON_WRAPPER_H_
#define _APRON_WRAPPER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ap_coeff.h"
#include "ap_dimension.h"
#include "ap_expr0.h"
#include "ap_manager.h"
#include "ap_global0.h"
#include "ap_generic.h"

#define elina_scalar_sgn ap_scalar_sgn
#define elina_scalar_t ap_scalar_t
#define ELINA_SCALAR_MPQ AP_SCALAR_MPQ
#define ELINA_SCALAR_DOUBLE AP_SCALAR_DOUBLE
#define ELINA_SCALAR_MPFR AP_SCALAR_MPFR
#define elina_scalar_free ap_scalar_free
#define elina_scalar_alloc ap_scalar_alloc
#define elina_scalar_init ap_scalar_init
#define elina_scalar_reinit ap_scalar_reinit
#define elina_scalar_infty ap_scalar_infty
#define elina_scalar_set ap_scalar_set
#define elina_scalar_set_int ap_scalar_set_int
#define elina_scalar_swap ap_scalar_swap
#define elina_scalar_neg ap_scalar_neg
#define elina_scalar_set_infty ap_scalar_set_infty
#define elina_scalar_equal ap_scalar_equal
#define elina_scalar_cmp ap_scalar_cmp
#define elina_scalar_cmp_int ap_scalar_cmp_int
#define elina_scalar_discr_t ap_scalar_discr_t
#define elina_scalar_fprint ap_scalar_fprint
#define elina_double_set_scalar ap_double_set_scalar
#define elina_scalar_set_double ap_scalar_set_double

#define elina_interval_t ap_interval_t
#define elina_interval_neg ap_interval_neg
#define elina_coeff_swap ap_coeff_swap
#define elina_interval_set_bottom ap_interval_set_bottom
#define elina_interval_set_top ap_interval_set_top
#define elina_interval_alloc ap_interval_alloc
#define elina_interval_array_alloc ap_interval_array_alloc
#define elina_interval_reinit ap_interval_reinit
#define elina_interval_set ap_interval_set
#define elina_interval_set_int ap_interval_set_int
#define elina_interval_set_scalar ap_interval_set_scalar
#define elina_interval_is_bottom ap_interval_is_bottom
#define elina_interval_is_leq ap_interval_is_leq
#define elina_interval_free ap_interval_free
#define elina_interval_is_top ap_interval_is_top
#define elina_interval_fprint ap_interval_fprint
#define elina_interval_array_free ap_interval_array_free
#define elina_interval_print ap_interval_print
#define elina_interval_equal ap_interval_equal


#define elina_coeff_t ap_coeff_t
#define elina_coeff_alloc ap_coeff_alloc
#define elina_coeff_init ap_coeff_init
#define elina_coeff_zero ap_coeff_zero
#define elina_coeff_set ap_coeff_set
#define elina_coeff_neg ap_coeff_neg
#define elina_coeff_reinit ap_coeff_reinit
#define elina_coeff_reduce ap_coeff_reduce
#define elina_coeff_set_scalar ap_coeff_set_scalar
#define elina_coeff_set_scalar_int ap_coeff_set_scalar_int
#define elina_coeff_set_interval ap_coeff_set_interval
#define elina_coeff_set_interval_scalar ap_coeff_set_interval_scalar
#define elina_coeff_clear ap_coeff_clear
#define elina_coeff_free ap_coeff_free
#define elina_coeff_set_interval_int ap_coeff_set_interval_int
#define elina_coeff_print ap_coeff_print
#define ELINA_COEFF_SCALAR AP_COEFF_SCALAR
#define ELINA_COEFF_INTERVAL AP_COEFF_INTERVAL

#define elina_dim_t ap_dim_t
#define elina_dimension_t ap_dimension_t
#define elina_dimchange_t ap_dimchange_t
#define elina_dimperm_t ap_dimperm_t
#define elina_dimperm_init ap_dimperm_init
#define elina_dimperm_clear ap_dimperm_clear
#define elina_dimperm_set_id ap_dimperm_set_id
#define elina_dimchange_init ap_dimchange_init
#define elina_dimchange_add_invert ap_dimchange_add_invert
#define elina_dimchange_clear ap_dimchange_clear
#define elina_dimchange_free ap_dimchange_free
#define elina_dimchange_alloc ap_dimchange_alloc
#define ELINA_DIM_MAX AP_DIM_MAX

#define elina_funid_t ap_funid_t
#define elina_funopt_t ap_funopt_t
#define elina_exc_t ap_exc_t
#define elina_manager_t ap_manager_t
#define elina_manager_alloc ap_manager_alloc
#define elina_manager_copy ap_manager_copy
#define elina_manager_free ap_manager_free
#define elina_fpu_init ap_fpu_init
#define elina_manager_set_abort_if_exception ap_manager_set_abort_if_exception
#define elina_manager_raise_exception ap_manager_raise_exception
#define ELINA_EXC_TIMEOUT AP_EXC_TIMEOUT
#define ELINA_EXC_OUT_OF_SPACE AP_EXC_OUT_OF_SPACE
#define ELINA_EXC_OVERFLOW AP_EXC_OVERFLOW
#define ELINA_EXC_NONE AP_EXC_NONE
#define ELINA_EXC_NOT_IMPLEMENTED AP_EXC_NOT_IMPLEMENTED
#define ELINA_EXC_SIZE AP_EXC_SIZE


#define ELINA_FUNID_COPY AP_FUNID_COPY
#define ELINA_FUNID_FREE AP_FUNID_FREE
#define ELINA_FUNID_CLOSURE AP_FUNID_CLOSURE
#define ELINA_FUNID_ASIZE AP_FUNID_ASIZE
#define ELINA_FUNID_MINIMIZE AP_FUNID_MINIMIZE
#define ELINA_FUNID_CANONICALIZE  AP_FUNID_CANONICALIZE
#define ELINA_FUNID_HASH  AP_FUNID_HASH
#define ELINA_FUNID_APPROXIMATE AP_FUNID_APPROXIMATE
#define ELINA_FUNID_FPRINT AP_FUNID_FPRINT
#define ELINA_FUNID_BOTTOM AP_FUNID_BOTTOM
#define ELINA_FUNID_TOP AP_FUNID_TOP
#define ELINA_FUNID_OF_BOX AP_FUNID_OF_BOX
#define ELINA_FUNID_DIMENSION AP_FUNID_DIMENSION
#define ELINA_FUNID_IS_BOTTOM AP_FUNID_IS_BOTTOM
#define ELINA_FUNID_IS_TOP AP_FUNID_IS_TOP
#define ELINA_FUNID_IS_LEQ AP_FUNID_IS_LEQ
#define ELINA_FUNID_IS_EQ AP_FUNID_IS_EQ
#define ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED AP_FUNID_IS_DIMENSION_UNCONSTRAINED
#define ELINA_FUNID_SAT_INTERVAL AP_FUNID_SAT_INTERVAL
#define ELINA_FUNID_SAT_LINCONS AP_FUNID_SAT_LINCONS
#define ELINA_FUNID_SAT_TCONS AP_FUNID_SAT_TCONS
#define ELINA_FUNID_BOUND_DIMENSION AP_FUNID_BOUND_DIMENSION
#define ELINA_FUNID_BOUND_LINEXPR AP_FUNID_BOUND_LINEXPR
#define ELINA_FUNID_BOUND_TEXPR AP_FUNID_BOUND_TEXPR
#define ELINA_FUNID_TO_BOX AP_FUNID_TO_BOX
#define ELINA_FUNID_TO_LINCONS_ARRAY AP_FUNID_TO_LINCONS_ARRAY
#define ELINA_FUNID_TO_TCONS_ARRAY AP_FUNID_TO_TCONS_ARRAY
#define ELINA_FUNID_MEET AP_FUNID_MEET
#define ELINA_FUNID_MEET_ARRAY AP_FUNID_MEET_ARRAY
#define ELINA_FUNID_MEET_LINCONS_ARRAY AP_FUNID_MEET_LINCONS_ARRAY
#define ELINA_FUNID_MEET_TCONS_ARRAY AP_FUNID_MEET_TCONS_ARRAY
#define ELINA_FUNID_JOIN AP_FUNID_JOIN
#define ELINA_FUNID_JOIN_ARRAY AP_FUNID_JOIN_ARRAY
#define ELINA_FUNID_ASSIGN_LINEXPR_ARRAY AP_FUNID_ASSIGN_LINEXPR_ARRAY
#define ELINA_FUNID_SUBSTITUTE_LINEXPR_ARRAY AP_FUNID_SUBSTITUTE_LINEXPR_ARRAY
#define ELINA_FUNID_ASSIGN_TEXPR_ARRAY AP_FUNID_ASSIGN_TEXPR_ARRAY
#define ELINA_FUNID_ADD_DIMENSIONS AP_FUNID_ADD_DIMENSIONS
#define ELINA_FUNID_REMOVE_DIMENSIONS AP_FUNID_REMOVE_DIMENSIONS
#define ELINA_FUNID_PERMUTE_DIMENSIONS AP_FUNID_PERMUTE_DIMENSIONS
#define ELINA_FUNID_FORGET_ARRAY AP_FUNID_FORGET_ARRAY
#define ELINA_FUNID_EXPAND AP_FUNID_EXPAND
#define ELINA_FUNID_FOLD AP_FUNID_FOLD
#define ELINA_FUNID_WIDENING AP_FUNID_WIDENING
#define ELINA_FUNID_SIZE AP_FUNID_SIZE

#define elina_abstract0_t ap_abstract0_t
#define elina_abstract0_is_bottom ap_abstract0_is_bottom
#define elina_abstract0_meet_lincons_array ap_abstract0_meet_lincons_array 
#define elina_abstract0_assign_linexpr_array ap_abstract0_assign_linexpr_array
#define elina_abstract0_top ap_abstract0_top
#define elina_abstract0_bottom ap_abstract0_bottom
#define elina_abstract0_is_top ap_abstract0_is_top
#define elina_abstract0_is_bottom ap_abstract0_is_bottom
#define elina_abstract0_is_eq ap_abstract0_is_eq
#define elina_abstract0_is_leq ap_abstract0_is_leq
#define elina_abstract0_meet ap_abstract0_meet
#define elina_abstract0_join ap_abstract0_join
#define elina_abstract0_widening ap_abstract0_widening
#define elina_abstract0_to_box ap_abstract0_to_box
#define elina_abstract0_of_box ap_abstract0_of_box
#define elina_abstract0_dimension ap_abstract0_dimension

#define elina_constyp_t ap_constyp_t
#define elina_linexpr_discr_t ap_linexpr_discr_t
#define elina_linterm_t ap_linterm_t
#define elina_linexpr_type_t ap_linexpr_type_t
#define elina_linexpr0_t ap_linexpr0_t
#define elina_linexpr0_alloc ap_linexpr0_alloc
#define elina_linexpr0_copy ap_linexpr0_copy
#define elina_linexpr0_free ap_linexpr0_free
#define elina_linexpr0_is_linear ap_linexpr0_is_linear
#define elina_linexpr0_is_quasilinear ap_linexpr0_is_quasilinear
#define elina_linexpr0_is_real ap_linexpr0_is_real
#define elina_linexpr0_set_cst_scalar_int ap_linexpr0_set_cst_scalar_int
#define elina_linexpr0_set_cst_interval_int ap_linexpr0_set_cst_interval_int
#define elina_linexpr0_set_coeff_scalar_double ap_linexpr0_set_coeff_scalar_double
#define elina_linexpr0_array_is_linear ap_linexpr0_array_is_linear
#define elina_linexpr0_array_free ap_linexpr0_array_free
#define elina_linexpr0_realloc ap_linexpr0_realloc
#define elina_linexpr0_add_dimensions ap_linexpr0_add_dimensions
#define elina_linexpr0_is_integer ap_linexpr0_is_integer
#define elina_linexpr0_fprint ap_linexpr0_fprint
#define ELINA_LINEXPR_DENSE AP_LINEXPR_DENSE
#define ELINA_LINEXPR_SPARSE AP_LINEXPR_SPARSE
#define ELINA_LINEXPR_LINEAR AP_LINEXPR_LINEAR
#define ELINA_LINEXPR_INTLINEAR AP_LINEXPR_INTLINEAR
#define ELINA_LINEXPR_QUASILINEAR AP_LINEXPR_QUASILINEAR
#define elina_linexpr0_ForeachLinterm ap_linexpr0_ForeachLinterm

#define elina_lincons0_t ap_lincons0_t
#define elina_lincons0_copy ap_lincons0_copy
#define elina_lincons0_array_t ap_lincons0_array_t
#define elina_lincons0_make ap_lincons0_make
#define elina_lincons0_array_make ap_lincons0_array_make
#define elina_lincons0_make_unsat ap_lincons0_make_unsat
#define elina_lincons0_is_unsat ap_lincons0_is_unsat
#define elina_lincons0_array_fprint ap_lincons0_array_fprint
#define elina_lincons0_array_clear ap_lincons0_array_clear
#define elina_lincons0_clear ap_lincons0_clear
#define elina_lincons0_array_is_quasilinear ap_lincons0_array_is_quasilinear
#define elina_lincons0_fprint ap_lincons0_fprint
#define elina_lincons0_is_sat ap_lincons0_is_sat
#define ELINA_CONS_DISEQ AP_CONS_DISEQ
#define ELINA_CONS_EQ AP_CONS_EQ
#define ELINA_CONS_SUPEQ AP_CONS_SUPEQ
#define ELINA_CONS_SUP AP_CONS_SUP
#define ELINA_CONS_EQMOD AP_CONS_EQMOD

#define elina_texpr0_t ap_texpr0_t
#define elina_texpr0_from_linexpr0 ap_texpr0_from_linexpr0
#define elina_texpr_rtype_t ap_texpr_rtype_t
#define elina_texpr_rdir_t ap_texpr_rdir_t
#define elina_texpr0_node_t ap_texpr0_node_t
#define ELINA_TEXPR_CST AP_TEXPR_CST
#define ELINA_TEXPR_DIM AP_TEXPR_DIM
#define ELINA_TEXPR_NODE AP_TEXPR_NODE
#define ELINA_RTYPE_INT AP_RTYPE_INT
#define ELINA_RDIR_ZERO AP_RDIR_ZERO
#define ELINA_RDIR_UP AP_RDIR_UP
#define ELINA_RDIR_DOWN AP_RDIR_DOWN
#define ELINA_RDIR_NEAREST AP_RDIR_NEAREST
#define ELINA_RTYPE_SINGLE AP_RTYPE_SINGLE
#define ELINA_RTYPE_QUAD AP_RTYPE_QUAD
#define ELINA_RTYPE_EXTENDED AP_RTYPE_EXTENDED
#define ELINA_RTYPE_DOUBLE AP_RTYPE_DOUBLE
#define ELINA_TEXPR_NEG AP_TEXPR_NEG
#define ELINA_TEXPR_CAST AP_TEXPR_CAST
#define ELINA_TEXPR_SQRT AP_TEXPR_SQRT
#define ELINA_TEXPR_ADD AP_TEXPR_ADD
#define ELINA_TEXPR_MUL AP_TEXPR_MUL
#define ELINA_TEXPR_DIV AP_TEXPR_DIV
#define ELINA_TEXPR_POW AP_TEXPR_POW
#define ELINA_TEXPR_MOD AP_TEXPR_MOD
#define ELINA_RTYPE_INT AP_RTYPE_INT
#define elina_texpr0_dim ap_texpr0_dim
#define elina_texpr0_binop ap_texpr0_binop
#define elina_texpr0_is_interval_cst ap_texpr0_is_interval_cst
#define elina_texpr0_is_interval_linear ap_texpr0_is_interval_linear
#define elina_texpr0_array_is_scalar ap_texpr0_array_is_scalar
#define elina_texpr0_array_is_interval_linear ap_texpr0_array_is_interval_linear
#define elina_texpr0_array_is_interval_linear ap_texpr0_array_is_interval_linear
#define elina_texpr0_add_dimensions ap_texpr0_add_dimensions
#define ELINA_RDIR_RND AP_RDIR_RND
#define ELINA_RTYPE_REAL AP_RTYPE_REAL
#define ELINA_TEXPR_SUB AP_TEXPR_SUB


#define elina_tcons0_t ap_tcons0_t
#define elina_tcons0_make ap_tcons0_make
#define elina_tcons0_fprint ap_tcons0_fprint
#define elina_tcons0_is_interval_linear ap_tcons0_is_interval_linear
#define elina_tcons0_array_t ap_tcons0_array_t
#define elina_tcons0_array_make ap_tcons0_array_make
#define elina_tcons0_array_clear ap_tcons0_array_clear
#define elina_tcons0_from_lincons0 ap_tcons0_from_lincons0



#ifdef __cplusplus
}
#endif

#endif

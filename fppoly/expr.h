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



#ifndef __EXPR_H_INCLUDED__
#define __EXPR_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif

#include "fppoly.h"

void elina_double_interval_add_expr_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup, double inf, double sup, double inf_expr, double sup_expr);

void elina_double_interval_add_cst_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup, double inf, double sup, double inf_expr, double sup_expr);

void elina_double_interval_mul_expr_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup, double inf, double sup, double inf_expr, double sup_expr);

void elina_double_interval_mul_cst_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup, double inf, double sup, double inf_expr, double sup_expr);

void expr_fprint(FILE * stream, expr_t *expr);

void expr_print(expr_t * expr);

expr_t * alloc_expr(void);

expr_t * create_dense_expr(double *coeff, double cst, size_t size);

expr_t * create_cst_expr(double l, double u);

expr_t * create_sparse_expr(double *coeff, double cst, size_t *dim, size_t size);

void free_expr(expr_t *expr);

expr_t * copy_cst_expr(expr_t *src);

expr_t * copy_expr(expr_t *src);

expr_t* concretize_dense_sub_expr(fppoly_internal_t *pr, expr_t * expr, double *inf, double *sup, size_t start, size_t size);

void merge_sparse_expr(expr_t *expr, size_t l, size_t m, size_t r);

void merge_sort_sparse_expr(expr_t *expr, size_t l, size_t r);

void sort_sparse_expr(expr_t *expr);

expr_t * multiply_expr(fppoly_internal_t *pr, expr_t *expr, double mul_inf, double mul_sup);

expr_t * multiply_cst_expr(fppoly_internal_t *pr, expr_t *expr, double mul_inf, double mul_sup);

void add_cst_expr(fppoly_internal_t *pr, expr_t * exprA, expr_t *exprB);

void add_expr(fppoly_internal_t *pr,expr_t * exprA, expr_t * exprB);

expr_t *extract_subexpr(expr_t *expr, size_t index_start, size_t num_neurons);

expr_t * lexpr_replace_bounds(fppoly_internal_t * pr, expr_t * expr, neuron_t ** neurons, bool is_activation);

expr_t * uexpr_replace_bounds(fppoly_internal_t * pr, expr_t * expr, neuron_t ** neurons, bool is_activation);

#ifdef __cplusplus
 }
#endif

#endif

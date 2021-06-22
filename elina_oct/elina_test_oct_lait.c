/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
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

#include "opt_oct.h"
#include "opt_oct_internal.h"
#include "opt_oct_hmat.h"
#include "opt_oct_lait.h"

int mat_head[6][7] = {
    {1, 1, 0, 0, 1},
    {1, 1, 1, 0, 1},
    {1, 1, 4, 0, 1},
    {1, -1, 4, 0, 1},
    {1, 1, 5, 0, 1},
    {1, -1, 5, 0, 1}
};

int mat_oa[18][7] = {
    {1, 1, 1, 0, 1},
    {1, 1, 2, 0, 1},
    {1, 1, 5, 0, 1},
    {1, -1, 5, 0, 1},
    {1, 1, 6, 0, 1},
    {1, -1, 6, 0, 1},
    {2, -1, 6, 1, 8, 0, 1},
    {2, 1, 6, 1, 8, 0, 1},
    {1, 1, 8, 0, 1},
    {2, -1, 6, -1, 8, 1, 1},
    {2, 1, 6, -1, 8, 1, 1},
    {1, -1, 8, 1, 1},
    {2, -1, 5, 1, 9, 0, 1},
    {2, 1, 5, 1, 9, 0, 1},
    {1, 1, 9, 0, 1},
    {2, -1, 5, -1, 9, 0, 1},
    {2, 1, 5, -1, 9, 0 ,1},
    {1, -1, 9, 0, 1}
};

int mat_ob[16][7] = {
    {1, 1, 1, 0, 1},
    {1, 1, 2, 0, 1},
    {1, 1, 5, 0, 1},
    {1, -1, 5, 0, 1},
    {1, 1, 6, 0, 1},
    {1, -1, 6, 0, 1},
    {1, 1, 7, -8, 1},
    {1, -1, 7, 8, 1},
    {2, -1, 6, 1, 8, 0, 1},
    {2, 1, 6, 1, 8, 0, 1},
    {1, 1, 8, 0, 1},
    {2, -1, 6, -1, 8, 0, 1},
    {2, 1, 6, -1, 8, 0, 1},
    {1, -1, 8, 0, 1},
    {1, 1, 9, -1, 1},
    {1, -1, 9, 1, 1},
};

opt_oct_t* generate_oct(elina_manager_t * man, int mat[][7], int num_cons, int num_vars) {
    elina_lincons0_array_t lincons0 = elina_lincons0_array_make(num_cons);

    for (int i=0; i<num_cons; i++) {
        int size = mat[i][0];
		elina_coeff_t *coeff;
		elina_linexpr0_t *linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, size);
        for (int j=0; j<size; j++) {
            int cf = mat[i][1+2*j];
            int dim = mat[i][2+2*j];
			elina_linterm_t *linterm = &linexpr0->p.linterm[j];
			linterm->dim = dim;
			coeff = &linterm->coeff;
			elina_scalar_set_to_int(coeff->val.scalar, cf, ELINA_SCALAR_MPQ);
        }

        int cst = mat[i][1+2*size];
        int cst_type = mat[i][2+2*size];
		elina_scalar_set_to_int((&linexpr0->cst)->val.scalar, cst, ELINA_SCALAR_MPQ);
		if (cst_type == 0) {
			lincons0.p[i].constyp = ELINA_CONS_EQ;
		} else if (cst_type == 1) {
			lincons0.p[i].constyp = ELINA_CONS_SUPEQ;
		} else if (cst_type == 2) {
			lincons0.p[i].constyp = ELINA_CONS_SUP;
		} else if (cst_type == 3) {
			lincons0.p[i].constyp = ELINA_CONS_DISEQ;
		}
		lincons0.p[i].linexpr0 = linexpr0;
    }

    opt_oct_t* op = opt_oct_top(man, num_vars, 0);
    op = opt_oct_meet_lincons_array(man, true, op, &lincons0);
    elina_lincons0_array_clear(&lincons0);

    return op;
}

int main(int argc, char** argv) {
    elina_manager_t * man = opt_oct_manager_alloc();

    // call opt_oct_lait_init at the start of the analysis to initialize Lait
    opt_oct_lait_init("/home/hej/ELINA/elina_oct", "/home/hej/ELINA/elina_oct/opt_oct_lait.pt");

    // the first join input
    opt_oct_t* oa = generate_oct(man, mat_oa, 18, 56);
    elina_lincons0_array_t arr_a = opt_oct_to_lincons_array(man, oa);
    elina_lincons0_array_fprint(stdout, &arr_a, NULL);
    elina_lincons0_array_clear(&arr_a);

    // the second join input
    opt_oct_t* ob = generate_oct(man, mat_ob, 16, 56);
    elina_lincons0_array_t arr_b = opt_oct_to_lincons_array(man, ob);
    elina_lincons0_array_fprint(stdout, &arr_b, NULL);
    elina_lincons0_array_clear(&arr_b);

    // the join output
    opt_oct_t* res = opt_oct_join(man, false, oa, ob);
    elina_lincons0_array_t arr_res = opt_oct_to_lincons_array(man, res);
    elina_lincons0_array_fprint(stdout, &arr_res, NULL);
    elina_lincons0_array_clear(&arr_res);

    // the abstract element at the loop head
    opt_oct_t* head = generate_oct(man, mat_head, 6, 6);
    elina_lincons0_array_t arr_head = opt_oct_to_lincons_array(man, head);
    elina_lincons0_array_fprint(stdout, &arr_head, NULL);
    elina_lincons0_array_clear(&arr_head);

    // applying Lait
    opt_oct_t* lait_res = opt_oct_lait(man, false, oa, ob, res, head, 0);
    elina_lincons0_array_t arr_lait = opt_oct_to_lincons_array(man, lait_res);
    elina_lincons0_array_fprint(stdout, &arr_lait, NULL);
    elina_lincons0_array_clear(&arr_lait);

    return 0;
}
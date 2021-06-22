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

#include "opt_pk.h"
#include "opt_pk_lait.h"

int mat_oa[18][17] = {
    {4, -2, 2, -2, 3, 1, 5, 2, 6, -1, 0},
    {6, -2, 3, 2, 4, -1, 5, 2, 9, -2, 10, 1, 11, 0, 1},
    {2, -1, 5, 1, 11, 0, 1},
    {5, 3, 5, -10, 6, -2, 7, -3, 11, 10, 12, 10, 1},
    {5, 3, 5, 10, 6, -2, 7, -3, 11, -10, 12, 10, 1},
    {2, 1, 3, -1, 9, 0, 1},
    {6, 2, 2, -1, 5, -2, 6, -2, 8, 1, 11, 2, 12, 0, 1},
    {1, -1, 11, 3, 1},
    {7, 20, 3, -10, 4, 7, 5, 2, 7, -20, 9, 10, 10, -7, 11, 0, 1},
    {3, 1, 5, 1, 7, -1, 11, 0, 1},
    {7, -2, 2, -2, 3, 1, 5, 4, 6, 2, 8, -1, 11, -2, 12, 0, 1},
    {6, -2, 2, -2, 3, 1, 5, 2, 6, 2, 8, -1, 11, 0, 1},
    {5, 2, 5, 2, 7, -10, 9, 10, 10, -7, 11, 5, 1},
    {3, -2, 3, 2, 4, -1, 5, 1, 1},
    {3, 4, 9, -2, 10, 1, 11, -1, 1},
    {3, -1, 9, 2, 10, -1, 11, 1, 1},
    {6, 4, 3, -2, 4, 1, 5, -4, 9, 4, 10, -2, 11, 1, 1},
    {3, 4, 3, -2, 4, 1, 5, -1, 1}
};

int mat_ob[15][17] = {
    {1, -1, 73, 2, 0},
    {1, 1, 7, -5, 0},
    {1, -1, 75, 3, 1},
    {1, 1, 75, -2, 1},
    {2, -1, 10, 1, 75, 0, 0},
    {4, -2, 2, -2, 3, 1, 5, 2, 6, -1, 0},
    {1, -1, 4, 2, 0},
    {2, 1, 2, -1, 8, 0, 0},
    {2, 1, 3, -1, 9, 0, 0},
    {2, 1, 5, -1, 11, 0, 0},
    {2, 4, 3, 1, 5, -5, 1},
    {1, -1, 5, 3, 1},
    {2, -2, 3, -1, 5, 5, 1},
    {2, -1, 3, 1, 6, 0, 1},
    {2, 1, 6, -1, 12, 0, 0}
};

int mat_head[5][17] = {
    {1, -1, 5, 1, 0},
    {1, -1 ,3, 0, 0},
    {1, -1, 2, 0, 0},
    {1, -1, 4, 0, 0},
    {1, -1, 6, 0, 0}
};

opt_pk_array_t* generate_poly(elina_manager_t * man, int mat[][17], int num_cons, int num_vars) {
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

	opt_pk_array_t* op = opt_pk_top(man, num_vars, 0);
	op = opt_pk_meet_lincons_array(man, true, op, &lincons0);
	elina_lincons0_array_clear(&lincons0);

	return op;
}

int main(int argc, char** argv) {
    elina_manager_t * man = opt_pk_manager_alloc(false);


    // call opt_pk_lait_init at the start of the analysis to initialize Lait
    opt_pk_lait_init("/home/hej/ELINA/elina_poly", "/home/hej/ELINA/elina_poly/opt_pk_lait.pt");

    // the first join input
    opt_pk_array_t* oa = generate_poly(man, mat_oa, 18, 76);
    elina_lincons0_array_t arr_a = opt_pk_to_lincons_array(man, oa);
    elina_lincons0_array_fprint(stdout, &arr_a, NULL);
    elina_lincons0_array_clear(&arr_a);

    // the second join input
    opt_pk_array_t* ob = generate_poly(man, mat_ob, 15, 76);
    elina_lincons0_array_t arr_b = opt_pk_to_lincons_array(man, ob);
    elina_lincons0_array_fprint(stdout, &arr_b, NULL);
    elina_lincons0_array_clear(&arr_b);

    // the join output
    opt_pk_array_t* res = opt_pk_join(man, false, oa, ob);
    elina_lincons0_array_t arr_res = opt_pk_to_lincons_array(man, res);
    elina_lincons0_array_fprint(stdout, &arr_res, NULL);
    elina_lincons0_array_clear(&arr_res);

    // the abstract element at the loop head
    opt_pk_array_t* head = generate_poly(man, mat_head, 5, 15);
    elina_lincons0_array_t arr_head = opt_pk_to_lincons_array(man, head);
    elina_lincons0_array_fprint(stdout, &arr_head, NULL);
    elina_lincons0_array_clear(&arr_head);

    // applying Lait
    opt_pk_array_t* lait_res = opt_pk_lait(man, false, oa, ob, res, head, 0);
    elina_lincons0_array_t arr_lait = opt_pk_to_lincons_array(man, lait_res);
    elina_lincons0_array_fprint(stdout, &arr_lait, NULL);
    elina_lincons0_array_clear(&arr_lait);

    return 0;
}
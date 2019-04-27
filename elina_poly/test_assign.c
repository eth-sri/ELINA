#include "opt_pk.h"

void test_assign(){
	elina_lincons0_array_t lincons = elina_lincons0_array_make(3);
	lincons.p[0].constyp = ELINA_CONS_SUPEQ;
	lincons.p[1].constyp = ELINA_CONS_SUPEQ;
	lincons.p[2].constyp = ELINA_CONS_SUPEQ;

	lincons.p[0].linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2);
	elina_linexpr0_set_cst_scalar_int(lincons.p[0].linexpr0,1);
        elina_linexpr0_set_coeff_scalar_int(lincons.p[0].linexpr0,0,1);
	elina_linexpr0_set_coeff_scalar_int(lincons.p[0].linexpr0,1,1);

	lincons.p[1].linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2);
	elina_linexpr0_set_cst_scalar_int(lincons.p[1].linexpr0,1);
        elina_linexpr0_set_coeff_scalar_int(lincons.p[1].linexpr0,2,-1);
 	elina_linexpr0_set_coeff_scalar_int(lincons.p[1].linexpr0,3,-1);
	
	lincons.p[2].linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
	elina_linexpr0_set_cst_scalar_int(lincons.p[2].linexpr0,1);
	elina_linexpr0_set_coeff_scalar_int(lincons.p[2].linexpr0,4,1);
	elina_linexpr0_set_coeff_scalar_int(lincons.p[2].linexpr0,5,1);
	
	elina_manager_t * man = opt_pk_manager_alloc(false);
	//generate first input
	opt_pk_array_t * oa = opt_pk_top(man, 6,0);
	
	//meet with constraints
	oa = opt_pk_meet_lincons_array(man,true,oa,&lincons);
	
	elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2);
	elina_linexpr0_set_cst_scalar_int(linexpr0,1);
        elina_linexpr0_set_coeff_scalar_int(linexpr0,2,-1);
 	elina_linexpr0_set_coeff_scalar_int(linexpr0,4,-1);
	
	elina_dim_t * tdim = (elina_dim_t *)malloc(sizeof(elina_dim_t));
	tdim[0] = 0;
	elina_linexpr0_t ** expr_array = (elina_linexpr0_t**)malloc(sizeof(elina_linexpr0_t*));
	expr_array[0] = linexpr0;

	printf("ELINA Input Polyhedra\n");
	elina_lincons0_array_t arr = opt_pk_to_lincons_array(man,oa);
  	elina_lincons0_array_fprint(stdout,&arr,NULL);
	printf("Assignment statement\n");
	printf("x%d = ",tdim[0]);
	elina_linexpr0_fprint(stdout,linexpr0,NULL);
	printf("\n");
  	fflush(stdout);

	oa = opt_pk_assign_linexpr_array(man,true,oa,tdim, expr_array,1,NULL);
	elina_linexpr0_free(linexpr0);
	free(expr_array);
	free(tdim);

	//meet with -x1 + 43 >= 0; 	

	// Print the result
	printf("ELINA Output Polyhedron\n");
	arr = opt_pk_to_lincons_array(man,oa);
  	elina_lincons0_array_fprint(stdout,&arr,NULL);
	printf("\n");
  	fflush(stdout);
}

int main(){
	test_assign();

}

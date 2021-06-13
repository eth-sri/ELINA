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

#include <time.h>
#include "zonoml.h"
#include "zonoml_internal.h"
#include "zonoml_fun.h"


elina_interval_t ** generate_random_box(size_t dim){
	elina_interval_t ** interval = (elina_interval_t **)malloc(dim*sizeof(elina_interval_t *));
	unsigned short int i;
	for(i=0; i < dim; i++){
		interval[i] = elina_interval_alloc();
		elina_interval_set_top(interval[i]);
	}
	for(i = 0; i < dim; i++){
		int inf = rand()%10;
		int sup = rand()%10 + inf;
		elina_scalar_set_to_int(interval[i]->inf,inf,ELINA_SCALAR_DOUBLE);
		elina_scalar_set_to_int(interval[i]->sup,sup,ELINA_SCALAR_DOUBLE);
	}
	return interval;
}


elina_linexpr0_t * generate_random_linexpr0(size_t dim){
	elina_coeff_t *cst, *coeff;
	unsigned short int j, k;
	elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,dim-1);
	cst = &linexpr0->cst;
	elina_scalar_set_to_int(cst->val.scalar,rand()%10,ELINA_SCALAR_DOUBLE);
        for(j=0; j < dim - 1; j++){
		elina_linterm_t * linterm = &linexpr0->p.linterm[j];
		linterm->dim = j;
		coeff = &linterm->coeff;
		elina_scalar_set_to_int(coeff->val.scalar,rand()%10+1,ELINA_SCALAR_DOUBLE);
	}
	return linexpr0;
}


elina_lincons0_array_t generate_random_lincons0_array(size_t dim, size_t nbcons){
	size_t i;
	unsigned short int j, k; 
	elina_coeff_t *cst, *coeff;
	elina_lincons0_array_t  lincons0 = elina_lincons0_array_make(nbcons);
	for(i=0; i < nbcons; i++){
		lincons0.p[i].constyp =  ELINA_CONS_SUPEQ;
		int v1 = rand()%dim;
		int v2;
		while(1){
			v2 = rand()%dim;
			if(v1!=v2){
				break;
			}
		}
		int coeff1 = rand()%2==0? -1 :1;
		int coeff2 = rand()%2==0? -1 :1;
		elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2);
		cst = &linexpr0->cst;
		elina_scalar_set_to_int(cst->val.scalar,rand()%10,ELINA_SCALAR_DOUBLE);
		elina_linterm_t * linterm = &linexpr0->p.linterm[0];
		linterm->dim = v1;
		coeff = &linterm->coeff;
		elina_scalar_set_to_int(coeff->val.scalar,coeff1,ELINA_SCALAR_DOUBLE);
		linterm = &linexpr0->p.linterm[1];
		linterm->dim = v2;
		coeff = &linterm->coeff;
		elina_scalar_set_to_int(coeff->val.scalar,coeff2,ELINA_SCALAR_DOUBLE);
		lincons0.p[i].linexpr0 = linexpr0;
	}
	return lincons0;
}


void test_zonoml(size_t dim, size_t nbcons){
	elina_manager_t * man = zonotope_manager_alloc();
	size_t i;
	//generate random box
	elina_interval_t ** interval = generate_random_box(dim);
	zonotope_t * zo = zonotope_of_box(man,dim,0,interval);
	printf("Input Interval\n");
	for(i = 0; i < dim; i++){
		printf("x%lu: ",i);
		elina_interval_print(interval[i]);
		printf(" ");
		elina_interval_free(interval[i]);
	}
	printf("\n");
	printf("ELINA Output zonotope\n");
	elina_lincons0_array_t arr = zonotope_to_lincons_array(man,zo);
  	elina_lincons0_array_fprint(stdout,&arr,NULL);
	printf("\n");
  	fflush(stdout);
 	elina_lincons0_array_clear(&arr);
	free(interval);

	//random assignments
	elina_dim_t * tdim = (elina_dim_t *)malloc(sizeof(elina_dim_t));
	tdim[0] = dim-1;
	elina_linexpr0_t ** expr_array = (elina_linexpr0_t**)malloc(sizeof(elina_linexpr0_t*));
	elina_linexpr0_t * linexpr0 = generate_random_linexpr0(dim);
	expr_array[0] = linexpr0;
	printf("Assignment statement\n");
	printf("x%d = ",tdim[0]);
	elina_linexpr0_fprint(stdout,linexpr0,NULL);
	printf("\n");
  	fflush(stdout);
	elina_lincons0_array_clear(&arr);
	//assign;
	zo = zonotope_assign_linexpr_array(man,true,zo,tdim, expr_array,1,NULL);
	elina_linexpr0_free(linexpr0);
	free(expr_array);
	free(tdim);

	//meet with -x1 + 43 >= 0; 	

	// Print the result
	printf("ELINA Output zonotope after assignment\n");
	arr = zonotope_to_lincons_array(man,zo);
  	elina_lincons0_array_fprint(stdout,&arr,NULL);
	printf("\n");
  	fflush(stdout);
 	elina_lincons0_array_clear(&arr);


	//random meet lincons
	elina_lincons0_array_t lincons0 = generate_random_lincons0_array(dim,nbcons);
	printf("input Lincons\n");
	elina_lincons0_array_fprint(stdout,&lincons0,NULL);
	fflush(stdout);
	zo = zonotope_meet_lincons_array(man,true,zo,&lincons0);

	printf("ELINA Output zonotope after meet lincons\n");
	arr = zonotope_to_lincons_array(man,zo);
  	elina_lincons0_array_fprint(stdout,&arr,NULL);
	printf("\n");
  	fflush(stdout);
 	elina_lincons0_array_clear(&arr);


	zonotope_free(man,zo);
	elina_manager_free(man);
	elina_lincons0_array_clear(&lincons0);
}




int main(int argc, char **argv){
	if(argc < 3){
		printf("The test requires two positive integers: (a) Number of variables and (b) Number of constraints");
		return 0;
	}
	size_t dim = atoi(argv[1]);
	size_t nbcons = atoi(argv[2]);
	if(dim <=0 || nbcons <=0){
		printf("The Input parameters should be positive\n");
		return 0;
	}
	printf("Testing zonoml\n");
	test_zonoml(dim,nbcons);
}

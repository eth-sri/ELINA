/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
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
#include "zonotope.h"
#include "zonotope_internal.h"

elina_manager_t * man;

elina_interval_t ** generate_random_box(unsigned short int dim){
	elina_interval_t ** interval = (elina_interval_t **)malloc(dim*sizeof(elina_interval_t *));
	unsigned short int i;
	for(i = 0; i < dim; i++){
		interval[i] = elina_interval_alloc();
		double inf = (double)rand()/RAND_MAX*2.0-1.0;
		double sup = inf + ((double)rand()/RAND_MAX*2.0);
		elina_interval_set_double(interval[i],inf,sup);
	}
	return interval;
}


elina_linexpr0_t * generate_random_linexpr0(unsigned short int dim){
	elina_coeff_t *cst, *coeff;
	unsigned short int j, k;
	elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2);
	cst = &linexpr0->cst;
	elina_scalar_set_to_int(cst->val.scalar,0,ELINA_SCALAR_DOUBLE);
	elina_linterm_t * linterm = &linexpr0->p.linterm[0];
	int v1 = rand()%dim;
	linterm->dim = v1;
	coeff = &linterm->coeff;
	elina_scalar_set_to_int(coeff->val.scalar,1,ELINA_SCALAR_DOUBLE);
	linterm = &linexpr0->p.linterm[1];
	int v2 = rand()%dim;
	while(v2==v1){
		v2 = rand()%dim;
	}
	linterm->dim = v2;
	coeff = &linterm->coeff;
	elina_scalar_set_to_int(coeff->val.scalar,-1,ELINA_SCALAR_DOUBLE);
	return linexpr0;
}


elina_lincons0_array_t generate_random_lincons0_array(unsigned short int dim, size_t nbcons){
	size_t i;
	unsigned short int j, k; 
	elina_coeff_t *cst, *coeff;
	elina_lincons0_array_t  lincons0 = elina_lincons0_array_make(nbcons);
	for(i=0; i < nbcons; i++){
		//lincons0.p[i].constyp = rand() %2 ? ELINA_CONS_SUPEQ : ELINA_CONS_EQ;
		lincons0.p[i].constyp = ELINA_CONS_SUPEQ;
		lincons0.p[i].linexpr0 = generate_random_linexpr0(dim);
	}
	return lincons0;
}

zonotope_t *test_meet_lincons(unsigned short int dim, size_t nbcons){
	elina_interval_t **interval_array = generate_random_box(dim);
	zonotope_t * z1 = zonotope_of_box(man, 0,dim,interval_array);
	//generate random constraints
	elina_lincons0_array_t lincons0 = generate_random_lincons0_array(dim,nbcons);
	//meet lincons
	//printf("ELINA Meet Lincons Inputs\n");
	//elina_lincons0_array_t arr1 = zonotope_to_lincons_array(man,z1);
  	//elina_lincons0_array_fprint(stdout,&arr1,NULL);
	for(unsigned short int i=0; i < dim; i++){
		elina_interval_fprint(stdout,interval_array[i]);
	}
	printf("\n");
	//zonotope_fprint(stdout,man,z1,NULL);
	elina_lincons0_array_fprint(stdout,&lincons0,NULL);
	fflush(stdout);
	//elina_lincons0_array_clear(&arr1);

	z1 = zonotope_meet_lincons_array(man,true,z1,&lincons0);
	
	//printf("ELINA Output Zonotope\n");
	//zonotope_fprint(stdout,man,z1,NULL);
	//elina_lincons0_array_t arr2 = zonotope_to_lincons_array(man,z1);
  	//elina_lincons0_array_fprint(stdout,&arr2,NULL);
  	//fflush(stdout);
	//elina_lincons0_array_clear(&arr2);

	for(unsigned short int i = 0; i < dim; i++){
		elina_interval_free(interval_array[i]);
	}
	
	elina_lincons0_array_clear(&lincons0);
	free(interval_array);
	return z1;
}

void test_join(unsigned short int dim, size_t nbcons){
	zonotope_t *z1 = test_meet_lincons(dim,nbcons);
	zonotope_t *z2 = test_meet_lincons(dim,nbcons);
	printf("Zonotope join inputs\n");
	zonotope_fprint(stdout,man,z1,NULL);
	zonotope_fprint(stdout,man,z2,NULL);
	fflush(stdout);
	zonotope_t *z3 = zonotope_join(man,false,z1,z2);
	printf("Zonotope join output\n");
	zonotope_fprint(stdout,man,z3,NULL);
	fflush(stdout);
	zonotope_free(man,z1);
	zonotope_free(man,z2);
	zonotope_free(man,z3);
}

void test_dana(){
    elina_interval_t **interval = elina_interval_array_alloc(2);
    elina_interval_set_top(interval[0]);
    elina_interval_set_top(interval[1]);
    elina_scalar_set_int(interval[0]->inf,0);
    elina_scalar_set_int(interval[1]->inf,1);
    elina_scalar_set_int(interval[1]->sup,3);
    zonotope_t * z = zonotope_of_box(man,0,2,interval);
    zonotope_fprint(stdout,man,z,NULL);
}

int main(int argc, char **argv){
	srand(time(NULL));
	unsigned short int dim = atoi(argv[1]);
	size_t nbcons = atoi(argv[2]);
	man = zonotope_manager_alloc();
	//test_join(dim,nbcons);
    test_of_box();
	elina_manager_free(man);
}

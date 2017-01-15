/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
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
#include "opt_pk.h"

elina_linexpr0_t * generate_random_linexpr0(unsigned short int dim){
	elina_coeff_t *cst, *coeff;
	unsigned short int j, k;
	elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,dim);
	cst = &linexpr0->cst;
	int r = rand()%10;
	elina_scalar_set_to_int(cst->val.scalar,r,ELINA_SCALAR_MPQ);
	k = 0;
	for(j=0; j < dim/3; j++, k++){
		elina_linterm_t * linterm = &linexpr0->p.linterm[k];
		linterm->dim = j;
		coeff = &linterm->coeff;
		int r = rand()%5;
		elina_scalar_set_to_int(coeff->val.scalar,r,ELINA_SCALAR_MPQ);
	}
	elina_linexpr0_reinit(linexpr0,k);
	return linexpr0;
}

elina_lincons0_array_t generate_random_lincons0_array(unsigned short int dim, size_t nbcons){
	size_t i;
	unsigned short int j, k; 
	elina_coeff_t *cst, *coeff;
	elina_lincons0_array_t  lincons0 = elina_lincons0_array_make(nbcons);
	for(i=0; i < nbcons/3; i++){
		lincons0.p[i].constyp = rand() %2 ? ELINA_CONS_SUPEQ : ELINA_CONS_EQ;
		int r = rand()%10;
		elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,dim);
		cst = &linexpr0->cst;
		elina_scalar_set_to_int(cst->val.scalar,r,ELINA_SCALAR_MPQ);
		k = 0;
		for(j=0; j < dim/3+1; j++, k++){
			elina_linterm_t * linterm = &linexpr0->p.linterm[k];
			linterm->dim = j;
			coeff = &linterm->coeff;
			int r = rand()%5;
			elina_scalar_set_to_int(coeff->val.scalar,r,ELINA_SCALAR_MPQ);
		}
		elina_linexpr0_reinit(linexpr0,k);
		lincons0.p[i].linexpr0 = linexpr0;
	}
	
	for(i=nbcons/3; i < (2*nbcons)/3; i++){
		lincons0.p[i].constyp = ELINA_CONS_SUPEQ;
		int r = rand()%10;
		elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,dim);
		cst = &linexpr0->cst;
		elina_scalar_set_to_int(cst->val.scalar,r,ELINA_SCALAR_MPQ);
		k = 0;
		for(j=dim/3+1; j < (2*dim)/3; j++, k++){
			elina_linterm_t * linterm = &linexpr0->p.linterm[k];
			linterm->dim = j;
			coeff = &linterm->coeff;
			int r = rand()%3;
			elina_scalar_set_to_int(coeff->val.scalar,r,ELINA_SCALAR_MPQ);
		}
		elina_linexpr0_reinit(linexpr0,k);
		lincons0.p[i].linexpr0 = linexpr0;
	}
	
	for(i=(2*nbcons)/3; i < nbcons; i++){
		lincons0.p[i].constyp = ELINA_CONS_SUPEQ;
		int r = rand()%10;
		elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,dim);
		cst = &linexpr0->cst;
		elina_scalar_set_to_int(cst->val.scalar,r,ELINA_SCALAR_MPQ);
		k = 0;
		for(j=(2*dim)/3; j < dim; j++, k++){
			elina_linterm_t * linterm = &linexpr0->p.linterm[k];
			linterm->dim = j;
			coeff = &linterm->coeff;
			int r = rand()%3;
			elina_scalar_set_to_int(coeff->val.scalar,r,ELINA_SCALAR_MPQ);
		}
		elina_linexpr0_reinit(linexpr0,k);
		lincons0.p[i].linexpr0 = linexpr0;
	}
	return lincons0;
}

void test_meetjoin(unsigned short int dim, size_t nbcons, bool meet){
	unsigned short int j,l=1;
	//generate random cosntraints	
	elina_lincons0_array_t lincons1 = generate_random_lincons0_array(dim,nbcons);
	elina_lincons0_array_t lincons2 = generate_random_lincons0_array(dim,nbcons);

	elina_manager_t * man = opt_pk_manager_alloc(false);
	//generate first input
	opt_pk_array_t * oa1 = opt_pk_top(man, dim,0);
	
	//meet with constraints
	opt_pk_array_t * oa2 = opt_pk_meet_lincons_array(man,false,oa1,&lincons1);


	//generate second input
	opt_pk_array_t * oa3 = opt_pk_top(man, dim,0);
	
	//meet with constraints
	opt_pk_array_t * oa4 = opt_pk_meet_lincons_array(man,false,oa3,&lincons2);

	// Print the ELINA result
	printf("ELINA Input Polyhedra\n");
	elina_lincons0_array_t arr1 = opt_pk_to_lincons_array(man,oa2);
  	elina_lincons0_array_fprint(stdout,&arr1,NULL);
	elina_lincons0_array_t arr2 = opt_pk_to_lincons_array(man,oa4);
  	elina_lincons0_array_fprint(stdout,&arr2,NULL);
  	fflush(stdout);
	elina_lincons0_array_clear(&arr1);
	elina_lincons0_array_clear(&arr2);
	// apply fold operation
	opt_pk_array_t * oa5 = meet ? opt_pk_meet(man,false,oa2,oa4) : opt_pk_join(man,false,oa2,oa4);	
	

	
	printf("ELINA Output Polyhedron\n");
	elina_lincons0_array_t arr3 = opt_pk_to_lincons_array(man,oa5);
  	elina_lincons0_array_fprint(stdout,&arr3,NULL);
  	fflush(stdout);
	elina_lincons0_array_clear(&arr3);

 	
	opt_pk_free(man,oa1);
	opt_pk_free(man,oa2);
	opt_pk_free(man,oa3);
	opt_pk_free(man,oa4);
	opt_pk_free(man,oa5);

	elina_manager_free(man);
	
	elina_lincons0_array_clear(&lincons1);
	elina_lincons0_array_clear(&lincons2);
}



void test_fold(unsigned short int dim, size_t nbcons){
	unsigned short int j,l=1;
	//generate random cosntraints	
	elina_lincons0_array_t lincons0 = generate_random_lincons0_array(dim,nbcons);
	//generate tdim
	unsigned short int size = dim/2;
	elina_dim_t * tdim = (elina_dim_t *)malloc(size*sizeof(elina_dim_t));
	tdim[0] = dim/2;
	for(j=dim/2+1; j < dim; j++){
		tdim[l] = j;
		l++;
	}

	//run with ELINA
	elina_manager_t * man = opt_pk_manager_alloc(false);
	opt_pk_array_t * oa1 = opt_pk_top(man, dim,0);
	
	//meet with constraints
	opt_pk_array_t * oa2 = opt_pk_meet_lincons_array(man,false,oa1,&lincons0);


	// Print the ELINA result
	printf("ELINA Input Polyhedron\n");
	elina_lincons0_array_t arr3 = opt_pk_to_lincons_array(man,oa2);
  	elina_lincons0_array_fprint(stdout,&arr3,NULL);
	printf("Dimensions: ");
	for(l=0; l < size; l++){
		printf("%d ",tdim[l]);
	}
	printf("\n");
  	fflush(stdout);
	// apply fold operation
	opt_pk_array_t * oa3 = opt_pk_fold(man,false,oa2,tdim,size);	
	

	
	printf("ELINA Output Polyhedron\n");
	elina_lincons0_array_t arr4 = opt_pk_to_lincons_array(man,oa3);
  	elina_lincons0_array_fprint(stdout,&arr4,NULL);
	printf("\n");
  	fflush(stdout);

 	elina_lincons0_array_clear(&arr3);
	elina_lincons0_array_clear(&arr4);
	opt_pk_free(man,oa1);
	opt_pk_free(man,oa2);
	opt_pk_free(man,oa3);
	elina_manager_free(man);
	
	free(tdim);	
	elina_lincons0_array_clear(&lincons0);
	
}


void test_expand(unsigned short int dim, size_t nbcons){
	unsigned short int j,l=1;
	//generate random cosntraints	
	elina_lincons0_array_t lincons0 = generate_random_lincons0_array(dim,nbcons);
	//generate tdim
	unsigned short int tdim = rand()%dim;
	unsigned short int dimsup = dim/3;
	

	//run with ELINA
	elina_manager_t * man = opt_pk_manager_alloc(false);
	opt_pk_array_t * oa1 = opt_pk_top(man, dim,0);
	
	
	
	//meet with constraints
	opt_pk_array_t * oa2 = opt_pk_meet_lincons_array(man,false,oa1,&lincons0);
	
	printf("ELINA Input Polyhedron\n");
	elina_lincons0_array_t arr3 = opt_pk_to_lincons_array(man,oa2);
  	elina_lincons0_array_fprint(stdout,&arr3,NULL);
	printf("tdim: %d dimsup: %d\n",tdim,dimsup);
  	fflush(stdout);
	// apply fold operation
	opt_pk_array_t * oa3 = opt_pk_expand(man,false,oa2,tdim,dimsup);	
	

	// Print the ELINA result
	printf("ELINA Output Polyhedron\n");
	elina_lincons0_array_t arr4 = opt_pk_to_lincons_array(man,oa3);
  	elina_lincons0_array_fprint(stdout,&arr4,NULL);
	printf("\n");
  	fflush(stdout);

 	elina_lincons0_array_clear(&arr3);
	elina_lincons0_array_clear(&arr4);
	opt_pk_free(man,oa1);
	opt_pk_free(man,oa2);
	opt_pk_free(man,oa3);
	elina_manager_free(man);
		
	elina_lincons0_array_clear(&lincons0);
}

void test_assign(unsigned short int dim, size_t nbcons){
	elina_manager_t * man = opt_pk_manager_alloc(false);
	opt_pk_array_t * oa1 = opt_pk_top(man, dim,0);
	//generate random constraints
	elina_lincons0_array_t lincons0 = generate_random_lincons0_array(dim,nbcons);
	opt_pk_array_t * oa2 = opt_pk_meet_lincons_array(man,false,oa1,&lincons0);

	
	elina_dim_t * tdim = (elina_dim_t *)malloc(sizeof(elina_dim_t));
	tdim[0] = rand()%dim;
	elina_linexpr0_t ** expr_array = (elina_linexpr0_t**)malloc(sizeof(elina_linexpr0_t*));
	elina_linexpr0_t * linexpr0 = generate_random_linexpr0(dim);
	expr_array[0] = linexpr0;
	printf("ELINA Input Polyhedron\n");
	elina_lincons0_array_t arr1 = opt_pk_to_lincons_array(man,oa2);
  	elina_lincons0_array_fprint(stdout,&arr1,NULL);
	printf("Assignment statement\n");
	printf("x%d = ",tdim[0]);
	elina_linexpr0_fprint(stdout,linexpr0,NULL);
	printf("\n");
  	fflush(stdout);
	elina_lincons0_array_clear(&arr1);
	//assign;
	opt_pk_array_t * oa3 = opt_pk_assign_linexpr_array(man,false,oa2,tdim, expr_array,1,NULL);
	elina_linexpr0_free(linexpr0);
	free(expr_array);
	free(tdim);

	//meet with -x1 + 43 >= 0; 	

	// Print the result
	printf("ELINA Output Polyhedron\n");
	elina_lincons0_array_t arr = opt_pk_to_lincons_array(man,oa3);
  	elina_lincons0_array_fprint(stdout,&arr,NULL);
	printf("\n");
  	fflush(stdout);
 	elina_lincons0_array_clear(&arr);


	opt_pk_free(man,oa1);
	opt_pk_free(man,oa2);
	opt_pk_free(man,oa3);
	elina_manager_free(man);
	elina_lincons0_array_clear(&lincons0);
}


int main(int argc, char **argv){
	if(argc < 3){
		printf("The test requires two positive integers: (a) Number of variables and (b) Number of constraints");
		return 0;
	}

	unsigned short int dim = atoi(argv[1]);
	size_t nbcons = atoi(argv[2]);
	if(dim <=0 || nbcons <=0){
		printf("The Input parameters should be positive\n");
		return 0;
	}
	printf("Testing Meet\n");
	test_meetjoin(dim,nbcons,true);
	printf("Testing Join\n");
	test_meetjoin(dim,nbcons,false);
	printf("Testing Assign\n");
	test_assign(dim,nbcons);
	printf("Testing Fold\n");
 	test_fold(6,6);
	printf("Testing Expand\n");
	test_expand(dim,nbcons);
}

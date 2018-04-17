/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "elina_coeff.h"


void test_set_scalar_int(elina_coeff_t *coeff){
	long int num = rand()%100;
	elina_coeff_set_scalar_int(coeff,num);
	printf("set scalar int : %ld coeff: ",num);
	elina_coeff_print(coeff);
	printf(" is zero: %d\n", elina_coeff_zero(coeff));
}


void test_set_scalar_mpq(elina_coeff_t *coeff){
	mpq_t mpq;
	mpq_init(mpq);
	long int p = rand()%1000;
	unsigned long int q= rand()%20+1;
	mpq_set_si(mpq,p,q);
	elina_coeff_set_scalar_mpq(coeff,mpq);
	printf("set scalar mpq : ");
	mpq_out_str(stdout,10,mpq);
	printf(" coeff: ");
	elina_coeff_print(coeff);
	printf(" is zero: %d\n", elina_coeff_zero(coeff));
	mpq_clear(mpq);
}

void test_set_scalar_frac(elina_coeff_t *coeff){
	long int p = rand()%100;
	unsigned long int q = rand()%20+1;
	elina_coeff_set_scalar_frac(coeff,p,q);
	printf("set scalar frac p : %ld q: %ld coeff: ",p,q);
	elina_coeff_print(coeff);
	printf("\n");
}


void test_set_scalar_double(elina_coeff_t *coeff){
	double d = (double)rand()/RAND_MAX*2.0-1.0;
	elina_coeff_set_scalar_double(coeff,d);
	printf("set scalar double : %g coeff: ",d);
	elina_coeff_print(coeff);
	printf(" is zero: %d\n", elina_coeff_zero(coeff));
}


void test_set_scalar_mpfr(elina_coeff_t *coeff){
	mpfr_t mpfr;
	mpfr_init(mpfr);
	double d = (double)rand()/RAND_MAX*2.0-1.0;
	mpfr_set_d(mpfr,d,GMP_RNDU);
	elina_coeff_set_scalar_mpfr(coeff,mpfr);
	printf("set scalar mpfr: ");
	mpfr_out_str(stdout,10,elina_scalar_print_prec,mpfr,GMP_RNDU);
	printf(" coeff: ");
	elina_coeff_print(coeff);
	printf(" is zero: %d\n", elina_coeff_zero(coeff));
	mpfr_clear(mpfr);
}


void test_set_interval_int(elina_coeff_t *coeff){
	long int inf = rand()%100;
	long int sup = rand()%100;
	elina_coeff_set_interval_int(coeff,inf,sup);
	printf("set interval int inf: %ld sup: %ld coeff: ",inf,sup);
	elina_coeff_print(coeff);
	printf(" is zero: %d\n", elina_coeff_zero(coeff));
}


void test_set_interval_mpq(elina_coeff_t *coeff){
	mpq_t inf1,sup1;
	mpq_init(inf1);
	mpq_init(sup1);
	long int p = rand()%1000;
	unsigned long int q= rand()%20+1;
	mpq_set_si(inf1,p,q);
	p = p + rand()%1000;
	mpq_set_si(sup1,p,q);
	elina_coeff_set_interval_mpq(coeff,inf1,sup1);
	printf("set interval mpq inf: ");
	mpq_out_str(stdout,10,inf1);
	printf(" sup: ");
	mpq_out_str(stdout,10,sup1);
	printf(" coeff: ");
	elina_coeff_print(coeff);
	printf(" is zero: %d\n", elina_coeff_zero(coeff));
	mpq_clear(inf1);
	mpq_clear(sup1);
}

void test_set_interval_frac(elina_coeff_t * coeff){
	long int p1 = rand()%100;
	unsigned long int q1 = rand()%20+1;
	long int p2 = rand()%100;
	unsigned long int q2 = rand()%20+1;
	elina_coeff_set_interval_frac(coeff,p1,q1,p2,q2);
	printf("set interval frac p1: %ld q1: %ld p2: %ld q2: %ld coeff: ",p1,q1,p2,q2);
	elina_coeff_print(coeff);
	printf(" is zero: %d\n", elina_coeff_zero(coeff));
}


void test_set_interval_double(elina_coeff_t *coeff){
	double inf = (double)rand()/RAND_MAX*2.0-1.0;
	double sup = (double)rand()/RAND_MAX*2.0-1.0;
	elina_coeff_set_interval_double(coeff,inf,sup);
	printf("set interval double inf: %g sup: %g coeff: ",inf,sup);
	elina_coeff_print(coeff);
	printf(" is zero: %d\n", elina_coeff_zero(coeff));
}


void test_set_interval_mpfr(elina_coeff_t * coeff){
	mpfr_t inf,sup;
	mpfr_init(inf);
	mpfr_init(sup);
	double d = (double)rand()/RAND_MAX*2.0-1.0;
	mpfr_set_d(inf,d,GMP_RNDU);
	
	d = (double)rand()/RAND_MAX*2.0-1.0;
	mpfr_set_d(sup,d,GMP_RNDU);

	elina_coeff_set_interval_mpfr(coeff,inf,sup);
	printf("set interval mpfr: ");
	mpfr_out_str(stdout,10,elina_scalar_print_prec,inf,GMP_RNDU);
	printf(" ");
	mpfr_out_str(stdout,10,elina_scalar_print_prec,sup,GMP_RNDU);
	printf(" coeff: ");
	elina_coeff_print(coeff);
	printf(" is zero: %d\n", elina_coeff_zero(coeff));
	mpfr_clear(inf);
	mpfr_clear(sup);
}	



void test_coeff_cmp(elina_coeff_t * coeff1, elina_coeff_t *coeff2){
	elina_coeff_t *coeff3 = elina_coeff_alloc(ELINA_COEFF_SCALAR);
	elina_coeff_t *coeff4 = elina_coeff_alloc(ELINA_COEFF_INTERVAL);
	int c = rand()%5;
	switch(c){
		case 0:
			test_set_scalar_int(coeff3);
			break;
		case 1:
			test_set_scalar_mpq(coeff3);
			break;
		case 2:
			test_set_scalar_frac(coeff3);
			break;
		case 3:
			test_set_scalar_double(coeff3);
			break;
		default:
			test_set_scalar_mpfr(coeff3);
			break;
	}

	c = rand()%5;
	switch(c){
		case 0:
			test_set_interval_int(coeff4);
			break;
		case 1:
			test_set_interval_mpq(coeff4);
			break;
		case 2:
			test_set_interval_frac(coeff4);
			break;
		case 3:
			test_set_interval_double(coeff4);
			break;
		default:
			test_set_interval_mpfr(coeff4);
			break;
	}

	printf("cmp scalar vs interval coeff: ");
	elina_coeff_print(coeff1);
	printf(" ");
	elina_coeff_print(coeff2);
	printf(" coeff1<=coeff2: %d coeff2<=coeff1: %d\n",elina_coeff_cmp(coeff1,coeff2), elina_coeff_cmp(coeff2,coeff1));


	printf("cmp scalar coeff: ");
	elina_coeff_print(coeff1);
	printf(" ");
	elina_coeff_print(coeff3);
	printf(" coeff1<=coeff2: %d coeff2<=coeff1: %d\n",elina_coeff_cmp(coeff1,coeff3), elina_coeff_cmp(coeff3,coeff1));

	printf("cmp interval coeff: ");
	elina_coeff_print(coeff4);
	printf(" ");
	elina_coeff_print(coeff2);
	printf(" coeff1<=coeff2: %d coeff2<=coeff1: %d\n",elina_coeff_cmp(coeff4,coeff2), elina_coeff_cmp(coeff2,coeff4));

	elina_coeff_free(coeff3);
	elina_coeff_free(coeff4);
}


void test_coeff_equality(elina_coeff_t * coeff1, elina_coeff_t *coeff2){
	
	
	elina_coeff_t *coeff3 = elina_coeff_alloc(ELINA_COEFF_SCALAR);
	elina_coeff_t *coeff4 = elina_coeff_alloc(ELINA_COEFF_INTERVAL);
	int c = rand()%5;
	switch(c){
		case 0:
			test_set_scalar_int(coeff3);
			break;
		case 1:
			test_set_scalar_mpq(coeff3);
			break;
		case 2:
			test_set_scalar_frac(coeff3);
			break;
		case 3:
			test_set_scalar_double(coeff3);
			break;
		default:
			test_set_scalar_mpfr(coeff3);
			break;
	}

	c = rand()%5;
	switch(c){
		case 0:
			test_set_interval_int(coeff4);
			break;
		case 1:
			test_set_interval_mpq(coeff4);
			break;
		case 2:
			test_set_interval_frac(coeff4);
			break;
		case 3:
			test_set_interval_double(coeff4);
			break;
		default:
			test_set_interval_mpfr(coeff4);
			break;
	}

	printf("equal scalar vs interval coeff: ");
	elina_coeff_print(coeff1);
	printf(" ");
	elina_coeff_print(coeff2);
	printf(" coeff1==coeff2: %d\n",elina_coeff_equal(coeff1,coeff2));

	printf("equal scalar coeff: ");
	elina_coeff_print(coeff1);
	printf(" ");
	elina_coeff_print(coeff3);
	printf(" coeff1==coeff2: %d\n",elina_coeff_equal(coeff1,coeff3));

	printf("equal interval coeff: ");
	elina_coeff_print(coeff4);
	printf(" ");
	elina_coeff_print(coeff2);
	printf(" coeff1==coeff2: %d\n",elina_coeff_equal(coeff4,coeff2));

	elina_coeff_free(coeff3);
	elina_coeff_free(coeff4);
}

void test_coeff_neg(elina_coeff_t * coeff1, elina_coeff_t *coeff2){
	elina_coeff_t *coeff3 = elina_coeff_alloc(ELINA_COEFF_SCALAR);
	elina_coeff_t *coeff4 = elina_coeff_alloc(ELINA_COEFF_INTERVAL);
	elina_coeff_neg(coeff3,coeff1);
	printf("scalar coeff: ");
	elina_coeff_print(coeff1);
	printf("neg coeff: ");
	elina_coeff_print(coeff3);
	printf("\n");

	elina_coeff_neg(coeff4,coeff2);
	printf("interval coeff: ");
	elina_coeff_print(coeff2);
	printf("neg coeff: ");
	elina_coeff_print(coeff4);
	printf("\n");

	elina_coeff_free(coeff3);
	elina_coeff_free(coeff4);
}

void test_coeff_reduce(){
	elina_coeff_t * coeff = elina_coeff_alloc(ELINA_COEFF_INTERVAL);
	long int num = rand()%100;
	elina_coeff_set_interval_int(coeff,num,num);
	printf("Before reduce: ");
	elina_coeff_print(coeff);
	printf("\n");
	elina_coeff_reduce(coeff);
	printf("After reduce: ");
	elina_coeff_print(coeff);
	printf("\n");
	elina_coeff_free(coeff);
}





int main(){
	srand (time(NULL));

	elina_coeff_t * coeff1 = elina_coeff_alloc(ELINA_COEFF_SCALAR);

	elina_coeff_t * coeff2 = elina_coeff_alloc(ELINA_COEFF_INTERVAL);

	// coeff is set to scalar of type int
	test_set_scalar_int(coeff1);
	
	//coeff is set to scalar of type mpq
	test_set_scalar_mpq(coeff1);
	
	//coeff is set to scalar of type fraction of two integers
	test_set_scalar_frac(coeff1);
	
	//coeff is set to scalar of type double
	test_set_scalar_double(coeff1);
	
	//coeff is set to scalar of type mpfr
	test_set_scalar_mpfr(coeff1);	
	
	// coeff is set to interval of type int
	test_set_interval_int(coeff2);
	
	//coeff is set to interval of type mpq
	test_set_interval_mpq(coeff2);
	
	//coeff is set to interval of type fraction of two integers
	test_set_interval_frac(coeff2);
	
	//coeff is set to interval of type double
	test_set_interval_double(coeff2);
	
	//coeff is set to interval of type mpfr
	test_set_interval_mpfr(coeff2);

	//set to coeff
	//test_set_coeff(coeff1, coeff2);	
	
	//test for comparison and equality 
	test_coeff_cmp(coeff1, coeff2);
	test_coeff_equality(coeff1, coeff2);
	
	
	//test for negation and reduce
	test_coeff_neg(coeff1, coeff2);
	test_coeff_reduce();


	elina_coeff_free(coeff1);
	elina_coeff_free(coeff2);
}

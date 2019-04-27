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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "elina_scalar.h"

elina_scalar_t * scalar1=NULL, *scalar2=NULL;

void test_set_int(){
	long int num = rand()%100;
	elina_scalar_set_int(scalar1,num);
	printf("set int : %ld scalar: ",num);
	elina_scalar_print(scalar1);
	printf("\n");
}

void test_cmp_int(){
	int r = rand()%100;
	printf("cmp int scalar: ");
	elina_scalar_print(scalar1);
	printf("int: %d scalar<=int: %d scalar==int: %d \n",r,elina_scalar_cmp_int(scalar1,r), elina_scalar_equal_int(scalar1,r));
}

void test_set_mpq(){
	mpq_t mpq;
	mpq_init(mpq);
	long int p = rand()%1000;
	unsigned long int q= rand()%20;
	mpq_set_si(mpq,p,q);
	elina_scalar_set_mpq(scalar1,mpq);
	printf("set mpq : ");
	mpq_out_str(stdout,10,mpq);
	printf(" scalar: ");
	elina_scalar_print(scalar1);
	printf("\n");
	mpq_clear(mpq);
}

void test_set_frac(){
	long int p = rand()%100;
	unsigned long int q = rand()%20;
	elina_scalar_set_frac(scalar1,p,q);
	printf("set frac p : %ld q: %ld scalar: ",p,q);
	elina_scalar_print(scalar1);
	printf("\n");
}


void test_set_double(){
	double d = (double)rand()/RAND_MAX*2.0-1.0;
	elina_scalar_set_double(scalar1,d);
	printf("set double : %g scalar: ",d);
	elina_scalar_print(scalar1);
	printf("\n");
}


void test_set_mpfr(){
	mpfr_t mpfr;
	mpfr_init(mpfr);
	double d = (double)rand()/RAND_MAX*2.0-1.0;
	mpfr_set_d(mpfr,d,GMP_RNDU);
	elina_scalar_set_mpfr(scalar1,mpfr);
	printf("set mpfr: ");
	mpfr_out_str(stdout,10,elina_scalar_print_prec,mpfr,GMP_RNDU);
	printf(" scalar: ");
	elina_scalar_print(scalar1);
	printf("\n");
	mpfr_clear(mpfr);
}	

void test_mpq_set(){
	mpq_t mpq;
	mpq_init(mpq);
	elina_mpq_set_scalar(mpq,scalar2,GMP_RNDU);
	printf("mpq_set scalar: ");
	elina_scalar_print(scalar2);
	printf(" mpq: ");
	mpq_out_str(stdout,10,mpq);
	printf("\n");
	mpq_clear(mpq);
}

void test_double_set(){
	elina_scalar_inv(scalar2,scalar1);
	double d;
	elina_double_set_scalar(&d,scalar2,GMP_RNDU);
	printf("double_set scalar: ");
	elina_scalar_print(scalar2);
	printf(" double: %g\n",d);
}

void test_mpfr_set(){
	elina_scalar_inv(scalar2,scalar1);
	mpfr_t mpfr;
	mpfr_init(mpfr);
	elina_mpfr_set_scalar(mpfr,scalar2,GMP_RNDU);
	printf("mpfr_set scalar: ");
	elina_scalar_print(scalar2);
	printf(" mpfr: ");
	mpfr_out_str(stdout,10,elina_scalar_print_prec,mpfr,GMP_RNDU);
	printf("\n");
	mpfr_clear(mpfr);
}

void test_set_scalar(){
	elina_scalar_set(scalar2,scalar1);
	printf("set scalar1: ");
	elina_scalar_print(scalar1);
	printf(" scalar2: ");
	elina_scalar_print(scalar2);
	printf(" scalar1==scalar2: %d\n",elina_scalar_equal(scalar1,scalar2));
}

void test_inv(){
	scalar2 = elina_scalar_alloc_set(scalar1);
	elina_scalar_inv(scalar2,scalar1);
	printf("inversion scalar1: ");
	elina_scalar_print(scalar1);
	printf("scalar2: ");
	elina_scalar_print(scalar2);
	printf(" scalar1<= scalar2: %d \n", elina_scalar_cmp(scalar1,scalar2));
}

void test_equality(){
	printf("equality scalar1: ");
	elina_scalar_print(scalar1);
	printf(" scalar2: ");
	elina_scalar_print(scalar2);
	printf(" scalar1==scalar2: %d\n", elina_scalar_equal(scalar1,scalar2));
}

void test_sgn(){
	elina_scalar_neg(scalar2,scalar1);
	printf("scalar1: ");
	elina_scalar_print(scalar1);
	printf(" scalar2: ");
	elina_scalar_print(scalar2);
	printf(" sgn(scalar1): %d sgn(scalar2): %d \n", elina_scalar_sgn(scalar1),elina_scalar_sgn(scalar2));
}

void test_swap(){
	elina_scalar_swap(scalar1,scalar2);
	printf("swap scalar1: ");
	elina_scalar_print(scalar1);
	printf("scalar2: ");
	elina_scalar_print(scalar2);
	printf("\n");
}

void test_infty(){
	elina_scalar_set_infty(scalar1,1);
	elina_scalar_set_infty(scalar2,-1);
	printf("infty scalar1: ");
	elina_scalar_print(scalar1);
	printf(" scalar2: ");
	elina_scalar_print(scalar2);
	printf(" isinfty(scalar1): %d isinfty(scalar2): %d\n", elina_scalar_infty(scalar1), elina_scalar_infty(scalar2));
	printf(" scalar1<=scalar2: %d \n", elina_scalar_cmp(scalar1,scalar2));
}

int main(){
	time(NULL);
	scalar1 = elina_scalar_alloc();

	// scalar is set to int
	test_set_int();
	
	//compare with an integer
	test_cmp_int();

	//scalar is set to mpq
	test_set_mpq();
	test_inv();	
	test_mpq_set();
	
	//scalar is set to fractions
	test_set_frac();
	
	//scalar is set to double
	test_set_double();
	test_double_set();

	//scalar is set to mpfr
	test_set_mpfr();	
	test_mpfr_set();

	//set to scalar
	test_set_scalar();	
	
	//test for equality and checking the sign of scalar
	test_equality();
	test_sgn();
	
	//test for swapping
	test_swap();
	
	//tests with infinities
	test_infty();	

	elina_scalar_free(scalar2);
	elina_scalar_free(scalar1);
}

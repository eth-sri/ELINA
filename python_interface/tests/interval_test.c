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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "elina_interval.h"





void test_set_int(elina_interval_t *interval1){
	long int inf = rand()%100;
	long int sup = rand()%100;
	elina_interval_set_int(interval1,inf,sup);
	printf("set int inf: %ld sup: %ld interval: ",inf,sup);
	elina_interval_print(interval1);
	printf(" is bottom: %d is top: %d\n", elina_interval_is_bottom(interval1),elina_interval_is_top(interval1));
}


void test_set_mpq(elina_interval_t *interval1){
	mpq_t inf1,sup1;
	mpq_init(inf1);
	mpq_init(sup1);
	long int p = rand()%1000;
	unsigned long int q= rand()%20+1;
	mpq_set_si(inf1,p,q);
	p = p + rand()%1000;
	mpq_set_si(sup1,p,q);
	elina_interval_set_mpq(interval1,inf1,sup1);
	printf("set mpq inf: ");
	mpq_out_str(stdout,10,inf1);
	printf(" sup: ");
	mpq_out_str(stdout,10,sup1);
	printf(" interval: ");
	elina_interval_print(interval1);
	printf(" is bottom: %d is top: %d\n", elina_interval_is_bottom(interval1),elina_interval_is_top(interval1));
	mpq_clear(inf1);
	mpq_clear(sup1);
}

void test_set_frac(elina_interval_t * interval1){
	long int p1 = rand()%100;
	unsigned long int q1 = rand()%20+1;
	long int p2 = rand()%100;
	unsigned long int q2 = rand()%20+1;
	elina_interval_set_frac(interval1,p1,q1,p2,q2);
	printf("set frac p1: %ld q1: %ld p2: %ld q2: %ld interval: ",p1,q1,p2,q2);
	elina_interval_print(interval1);
	printf(" is bottom: %d is top: %d\n", elina_interval_is_bottom(interval1),elina_interval_is_top(interval1));
}


void test_set_double(elina_interval_t *interval1){
	double inf = (double)rand()/RAND_MAX*2.0-1.0;
	double sup = (double)rand()/RAND_MAX*2.0-1.0;
	elina_interval_set_double(interval1,inf,sup);
	printf("set double inf: %g sup: %g interval: ",inf,sup);
	elina_interval_print(interval1);
	printf(" is bottom: %d is top: %d\n", elina_interval_is_bottom(interval1),elina_interval_is_top(interval1));
}


void test_set_mpfr(elina_interval_t * interval1){
	mpfr_t inf,sup;
	mpfr_init(inf);
	mpfr_init(sup);
	double d = (double)rand()/RAND_MAX*2.0-1.0;
	mpfr_set_d(inf,d,GMP_RNDU);
	
	d = (double)rand()/RAND_MAX*2.0-1.0;
	mpfr_set_d(sup,d,GMP_RNDU);

	elina_interval_set_mpfr(interval1,inf,sup);
	printf("set mpfr: ");
	mpfr_out_str(stdout,10,elina_scalar_print_prec,inf,GMP_RNDU);
	printf(" ");
	mpfr_out_str(stdout,10,elina_scalar_print_prec,sup,GMP_RNDU);
	printf(" interval: ");
	elina_interval_print(interval1);
	printf(" is bottom: %d is top: %d\n", elina_interval_is_bottom(interval1),elina_interval_is_top(interval1));
	mpfr_clear(inf);
	mpfr_clear(sup);
}	



void test_cmp(elina_interval_t * interval1, elina_interval_t *interval2){
	long int inf2 = rand()%100;
	long int sup2 = inf2+ rand()%100; 
	elina_interval_set_int(interval2,inf2,sup2);
	printf("cmp interval: ");
	elina_interval_print(interval1);
	printf(" ");
	elina_interval_print(interval2);
	printf(" interval1<=interval2: %d interval1==interval2: %d \n",elina_interval_cmp(interval1,interval2), elina_interval_equal(interval1,interval2));
}


void test_set_interval(elina_interval_t * interval1, elina_interval_t *interval2){
	elina_interval_set(interval2,interval1);
	printf("set interval1: ");
	elina_interval_print(interval1);
	printf(" interval2: ");
	elina_interval_print(interval2);
	printf(" interval1==interval2: %d\n",elina_interval_equal(interval1,interval2));
}



void test_equality(elina_interval_t * interval1, elina_interval_t *interval2){
	elina_interval_set_bottom(interval1);
	printf("equality interval1: ");
	elina_interval_print(interval1);
	elina_interval_set_top(interval2);
	printf(" interval2: ");
	elina_interval_print(interval2);
	printf(" interval1==interval2: %d\n", elina_interval_equal(interval1,interval2));
}



void test_swap(elina_interval_t * interval1, elina_interval_t *interval2){
	elina_interval_swap(interval1,interval2);
	printf("swap interval1: ");
	elina_interval_print(interval1);
	printf("interval2: ");
	elina_interval_print(interval2);
	printf("interval1 is bottom: %d is top: %d\n", elina_interval_is_bottom(interval1),elina_interval_is_top(interval1));
	printf("interval2 is bottom: %d is top: %d\n", elina_interval_is_bottom(interval2),elina_interval_is_top(interval2));
}

void test_neg(elina_interval_t * interval1, elina_interval_t *interval2){
	elina_interval_neg(interval1,interval2);
	printf("neg interval1: ");
	elina_interval_print(interval1);
	printf("interval2: ");
	elina_interval_print(interval2);
	printf("interval1 is bottom: %d is top: %d\n", elina_interval_is_bottom(interval1),elina_interval_is_top(interval1));
	printf("interval2 is bottom: %d is top: %d\n", elina_interval_is_bottom(interval2),elina_interval_is_top(interval2));
}


void test_interval_array(){
	
	size_t size= rand() %100;
	elina_interval_t ** interval_array = elina_interval_array_alloc(size);
	for(size_t i =0; i < size; i++){
		int c = rand()%5;
		switch(c){
			case 0:
				test_set_int(interval_array[i]);
				break;
			case 1:
				test_set_mpq(interval_array[i]);
				break;
			case 2:
				test_set_frac(interval_array[i]);
				break;
			case 3:
				test_set_double(interval_array[i]);
				break;
			default:
				test_set_mpfr(interval_array[i]);
				break;
		}
	}
	
	printf("Interval array of size: %lu\n",size);
	for(size_t i=0; i < size; i++){
		elina_interval_print(interval_array[i]);
		printf(" ");
	}
	printf("\n");
	elina_interval_array_free(interval_array,size);
}


int main(){
	srand (time(NULL));

	elina_interval_t * interval1 = elina_interval_alloc();

	elina_interval_t * interval2 = elina_interval_alloc();

	// interval is set to int
	test_set_int(interval1);
	
	

	//interval is set to mpq
	test_set_mpq(interval1);
	
	
	
	//interval is set to fractions
	test_set_frac(interval1);
	
	//interval is set to double
	test_set_double(interval1);
	

	//interval is set to mpfr
	test_set_mpfr(interval1);	
	

	//set to interval
	test_set_interval(interval1, interval2);	
	
	//test for comparison and equality 
	test_cmp(interval1, interval2);
	test_equality(interval1, interval2);
	
	
	//test for swapping and negation
	test_swap(interval1, interval2);
	test_neg(interval1, interval2);
	
	test_interval_array();

	elina_interval_free(interval2);
	elina_interval_free(interval1);
}

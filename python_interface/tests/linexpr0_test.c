#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "elina_linexpr0.h"


void test_set_linexpr_scalar_int(elina_linexpr0_t *linexpr){
	size_t size = elina_linexpr0_size(linexpr);
	size_t new_size = size/2;
	if(!new_size){
		new_size = 1;
	}
	elina_linexpr0_realloc(linexpr,new_size);
	size_t i;
	int num = rand()%100;
	elina_linexpr0_set_cst_scalar_int(linexpr,num);
	for(i=0; i < new_size; i++){
		num = rand()%100;
		elina_linexpr0_set_coeff_scalar_int(linexpr,i,num);
	}
	printf("set scalar int size : %ld linexpr: ",new_size);
	elina_linexpr0_print(linexpr,NULL);
	printf(" is linear: %d is quasilinear: %d\n", elina_linexpr0_is_linear(linexpr), elina_linexpr0_is_quasilinear(linexpr));
	size_t intdim = rand()%new_size;
	size_t realdim = rand()%new_size;
	printf("dim: %zu is integer: %d dim: %zu is real: %d\n",intdim,elina_linexpr0_is_integer(linexpr,intdim),realdim,elina_linexpr0_is_real(linexpr,realdim));
}


void test_set_linexpr_scalar_frac(elina_linexpr0_t *linexpr){
	size_t size = elina_linexpr0_size(linexpr);
	size_t new_size = size/2;
	if(!new_size){
		new_size = 1;
	}
	elina_linexpr0_realloc(linexpr,new_size);
	size_t i;
	
	long int p = rand()%1000;
	unsigned long int q= rand()%20+1;
	elina_linexpr0_set_cst_scalar_frac(linexpr,p,q);
	for(i=0; i < new_size; i++){
		p = rand()%1000;
		q= rand()%20+1;
		elina_linexpr0_set_coeff_scalar_frac(linexpr,i,p,q);
	}
	printf("set scalar frac size : %ld linexpr: ",new_size);
	elina_linexpr0_print(linexpr,NULL);
	printf(" is linear: %d is quasilinear: %d\n", elina_linexpr0_is_linear(linexpr), elina_linexpr0_is_quasilinear(linexpr));
	size_t intdim = rand()%(new_size+1);
	size_t realdim = rand()%(new_size+1);
	printf("dim: %zu is integer: %d dim: %zu is real: %d\n",intdim,elina_linexpr0_is_integer(linexpr,intdim),realdim,elina_linexpr0_is_real(linexpr,realdim));
}


void test_set_linexpr_scalar_double(elina_linexpr0_t *linexpr){
	size_t size = elina_linexpr0_size(linexpr);
	size_t new_size = size/2;
	if(!new_size){
		new_size = 1;
	}
	elina_linexpr0_realloc(linexpr,new_size);
	size_t i;
	double d = (double)rand()/RAND_MAX*2.0-1.0;
	elina_linexpr0_set_cst_scalar_double(linexpr,d);
	for(i=0; i < new_size; i++){
		d = (double)rand()/RAND_MAX*2.0-1.0;
		elina_linexpr0_set_coeff_scalar_double(linexpr,i,d);
	}
	printf("set scalar double size : %ld linexpr: ",new_size);
	elina_linexpr0_print(linexpr,NULL);
	printf(" is linear: %d is quasilinear: %d\n", elina_linexpr0_is_linear(linexpr), elina_linexpr0_is_quasilinear(linexpr));
	size_t intdim = rand()%new_size;
	size_t realdim = rand()%new_size;
	printf("dim: %zu is integer: %d dim: %zu is real: %d\n",intdim,elina_linexpr0_is_integer(linexpr,intdim),realdim,elina_linexpr0_is_real(linexpr,realdim));
}



void test_set_linexpr_interval_int(elina_linexpr0_t *linexpr){
	size_t size = elina_linexpr0_size(linexpr);
	size_t new_size = size/2;
	if(!new_size){
		new_size = 1;
	}
	elina_linexpr0_realloc(linexpr,new_size);
	size_t i;
	int inf = rand()%100;
	int sup = inf + rand()%100;
	elina_linexpr0_set_cst_interval_int(linexpr,inf,sup);
	for(i=0; i < new_size; i++){
		inf = rand()%100;
		sup = inf + rand()%100;
		elina_linexpr0_set_coeff_interval_int(linexpr,i,inf,sup);
	}

	printf("set interval int size : %ld linexpr: ",new_size);
	elina_linexpr0_print(linexpr,NULL);
	printf(" is linear: %d is quasilinear: %d\n", elina_linexpr0_is_linear(linexpr), elina_linexpr0_is_quasilinear(linexpr));
	size_t intdim = rand()%new_size;
	size_t realdim = rand()%new_size;
	printf("dim: %zu is integer: %d dim: %zu is real: %d\n",intdim,elina_linexpr0_is_integer(linexpr,intdim),realdim,elina_linexpr0_is_real(linexpr,realdim));
}



void test_set_linexpr_interval_frac(elina_linexpr0_t *linexpr){
	size_t size = elina_linexpr0_size(linexpr);
	size_t new_size = size/2;
	if(!new_size){
		new_size = 1;
	}
	elina_linexpr0_realloc(linexpr,new_size);
	size_t i;
	int numinf = rand()%100;
	unsigned short int deninf = rand()%20+1;	

	int numsup = numinf + rand()%100;
	unsigned short int densup = deninf;
	elina_linexpr0_set_cst_interval_frac(linexpr,numinf,deninf,numsup,densup);
	for(i=0; i < new_size; i++){
		numinf = rand()%100;
		deninf = rand()%20+1;
		
		numsup = numinf + rand()%100;
		densup = deninf;

		elina_linexpr0_set_coeff_interval_frac(linexpr,i,numinf,deninf,numsup,densup);
	}

	printf("set interval frac size : %ld linexpr: ",new_size);
	elina_linexpr0_print(linexpr,NULL);
	printf(" is linear: %d is quasilinear: %d\n", elina_linexpr0_is_linear(linexpr), elina_linexpr0_is_quasilinear(linexpr));
	size_t intdim = rand()%new_size;
	size_t realdim = rand()%new_size;
	printf("dim: %zu is integer: %d dim: %zu is real: %d\n",intdim,elina_linexpr0_is_integer(linexpr,intdim),realdim,elina_linexpr0_is_real(linexpr,realdim));
}


void test_set_linexpr_interval_double(elina_linexpr0_t *linexpr){
	size_t size = elina_linexpr0_size(linexpr);
	size_t new_size = size/2;
	if(!new_size){
		new_size = 1;
	}
	elina_linexpr0_realloc(linexpr,new_size);
	size_t i;
	double inf = (double)rand()/RAND_MAX*2.0-1.0;
	double sup = inf + (double)rand()/RAND_MAX*2.0-1.0;
	elina_linexpr0_set_cst_interval_double(linexpr,inf,sup);
	for(i=0; i < new_size; i++){
		inf = (double)rand()/RAND_MAX*2.0-1.0;
		sup = inf + (double)rand()/RAND_MAX*2.0-1.0;
		elina_linexpr0_set_coeff_interval_double(linexpr,i,inf,sup);
	}

	printf("set interval double size : %ld linexpr: ",new_size);
	elina_linexpr0_print(linexpr,NULL);
	printf(" is linear: %d is quasilinear: %d\n", elina_linexpr0_is_linear(linexpr), elina_linexpr0_is_quasilinear(linexpr));
	size_t intdim = rand()%new_size;
	size_t realdim = rand()%new_size;
	printf("dim: %zu is integer: %d dim: %zu is real: %d\n",intdim,elina_linexpr0_is_integer(linexpr,intdim),realdim,elina_linexpr0_is_real(linexpr,realdim));
}




void test_linexpr_compare(elina_linexpr0_t *linexpr1, elina_linexpr0_t *linexpr2,
			  elina_linexpr0_t *linexpr3, elina_linexpr0_t *linexpr4){

	printf("linexpr1: ");
	elina_linexpr0_print(linexpr1,NULL);
	printf("\n");

	printf("linexpr2: ");
	elina_linexpr0_print(linexpr2,NULL);
	printf("\n");
	
	printf("linexpr3: ");
	elina_linexpr0_print(linexpr3,NULL);
	printf("\n");

	printf("linexpr4: ");
	elina_linexpr0_print(linexpr4,NULL);
	printf("\n");

	printf(" linexpr1<=linexpr2: %d linexpr2<=linexpr1: %d\n",elina_linexpr0_compare(linexpr1,linexpr2), elina_linexpr0_compare(linexpr2,linexpr1));
	printf(" linexpr1<=linexpr3: %d linexpr3<=linexpr1: %d\n",elina_linexpr0_compare(linexpr1,linexpr3), elina_linexpr0_compare(linexpr3,linexpr1));
	printf(" linexpr1<=linexpr4: %d linexpr4<=linexpr1: %d\n",elina_linexpr0_compare(linexpr1,linexpr4), elina_linexpr0_compare(linexpr4,linexpr1));


	printf(" linexpr2<=linexpr3: %d linexpr3<=linexpr2: %d\n",elina_linexpr0_compare(linexpr2,linexpr3), elina_linexpr0_compare(linexpr3,linexpr3));
	printf(" linexpr2<=linexpr4: %d linexpr4<=linexpr2: %d\n",elina_linexpr0_compare(linexpr2,linexpr4), elina_linexpr0_compare(linexpr4,linexpr2));
	
	printf(" linexpr3<=linexpr4: %d linexpr4<=linexpr3: %d\n",elina_linexpr0_compare(linexpr3,linexpr4), elina_linexpr0_compare(linexpr4,linexpr3));
	
}


void test_linexpr_equality(elina_linexpr0_t *linexpr1, elina_linexpr0_t *linexpr2,
			  elina_linexpr0_t *linexpr3, elina_linexpr0_t *linexpr4){

	printf("linexpr1: ");
	elina_linexpr0_print(linexpr1,NULL);
	printf("\n");

	printf("linexpr2: ");
	elina_linexpr0_print(linexpr2,NULL);
	printf("\n");
	
	printf("linexpr3: ");
	elina_linexpr0_print(linexpr3,NULL);
	printf("\n");

	printf("linexpr4: ");
	elina_linexpr0_print(linexpr4,NULL);
	printf("\n");

	printf(" linexpr1==linexpr2: %d\n",elina_linexpr0_equal(linexpr1,linexpr2));
	printf(" linexpr1==linexpr3: %d\n",elina_linexpr0_equal(linexpr1,linexpr3));
	printf(" linexpr1==linexpr4: %d\n",elina_linexpr0_equal(linexpr1,linexpr4));


	printf(" linexpr2==linexpr3: %d \n",elina_linexpr0_equal(linexpr2,linexpr3));
	printf(" linexpr2==linexpr4: %d \n",elina_linexpr0_equal(linexpr2,linexpr4));
	
	printf(" linexpr3==linexpr4: %d \n",elina_linexpr0_equal(linexpr3,linexpr4));
	
}



int main(){
	srand (time(NULL));
	size_t size = rand()%20+3;

	elina_linexpr0_t * linexpr1 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,size);
	elina_linexpr0_t * linexpr2 = elina_linexpr0_alloc(ELINA_LINEXPR_DENSE,size);
	elina_linexpr0_t * linexpr3 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,size);
	elina_linexpr0_t * linexpr4 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,size);

	// create sparse and dense linexpr of type scalar int
	test_set_linexpr_scalar_int(linexpr1);
	fflush(stdout);
	test_set_linexpr_scalar_int(linexpr2);
	fflush(stdout);
	
	//create sparse and dense linexpr of type scalar fraction of two ints
	test_set_linexpr_scalar_frac(linexpr1);
	fflush(stdout);
	test_set_linexpr_scalar_frac(linexpr2);
	fflush(stdout);
	
	//create sparse and dense linexpr of type scalar double
	test_set_linexpr_scalar_double(linexpr1);
	fflush(stdout);
	test_set_linexpr_scalar_double(linexpr2);
	fflush(stdout);

	// create sparse and dense linexpr of type interval int
	test_set_linexpr_interval_int(linexpr3);
	fflush(stdout);
	test_set_linexpr_interval_int(linexpr4);
	fflush(stdout);

	
	//create sparse and dense linexpr of type interval fraction of two ints
	test_set_linexpr_interval_frac(linexpr3);
	fflush(stdout);
	test_set_linexpr_interval_frac(linexpr4);
	fflush(stdout);
	
	//create sparse and dense linexpr of type interval double
	test_set_linexpr_interval_double(linexpr3);
	fflush(stdout);
	test_set_linexpr_interval_double(linexpr4);
	fflush(stdout);

	
	//test for comparison and equality 
	test_linexpr_compare(linexpr1, linexpr2, linexpr3, linexpr4);
	fflush(stdout);
	test_linexpr_equality(linexpr1, linexpr2, linexpr3, linexpr4);
	fflush(stdout);
	

	elina_linexpr0_free(linexpr1);
	elina_linexpr0_free(linexpr2);
	elina_linexpr0_free(linexpr3);
	elina_linexpr0_free(linexpr4);
}

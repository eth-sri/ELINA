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
#include "elina_lincons0.h"


elina_lincons0_array_t generate_random_lincons0_array(unsigned short int dim, size_t nbcons){
	size_t i;
	unsigned short int j, k; 
	elina_coeff_t *cst, *coeff;
	elina_lincons0_array_t  lincons0 = elina_lincons0_array_make(nbcons);
//	elina_lincons0_array_fprint(stdout, &lincons0, NULL);

	for(i=0; i < nbcons/3; i++){
		lincons0.p[i].constyp = rand() %2 ? ELINA_CONS_SUPEQ : ELINA_CONS_EQ;
		double d = (double)rand()/RAND_MAX*2.0-1.0;
		elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,dim);
		cst = &linexpr0->cst;
		elina_scalar_set_double(cst->val.scalar,d);
		k = 0;
		for(j=0; j < dim/3+1; j++, k++){
			elina_linterm_t * linterm = &linexpr0->p.linterm[k];
			linterm->dim = j;
			coeff = &linterm->coeff;
			double d = (double)rand()/RAND_MAX*2.0-1.0;
			elina_scalar_set_double(coeff->val.scalar,d);
		}
		elina_linexpr0_realloc(linexpr0,k);
		lincons0.p[i].linexpr0 = linexpr0;
	}
	
	for(i=nbcons/3; i < (2*nbcons)/3; i++){
		lincons0.p[i].constyp = ELINA_CONS_SUPEQ;
		double d = (double)rand()/RAND_MAX*2.0-1.0;
		elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,dim);
		cst = &linexpr0->cst;
		elina_scalar_set_double(cst->val.scalar,d);
		k = 0;
		for(j=dim/3+1; j < (2*dim)/3; j++, k++){
			elina_linterm_t * linterm = &linexpr0->p.linterm[k];
			linterm->dim = j;
			coeff = &linterm->coeff;
			double d = (double)rand()/RAND_MAX*2.0-1.0;
			elina_scalar_set_double(coeff->val.scalar,d);
		}
		elina_linexpr0_realloc(linexpr0,k);
		lincons0.p[i].linexpr0 = linexpr0;
	}
	
	for(i=(2*nbcons)/3; i < nbcons; i++){
		lincons0.p[i].constyp = ELINA_CONS_SUPEQ;
		double d = (double)rand()/RAND_MAX*2.0-1.0;
		elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,dim);
		cst = &linexpr0->cst;
		elina_scalar_set_double(cst->val.scalar,d);
		k = 0;
		for(j=(2*dim)/3; j < dim; j++, k++){
			elina_linterm_t * linterm = &linexpr0->p.linterm[k];
			linterm->dim = j;
			coeff = &linterm->coeff;
			double d = (double)rand()/RAND_MAX*2.0-1.0;
			elina_scalar_set_double(coeff->val.scalar,d);
		}
		elina_linexpr0_realloc(linexpr0,k);
		lincons0.p[i].linexpr0 = linexpr0;
	}
	return lincons0;
}


int main(){
	srand (time(NULL));
	unsigned short int dim = rand()%50+3;
	size_t nbcons = rand()%50+3;
	
	elina_lincons0_array_t array = generate_random_lincons0_array(dim,nbcons);
	printf("Lincons array\n");
	elina_lincons0_array_fprint(stdout,&array,NULL);
	printf("is linear: %d is quasilinear: %d type of array: %d\n",elina_lincons0_array_is_linear(&array),elina_lincons0_array_is_quasilinear(&array),
		elina_lincons0_array_type(&array));
	fflush(stdout);
	
	elina_dimperm_t* perm = elina_dimperm_alloc(dim);
	for(unsigned short int i =0; i < dim; i++){
		perm->dim[i] = (i+4)%dim;
	}
	
	elina_lincons0_array_t perm_array = elina_lincons0_array_permute_dimensions(&array,perm);
	printf("permutation\n");
	elina_dimperm_fprint(stdout,perm);
	printf("permuted array\n");
	elina_lincons0_array_fprint(stdout,&perm_array,NULL);
	fflush(stdout);

	size_t intdim = rand()%5+1;
	size_t realdim = rand()%5+1;
	elina_dimchange_t * dimchange = elina_dimchange_alloc(intdim, realdim);
	for(size_t i =0; i < intdim+realdim; i++){
		dimchange->dim[i] = rand()%dim;
	}
	elina_lincons0_array_t add_array = elina_lincons0_array_add_dimensions(&array,dimchange);
	printf("dimension add array\n");
	elina_dimchange_fprint(stdout,dimchange);
	printf("array after adding dimension\n");
	elina_lincons0_array_fprint(stdout,&add_array,NULL);
	fflush(stdout);


	elina_dimperm_free(perm);
	elina_dimchange_free(dimchange);
	elina_lincons0_array_clear(&array);
	elina_lincons0_array_clear(&perm_array);
	elina_lincons0_array_clear(&add_array);
}


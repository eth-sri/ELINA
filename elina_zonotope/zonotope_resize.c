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

#include "zonotope_resize.h"


/*********************/
/* Add dimensions */
/*********************/
zonotope_t* zonotope_add_dimensions(elina_manager_t* man, bool destructive, zonotope_t* z, elina_dimchange_t* dimchange, bool project)
{
  //printf("add input %d\n",destructive);
	//zonotope_fprint(stdout,man,z,NULL);
	//fflush(stdout);
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ADD_DIMENSIONS);

    zonotope_t* res = destructive ? z : zonotope_copy(man, z);
    size_t intdim = z->intdim + dimchange->intdim;
    size_t dims = z->dims + (dimchange->intdim + dimchange->realdim);
    res->box = realloc(res->box,dims*sizeof(elina_interval_t *));
    res->paf = realloc(res->paf,dims*sizeof(zonotope_aff_t*));
    size_t i = 0;
    int j = 0;
    for (i=0; i<dimchange->intdim + dimchange->realdim; i++) {
	if (res->dims == dimchange->dim[i]) {
	    /* add in the last place */
	    res->box[res->dims] = elina_interval_alloc();
	} else {
	    /* increment */
	    for (j=(int)-1+res->dims;j>=(int)dimchange->dim[i];j--) {
		res->box[j+1] = elina_interval_alloc();
		res->paf[j+1] = res->paf[j];
		elina_interval_set(res->box[j+1],res->box[j]);
	    }
	}
	res->paf[dimchange->dim[i]] = project ? zonotope_aff_alloc_init(pr) : pr->top;
	res->paf[dimchange->dim[i]]->pby++;
	if (project) elina_interval_set_int(res->box[dimchange->dim[i]], 0,0);
	else elina_interval_set_top(res->box[dimchange->dim[i]]);
	res->dims++;
    }
    res->intdim = intdim;
	//printf("add output\n");
	//zonotope_fprint(stdout,man,res,NULL);
	//fflush(stdout);
    return res;
}


 
/*********************/
/* Remove dimensions, Optimized for AI^2, Linear cost compared with quadratic with APRON */
/*********************/
zonotope_t* zonotope_remove_dimensions(elina_manager_t* man, bool destructive, zonotope_t* z, elina_dimchange_t* dimchange)
{
   
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_REMOVE_DIMENSIONS);
     //printf("remove input %d\n",destructive);
	//zonotope_fprint(stdout,man,z,NULL);
	//fflush(stdout);
    //zonotope_t* res = destructive ? z : zonotope_copy(man, z);
    size_t i = 0;
    size_t j = 0;
    size_t num_remove = dimchange->intdim + dimchange->realdim;
    size_t num_rem = z->dims - num_remove;
    char * map = (char *)calloc(z->dims,sizeof(char));
    elina_dim_t * var_rem = (elina_dim_t*)malloc(num_rem*sizeof(elina_dim_t));
    size_t count = 0;
    
    for(;i < num_remove; i++){
	map[dimchange->dim[i]] = 1;
    }
    
    for(i=0; i < z->dims; i++){
	if(!map[i]){
		//printf("AI2 %d %d\n",i,count);
		//fflush(stdout);
		var_rem[count] = i;
		count++;
	}
    }
    zonotope_t* res = (zonotope_t *)malloc(sizeof(zonotope_t));
    res->intdim = z->intdim -dimchange->intdim;
    res->dims = num_rem;
    res->size = z->size;
    res->nsymcons = (elina_dim_t*)calloc(res->size, sizeof(elina_dim_t));
    memcpy((void *)res->nsymcons, (void *)z->nsymcons, (res->size)*sizeof(elina_dim_t));
    res->gamma = (elina_interval_t**)calloc(res->size, sizeof(elina_interval_t*));
    res->abs = elina_abstract0_copy(pr->manNS, z->abs);
    res->hypercube = z->hypercube;
    res->box = elina_interval_array_alloc(res->dims);
    res->gn = 0;
    res->g = NULL;
    res->paf = (zonotope_aff_t**)malloc(res->dims*sizeof(zonotope_aff_t*));
    for (i=0; i< num_rem; i++) {
	j = var_rem[i];
	res->paf[i] = zonotope_aff_alloc_init(pr);
	//printf("coming here: %d %d %d %d %p %d\n",i,j,num_rem,num_remove,res->paf[i],res->dims);
	//fflush(stdout);
	memcpy((void *)res->paf[i], (void *)z->paf[j], sizeof(zonotope_aff_t));
	//zonotope_aff_check_free(pr,res->paf[dimchange->dim[i]]);
	//res->paf[dimchange->dim[i]] = NULL;
	//for (j=dimchange->dim[i];j<-1+res->dims;j++) {
	  //  res->paf[j] = res->paf[j+1];
	elina_interval_set(res->box[i],z->box[j]);
	
	res->paf[i]->pby++;
	//}
    }
    
   
    size_t nsymcons_size = zonotope_noise_symbol_cons_get_dimension(pr, z);
    for (i=0; i<nsymcons_size; i++) {
	if (z->gamma[i]) {
	    if (z->gamma[i] != pr->ap_muu) res->gamma[i] = elina_interval_alloc_set(z->gamma[i]);
	    else res->gamma[i] = pr->ap_muu;
	} else {
	    printf("zonotope_copy, unconsistent gamma for Zonotope abstract object\n");
	}
    }
    
    //res->intdim = z->intdim - dimchange->intdim;
    //res->dims = z->dims - (dimchange->intdim + dimchange->realdim);
    //res->box = realloc(res->box,res->dims*sizeof(elina_interval_t *));
    //res->paf = realloc(res->paf,res->dims*sizeof(zonotope_aff_t*));
    //printf("remove output %p\n",res);
	//zonotope_fprint(stdout,man,res,NULL);
	//fflush(stdout);
    free(map);
    free(var_rem);
    zonotope_t * tmp = z; 
    if(destructive){
	z = res;
	//zonotope_free(man,tmp);
	return z;
    }
    return res;
}


zonotope_t* zonotope_permute_dimensions(elina_manager_t* man, bool destructive, zonotope_t* z, elina_dimperm_t* permutation)
{
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_PERMUTE_DIMENSIONS);
    zonotope_aff_t* tmp = NULL;
    zonotope_t* res = destructive ? z : zonotope_copy(man, z);
    size_t i = 0;
    for (i=0; i<permutation->size; i++) {
	tmp = res->paf[i];
	res->paf[i] = res->paf[permutation->dim[i]];
	res->paf[permutation->dim[i]] = tmp;
	elina_interval_swap(res->box[i],res->box[permutation->dim[i]]);
    }
    return res;
}

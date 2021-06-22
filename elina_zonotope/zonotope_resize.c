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

#include "zonotope_resize.h"


/*********************/
/* Add dimensions */
/*********************/
zonotope_t* zonotope_add_dimensions(elina_manager_t* man, bool destructive, zonotope_t* z, elina_dimchange_t* dimchange, bool project)
{
    //printf("add input %d\n",destructive);
    //zonotope_fprint(stdout,man,z,NULL);
    //fflush(stdout);
    start_timing();
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ADD_DIMENSIONS);
    
    zonotope_t* res = destructive ? z : zonotope_copy(man, z);
    size_t intdim = z->intdim + dimchange->intdim;
    size_t dims = z->dims + (dimchange->intdim + dimchange->realdim);
    res->box_inf = realloc(res->box_inf,dims*sizeof(double));
    res->box_sup = realloc(res->box_sup,dims*sizeof(double));
    res->paf = realloc(res->paf,dims*sizeof(zonotope_aff_t*));
    size_t i = 0;
    int j = 0;
    for (i=0; i<dimchange->intdim + dimchange->realdim; i++) {
        if (res->dims == dimchange->dim[i]) {
            /* add in the last place */
            res->box_inf[res->dims] = INFINITY;
	    res->box_sup[res->dims] = INFINITY;
        } else {
            /* increment */
            for (j=(int)-1+res->dims;j>=(int)dimchange->dim[i];j--) {
                res->paf[j+1] = res->paf[j];
		res->box_inf[j+1] = res->box_inf[j];
		res->box_sup[j+1] = res->box_sup[j];
            }
        }
        res->paf[dimchange->dim[i]] = project ? zonotope_aff_alloc_init(pr) : pr->top;
        res->paf[dimchange->dim[i]]->pby++;
        if (project){
	    res->box_inf[dimchange->dim[i]]= 0.0;
	    res->box_sup[dimchange->dim[i]]= 0.0;
        }
        else {
	    res->box_inf[dimchange->dim[i]] = INFINITY;
	    res->box_sup[dimchange->dim[i]] = INFINITY;
	}
        res->dims++;
    }
    res->intdim = intdim;
    //printf("add output\n");
    //zonotope_fprint(stdout,man,res,NULL);
    //fflush(stdout);
    record_timing(zonotope_add_dimension_time);
    return res;
}



/*********************/
/* Remove dimensions, Optimized for AI^2, Linear cost compared with quadratic with APRON */
/*********************/
zonotope_t* zonotope_remove_dimensions(elina_manager_t* man, bool destructive, zonotope_t* z, elina_dimchange_t* dimchange)
{
    start_timing();
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_REMOVE_DIMENSIONS);
    //printf("remove input %d\n",destructive);
    //zonotope_fprint(stdout,man,z,NULL);
    //for(int i=0; i < dimchange->intdim+dimchange->realdim;i++){
      //printf("%d ",dimchange->dim[i]);
    //}
    //printf("\n");
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
   
    res->box_inf = (double*)malloc(res->dims*sizeof(double));
    res->box_sup = (double*)malloc(res->dims*sizeof(double));
    res->gn = 0;
    //res->g = NULL;
    res->paf = (zonotope_aff_t**)malloc(res->dims*sizeof(zonotope_aff_t*));
    for (i=0; i< num_rem; i++) {
        j = var_rem[i];
        res->paf[i] = zonotope_aff_alloc_init(pr);
        //printf("coming here: %d %d %d %d %p %d\n",i,j,num_rem,num_remove,res->paf[i],res->dims);
        //fflush(stdout);
        //res->paf[i] = (zonotope_aff_t *)malloc(sizeof(zonotope_aff_t*));
        res->paf[i] = z->paf[j];
        //memcpy((void *)res->paf[i], (void *)z->paf[j], sizeof(zonotope_aff_t));
        //zonotope_aff_check_free(pr,res->paf[dimchange->dim[i]]);
        //res->paf[dimchange->dim[i]] = NULL;
        //for (j=dimchange->dim[i];j<-1+res->dims;j++) {
        //  res->paf[j] = res->paf[j+1];
	res->box_inf[i] = z->box_inf[j];
	res->box_sup[i] = z->box_sup[j];
        if(!destructive)
        res->paf[i]->pby++;
        //}
    }
    
    
    size_t nsymcons_size = zonotope_noise_symbol_cons_get_dimension(pr, z);
    for (i=0; i<nsymcons_size; i++) {
        if (z->gamma[i]) {
            if (z->gamma[i] != pr->ap_muu) res->gamma[i] = elina_interval_alloc_set(z->gamma[i]);
            else res->gamma[i] = pr->ap_muu;
        } else {
	    fprintf(stderr,"zonotope_remove, unconsistent gamma for Zonotope abstract object\n");
        }
    }
    
    //res->intdim = z->intdim - dimchange->intdim;
    //res->dims = z->dims - (dimchange->intdim + dimchange->realdim);
    
    //res->paf = realloc(res->paf,res->dims*sizeof(zonotope_aff_t*));
    //printf("remove output %p\n",res);
    //zonotope_fprint(stdout,man,res,NULL);
    //fflush(stdout);
    
    zonotope_t * tmp = z;
    if(destructive){
        //z = res;
        //zonotope_free(man,tmp);
        //return z;
        
        //for(i=0; i < z->dims; z++){
          //  if(map[i]){
            //    zonotope_aff_check_free(pr,z->paf[i]);
            //}
        //}
        //free(z->paf);
        //z->paf=NULL;
       
        //for (i=0; i<nsymcons_size; i++) {
          //  if (z->gamma[i]) {
            //    if (z->gamma[i] != pr->ap_muu) elina_interval_free(z->gamma[i]);
            //}
        //}
        //free(z->gamma);
        //z->gamma = NULL;
        //free(z->nsymcons);
        //z->nsymcons = NULL;
        //elina_abstract0_free(pr->manNS, z->abs);
        //free(z);
    }
    free(map);
    free(var_rem);
    record_timing(zonotope_remove_dimension_time);
    return res;
}


zonotope_t* zonotope_permute_dimensions(elina_manager_t* man, bool destructive, zonotope_t* z, elina_dimperm_t* permutation)
{
	
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_PERMUTE_DIMENSIONS);
    //zonotope_aff_t* tmp = NULL;
   // printf("permutation input %d\n",destructive);
    //zonotope_fprint(stdout,man,z,NULL);
    //for(int i=0; i < permutation->size;i++){
      //printf("%d ",permutation->dim[i]);
    //}
    //printf("\n");
    //fflush(stdout);
    start_timing();
    zonotope_t* res = (zonotope_t *)malloc(sizeof(zonotope_t));
    res->intdim = z->intdim;
    res->dims = z->dims;
    res->size = z->size;
    res->nsymcons = (elina_dim_t*)calloc(res->size, sizeof(elina_dim_t));
    memcpy((void *)res->nsymcons, (void *)z->nsymcons, (res->size)*sizeof(elina_dim_t));
    res->gamma = (elina_interval_t**)calloc(res->size, sizeof(elina_interval_t*));
    res->abs = elina_abstract0_copy(pr->manNS, z->abs);
    res->hypercube = z->hypercube;
    res->box_inf = (double*)malloc(res->dims*sizeof(double));
    res->box_sup = (double*)malloc(res->dims*sizeof(double));
    res->gn = 0;
    //res->g = NULL;
    res->paf = (zonotope_aff_t**)malloc(res->dims*sizeof(zonotope_aff_t*));
    size_t i = 0;
    size_t j = 0;
    for (i=0; i< permutation->size; i++) {
        j = permutation->dim[i];
        res->paf[j] = zonotope_aff_alloc_init(pr);
        //printf("coming here: %d %d %d %d %p %d\n",i,j,num_rem,num_remove,res->paf[i],res->dims);
        //fflush(stdout);
        memcpy((void *)res->paf[j], (void *)z->paf[i], sizeof(zonotope_aff_t));
        //zonotope_aff_check_free(pr,res->paf[dimchange->dim[i]]);
        //res->paf[dimchange->dim[i]] = NULL;
        //for (j=dimchange->dim[i];j<-1+res->dims;j++) {
        //  res->paf[j] = res->paf[j+1];
	res->box_inf[j] = z->box_inf[i];
        res->box_sup[j] = z->box_sup[i];
        if(!destructive)
        res->paf[i]->pby++;
        //}
    }
    
    
    size_t nsymcons_size = zonotope_noise_symbol_cons_get_dimension(pr, z);
    for (i=0; i<nsymcons_size; i++) {
        if (z->gamma[i]) {
            if (z->gamma[i] != pr->ap_muu) res->gamma[i] = elina_interval_alloc_set(z->gamma[i]);
            else res->gamma[i] = pr->ap_muu;
        } else {
	    fprintf(stderr,"zonotope_permute, unconsistent gamma for Zonotope abstract object\n");
        }
    }
    //if(destructive){
      //  z = res;
        //zonotope_free(man,tmp);
        //return z;
    //}
    
    //char *map = (char *)calloc(permutation->size,sizeof(char));
    //for (i=0; i<permutation->size; i++) {
    //if(!map[i] && !map[permutation->dim[i]]){
    //      tmp = res->paf[i];
    //    res->paf[i] = res->paf[z->dim[i]];
    //res->paf[permutation->dim[i]] = tmp;
    
    
    //map[i] = 1;
    //map[permutation->dim[i]] = 1;
    //}
    //} 
    //printf("permutation output %p\n",res);
    //zonotope_fprint(stdout,man,res,NULL);
    //fflush(stdout);
    //free(map);
    record_timing(zonotope_permute_dimension_time);
    return res;
}

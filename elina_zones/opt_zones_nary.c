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

#include "opt_mat.h"

/*******
Standard Meet Operator
****/

opt_zones_t* opt_zones_meet(elina_manager_t* man, bool destructive, opt_zones_t* o1, opt_zones_t* o2)
{
 
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_MEET,0);
  opt_zones_mat_t* oz;
  
 
   if((o1->dim != o2->dim) || (o1->intdim != o2->intdim))return NULL;
  if ((!o1->closed && !o1->m) || (!o2->closed && !o2->m))
    /* one argument is empty */
    return opt_zones_set_mat(pr,o1,NULL,NULL,destructive);
  else {
    opt_zones_mat_t * oz1 = o1->closed ? o1->closed : o1->m;
    opt_zones_mat_t * oz2 = o2->closed ? o2->closed : o2->m;
    oz = destructive ? oz1 : opt_zones_mat_alloc(o1->dim);
    #if defined(TIMING)
	start_timing();
    #endif
    meet_zones_mat(oz,oz1,oz2,o1->dim,destructive);
    /* optimal, but not closed */
    #if defined(TIMING)
	record_timing(zones_meet_time);
    #endif
    opt_zones_t * r = opt_zones_set_mat(pr,o1,oz,NULL,destructive);
    
   
    return r;
  }
}
array_comp_list_t * compute_finest_zones(opt_zones_mat_t * oz, unsigned short int dim){
	array_comp_list_t *res = create_array_comp_list();
	comp_list_t * cl = oz->acl->head;
	double * m = oz->mat;
	unsigned short int n = dim+1;
	while(cl!=NULL){
		unsigned short int comp_size = cl->size;
		unsigned short int *ca = to_sorted_array(cl,dim);
		unsigned short int i,j;
		for(i=0; i < comp_size; i++){
			unsigned short int i1 = ca[i];
			if(m[n*(i1+1)]!=INFINITY || m[i1+1]!=INFINITY){
				comp_list_t * cl = create_comp_list();
				insert_comp(cl,i1);
				insert_comp_list(res,cl);
			}
		}
		for(i=0; i < comp_size; i++){
			unsigned short int i1 = ca[i];
			for(j=0; j < i; j++){
				unsigned short int j1 = ca[j];
				if(m[n*(i1+1)+(j1+1)]!=INFINITY || m[n*(j1+1)+i1+1]!=INFINITY){
					comp_list_t * cl = create_comp_list();
					insert_comp(cl,i1);
					insert_comp(cl,j1);
					insert_comp_list_with_union(res,cl,dim);
				}
			}
		}
		free(ca);
		cl = cl->next;
	}
	return res;
}
/*****
Standard Join Operator
***/
opt_zones_t* opt_zones_join(elina_manager_t* man, bool destructive, opt_zones_t* o1, opt_zones_t* o2)
{

 
 opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_JOIN,0);
  
 if((o1->dim != o2->dim) || (o1->intdim != o2->intdim))return NULL;
 if (pr->funopt->algorithm>=0) {
   opt_zones_mat_t* oz1 = o1->closed ? o1->closed : o1->m;
   if(oz1 && !oz1->is_dense){
	opt_zones_sparse_weak_closure(pr,o1);
   }
   else{
	opt_zones_cache_closure(pr,o1);
   }
   opt_zones_mat_t* oz2 = o2->closed ? o2->closed : o2->m; 
   if(oz2 && !oz2->is_dense){
	 opt_zones_sparse_weak_closure(pr,o2);
   }
   else{
	 opt_zones_cache_closure(pr,o2);
   }
  
   
 }
 opt_zones_t * r;
 if (!o1->closed && !o1->m) {
   if (!o2->closed && !o2->m){
     /* both empty */
	
     r = opt_zones_set_mat(pr,o1,NULL,NULL,destructive);
	
	return r;
   }
   else{
     
     /* a1 empty, a2 not empty */
     r = opt_zones_set_mat(pr,o1,opt_zones_mat_copy(o2->m,o1->dim),
			opt_zones_mat_copy(o2->closed,o1->dim),destructive);
	
	return r;
     }
 }
 else if (!o2->closed && !o2->m){
   /* a1 not empty, a2 empty */
  
   r = opt_zones_set_mat(pr,o1,o1->m,o1->closed,destructive);
	
	return r;
  }
 else {
   /* not empty */
  
   opt_zones_mat_t* oz1 = o1->closed ? o1->closed : o1->m;
   opt_zones_mat_t* oz2 = o2->closed ? o2->closed : o2->m;
   opt_zones_mat_t * oz = destructive ? oz1 : opt_zones_mat_alloc(o1->dim);
   size_t i;
   man->result.flag_exact = false;
    #if defined(TIMING)
	start_timing();
    #endif
   if(!oz1->is_dense && !oz2->is_dense){
	sparse_join_zones_mat(oz,oz1,oz2,o1->dim,destructive);
   }
   else{
   	join_zones_mat(oz,oz1,oz2,o1->dim,destructive);
   }
    #if defined(TIMING)
	record_timing(zones_join_time);
    #endif
   if (o1->closed && o2->closed) {
     /* result is closed and optimal on Q */
     if (zone_num_incomplete || o1->intdim) zones_flag_incomplete;
     r = opt_zones_set_mat(pr,o1,NULL,oz,destructive);
   }
   else {
     /* not optimal, not closed */
     zone_flag_algo;
     r =  opt_zones_set_mat(pr,o1,oz,NULL,destructive);
   }
	
   return r;
 }
}

/****
	Standard Widening Operator
***/

opt_zones_t* opt_zones_widening(elina_manager_t* man, opt_zones_t* o1, opt_zones_t* o2)
{
  //printf("widening start %p %p\n",o1,o2);
  //fflush(stdout);
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_WIDENING,0);
  int algo = pr->funopt->algorithm;
  opt_zones_t* r;
  
  
  if((o1->dim != o2->dim) || (o1->intdim != o2->intdim))return NULL;
  if (algo>=0) {
	opt_zones_mat_t * oz2 = o2->closed ? o2->closed : o2->m;
	if(oz2 && !oz2->is_dense){
		opt_zones_sparse_weak_closure(pr,o2);
	}
	else{
		opt_zones_cache_closure(pr,o2);
	}
  }
  if (!o1->closed && !o1->m){
    /* o1 definitively closed */
    r = opt_zones_copy_internal(pr,o2);
    
  }
  else if (!o2->closed && !o2->m){
   /* o2 definitively closed */
    r = opt_zones_copy_internal(pr,o1);
    
  }
  else {
    /* work on the origial left matrix, not the closed cache! */
    opt_zones_mat_t * oz1 = o1->m ? o1->m : o1->closed;
    opt_zones_mat_t * oz2 = o2->closed ? o2->closed : o2->m;
    size_t i;
    r = opt_zones_alloc_internal(pr,o1->dim,o1->intdim);
    r->m = opt_zones_mat_alloc(o1->dim);
    //posix_memalign((void **)&(r->m),32,size*sizeof(double));
    if (algo==opt_zones_pre_widening || algo==-opt_zones_pre_widening) {
      /* degenerate hull: NOT A PROPER WIDENING, use with care */
	join_zones_mat(r->m,oz1,oz2,o1->dim,false);
    }
    else {
      /* standard widening */
	#if defined (TIMING)
 		start_timing();
	#endif
	 if(!oz1->is_dense && !oz2->is_dense){
		
		sparse_widening_zones_mat(r->m,oz1,oz2,o1->dim);
		
  	 }
	else{
        	widening_zones_mat(r->m,oz1,oz2,o1->dim);
	}
	 #if defined(TIMING)
		record_timing(zones_widening_time);
    	#endif
     
    }
  }
  return r;
}

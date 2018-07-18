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


#include "opt_mat.h"

/****************
Projection:
*****************/

opt_zones_t* opt_zones_forget_array(elina_manager_t* man,
			bool destructive, opt_zones_t* o,
			elina_dim_t* tdim, int size,
			bool project)
{
   
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_FORGET_ARRAY,0);
  if (pr->funopt->algorithm>=0){
	opt_zones_mat_t* oz = o->closed ? o->closed : o->m;
	if(oz && !oz->is_dense){
		opt_zones_sparse_weak_closure(pr,o);
	}
	else{
		opt_zones_cache_closure(pr,o);
	}
  }
  if (!o->closed && !o->m)
    /* definitively empty */
    return opt_zones_set_mat(pr,o,NULL,NULL,destructive);
  else {
    opt_zones_mat_t* oz = o->closed ? o->closed : o->m;
    if (!destructive) oz = opt_zones_mat_copy(oz,o->dim);
    #if defined(TIMING)
  	start_timing();
    #endif
    forget_array_zones_mat(oz,tdim,o->dim,size,project);
    #if defined(TIMING)
  	record_timing(zones_forget_array_time);
    #endif
    opt_zones_t *r;
    if (o->closed) {
      /* result is exact on Q, and closed if forget, not project */
      if (zone_num_incomplete || o->intdim) zones_flag_incomplete;
      if (project) r = opt_zones_set_mat(pr,o,oz,NULL,destructive);
      else r = opt_zones_set_mat(pr,o,NULL,oz,destructive);
    }
    else {
      /* not exact, not closed */
      zone_flag_algo;
      r = opt_zones_set_mat(pr,o,oz,NULL,destructive);
    }
    
    return r;
  }
}



/******
Add Dimensions
*****/

opt_zones_t* opt_zones_add_dimensions(elina_manager_t* man,
			  bool destructive, opt_zones_t* o,
			  elina_dimchange_t* dimchange,
			  bool project)
{
   
   
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_ADD_DIMENSIONS,0);
  opt_zones_mat_t* src = o->closed ? o->closed : o->m;
  //print_mat(src,o->dim);
  opt_zones_mat_t* dst;
  size_t i, nb = dimchange->intdim+dimchange->realdim;
  opt_zones_t* r;
  if (!src) dst = NULL;
  else {
    /* check */
    
  for (i=0;i<nb;i++) {
    if(dimchange->dim[i] > o->dim) return NULL;
    if(i &&(dimchange->dim[i-1]> dimchange->dim[i])) return NULL;
   }
    /* insert variables */
    int dim = o->dim + nb;
    dst = opt_zones_mat_alloc_top(dim);
     
    #if defined(TIMING)
	start_timing();
    #endif
    
    opt_zones_mat_addrem_dimensions(dst,src,dimchange->dim,
			   nb,1,o->dim,true);
    
    if(dst->is_dense){
	free_array_comp_list(dst->acl);
    }
    int count = dst->nni;
    /* set new variables to 0, if necessary */
    if (project) {
      unsigned short int n = dim + 1;
      for (i=0;i<nb;i++) {
	double *mm = dst->mat;
	unsigned short int v = i+dimchange->dim[i]+1;
	mm[v] = 0;
	mm[n*v]  = 0;
	count = count + 2;
      }
    }
    dst->nni = count;
     #if defined(TIMING)
	record_timing(zones_add_dimension_time);
    #endif
  }
  /* always exact, respect closure if embedding, not projecting */
  if (o->closed && !project) r = opt_zones_set_mat(pr,o,NULL,dst,destructive);
  else r = opt_zones_set_mat(pr,o,dst,NULL,destructive);
  r->dim += nb;
  r->intdim += dimchange->intdim;
  
  return r;
}


/********
Remove dimensions
********/

opt_zones_t* opt_zones_remove_dimensions(elina_manager_t* man,
			     bool destructive, opt_zones_t* o,
			     elina_dimchange_t* dimchange)
{
 
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_REMOVE_DIMENSIONS,0);
  opt_zones_mat_t *src, *dst;
  unsigned short int i, nb = dimchange->intdim+dimchange->realdim;
  opt_zones_t* r;
  if (pr->funopt->algorithm>=0){
	opt_zones_mat_t* oz = o->closed ? o->closed : o->m;
	if(oz && !oz->is_dense){
		opt_zones_sparse_weak_closure(pr,o);
	}
	else{
		opt_zones_cache_closure(pr,o);
	}
  }
  src = o->closed ? o->closed : o->m;
  //print_mat(src,o->dim);
  if (!src) dst = NULL;
  else {
    /* check */
     for (i=0;i<nb;i++) {
		if(dimchange->dim[i] >= o->dim) return NULL;
		if(i && (dimchange->dim[i-1]>=dimchange->dim[i]))return NULL;
     }
    /* remove variables */
    int dim = o->dim - nb;
    
    dst = opt_zones_mat_alloc(dim);
     #if defined(TIMING)
	start_timing();
    #endif
    opt_zones_mat_addrem_dimensions(dst,src,dimchange->dim,
			   nb,1,o->dim,false);

    if(dst->is_dense){
	free_array_comp_list(dst->acl);
    }
     #if defined(TIMING)
	record_timing(zones_add_dimension_time);
    #endif
  }

  if (o->closed) {
    /* result is exact on Q, and closed */
    if (zone_num_incomplete || o->intdim) zones_flag_incomplete;
    r = opt_zones_set_mat(pr,o,NULL,dst,destructive);
  }
  else {
    /* not exact, not closed */
    zone_flag_algo;
    r = opt_zones_set_mat(pr,o,dst,NULL,destructive);
  }
  r->dim -= nb;
  r->intdim -= dimchange->intdim;
  
  return r;
}


/***********
Permute Dimensions
***********/

opt_zones_t* opt_zones_permute_dimensions(elina_manager_t* man,
			      bool destructive, opt_zones_t* o,
			      elina_dimperm_t* permutation)
{
  
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_ADD_DIMENSIONS,0);
  opt_zones_mat_t* src = o->closed ? o->closed : o->m;
  //print_mat(src,o->dim);
  opt_zones_mat_t* dst;
  if(permutation->size!=o->dim)return NULL;
  if (!src) dst = NULL;
  else {
    /* check (only bounds, not injectivity) */
    int i,j;    
    /* permuted copy */
     for (i=0;i<o->dim;i++){
          if(permutation->dim[i]>=o->dim)return NULL;
     }

    dst = opt_zones_mat_alloc(o->dim);
     #if defined(TIMING)
	start_timing();
    #endif
    
    opt_zones_mat_permute(dst,src,o->dim,o->dim,permutation->dim);
     
    if(dst->is_dense){
	free_array_comp_list(dst->acl);
    }
     #if defined(TIMING)
	record_timing(zones_permute_dimension_time);
    #endif
  }
  /* always exact, respects closure */
  opt_zones_t * r;
  if (o->closed) r= opt_zones_set_mat(pr,o,NULL,dst,destructive);
  else r = opt_zones_set_mat(pr,o,dst,NULL,destructive);
  
  return r;
}

/******************
	Expand Operator
*******************/

opt_zones_t* opt_zones_expand(elina_manager_t* man,
		  bool destructive, opt_zones_t* o,
		  elina_dim_t dim,
		  size_t n)
{
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_EXPAND,0);
  opt_zones_mat_t* src = o->closed ? o->closed : o->m;
  size_t i, j;
  elina_dim_t pos = (dim < o->intdim) ? o->intdim : o->dim;
  opt_zones_mat_t* dst;
  opt_zones_t* r;
  if (!src) dst = NULL;
  else {
	
    /* insert n variables at pos */
    dst = opt_zones_mat_alloc_top(o->dim+n);
    opt_zones_mat_addrem_dimensions(dst,src,&pos,1,n,o->dim,true);
    //#if defined(TIMING)
  	//start_timing();
    //#endif
    double *mm = dst->mat;
    unsigned short int sn = o->dim + 1;
    unsigned short int dn = sn + n; 
    if(!src->is_dense){
	    for (i=0;i<n;i++) {

	      /* copy binary constraints */
	      comp_list_t * cl = create_comp_list();
	      comp_list_t * cj = find(dst->acl,dim);	
	      //for (j=0;j<2*dim;j++) {
	      
	      if(cj != NULL){
		comp_t *c = cj->head;
		while(c != NULL){
			unsigned short int  j = c->num;
			ini_relation_zones(mm,pos+i+1,j+1,o->dim+n);
			if(j==dim){
				c = c->next;
				continue;
			}
		
			mm[dn*(pos+i+1) + j+1] = mm[dn*(dim+1) + j + 1];
			mm[dn*(j+1) + pos + i +1] = mm[dn*(j+1) + dim + 1];
			insert_comp(cl,j);
			c = c->next;
		}
	      }
	
	      int flag = 0;
	      comp_list_t * cp = find(dst->acl,pos+i);
	      if(cp!=NULL){
		if(cp==cj){
			mm[dn*(pos+i+1) + pos+i+1] = mm[dn*(dim+1) + pos+i+1];
		}
		else{
		
			mm[dn*(pos+i+1) + pos+i+1] = INFINITY;
			//flag = 1;
		}
	      }
	      else{
		
			mm[dn*(pos+i+1) + pos+i+1] = INFINITY;
			//flag = 1;
	      }
	      /* copy unary constraints */
	      if(cj!=NULL){
		mm[dn*(pos+i+1)] = mm[dn*(dim+1)];
		mm[pos+i+1] = mm[dim+1];
	      }
	      else{
			mm[dn*(pos+i+1)] = INFINITY;
			mm[pos+i+1] = INFINITY;
	      }
		insert_comp(cl,pos+i);
		insert_comp_list_with_union(dst->acl,cl,o->dim+n);
	    }
		 
  }
  else{
	free_array_comp_list(dst->acl);
	 /* copy binary and unary constraints */
	for (i=0;i<n;i++) {
		int src_ind,dest_ind;
		 #if defined(VECTOR)
		      v_double_type src;
		      double * pd = mm + dn*(dim+1);
		      double * pp = mm + dn*(pos+i+1);
		      for (j=0;j<(dim)/v_length;j++) {
			src = v_load_double(pd+j*v_length);
			v_store_double(pp+j*v_length,src); 
		      }
		      for(j = (dim/v_length)*v_length;j < dim; j++){
			 mm[dn*(pos+i+1) + j] = mm[dn*(dim+1) + j];
	      	      }	
		      for (j=dim+1;j<sn/v_length;j++) {
			src = v_load_double(pd+j*v_length);
			v_store_double(pp+j*v_length,src); 
		      }
		      for(j = (sn/v_length)*v_length;j < sn; j++){
			 mm[dn*(pos+i+1) + j] = mm[dn*(dim+1) + j];
	      	      }	
	      #else
		      for(j= 0; j < dim;j++ ){
			  mm[dn*(pos+i+1) + j] = mm[dn*(dim+1) + j];
		      }
		      for(j=dim+1; j < sn; j++){
			  mm[dn*(pos+i+1) + j] = mm[dn*(dim+1) + j];
		      }
		   
	      #endif
	      
	      for(j= 0; j < dim;j++ ){
		  mm[dn*j + pos + i+1] = mm[dn*(j+1) + dim+1];
	      }
	      for(j=dim+1; j < sn; j++){
		  mm[dn*j + pos + i+1] = mm[dn*(j+1) + dim+1];
	      }
     
     	}
     //#if defined(TIMING)
  	//record_timing(expand_time);
     //#endif
   }
   
    int dst_size = dn*dn;
    int src_size = sn*sn;
    /*  exact, generally not closed */
    dst->nni = min(dst_size,src->nni + dst_size-src_size);
  }
  r = opt_zones_set_mat(pr,o,dst,NULL,destructive);
  r->dim += n;
  if (dim<o->intdim) r->intdim += n;
  
  return r;
}


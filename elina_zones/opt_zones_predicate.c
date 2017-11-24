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


#include "opt_mat.h"

bool opt_zones_is_bottom(elina_manager_t* man, opt_zones_t* o)
{
  //printf("is bottom: start\n");
  //fflush(stdout);
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_IS_BOTTOM,0);
  if (pr->funopt->algorithm>=0){ 
	opt_zones_mat_t * oz = o->closed ? o->closed : o->m;
	if(oz && !oz->is_dense){
		opt_zones_sparse_weak_closure(pr,o);	
	}
	else{
		opt_zones_cache_closure(pr,o);
	}
  }
  //m = o->closed ? o->closed : o->m;
  
  if (o->closed) {
    /* definitively non empty on Q */
   
    if (zone_num_incomplete || o->intdim) { zones_flag_incomplete; }
    
    return false;
  }
  else if (!o->m){
    /* definitively empty */
    
    return true;
  }
  else {
    /* no closure => we don't know */
    zone_flag_algo;
   
    return false;
  }
}

bool opt_zones_is_top(elina_manager_t* man, opt_zones_t* o)
{
  
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_IS_TOP,0);
  int i,j;
  opt_zones_mat_t* oz = o->m ? o->m : o->closed;
  
  if (!oz) return false;
  #if defined(TIMING)
	start_timing();
  #endif
  bool res = is_top_zones_mat(oz,o->dim);
  #if defined(TIMING)
        record_timing(zones_is_top_time);
  #endif
  
  return res;
}


bool opt_zones_is_leq(elina_manager_t* man, opt_zones_t* o1, opt_zones_t* o2)
{
  
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_IS_LEQ,0);
   if((o1->dim != o2->dim) || (o1->intdim != o2->intdim))return false;
   if (pr->funopt->algorithm>=0){ 
	opt_zones_mat_t * oz1 = o1->closed ? o1->closed : o1->m;
    	if(oz1 && !oz1->is_dense){
		opt_zones_sparse_weak_closure(pr,o1);	
    	}
    	else{
		opt_zones_cache_closure(pr,o1);
    	}
  }
  if (!o1->closed && !o1->m) {
    /* a1 definitively empty */
    return true;
  }
  else if (!o2->closed && !o2->m) {
    /* a2 definitively empty */
    if (o1->closed) {
      /* a1 not empty on Q */
      if (zone_num_incomplete || o1->intdim) { zones_flag_incomplete; }
      return false;
    }
    else { zone_flag_algo; return false; }
  }
  else {
    opt_zones_mat_t *oz1 = o1->closed ? o1->closed : o1->m;
    opt_zones_mat_t *oz2 = o2->closed ? o2->closed : o2->m;
    #if defined(TIMING)
	start_timing();
    #endif
    bool res= is_lequal_zones_mat(oz1, oz2, o1->dim);
    #if defined(TIMING)
	record_timing(zones_is_lequal_time);
    #endif
   
    return res;
  }
}


bool opt_zones_is_eq(elina_manager_t* man, opt_zones_t* o1, opt_zones_t* o2)
{
 
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_IS_EQ,0);
   if((o1->dim != o2->dim) || (o1->intdim != o2->intdim))return false;
  if (pr->funopt->algorithm>=0) {
     
    opt_zones_mat_t * oz1 = o1->closed ? o1->closed : o1->m;
    if(oz1 && !oz1->is_dense){
	opt_zones_sparse_weak_closure(pr,o1);	
    }
    else{
	opt_zones_cache_closure(pr,o1);
    }
    opt_zones_mat_t * oz2 = o2->closed ? o2->closed : o2->m;
    if(oz2 && !oz2->is_dense){
	opt_zones_sparse_weak_closure(pr,o2);	
    }
    else{
	opt_zones_cache_closure(pr,o2);
    }
  }
  if (!o1->closed && !o1->m) {
    if (!o2->closed && !o2->m) {
      /* both are empty */
      return true;
    }
    else if (o2->closed) {
      /* a1 empty, e2 not empty on Q */
      if (zone_num_incomplete || o1->intdim) { zones_flag_incomplete; }
      return false;
    }
    else { zone_flag_algo; return false; }
  }
  else if (!o2->closed && !o2->m) {
    if (o1->closed) {
      /* a2 empty, e1 not empty on Q */
      if (zone_num_incomplete || o1->intdim) { zones_flag_incomplete; }
      return false;
    }
    else { zone_flag_algo; return false; }
  }
  else {
    
    opt_zones_mat_t *oz1 = o1->closed ? o1->closed : o1->m;
    opt_zones_mat_t *oz2 = o2->closed ? o2->closed : o2->m;
	
    #if defined(TIMING)
	start_timing();
    #endif
    bool res = is_equal_zones_mat(oz1,oz2,o1->dim);
     #if defined(TIMING)
	record_timing(zones_is_equal_time);
    #endif
    
    return res;
  }
}

elina_tcons0_array_t opt_zones_to_tcons_array(elina_manager_t* man, opt_zones_t* o)
{
  return elina_generic_to_tcons_array(man,o);
}


elina_interval_t** opt_zones_to_box(elina_manager_t* man, opt_zones_t* o)
{
 
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_TO_BOX,0);
  elina_interval_t** in = elina_interval_array_alloc(o->dim);
  unsigned short int i;
  if (pr->funopt->algorithm>=0){ 
	opt_zones_mat_t * oz = o->closed ? o->closed : o->m;
	if(oz && !oz->is_dense){
		opt_zones_sparse_weak_closure(pr,o);	
	}
	else{
		opt_zones_cache_closure(pr,o);
	}
  }
   #if defined(TIMING)
	start_timing();
    #endif
  if (!o->closed && !o->m) {
    /* definitively empty */
    for (i=0;i<o->dim;i++)
      elina_interval_set_bottom(in[i]);
  }
  else {
    /* put variable bounds */
    opt_zones_mat_t* oz = o->closed ? o->closed : o->m;
    double *m = oz->mat;
    unsigned short int n = o->dim + 1;
    if(!oz->is_dense){
	    array_comp_list_t * acl = oz->acl;
	     for (i=0;i<o->dim;i++){
	     	elina_interval_set_top(in[i]);
	     }
	    comp_list_t * cl = acl->head;
	    while(cl!=NULL){
		    unsigned short int comp_size = cl->size;
		    comp_t * c = cl->head;
		    for (i=0;i<comp_size;i++){
		      int i1 = c->num;
		      zones_interval_of_bounds(in[i1],
					 m[n*(i1+1)],m[i1+1]);
			
			
			 c = c->next;
		    }
		   cl = cl->next;
	    }
    }
    else{
	
	 for (i=1;i<n;i++){
     		 zones_interval_of_bounds(in[i-1],
			 m[n*i],m[i]);
		
    	}
    }
    man->result.flag_exact = false;
    if (!o->closed) zone_flag_algo;
    else if (zone_num_incomplete || o->intdim) zones_flag_incomplete;
    else if (pr->conv) zone_flag_conv;
   
  }
   #if defined(TIMING)
	record_timing(zones_to_box_time);
    #endif
  
  return in;
}

elina_lincons0_array_t opt_zones_to_lincons_array(elina_manager_t* man, opt_zones_t* o)
{
 
  elina_lincons0_array_t ar;
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_TO_LINCONS_ARRAY,0);
  
  if (!o->closed && !o->m) {
    /* definitively empty */
    ar = elina_lincons0_array_make(1);
    ar.p[0] = elina_lincons0_make_unsat();
  }
  else {
    /* put non-oo constraint bounds only */
    opt_zones_mat_t* oz = o->closed ? o->closed : o->m;
    
    double *m = oz->mat;
    int s=0;
    unsigned short int i,j;
    unsigned short int n = o->dim+1;
    int size = n*n;
   
    ar = elina_lincons0_array_make(size);
    
    if(!oz->is_dense){
	    array_comp_list_t * acl = oz->acl;
	    
	    comp_list_t * cl = acl->head;
	    while(cl!=NULL){
		unsigned short int *ca = to_sorted_array(cl,n);
		unsigned short int comp_size = cl->size;
		for(i=0; i < comp_size; i++){
			unsigned short int i1 = ca[i]+1;
			
			for(j=0; j < comp_size; j++){
				unsigned short int j1 = ca[j]+1;
				if((i1==j1) ||(m[n*i1+j1]==INFINITY)){
					continue;
				}
				ar.p[s] = zones_lincons_of_bound(pr,i1,j1,m[n*i1+j1]);
				s++;
			}
			if(m[n*i1]!=INFINITY){
		      		ar.p[s] = zones_lincons_of_bound(pr,i1,0,m[n*i1]);
		     	 	s++;
	      		}
	      		if(m[i1]!=INFINITY){
	      			ar.p[s] = zones_lincons_of_bound(pr,0,i1,m[i1]);
	      			s++;
	      		}
		}
		free(ca);
		cl = cl->next;
	    }
	   
    }
    else{
	for (i=0;i<n;i++){
      		for (j=0;j<n;j++,m++) {
			if ((i==j) || (*m==INFINITY)) continue;

			ar.p[s] = zones_lincons_of_bound(pr,i,j,*m);
			s++;
      		}
	}
   }
    ar.size = s;
   // m = o->closed ? o->closed : o->m;
     
    if (pr->conv) zone_flag_conv;
  }
  
  return ar;
}

/**************************
	Bound Dimension
***************************/
elina_interval_t* opt_zones_bound_dimension(elina_manager_t* man,
				   opt_zones_t* o, elina_dim_t dim)
{
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_BOUND_DIMENSION,0);
  elina_interval_t* r = elina_interval_alloc();
  if(dim>=o->dim){
	elina_interval_free(r);
	return NULL;
  }
  unsigned short int n = o->dim+1;
  if (pr->funopt->algorithm>=0){ 
	opt_zones_mat_t * oz = o->closed ? o->closed : o->m;
	if(oz && !oz->is_dense){
		opt_zones_sparse_weak_closure(pr,o);	
	}
	else{
		opt_zones_cache_closure(pr,o);
	}
  }
  if (!o->closed && !o->m) {
    /* really empty */
    elina_interval_set_bottom(r);
  }
  else if (o->closed) {
    /* optimal in Q */
    opt_zones_mat_t * oz = o->closed;
    double *mm = oz->mat;
    if(!oz->is_dense){
	if(find(oz->acl,dim)==NULL){
		elina_interval_set_top(r);
	}
	else{
		zones_interval_of_bounds(r,mm[n*(dim+1)],mm[dim+1]);
	}
    }
    else{
    	zones_interval_of_bounds(r,mm[n*(dim+1)],mm[dim+1]);
    }
    if (zone_num_incomplete || o->intdim) zones_flag_incomplete;
    else if (pr->conv) zone_flag_conv;
  }
  else {
    /* not optimal */
    opt_zones_mat_t * oz = o->m;
    double *mm = oz->mat;
    if(!oz->is_dense){
	if(find(oz->acl,dim)==NULL){
		elina_interval_set_top(r);
	}
	else{
		zones_interval_of_bounds(r,mm[n*(dim+1)],mm[dim+1]);
	}
    }
    else{
    	zones_interval_of_bounds(r,mm[n*(dim+1)],mm[dim+1]);
    }
    zone_flag_algo;
  }
  return r;
}

/******************************
 Is dimension unconstrained
*****************************/
bool opt_zones_is_dimension_unconstrained(elina_manager_t* man, opt_zones_t* o,
				    elina_dim_t dim)
{
  opt_zones_internal_t* pr =
    opt_zones_init_from_manager(man,ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED,0);
  if(dim>=o->dim){
	return false;
  }
  
  if (!o->closed && !o->m)
    /* definitively empty */
    return false;
  else {
    unsigned short int n = o->dim+1;
    opt_zones_mat_t * oz = o->closed ? o->closed : o->m;
    double * m = oz->mat;
    size_t i;
    if(!oz->is_dense){
	
	array_comp_list_t *acl = oz->acl;
	comp_list_t * cl = find(acl,dim);
	if(cl==NULL){
		return true;
	}
	unsigned short int *ca = to_sorted_array(cl,o->dim);
	unsigned short int comp_size = cl->size;
	unsigned short int j;
	
	for(j=0; j< comp_size; j++){
		unsigned short int j1 = ca[j];
		if(j1==dim){
			if(m[j1+1]!=INFINITY || m[n*(j1+1)]!=INFINITY){
				return false;
			}
		}
		else{
			if(m[n*(j1+1) + dim+1]!=INFINITY){
				return false;
			}
			if(m[n*(dim+1) + j1+1]!=INFINITY){
				return false;
			}
		}
	} 
	free(ca);
		
    }
    else{
	
	if(m[dim+1]!=INFINITY || m[n*(dim+1)]!=INFINITY){
		return false;
	}
	for (i=0;i<o->dim;i++) {
		if(i==dim){
			continue;
		}
      		if (m[n*(i+1)+dim+1]!=INFINITY){
			 return false;
		}
      		if (m[n*(dim+1)+i+1]!=INFINITY){
			 return false;
		}
    	}
    }
    
    return true;
  }
}



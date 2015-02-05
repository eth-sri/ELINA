/*
	Copyright 2015 Software Reliability Lab, ETH Zurich

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

		http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
*/


#include "opt_oct_hmat.h"

/*******
Standard Meet Operator
****/

opt_oct_t* opt_oct_meet(ap_manager_t* man, bool destructive, opt_oct_t* o1, opt_oct_t* o2)
{
  
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_MEET,0);
  opt_oct_mat_t* oo;
  
   if((o1->dim != o2->dim) || (o1->intdim != o2->intdim))return NULL;
  if ((!o1->closed && !o1->m) || (!o2->closed && !o2->m))
    /* one argument is empty */
    return opt_oct_set_mat(pr,o1,NULL,NULL,destructive);
  else {
    int size = 2*(o1->dim)*(o1->dim + 1);
    opt_oct_mat_t * oo1 = o1->closed ? o1->closed : o1->m;
    opt_oct_mat_t * oo2 = o2->closed ? o2->closed : o2->m;
    oo = destructive ? oo1 : opt_hmat_alloc(size);
    meet_half(oo,oo1,oo2,o1->dim,destructive);
    /* optimal, but not closed */
    return opt_oct_set_mat(pr,o1,oo,NULL,destructive);
  }
}

/*****
Standard Join Operator
***/
opt_oct_t* opt_oct_join(ap_manager_t* man, bool destructive, opt_oct_t* o1, opt_oct_t* o2)
{
 opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_JOIN,0);
 
 if((o1->dim != o2->dim) || (o1->intdim != o2->intdim))return NULL;
 int size = 2*(o1->dim)*(o1->dim + 1);
 if (pr->funopt->algorithm>=0) {
   opt_oct_cache_closure(pr,o1);
   opt_oct_cache_closure(pr,o2);
 }
 if (!o1->closed && !o1->m) {
   if (!o2->closed && !o2->m)
     /* both empty */
     return opt_oct_set_mat(pr,o1,NULL,NULL,destructive);
   else{
     
     /* a1 empty, a2 not empty */
     return opt_oct_set_mat(pr,o1,opt_hmat_copy(o2->m,o1->dim),
			opt_hmat_copy(o2->closed,o1->dim),destructive);
     }
 }
 else if (!o2->closed && !o2->m){
   /* a1 not empty, a2 empty */
   return opt_oct_set_mat(pr,o1,o1->m,o1->closed,destructive);
  }
 else {
   /* not empty */
   opt_oct_mat_t* oo1 = o1->closed ? o1->closed : o1->m;
   opt_oct_mat_t* oo2 = o2->closed ? o2->closed : o2->m;
   opt_oct_mat_t * oo = destructive ? oo1 : opt_hmat_alloc(size);
   size_t i;
   man->result.flag_exact = false;
   join_half(oo,oo1,oo2,o1->dim,destructive);
   if (o1->closed && o2->closed) {
     /* result is closed and optimal on Q */
     if (num_incomplete || o1->intdim) flag_incomplete;
     return opt_oct_set_mat(pr,o1,NULL,oo,destructive);
   }
   else {
     /* not optimal, not closed */
     flag_algo;
     return opt_oct_set_mat(pr,o1,oo,NULL,destructive); 
   }
 }
}

/****
	Standard Widening Operator
***/

opt_oct_t* opt_oct_widening(ap_manager_t* man, opt_oct_t* o1, opt_oct_t* o2)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_WIDENING,0);
  int algo = pr->funopt->algorithm;
  opt_oct_t* r;
  
  if((o1->dim != o2->dim) || (o1->intdim != o2->intdim))return NULL;
  if (algo>=0) opt_oct_cache_closure(pr,o2);
  if (!o1->closed && !o1->m){
    /* o1 definitively closed */
    r = opt_oct_copy_internal(pr,o2);
    
  }
  else if (!o2->closed && !o2->m){
   /* o2 definitively closed */
    r = opt_oct_copy_internal(pr,o1);
    
  }
  else {
    /* work on the origial left matrix, not the closed cache! */
    opt_oct_mat_t * oo1 = o1->m ? o1->m : o1->closed;
    opt_oct_mat_t * oo2 = o2->closed ? o2->closed : o2->m;
    size_t i;
    r = opt_oct_alloc_internal(pr,o1->dim,o1->intdim);
    int size = 2*(r->dim)*(r->dim + 1);
    r->m = opt_hmat_alloc(size);
    //posix_memalign((void **)&(r->m),32,size*sizeof(double));
    if (algo==opt_oct_pre_widening || algo==-opt_oct_pre_widening) {
      /* degenerate hull: NOT A PROPER WIDENING, use with care */
	join_half(r->m,oo1,oo2,o1->dim,false);
    }
    else {
      /* standard widening */
        widening_half(r->m,oo1,oo2,o1->dim);
    }
  }
  return r;
}


/* ============================================================ */
/* Perturbation */
/* ============================================================ */

opt_oct_t* opt_oct_add_epsilon(ap_manager_t* man, opt_oct_t* o, ap_scalar_t* epsilon)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_WIDENING,2);
  opt_oct_t* r = opt_oct_alloc_internal(pr,o->dim,o->intdim);
  opt_oct_mat_t * oo;
  oo = o->m ? o->m : o->closed;
  if (oo) {
    double *m = oo->mat;
    size_t i;
    /* compute max of finite bounds */
    pr->tmp[0] = 0;
    int size = 2*(o->dim)*(o->dim+1);
    r->m = opt_hmat_alloc(size);
    free_array_comp_list(r->m->acl);
    double *mm = r->m->mat;
    if(!oo->is_dense){
	r->m->is_dense = false;
	r->m->ti = false;
	r->m->acl = copy_array_comp_list(oo->acl);
	comp_list_t * cl = oo->acl->head;
	while(cl != NULL){
		unsigned short int * ca = to_sorted_array(cl,o->dim);
		unsigned short int comp_size = cl->size;
		for(int i = 0; i < 2*comp_size; i++){
			unsigned short int i1 = (i%2==0)? 2*ca[i/2] : 2*ca[i/2]+1;
			for(int j = 0; j < 2*comp_size; j++){
				unsigned short int j1 = (j%2==0) ? 2*ca[j/2] : 2*ca[j/2] + 1;
				if(j1 > (i1|1)){
					break;
				} 
				int ind = j1 + (((i1 + 1)*(i1 + 1))/2);
				if(m[ind]==INFINITY){
					continue;
				} 
				if(m[ind] >= 0){
					pr->tmp[0] = max(pr->tmp[0],m[ind]);
				}
				else{
					pr->tmp[1] = -1*m[ind];
					pr->tmp[0] = max(pr->tmp[0],pr->tmp[1]);
				}
			}
		}
		cl =cl->next; 
		free(ca);
		
	}
	 /* multiply by epsilon */
	 opt_bound_of_scalar(pr,&pr->tmp[1],epsilon,false,false);
	 pr->tmp[0] = pr->tmp[0]*pr->tmp[1];
	/* enlarge bounds */
	cl = oo->acl->head;
	while(cl != NULL){
		unsigned short int * ca = to_sorted_array(cl,o->dim);
		unsigned short int comp_size = cl->size;
		for(int i = 0; i < 2*comp_size; i++){
			unsigned short int i1 = (i%2==0)? 2*ca[i/2] : 2*ca[i/2]+1;
			for(int j = 0; j < 2*comp_size; j++){
				unsigned short int j1 = (j%2==0) ? 2*ca[j/2] : 2*ca[j/2] + 1;
				if(j1 > (i1|1)){
					break;
				} 
				int ind = j1 + (((i1 + 1)*(i1 + 1))/2);
				mm[ind] = m[ind] + pr->tmp[0];
			}
		}
		cl =cl->next; 
		free(ca);		
	}
	/*for(int i = 0; i < o->dim; i++){
		//if(find(oo->acl,i)==NULL){
			int i1 = 2*i;
			int i2 = i1+1;
			int ind1 = i1 + (((i1+1)*(i1+1))/2);
			mm[ind1] = pr->tmp[0];
			int ind2 = i2 + (((i2+1)*(i2+1))/2);
			mm[ind2] = pr->tmp[0];
			comp_list_t * cl1 = create_comp_list();
			insert_comp(cl1,i);
			insert_comp_list(oo->acl,cl1);
		//}
	}*/
    }
    else{
	    r->m->is_dense = true;
	    r->m->ti = true;
	    for (i=0;i<size;i++) {
	      if (m[i]==INFINITY){
		 continue;
	      }
	      if (m[i]>=0){ 
		 pr->tmp[0] = max(pr->tmp[0],m[i]);
	      }
	      else {
		pr->tmp[1] = -1*m[i];
		pr->tmp[0] = max(pr->tmp[0],pr->tmp[1]);
	      }
	    }
	    /* multiply by epsilon */
	    opt_bound_of_scalar(pr,&pr->tmp[1],epsilon,false,false);
	    pr->tmp[0] = pr->tmp[0]*pr->tmp[1];
	    /* enlarge bounds */
	    #if defined(VECTOR)
	  	__m256d val = _mm256_set1_pd(pr->tmp[0]);
		for(i = 0; i < size/4; i++){
			__m256d op1 = _mm256_loadu_pd(m + i*4);
			__m256d op2 = _mm256_add_pd(val,op1);
			_mm256_storeu_pd(mm + i*4,op2);
		}  
	    #else
		for (i=0;i<(size/4)*4;i++){
	      		mm[i] = m[i] + pr->tmp[0];
	    	}
	    #endif
		for(i = (size/4)*4;i < size; i++){
			mm[i] = m[i] + pr->tmp[0];
		}
    }
     r->m->nni = oo->nni;
  }
 
  return r;
}

ap_abstract0_t* 
ap_abstract0_opt_oct_add_epsilon(ap_manager_t* man, 
			     ap_abstract0_t* a1, 
			     ap_scalar_t* epsilon)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_WIDENING,0);
  opt_oct_t* o = (opt_oct_t*) (a1->value);
  if(man->library!=a1->man->library){
	     return abstract0_of_opt_oct(man,opt_oct_alloc_top(pr,o->dim,o->intdim));
  }
  return abstract0_of_opt_oct(man,opt_oct_add_epsilon(man,o,epsilon));
}

opt_oct_t* opt_oct_add_epsilon_bin(ap_manager_t* man, opt_oct_t* o1, opt_oct_t* o2, 
			   ap_scalar_t* epsilon)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_WIDENING,2);
  opt_oct_t* r;
  if((o1->dim!=o2->dim) || (o1->intdim!=o2->intdim)){
	return NULL;
  }
  if (!o1->closed && !o1->m){
    /* a1 definitely empty */
    r = opt_oct_copy_internal(pr,o2);
  }
  else if (!o2->closed && !o2->m){
    /* a2 definitely empty */
    r = opt_oct_copy_internal(pr,o1);
  }
  else {
    opt_oct_mat_t * oo1 = o1->m ? o1->m : o1->closed;
    opt_oct_mat_t * oo2 = o2->m ? o2->m : o2->closed;
    double *m1 = oo1->mat;
    double *m2 = oo2->mat;
    int size = 2*(o1->dim)*(o1->dim+1);
    size_t i;
    r = opt_oct_alloc_internal(pr,o1->dim,o1->intdim);
    r->m = opt_hmat_alloc(size);
    /* get max abs of non +oo coefs in m2, times epsilon */
    pr->tmp[0] = 0;
    double *mm = r->m->mat;
    free_array_comp_list(r->m->acl);
    if(!oo1->is_dense){
	r->m->is_dense = false;
	r->m->ti = false;
	r->m->acl = copy_array_comp_list(oo1->acl);
	comp_list_t * cl = oo2->acl->head;
	while(cl != NULL){
		unsigned short int * ca = to_sorted_array(cl,o1->dim);
		unsigned short int comp_size = cl->size;
		
		for(int i = 0; i < 2*comp_size; i++){
			unsigned short int i1 = (i%2==0)? 2*ca[i/2] : 2*ca[i/2]+1;
			for(int j = 0; j < 2*comp_size; j++){
				unsigned short int j1 = (j%2==0) ? 2*ca[j/2] : 2*ca[j/2] + 1;
				if(j1 > (i1|1)){
					break;
				} 
				int ind = j1 + (((i1 + 1)*(i1 + 1))/2);
				if(m2[ind]==INFINITY){
					continue;
				} 
				if(m2[ind] >= 0){
					pr->tmp[0] = max(pr->tmp[0],m2[ind]);
				}
				else{
					pr->tmp[1] = -1*m2[ind];
					pr->tmp[0] = max(pr->tmp[0],pr->tmp[1]);
				}
			}
		}
		cl =cl->next; 
		free(ca);
		
	}
	 /* multiply by epsilon */
	 opt_bound_of_scalar(pr,&pr->tmp[1],epsilon,false,false);
	 pr->tmp[0] = pr->tmp[0] * pr->tmp[1];
	 /* enlarge unstable coefficients in o1 */
	 cl = oo1->acl->head;
	 while(cl!=NULL){
		unsigned short int * ca = to_sorted_array(cl,o1->dim);
		unsigned short int comp_size = cl->size;
		if(!oo2->ti){
			for(int i = 0; i < comp_size; i++){
				int i1 = ca[i];
				for(int j = 0; j <=i; j++){
					int j1 = ca[j];
					handle_binary_relation(m2,oo2->acl,i1,j1,o1->dim);
				}
			}
		}
		for(int i = 0; i < 2*comp_size; i++){
			int i1 = (i%2==0)? 2*ca[i/2] : 2*ca[i/2]+1;
			for(int j = 0; j < 2*comp_size; j++){
				int j1 = (j%2==0)? 2*ca[j/2] : 2*ca[j/2]+1;
				if(j1 > (i1 | 1)){
					break;
				}
				int ind = j1 + (((i1+1)*(i1+1))/2);
				if(m1[ind] < m2[ind]){
					mm[ind] = m2[ind] + pr->tmp[0];
				}
				else{
					mm[ind] = m1[ind];
				}
			}
		}
		cl = cl->next;
		free(ca);
	}
    }
    else{
	    r->m->is_dense = true;
	    r->m->ti = true;
	    for (i=0;i<size;i++) {
	      if (m2[i]==INFINITY){
		 continue;
	      }
	      if (m2[i]>=0){ 
		pr->tmp[0] = max(pr->tmp[0],m2[i]);
	      }
	      else {
		pr->tmp[1] = -1*m2[i];
		pr->tmp[0] = max(pr->tmp[0],pr->tmp[1]);
	      }
	    }
	    /* multiply by epsilon */
	    opt_bound_of_scalar(pr,&pr->tmp[1],epsilon,false,false);
	    pr->tmp[0] = pr->tmp[0] * pr->tmp[1];
	    /* enlarge unstable coefficients in o1 */
	    for (i=0;i<size;i++){
	      if (m1[i] < m2[i]) {
		mm[i] = m2[i] + pr->tmp[0];
	      }
	      else {
		mm[i] = m1[i];
	      }
	    }
    }
    r->m->nni = oo1->nni;
  }
  return r;
}

ap_abstract0_t* 
ap_abstract0_opt_oct_add_epsilon_bin(ap_manager_t* man, 
				 ap_abstract0_t* a1, 
				 ap_abstract0_t* a2, 
				 ap_scalar_t* epsilon)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_WIDENING,0);
  opt_oct_t* a = (opt_oct_t*) (a1->value);
  if((man->library!=a1->man->library) || (man->library!=a2->man->library)){
	return abstract0_of_opt_oct(man,opt_oct_alloc_top(pr,a->dim,a->intdim));
  }
  return abstract0_of_opt_oct(man,opt_oct_add_epsilon_bin(man,a1->value,a2->value,epsilon));
}



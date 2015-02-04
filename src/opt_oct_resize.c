#include "opt_oct_hmat.h"

/****************
Projection:
*****************/

opt_oct_t* opt_oct_forget_array(ap_manager_t* man,
			bool destructive, opt_oct_t* o,
			ap_dim_t* tdim, int size,
			bool project)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_FORGET_ARRAY,0);
  if (pr->funopt->algorithm>=0) opt_oct_cache_closure(pr,o);
  if (!o->closed && !o->m)
    /* definitively empty */
    return opt_oct_set_mat(pr,o,NULL,NULL,destructive);
  else {
    opt_oct_mat_t* oo = o->closed ? o->closed : o->m;
    int mat_size = 2*o->dim*(o->dim + 1);
    if (!destructive) oo = opt_hmat_copy(oo,o->dim);
    
    forget_array_avx_half(oo,tdim,o->dim,size,project);
    
    if (o->closed) {
      /* result is exact on Q, and closed if forget, not project */
      if (num_incomplete || o->intdim) flag_incomplete;
      if (project) return opt_oct_set_mat(pr,o,oo,NULL,destructive);
      else return opt_oct_set_mat(pr,o,NULL,oo,destructive);
    }
    else {
      /* not exact, not closed */
      flag_algo;
      return opt_oct_set_mat(pr,o,oo,NULL,destructive);
    }
  }
}



/******
Add Dimensions

*****/
opt_oct_t* opt_oct_add_dimensions(ap_manager_t* man,
			  bool destructive, opt_oct_t* o,
			  ap_dimchange_t* dimchange,
			  bool project)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_ADD_DIMENSIONS,0);
  opt_oct_mat_t* src = o->closed ? o->closed : o->m;
  opt_oct_mat_t* dst;
  size_t i, nb = dimchange->intdim+dimchange->realdim;
  opt_oct_t* r;
  if (!src) dst = NULL;
  else {
    /* check */
    
  for (i=0;i<nb;i++) {
    if(dimchange->dim[i] > o->dim) return NULL;
    if(i &&(dimchange->dim[i-1]> dimchange->dim[i])) return NULL;
   }
    /* insert variables */
    int dim = o->dim + nb;
    int size = 2*dim*(dim + 1); 
    dst = opt_hmat_alloc_top(dim);
    opt_hmat_addrem_dimensions(dst,src,dimchange->dim,
			   nb,1,o->dim,true);
    if(dst->is_dense){
	free_array_comp_list(dst->acl);
    }
    int count = dst->nni;
    /* set new variables to 0, if necessary */
    if (project) {
      for (i=0;i<nb;i++) {
	double *mm = dst->mat;
	size_t v = 2*(i+dimchange->dim[i]);
	mm[opt_matpos(v+1,v)] = 0;
	mm[opt_matpos(v,v+1)]  = 0;
	count = count + 2;
      }
    }
    dst->nni = count;
  }
  /* always exact, respect closure if embedding, not projecting */
  if (o->closed && !project) r = opt_oct_set_mat(pr,o,NULL,dst,destructive);
  else r = opt_oct_set_mat(pr,o,dst,NULL,destructive);
  r->dim += nb;
  r->intdim += dimchange->intdim;
  return r;
}


/********
Remove dimensions
********/

opt_oct_t* opt_oct_remove_dimensions(ap_manager_t* man,
			     bool destructive, opt_oct_t* o,
			     ap_dimchange_t* dimchange)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_REMOVE_DIMENSIONS,0);
  opt_oct_mat_t *src, *dst;
  size_t i, nb = dimchange->intdim+dimchange->realdim;
  opt_oct_t* r;
  if (pr->funopt->algorithm>=0) opt_oct_cache_closure(pr,o);
  src = o->closed ? o->closed : o->m;
  if (!src) dst = NULL;
  else {
    /* check */
     for (i=0;i<nb;i++) {
		if(dimchange->dim[i] >= o->dim) return NULL;
		if(i && (dimchange->dim[i-1]>=dimchange->dim[i]))return NULL;
     }
    /* remove variables */
    int dim = o->dim - nb;
    int size = 2*dim*(dim + 1);
    dst = opt_hmat_alloc(size);
    //posix_memalign((void **)&mm,32,size*sizeof(double));
    opt_hmat_addrem_dimensions(dst,src,dimchange->dim,
			   nb,1,o->dim,false);
    if(dst->is_dense){
	free_array_comp_list(dst->acl);
    }
  }

  if (o->closed) {
    /* result is exact on Q, and closed */
    if (num_incomplete || o->intdim) flag_incomplete;
    r = opt_oct_set_mat(pr,o,NULL,dst,destructive);
  }
  else {
    /* not exact, not closed */
    flag_algo;
    r = opt_oct_set_mat(pr,o,dst,NULL,destructive);
  }
  r->dim -= nb;
  r->intdim -= dimchange->intdim;
  return r;
}


/***********
Permute Dimensions
***********/

opt_oct_t* opt_oct_permute_dimensions(ap_manager_t* man,
			      bool destructive, opt_oct_t* o,
			      ap_dimperm_t* permutation)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_ADD_DIMENSIONS,0);
  opt_oct_mat_t* src = o->closed ? o->closed : o->m;
  opt_oct_mat_t* dst;
  if(permutation->size!=o->dim)return NULL;
  if (!src) dst = NULL;
  else {
    /* check (only bounds, not injectivity) */
    int i,j;    
    /* permuted copy */
     for (i=0;i<o->dim;i++){
          if(permutation->dim[i]>=o->dim)return NULL;
     }

    int size = 2*o->dim*(o->dim + 1);
    dst = opt_hmat_alloc(size);
    
    opt_hmat_permute(dst,src,o->dim,o->dim,permutation->dim);
    
    if(dst->is_dense){
	free_array_comp_list(dst->acl);
    }
  }
  /* always exact, respects closure */
  if (o->closed) return opt_oct_set_mat(pr,o,NULL,dst,destructive);
  else return opt_oct_set_mat(pr,o,dst,NULL,destructive);
}

opt_oct_t* opt_oct_expand(ap_manager_t* man,
		  bool destructive, opt_oct_t* o,
		  ap_dim_t dim,
		  size_t n)
{
  
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_EXPAND,0);
  opt_oct_mat_t* src = o->closed ? o->closed : o->m;
  size_t i, j;
  ap_dim_t pos = (dim < o->intdim) ? o->intdim : o->dim;
  opt_oct_mat_t* dst;
  opt_oct_t* r;
  if (!src) dst = NULL;
  else {
	
    /* insert n variables at pos */
    dst = opt_hmat_alloc_top(o->dim+n);
    opt_hmat_addrem_dimensions(dst,src,&pos,1,n,o->dim,true);
    
    double *mm = dst->mat;
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
			ini_relation(mm,pos+i,j,o->dim+n);
			if(j==dim){
				c = c->next;
				continue;
			}
		
			mm[opt_matpos2(2*(pos+i)  ,2*j)] = mm[opt_matpos2(2*dim  ,2*j)];
			mm[opt_matpos2(2*(pos+i)+1,2*j)] = mm[opt_matpos2(2*dim+1,2*j)];
			mm[opt_matpos2(2*(pos+i)  ,2*j+1)] = mm[opt_matpos2(2*dim  ,2*j+1)];
			mm[opt_matpos2(2*(pos+i)+1,2*j+1)] = mm[opt_matpos2(2*dim+1,2*j+1)];
			insert_comp(cl,j);
			c = c->next;
		}
		
	      }
	
	      int flag = 0;
	      comp_list_t * cp = find(dst->acl,pos+i);
	      if(cp!=NULL){
		if(cp==cj){
	      		mm[opt_matpos2(2*(pos+i)  ,2*(pos+i))] = mm[opt_matpos2(2*dim  ,2*(pos+i))];
	      		mm[opt_matpos2(2*(pos+i)+1,2*(pos+i)+1)] = mm[opt_matpos2(2*dim+1,2*(pos+i)+1)];
			//mm[opt_matpos2(2*(pos+i)  ,2*(pos+i))] = 0;
	      		//mm[opt_matpos2(2*(pos+i)+1,2*(pos+i)+1)] = 0;
			//flag = 1;
		}
		else{
		
			mm[opt_matpos2(2*(pos+i)  ,2*(pos+i))]  = INFINITY;
			mm[opt_matpos2(2*(pos+i)+1,2*(pos+i)+1)] = INFINITY;
			//flag = 1;
		}
	      }
	      else{
		
			mm[opt_matpos2(2*(pos+i)  ,2*(pos+i))]  = INFINITY;
			mm[opt_matpos2(2*(pos+i)+1,2*(pos+i)+1)] = INFINITY;
			//flag = 1;
	      }
	      /* copy unary constraints */
	      if(cj!=NULL){
	      	mm[opt_matpos2(2*(pos+i),2*(pos+i)+1)] = mm[opt_matpos2(2*dim,2*dim+1)];
		//flag = 1;
	      //}
	      //if(mm[opt_matpos2(2*dim+1,2*dim)] != INFINITY){
	      	mm[opt_matpos2(2*(pos+i)+1,2*(pos+i))] = mm[opt_matpos2(2*dim+1,2*dim)];
		//flag = 1;
	      }
	      else{
			mm[opt_matpos2(2*(pos+i),2*(pos+i)+1)] = INFINITY;
			mm[opt_matpos2(2*(pos+i)+1,2*(pos+i))] = INFINITY;
	      }
	      //if(cl->size > 0 || flag){
			insert_comp(cl,pos+i);
			insert_comp_list_with_union(dst->acl,cl,o->dim+n);
	      //}
	      //else{
		//	free_comp_list(cl);
	      //}
	    }
		 
  }
  else{
	free_array_comp_list(dst->acl);
	for (i=0;i<n;i++) {
		int src_ind,dest_ind;
		 #if defined(VECTOR)
			__m256d src;
	      /* copy binary constraints */
		      for (j=0;j<(2*dim)/4;j++) {
			//mm[opt_matpos(2*(pos+i),j)] = mm[opt_matpos(2*dim,j)];
			//int op = 2*(pos+i);
			//op = ((op+1)*(op+1))/2;
			dest_ind = opt_matpos(2*(pos+i),j*4);
			src_ind = opt_matpos(2*dim,j*4);
			src = _mm256_loadu_pd(mm+src_ind);
			_mm256_storeu_pd(mm+dest_ind,src); 
			dest_ind = opt_matpos(2*(pos+i)+1,j*4);
			src_ind = opt_matpos(2*dim+1,j*4);
			src = _mm256_loadu_pd(mm+src_ind);
			_mm256_storeu_pd(mm+dest_ind,src); 
			//count = count + 8;
			//mm[opt_matpos(2*(pos+i)+1,j)] = mm[opt_matpos(2*dim+1,j)];
		      }
	      #else
		   for(j= 0; j < ((2*dim)/4)*4;j++ ){
			mm[opt_matpos(2*(pos+i),j)] = mm[opt_matpos(2*dim,j)];
			mm[opt_matpos(2*(pos+i)+1,j)] = mm[opt_matpos(2*dim+1,j)]; 
			//count = count + 2;
		   }
	      #endif
	      for(j = ((2*dim)/4)*4;j < 2*dim; j++){
		mm[opt_matpos(2*(pos+i),j)] = mm[opt_matpos(2*dim,j)];
		mm[opt_matpos(2*(pos+i)+1,j)] = mm[opt_matpos(2*dim+1,j)];
		//count = count + 2;
	      }	

	      for (j=2*dim+2;j<2*(o->dim+n);j++) {
		mm[opt_matpos2(2*(pos+i),j)] = mm[opt_matpos(j^1,2*dim+1)];
		mm[opt_matpos2(2*(pos+i)+1,j)] = mm[opt_matpos(j^1,2*dim)];
		//count = count + 2;
	      }

	      /* copy unary constraints */
	      mm[opt_matpos2(2*(pos+i),2*(pos+i)+1)] = mm[opt_matpos2(2*dim,2*dim+1)];
	      mm[opt_matpos2(2*(pos+i)+1,2*(pos+i))] = mm[opt_matpos2(2*dim+1,2*dim)];
	      //count = count + 2;
	    }
     
     }
    
  }
  int dst_dim = o->dim+n;
  int dst_size = 2*dst_dim*(dst_dim+1);
  int src_size = 2*o->dim*(o->dim+1);
  /*  exact, generally not closed */
  dst->nni = min(dst_size,src->nni + dst_size-src_size);
  
  r = opt_oct_set_mat(pr,o,dst,NULL,destructive);
  r->dim += n;
  if (dim<o->intdim) r->intdim += n;
  
  return r;
}

opt_oct_t* opt_oct_fold(ap_manager_t* man,
		bool destructive, opt_oct_t* o,
		ap_dim_t* tdim,
		size_t size)
{
  
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_FOLD,0);
  opt_oct_mat_t* src;
  opt_oct_mat_t* dst;
  opt_oct_t* r;
  opt_oct_mat_t *oo;
  if (pr->funopt->algorithm>=0) opt_oct_cache_closure(pr,o);
  src = o->closed ? o->closed : o->m;
  if (!src) dst = NULL;
  else {
    /* check, assuming tdim[0..(size-1)] is strictly increasing */
    
    size_t i,j;
    if(size <= 0){
	return NULL;
    }
	
    for (i=1;i<size;i++) {
      if(tdim[i-1] >= tdim[i]){
	  return NULL;
      }
    }
	
    if(tdim[size-1] >= o->dim){
	return NULL;
    }
    double *m = src->mat;
    
    oo = opt_hmat_alloc(2*o->dim*(o->dim+1));
   
    opt_hmat_set_array(oo->mat,m,2*o->dim*(o->dim+1));
    oo->is_dense = src->is_dense;
    oo->nni = src->nni;
    if(!src->is_dense){
    	    oo->acl = copy_array_comp_list(src->acl);
    
	    /* merge binary constraints */
	    comp_list_t * cto = find(oo->acl,tdim[0]);
	    int count = 0;
	    for (j=0;j<tdim[0];j++) {
	      comp_list_t *cj = find(oo->acl,j);
	      if(cj!=cto){
		 continue;
	      }
	      double* mm1 = oo->mat + opt_matpos2(tdim[0]*2  ,2*j);
	      double* mm2 = oo->mat + opt_matpos2(tdim[0]*2+1,2*j);
	      double* mm3 = oo->mat + opt_matpos2(tdim[0]*2  ,2*j+1);
	      double* mm4 = oo->mat + opt_matpos2(tdim[0]*2+1,2*j+1);
	      for (i=1;i<size;i++) {
		comp_list_t * cti = find(oo->acl,tdim[i]);
		if(cj!=cti){
			*mm1 = INFINITY;
			*mm2 = INFINITY;
			*mm3 = INFINITY;
			*mm4 = INFINITY;
			/****
				Remove tdim[0] only when number of components is 2 with j as the other component
			*****/
			count++;
			break;
		}
		else{
			*mm1 = max(*mm1,m[opt_matpos2(tdim[i]*2  ,2*j)]);
			*mm2 = max(*mm2,m[opt_matpos2(tdim[i]*2+1,2*j)]);
			*mm3 = max(*mm3,m[opt_matpos2(tdim[i]*2  ,2*j+1)]);
			*mm4 = max(*mm4,m[opt_matpos2(tdim[i]*2+1,2*j+1)]);
		}
	      }
	    }
	    for (j=tdim[0]+1;j<o->dim;j++) {
	      comp_list_t *cj = find(oo->acl,j);
	      if(cj!=cto){
		 continue;
	      }
	      double * mm1 = oo->mat + opt_matpos2(tdim[0]*2  ,2*j);
	      double * mm2 = oo->mat + opt_matpos2(tdim[0]*2+1,2*j);
	      double* mm3 = oo->mat + opt_matpos2(tdim[0]*2  ,2*j+1);
	      double* mm4 = oo->mat + opt_matpos2(tdim[0]*2+1,2*j+1);
	      for (i=1;i<size;i++) {
		comp_list_t *cti = find(oo->acl,tdim[i]);
		if(cj!=cti){
			*mm1 = INFINITY;
			*mm2 = INFINITY;
			*mm3 = INFINITY;
			*mm4 = INFINITY;
			count++;
			break;
		}
		else{
			*mm1 = max(*mm1,m[opt_matpos2(tdim[i]*2  ,2*j)]);
			*mm2 = max(*mm2,m[opt_matpos2(tdim[i]*2+1,2*j)]);
			*mm3 = max(*mm3,m[opt_matpos2(tdim[i]*2  ,2*j+1)]);
			*mm4 = max(*mm4,m[opt_matpos2(tdim[i]*2+1,2*j+1)]);
		}
	      }
	    }
	
	    /* merge unary constraints */
	    {
	      double* mm1 = oo->mat + opt_matpos2(tdim[0]*2  ,tdim[0]*2+1);
	      double* mm2 = oo->mat + opt_matpos2(tdim[0]*2+1,tdim[0]*2  );
	      for (i=1;i<size;i++) {
		comp_list_t *cti = find(oo->acl,tdim[i]);
		if(cti==NULL){
			*mm1 = INFINITY;
			*mm2 = INFINITY;
			count++;
			break;
		}
		else{
			*mm1 = max(*mm1,m[opt_matpos2(tdim[i]*2,tdim[i]*2+1)]);
			*mm2 = max(*mm2,m[opt_matpos2(tdim[i]*2+1,tdim[i]*2)]);
		}
	      }
	    }
	    if((cto!=NULL) && (count==(cto->size))){
			remove_comp(cto,tdim[0]);
			if(cto->size==0){
				remove_comp_list(oo->acl,cto);
			}
	   }
     }
     else{
		
	    /* merge binary constraints */
	    for (j=0;j<2*tdim[0];j++) {
	      double* mm1 = oo->mat+opt_matpos2(tdim[0]*2,j);
	      double* mm2 = oo->mat+opt_matpos2(tdim[0]*2+1,j);
	      if((*mm1)!=INFINITY){
	      	for (i=1;i<size;i++) {
			*mm1 =  max(*mm1,m[opt_matpos2(tdim[i]*2,j)]);
			if((*mm1)==INFINITY){
				oo->nni--;
				break;
			}
	      	}
	      }
	
		if((*mm2)!=INFINITY){
			for (i=1;i<size;i++) {
				*mm2 = max(*mm2,m[opt_matpos2(tdim[i]*2+1,j)]);
				if((*mm2)==INFINITY){
					oo->nni--;
					break;
				}
			}
		 }
	       }


	     for (j=2*(tdim[0]+1);j<2*o->dim;j++) {
	       double* mm1 = oo->mat+opt_matpos2(tdim[0]*2,j);
	       double* mm2 = oo->mat+opt_matpos2(tdim[0]*2+1,j);
	       if((*mm1)!=INFINITY){
	      	  for (i=1;i<size;i++) {
			*mm1 = max(*mm1,m[opt_matpos2(tdim[i]*2,j)]);
			if((*mm1)==INFINITY){
				oo->nni--;
				break;
			}
		  }
	       	}
		if((*mm2)!=INFINITY){
		    for (i=1;i<size;i++) {
			*mm2 = max(*mm2,m[opt_matpos2(tdim[i]*2+1,j)]);
			if((*mm2)==INFINITY){
				oo->nni--;
				break;
			}
		    }
		 }
	      }
	    

	    /* merge unary constraints */
	    {
	      double* mm1 = oo->mat+opt_matpos2(tdim[0]*2,tdim[0]*2+1);
	      double* mm2 = oo->mat+opt_matpos2(tdim[0]*2+1,tdim[0]*2);
	      if((*mm1)!=INFINITY){
	      	for (i=1;i<size;i++) {
			*mm1 = max(*mm1,m[opt_matpos2(tdim[i]*2,tdim[i]*2+1)]);
			if((*mm1)==INFINITY){
				oo->nni--;
				break;
			}
		}
	      }
		if((*mm2)!=INFINITY){
		   for (i=1;i<size;i++) {
			*mm2 = max(*mm2,m[opt_matpos2(tdim[i]*2+1,tdim[i]*2)]);
			if((*mm2)==INFINITY){
				oo->nni--;
				break;
			}
		   }
		}
	      }
	
     }
    
    
    /* destroy all dimensions in tdim except the first one */
    dst = opt_hmat_alloc_top(o->dim-size+1);
    opt_hmat_addrem_dimensions(dst,oo,tdim+1,size-1,1,o->dim,false);
    if(dst->is_dense){
	free_array_comp_list(dst->acl);
    }
    double *mm = dst->mat;
    /* reset diagonal elements */
    mm[opt_matpos(tdim[0]*2,tdim[0]*2  )] = 0;
    mm[opt_matpos(tdim[0]*2+1,tdim[0]*2+1)] = 0;

    man->result.flag_exact = false;
  }
  
  if (o->closed) {
    /* result is optimal on Q, not closed */
    if (num_incomplete || o->intdim) flag_incomplete;
    r = opt_oct_set_mat(pr,o,dst,NULL,destructive);
  }
  else {
    /* not exact, not closed */
    flag_algo;
    r = opt_oct_set_mat(pr,o,dst,NULL,destructive);
  }
  r->dim -= size-1;
  if (tdim[0]<r->intdim) r->intdim -= size-1;
  opt_hmat_free(oo); 
   
  return r;
}


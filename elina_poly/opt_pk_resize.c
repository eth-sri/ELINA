/*
	Copyright 2016 Software Reliability Lab, ETH Zurich

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


/* ********************************************************************** */
/* opt_pk_resize.c: change and permutation of dimensions  */
/* ********************************************************************** */


#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_project.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"
#include "opt_pk_representation.h"
#include "opt_pk_user.h"
#include "opt_pk_constructor.h"



/* ********************************************************************** */
/* 		Add Dimensions with constraints*/
/* ********************************************************************** */

opt_pk_array_t* opt_pk_add_dimensions_cons(elina_manager_t* man,
			bool destructive,
			opt_pk_array_t* oa,
			elina_dimchange_t* dimchange,
			bool project)
{
  opt_pk_array_t* op;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_ADD_DIMENSIONS);
  opt_pk_internal_realloc_lazy(opk,oa->maxcols+dimchange->intdim+dimchange->realdim - 2);
  
  array_comp_list_t * acla = oa->acl;
  unsigned short int size = dimchange->intdim + dimchange->realdim;
  /* Return empty if empty */
  if(oa->is_bottom || !acla){
	 man->result.flag_best = man->result.flag_exact = true;
	 if (destructive){
	      	oa->maxcols += size;
	      	return oa;
	 }
	 else {
	      	return opt_pk_bottom(man,
			oa->maxcols + size - opk->dec,
			0);
	}
  }
  unsigned short int num_compa = acla->size;
  unsigned short int i, l,k;
  /*************************************
		Minimized the input
  *************************************/
  opt_pk_t ** poly_a = oa->poly;
  if(opk->funopt->algorithm>0){
	for(k=0; k < num_compa; k++){
	     opt_pk_t * oak = poly_a[k];
	     opt_poly_chernikova(man,oak,"add dimensions");
	     if (opk->exn){
	         opk->exn = ELINA_EXC_NONE;
		/* we can still operate on the available matrix */
    	     }
	     if(!oak->C && !oak->F){
		 man->result.flag_best = man->result.flag_exact = true;
		 if (destructive){
		     opt_poly_set_bottom(opk,oa);
		     oa->maxcols += dimchange->intdim + dimchange->realdim;
		     return oa;
		 }
		 else {
		      return opt_pk_bottom(man,
			     oa->maxcols + dimchange->intdim - 2,
			     0);
		 }
	     }
	}
  }
  /*******************************
		Create mapping for independent components
  ********************************/
  unsigned short int maxcols = oa->maxcols;
  unsigned short int * cmap = (unsigned short int *)calloc(maxcols+1,sizeof(unsigned short int));
  unsigned short int * ncmap = (unsigned short int *)calloc(size,sizeof(unsigned short int));
  l = 0;
  k = opk->dec;
  elina_dim_t * dim = dimchange->dim;
  for(i=k; i <=maxcols; i++){
	unsigned short int var = dim[l] + opk->dec;
	while((l < size) && (var==i)){
		ncmap[l] = k;
		l++;
		var = dim[l] + opk->dec;
		k++;
	}
	cmap[i] = k;
	//printf("cmap[%d]: %d\n",i,k);
	k++;
  }
  /*****************************************
	Handle independent components
  ******************************************/
  array_comp_list_t * acl = create_comp_list();
  if(project){
	for(l=0; l < size; l++){
		comp_list_t * cl = create_comp_list();
		insert_comp(cl,ncmap[l]);
		insert_comp_list(acl,cl);
	}
  }

  comp_list_t * cla = acla->head;
  while(cla != NULL){
	comp_list_t * cl = create_comp_list();
	comp_t * c = cla->head;
	while(c != NULL){
		unsigned short int numa = c->num;
		unsigned short int num = cmap[numa];
		//printf("numa: %d num: %d\n",numa,num);
		insert_comp(cl,num);
		c = c->next;
	}
	insert_comp_list(acl,cl);
  	cla = cla->next;
  }
  
  opt_pk_t ** poly = destructive ? poly_a : (opt_pk_t **)malloc(num_compa*sizeof(opt_pk_t *));
  if(!destructive){
	  for(k=0;  k< num_compa; k++){
		  unsigned short int k1 = num_compa - k - 1;
		  opt_pk_t * src = poly_a[k];
		  poly[k1] = opt_poly_alloc(src->intdim,src->realdim);
		  poly[k1]->nbeq = src->nbeq;
		  poly[k1]->C = src->C ? opt_matrix_copy(src->C) : NULL;
		  poly[k1]->F = src->F ? opt_matrix_copy(src->F) : NULL;
		  poly[k1]->satC = src->satC ? opt_satmat_copy(src->satC) : NULL;
		  poly[k1]->satF = src->satF ? opt_satmat_copy(src->satF) : NULL;
		  poly[k1]->nbline = src->nbline;
		  poly[k1]->status = src->status;
		  poly[k1]->is_minimized = src->is_minimized;
		  //op = opt_matrix_add_dimensions(opk, destructive, oa->C, dimchange, project);
	  }
  }
  
  if(project){
	unsigned short int num_comp = project ? num_compa + size : num_compa;
	if(destructive){
		 poly = (opt_pk_t **)realloc(poly_a, num_comp*sizeof(opt_pk_t*));
	}else{
		poly = (opt_pk_t **)realloc(poly, num_comp*sizeof(opt_pk_t*));
	}
	// Handle project with generators
	for(i = 0; i < size; i++){
		unsigned short int i1 = num_comp - i - 1;
		poly[i1] = opt_poly_alloc(1,0);
		poly[i1]->nbeq = 1;
		opt_matrix_t *mat = opt_matrix_alloc(1,maxcols + size,false);
		opt_numint_t **p = mat->p;
		opt_numint_t *pi = p[0];
		pi[0] = 0;
		pi[1] = 0;
		pi[ncmap[i]] = 1;
		poly[i1]->C = mat;
		poly[i1]->is_minimized = true;
	}
  }
  free(cmap);
  free(ncmap);
  if(destructive){
	op = oa;
	op->maxcols = maxcols + size;
	free(acla);
	op->acl = copy_array_comp_list(acl);
  }
  else{
  	op = opt_pk_array_alloc(poly,acl,maxcols+size);
	
  }
  return op;
}


/***********************************
	Add dimensions
************************************/
opt_pk_array_t* opt_pk_add_dimensions(elina_manager_t* man,
			bool destructive,
			opt_pk_array_t* oa,
			elina_dimchange_t* dimchange,
			bool project){
	#if defined (TIMING)
		start_timing();
	#endif

        opt_pk_array_t * op;
	op = opt_pk_add_dimensions_cons(man,destructive,oa,dimchange,project);
	#if defined (TIMING)
		record_timing(add_dimension_time);
	#endif
	return op;
}


/*******************************
	Remove Dimensions
********************************/
opt_pk_array_t* opt_pk_remove_dimensions(elina_manager_t* man,
			   bool destructive,
			   opt_pk_array_t* oa,
			   elina_dimchange_t* dimchange)
{
  //printf("remove start %p\n",oa);
  //elina_lincons0_array_t arr1 = opt_pk_to_lincons_array(man,oa);
  //elina_lincons0_array_fprint(stdout,&arr1,NULL);
  //fflush(stdout);
   #if defined(TIMING)
 	 start_timing();
  #endif  
  opt_pk_array_t* op;
  size_t dimsup;
  dimsup = dimchange->intdim+dimchange->realdim;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_REMOVE_DIMENSIONS);
  array_comp_list_t * acla = oa->acl;
  /* Return empty if empty */
  if(oa->is_bottom || !acla){
	man->result.flag_best = man->result.flag_exact = true;
      if (destructive){
          oa->maxcols -= dimsup;
	  #if defined(TIMING)
 	 		record_timing(remove_dimension_time);
  	  #endif
          return oa;
      }
      else {
          // Fix Me
	  #if defined(TIMING)
 	 		record_timing(remove_dimension_time);
  	  #endif
          return opt_pk_bottom(man,
                               oa->maxcols - dimsup - opk->dec,
                               dimchange->realdim);
      }	
    	
  }
  unsigned short int num_compa = acla->size;
  unsigned short int maxcols = oa->maxcols;
  unsigned short int k;
  opt_pk_t ** poly_a = oa->poly;
  for(k=0; k < num_compa; k++){
        opt_pk_t * oak = poly_a[k];
	if(opk->funopt->algorithm < 0){
		opt_poly_obtain_F(man,oak,"convert to gen");
	}
	else{
		opt_poly_chernikova(man,oak,"convert to gen");
	}
		
		//if overflow exception
	if(opk->exn){
	   opk->exn = ELINA_EXC_NONE;
	   if (!oak->F){
		man->result.flag_best = man->result.flag_exact = false;
		opt_pk_array_t * op = destructive ? oa : opt_pk_array_alloc(NULL,NULL,oa->maxcols);
		op->maxcols -= dimsup;
		opt_poly_set_top(opk,op);
		record_timing(remove_dimension_time);
		return op;
	   }
	}

	if(!oak->C && !oak->F){
	    man->result.flag_best = man->result.flag_exact = true;
   	    if (destructive){
		oa->maxcols -= dimsup;
	  	#if defined(TIMING)
 	 	     record_timing(remove_dimension_time);
  	        #endif
          	return oa;
    	    }
	    else{
		#if defined(TIMING)
 	 	     record_timing(remove_dimension_time);
  	  	#endif
          	return opt_pk_bottom(man,
                       oa->maxcols - dimsup - opk->dec,
                       dimchange->realdim);
	    }
	}
  }

  elina_dim_t * dima = dimchange->dim;
  
  
  /*********************************
	Handle independent components
  *********************************/
  unsigned short int * map = (unsigned short int *)calloc(maxcols, sizeof(unsigned short int));
 
  unsigned short int l = 0;
  for(k=opk->dec; k < maxcols; k++){
	unsigned short int var = dima[l] + opk->dec;
	if((l < dimsup) && (k==var)){
		map[k] = maxcols+1;
		l++;
	}
	else{
		map[k] = k - l;
	}
  }

  comp_list_t * cla = acla->head; 
  array_comp_list_t * acl = create_array_comp_list();
  unsigned short int num_comp = 0;
  while(cla!=NULL){
	comp_list_t * cl = create_comp_list();
	comp_t * c = cla->head;
	while(c!=NULL){
		unsigned short int numa = c->num;
		unsigned short int num = map[numa];
		if(num != (maxcols +1)){
			insert_comp(cl,num);
		}
		c = c->next;
	}
	if(cl->size){
		insert_comp_list(acl,cl);
		num_comp++;
	}
	else{
		free_comp_list(cl);
	}
	cla = cla->next;
  }
  /*********************************
	Remove variables from the blocks
  **********************************/
  cla = acla->head;
  if(destructive){
	k=0;
	while(k < num_compa){
		//printf("AA %d %d\n",k,num_compa);
		//fflush(stdout);
		unsigned short int comp_size = cla->size;
		unsigned short int * ca = to_sorted_array(cla,maxcols);
		opt_pk_t * oak = poly_a[k];
    		opt_pk_t * ot;
		/**************************
			Find the variables to remove for this component.
		****************************/
		elina_dim_t * ndim = (elina_dim_t *)calloc(comp_size, sizeof(elina_dim_t));
		unsigned short int size = 0;
		unsigned short int i,j;
        	for(i=0; i < comp_size; i++){
			unsigned short int num = ca[i];
			for(j=0; j < dimsup; j++){
				unsigned short int var = dima[j] + opk->dec;
				if(var==num){
					ndim[size] = i;
					size++;
				}
			}
		}
		
		if(size==comp_size){
			comp_list_t * tcla = cla->next;
			remove_comp_list(acla,cla);
			cla = tcla;
			free(ca);
			free(ndim);
			opt_pk_t * tpoly = poly_a[k];
			unsigned short int k1;
			for(k1=k; k1 < num_compa - 1; k1++){
				poly_a[k1] = poly_a[k1+1];
			}
			opt_poly_clear(tpoly);
			num_compa--;
			continue;
		}
		else if(size){
			//ot = opt_poly_alloc(oak->intdim, oak->realdim);
			if(oak->C){
			   opt_matrix_free(oak->C);
			   oak->C = NULL;
			}
			if(oak->satC){
			   opt_satmat_free(oak->satC);
			   oak->satC = NULL;
			}
			if(oak->satF){
			   opt_satmat_free(oak->satF);
			   oak->satF = NULL;
			}
			elina_dimchange_t dimchange1;
			dimchange1.dim = ndim;
			dimchange1.intdim = size;
			dimchange1.realdim = 0;
			poly_a[k]->F = opt_matrix_remove_dimensions(opk,false,oak->F,&dimchange1);
			poly_a[k]->nbline = oak->nbline;
			poly_a[k]->nbeq = 0;
			if (opk->funopt->algorithm>0){
			    opt_poly_chernikova(man,oak,"of the result");
			    if (opk->exn){
				opk->exn = ELINA_EXC_NONE;
			    }
		        }
			man->result.flag_best = man->result.flag_exact =
			     dimchange1.intdim==0;
		}
		cla = cla->next;
		k++;
	}
	
	array_comp_list_t * tmp = acl;
        oa->acl = copy_array_comp_list(acl);
	
	free_array_comp_list(tmp);
	oa->maxcols = oa->maxcols - dimsup;
	#if defined(TIMING)
 	 	record_timing(remove_dimension_time);
  	#endif
	return oa;
  }
  else{
	opt_pk_t ** poly = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t*));
	unsigned short int k1 = 0;
	for(k=0; k < num_compa; k++){
		unsigned short int comp_size = cla->size;
		unsigned short int * ca = to_sorted_array(cla,maxcols);
		opt_pk_t * oak = poly_a[k];
		
    		opt_pk_t * ot;
		/**************************
			Find the variables to remove for this component.
		****************************/
		elina_dim_t * ndim = (elina_dim_t *)calloc(comp_size, sizeof(elina_dim_t));
		unsigned short int size = 0;
		unsigned short int i,j;
        	for(i=0; i < comp_size; i++){
			unsigned short int num = ca[i];
			for(j=0; j < dimsup; j++){
				unsigned short int var = dima[j] + opk->dec;
				if(var==num){
					ndim[size] = i;
					size++;
				}
			}
		}
		if(size==comp_size){
			free(ca);
			free(ndim);
			cla = cla->next;
			continue;
		}
		else if(size){
		        elina_dimchange_t dimchange1;
			dimchange1.dim = ndim;
			dimchange1.intdim = size;
			dimchange1.realdim = 0;
			poly[k1] = opt_poly_alloc(oak->intdim - size,
				                      oak->realdim);
			poly[k1]->F = opt_matrix_remove_dimensions(opk,false,oak->F,&dimchange1);
			poly[k1]->nbeq = 0;
			poly[k1]->nbline = oak->nbline;
			if (opk->funopt->algorithm>0){
			    opt_poly_chernikova(man,poly[k1],"of the result");
			    if (opk->exn){
				opk->exn = ELINA_EXC_NONE;
			    }
			}
			man->result.flag_best = man->result.flag_exact =
				    dimchange1.intdim==0;
		}
		else{
			poly[k1] = opt_poly_alloc(oak->intdim,oak->realdim);
			opt_poly_copy(poly[k1],oak);
			poly[k1]->is_minimized = oak->is_minimized;
		}

		k1++;
		free(ca);
        	free(ndim);
		cla = cla->next;
	}
	poly = (opt_pk_t **)realloc(poly,k1*sizeof(opt_pk_t*));
	array_comp_list_t * res = copy_array_comp_list(acl);
	free_array_comp_list(acl);
	op = opt_pk_array_alloc(poly,res,maxcols - dimsup);
	#if defined(TIMING)
 	 	record_timing(remove_dimension_time);
   	#endif
	return op;
  }
  
}


/*******************************
	Permute Dimensions
********************************/
opt_pk_array_t* opt_pk_permute_dimensions(elina_manager_t* man,
			    bool destructive,
			    opt_pk_array_t* oa,
			    elina_dimperm_t* permutation)
{
  #if defined(TIMING)
 	 start_timing();
  #endif
  opt_pk_array_t* op;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_PERMUTE_DIMENSIONS);
  array_comp_list_t * acla = oa->acl;
  if(oa->is_bottom || !acla){
	if(destructive){
		#if defined(TIMING)
 	 		record_timing(permute_dimension_time);
  		#endif
		return oa;
	}
	else{
		#if defined(TIMING)
 	 		record_timing(permute_dimension_time);
  		#endif
		return opt_pk_bottom(man,oa->maxcols -2,0);
	}
  }
  unsigned short int num_comp = acla->size;
  opt_pk_t ** poly_a = oa->poly;
   unsigned short int k;
  /***************************************
	Minimize the input
  ***************************************/
  if(opk->funopt->algorithm>0){
	for(k=0; k < num_comp; k++){
	    if (opk->funopt->algorithm>0){
    		/* Minimize the argument */
    		opt_poly_chernikova(man,poly_a[k],"of the argument");
    		if (opk->exn){
     		    opk->exn = ELINA_EXC_NONE;
      		/* we can still operate on the available matrix */
    		}
  	    }
	}
  }
  unsigned short int **nca_arr = (unsigned short int **)calloc(num_comp,sizeof(unsigned short int **));
  unsigned short int maxcols = oa->maxcols;
  /**********************
	Handle the independent components
  ***********************/
  array_comp_list_t * acl;  
  acl = create_array_comp_list();
  comp_list_t * cla = acla->head;
  elina_dim_t * dima = permutation->dim;
  k = 0;
  while(cla!=NULL){
	comp_list_t *cl = create_comp_list();
	comp_t * c = cla->head;
	while(c!=NULL){
		unsigned short int numa = c->num - opk->dec;
		unsigned short int num = dima[numa] + opk->dec; 
		//printf("numa: %d num: %d\n",numa,dima[numa]);
		insert_comp(cl,num);
		c = c->next;
	}
	//print_comp_list(cl,oa->maxcols);
	insert_comp_list(acl,cl);
	nca_arr[k] = to_sorted_array(cl,maxcols);
	k++;
	cla = cla->next;
  }
  
  opt_pk_t ** poly;
  poly = destructive ? poly_a : (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
  if(!destructive){
	for(k=0; k < num_comp; k++){
		opt_pk_t * src = poly_a[k];
		unsigned short int k1 = num_comp - k - 1;
		poly[k1] = opt_poly_alloc(src->intdim,src->realdim);
		poly[k1]->nbeq = src->nbeq;
		poly[k1]->nbline = src->nbline;
		poly[k1]->satC = src->satC ? opt_satmat_copy(src->satC) : NULL;
		poly[k1]->satF = src->satF ? opt_satmat_copy(src->satF) : NULL;
		poly[k1]->status = src->status;
		poly[k1]->is_minimized = src->is_minimized;
	}
	op = opt_pk_array_alloc(poly,acl,maxcols);
  }
  else{
	//
        op = oa;
	op->acl = copy_array_comp_list(acl);
  }
  cla = acla->head;
  
  for(k=0; k < num_comp; k++){
      unsigned short int k1 = num_comp - k - 1;
      unsigned short int k2 = destructive ? k : k1;
	  unsigned short int comp_size = cla->size;
	  unsigned short int * ca = to_sorted_array(cla,maxcols);
	  unsigned short int * nca = nca_arr[k];
	  opt_pk_t *oak = poly_a[k];
	  man->result.flag_best = man->result.flag_exact = true;
	  elina_dim_t * dim = (elina_dim_t *)calloc(comp_size, sizeof(elina_dim_t));
	  unsigned short int i,j;
	  for(i=0; i < comp_size; i++){
		unsigned short int nvar = nca[i] - opk->dec;
		for(j=0; j < comp_size; j++){
			unsigned short int var = ca[j];
			if(dima[var-opk->dec]==nvar){
				//printf("nvar: %d %d %d");
				dim[j] = i;
				break;
			}
		}
		//dim[i] = permutation->dim[var-opk->dec];
		
	  }
	  if(oak->C){
	   	poly[k2]->C = opt_matrix_permute_dimensions(opk,destructive,oak->C,dim);
	  }
	  if(oak->F){
		 poly[k2]->F = opt_matrix_permute_dimensions(opk,destructive,oak->F,dim);
	  }
	  cla = cla->next;
	  free(ca);
	  free(dim);
	  free(nca);
  }
  if(destructive){
	free_array_comp_list(acla);
  }
  free(nca_arr);
  #if defined(TIMING)
 	 record_timing(permute_dimension_time);
  #endif
  return op;
}

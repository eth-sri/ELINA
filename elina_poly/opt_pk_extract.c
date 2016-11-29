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
/* opt_pk_extract.c: property extraction */
/* ********************************************************************** */

#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"
#include "opt_pk_user.h"
#include "opt_pk_representation.h"


/* Bounding the value of a dimension in a matrix of generators. */


/* ====================================================================== */
/* Bounding the value of a dimension in a polyhedra */
/* ====================================================================== */

elina_interval_t* opt_pk_bound_dimension_gen(elina_manager_t* man,
				  opt_pk_array_t* oa,
				  elina_dim_t dim){
  itv_t itv;
  elina_interval_t* interval;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_BOUND_DIMENSION);

  interval = elina_interval_alloc();
  elina_interval_reinit(interval,ELINA_SCALAR_MPQ);

  if(oa->is_bottom || !oa->acl){
	elina_interval_set_bottom(interval);
	man->result.flag_exact = man->result.flag_best = true;
        return interval;
  }
  elina_dim_t var = dim+2;
  comp_list_t * cl = find(oa->acl,var);
  if(cl==NULL){
	elina_interval_set_top(interval);
	return interval;
  }
  short int k = find_index(oa->acl,cl);
  opt_pk_t ** poly_a = oa->poly;
  opt_pk_t * oak = poly_a[k];
  if (opk->funopt->algorithm>0)
    opt_poly_chernikova(man,oak,NULL);
  else
    opt_poly_obtain_F(man,oak,NULL);

  if (opk->exn){
    opk->exn = ELINA_EXC_NONE;
    elina_interval_set_top(interval);
    return interval;
  }

  if (!oak->F){ /* po is empty */
    elina_interval_set_bottom(interval);
    man->result.flag_exact = man->result.flag_best = true;
    return interval;
  }

  itv_init(itv);
  unsigned short int * ca = to_sorted_array(cl,oa->maxcols);
  unsigned short int i, new_dim=0;
  unsigned short int comp_size = cl->size;
  for(i=0; i < comp_size; i++){
	if(ca[i]==var){
		new_dim = i;
		break;
	}
  }

  opt_generator_bound_dimension(opk,itv,new_dim,oak->F);
  free(ca);
  elina_interval_set_itv(opk->itv,interval, itv);
  itv_clear(itv);
  man->result.flag_exact = man->result.flag_best = 
    new_dim<oak->intdim ? false : true;
  return interval;

}


elina_interval_t* opt_pk_bound_dimension(elina_manager_t* man,
				  opt_pk_array_t* oa,
				  elina_dim_t dim){
	#if defined(TIMING)
		start_timing();
	#endif

	return opt_pk_bound_dimension_gen(man,oa,dim);
	
	#if defined(TIMING)
		record_timing(bound_dimension_time);
	#endif
}

/* ====================================================================== */
/* Bounding the value of a linear expression in a polyhedra */
/* ====================================================================== */



/* ====================================================================== */
/* Converting to a set of constraints */
/* ====================================================================== */

elina_lincons0_array_t opt_pk_to_lincons_array(elina_manager_t* man,
					opt_pk_array_t* oa)
{
  //printf("to lincons: %p",oa);
  //fflush(stdout);
  elina_lincons0_array_t array;
  size_t i,l;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_TO_LINCONS_ARRAY);

  man->result.flag_exact = man->result.flag_best = true;
  array_comp_list_t * acla = oa->acl;
  if(oa->is_bottom || !acla){
    array = elina_lincons0_array_make(1);
    array.p[0] = elina_lincons0_make_unsat();
    return array;
  }
  if (acla->size==0){
    //opk->exn = ELINA_EXC_NONE;
    man->result.flag_exact = man->result.flag_best = true;
    array = elina_lincons0_array_make(0);
    return array;
  }
  
  unsigned short int num_compa = acla->size;
 
  unsigned short int k;
  opt_pk_t ** poly_a = oa->poly;
  size_t nbcons = 0;
  for(k=0; k < num_compa; k++){
	opt_pk_t * oak = poly_a[k];
	opt_poly_chernikova(man,oak,"to lincons array");
	if(opk->exn){
		opk->exn = ELINA_EXC_NONE;
    		man->result.flag_exact = man->result.flag_best = false;
    		array = elina_lincons0_array_make(0);
    		return array;
  	}
  	if (!oak->C){ /* po is empty */
	    array = elina_lincons0_array_make(1);
	    array.p[0] = elina_lincons0_make_unsat();
	    return array;
 	}
	opt_matrix_t* oc = oak->C;
	nbcons = nbcons + oc->nbrows;
  }

  array = elina_lincons0_array_make(nbcons);
  comp_list_t * cla = acla->head;
  l=0; 
  for(k=0; k < num_compa; k++){
	opt_pk_t * oak = poly_a[k];
	opt_matrix_t* oc = oak->C;
	//opt_matrix_sort_rows(opk,oc);
	unsigned short int * ca = to_sorted_array(cla,oa->maxcols);
	for (i=0; i<oc->nbrows; i++){
    		if (! opt_vector_is_dummy_constraint(opk,
				     oc->p[i], oc->nbcolumns)){
			
      			array.p[l] = opt_lincons0_of_vector(opk, oc->p[i], ca, cla->size, oa->maxcols);
     			 l++;
    		}
  	}
	free(ca);
	cla = cla->next;
  } 
  
  array.size = l;
  return array;
}

//elina_tcons0_array_t pk_to_tcons_array(elina_manager_t* man,
//				    opt_pk_t* op)
//{
  //return elina_generic_to_tcons_array(man,po);
//}

/* ====================================================================== */
/* Converting to a box */
/* ====================================================================== */

elina_interval_t** opt_pk_to_box(elina_manager_t* man,
			  opt_pk_array_t* op)
{
  #if defined(TIMING)
       start_timing();
  #endif
  elina_interval_t** interval;
  itv_t* titv;
  unsigned short int i,dim;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_TO_BOX);
   unsigned short int maxcols = op->maxcols;
  dim = maxcols - 2;
  interval = elina_interval_array_alloc(dim);
  array_comp_list_t * acl = op->acl;
  if(op->is_bottom || !acl){
	for (i=0; i<dim; i++){
      		elina_interval_set_bottom(interval[i]);
    	}
	#if defined(TIMING)
 	 	record_timing(poly_to_box_time);
  	#endif
	return interval;
  }
  else{
	for (i=0; i<dim; i++){
     		     elina_interval_set_top(interval[i]);
    	}
  }

  
  unsigned short int num_comp = acl->size;
  unsigned short int k;
  opt_pk_t ** poly = op->poly;
  /***********************************
	Minimize the input
  ************************************/
  
  for(k=0; k < num_comp; k++){
      if(opk->funopt->algorithm>=0){
	 opt_poly_chernikova(man,poly[k],"to box");
      }
      else{
	  opt_poly_obtain_F(man,poly[k],NULL);
      }
      if (opk->exn){
    	  opk->exn = ELINA_EXC_NONE;
    	  man->result.flag_exact = man->result.flag_best = false;
	  #if defined(TIMING)
 	       record_timing(poly_to_box_time);
  	  #endif
	  return interval;
      }
		
      if(!poly[k]->F){
	 for (i=0; i<dim; i++){
      	      elina_interval_set_bottom(interval[i]);
    	 }
	 #if defined(TIMING)
 	      record_timing(poly_to_box_time);
  	 #endif
	 return interval;
      }
	
  }
  comp_list_t * cl = acl->head;
  for(k=0; k < num_comp; k++){
  	opt_pk_t * ok = poly[k];
	unsigned short int comp_size = cl->size;
	unsigned short int * ca = to_sorted_array(cl,maxcols);
	titv = opt_generator_to_box(opk,ok->F);	
	for (i=0; i<comp_size; i++){
	     unsigned short int var = ca[i] - opk->dec;
      	     elina_interval_set_itv(opk->itv,interval[var],titv[i]);
    	}
	itv_array_free(titv,comp_size);
	free(ca);
	cl = cl->next;
  }
  man->result.flag_exact = man->result.flag_best = true;

  #if defined(TIMING)
       record_timing(poly_to_box_time);
  #endif
  return interval;
}



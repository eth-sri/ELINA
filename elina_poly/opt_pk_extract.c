/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
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

/* ********************************************************************** */
/* opt_pk_extract.c: property extraction */
/* ********************************************************************** */

#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"
#include "opt_pk_user.h"
#include "opt_pk_meetjoin.h"
#include "opt_pk_representation.h"
#include "opt_pk_test.h"


/* Bounding the value of a dimension in a matrix of generators. */


/* ====================================================================== */
/* Bounding the value of a dimension in a polyhedra */
/* ====================================================================== */

elina_interval_t* opt_pk_bound_dimension_gen(elina_manager_t* man,
				  opt_pk_array_t* oa,
				  elina_dim_t dim){
 
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

  unsigned short int * ca = to_sorted_array(cl,oa->maxcols);
  unsigned short int i, new_dim=0;
  unsigned short int comp_size = cl->size;
  for(i=0; i < comp_size; i++){
	if(ca[i]==var){
		new_dim = i;
		break;
	}
  }

  opt_generator_bound_dimension(opk,interval,new_dim,oak->F);
  free(ca);
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
/* Bounding the value of a linear expression constrained under oa */
/* ====================================================================== */

elina_interval_t* opt_pk_bound_linexpr(elina_manager_t* man,
				opt_pk_array_t* oa,
				elina_linexpr0_t* expr)
{
  bool exact;
  elina_interval_t* interval;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_BOUND_LINEXPR);

  array_comp_list_t * acla = oa->acl;
 
  interval = elina_interval_alloc();
  elina_interval_reinit(interval,ELINA_SCALAR_MPQ);
  if(oa->is_bottom || !acla){
    elina_interval_set_bottom(interval);
    man->result.flag_exact = man->result.flag_best = true;
    return interval;
  }
  
  unsigned short int num_compa = acla->size;
  opt_pk_t **poly_a = oa->poly;
  unsigned short int k, maxcols = oa->maxcols;
  for(k=0; k < num_compa; k++){
      opt_pk_t * oak = poly_a[k];
      if (opk->funopt->algorithm>0){
	    opt_poly_chernikova(man,oak,"Bound Linexpr");
      }
      else{
	    opt_poly_obtain_F(man,oak,"Bound Linexpr");
      }
      if (opk->exn){
	    opk->exn = ELINA_EXC_NONE;
	    elina_interval_set_top(interval);
	    return interval;
      }
      if (!oak->F){ /* oa is empty */
	    elina_interval_set_bottom(interval);
	    man->result.flag_exact = man->result.flag_best = true;
	    return interval;
      }
  }

  if(expr->size==0){
	elina_coeff_t *cst = &expr->cst;
	if(cst->discr==ELINA_COEFF_SCALAR){
		elina_scalar_set(interval->inf,cst->val.scalar);
		elina_scalar_set(interval->sup,cst->val.scalar);
	}
	else{
		elina_scalar_set(interval->inf,cst->val.interval->inf);
		elina_scalar_set(interval->sup,cst->val.interval->sup);
	}
	return interval;
  }
  // compute the block of variables B in the linear expression
  comp_list_t * clb = linexpr0_to_comp_list(opk,expr);
  array_comp_list_t *aclb = create_array_comp_list();
  insert_comp_list(aclb,clb);
  // compute the fusion of B with acla
  array_comp_list_t * acl = union_array_comp_list(acla,aclb,maxcols);
 
  // compute the blocks of acla that intersect with B
  char * map = create_map(clb,maxcols);
  comp_list_t * cla = acla->head;
  char * intersect_map = (char *)calloc(num_compa,sizeof(char));
  size_t * num_vertex_a = (size_t *)calloc(num_compa,sizeof(size_t));
  size_t num_vertex=0;
  size_t nbgen=0;
  bool flag = false;
  for(k=0; k < num_compa; k++){
      if(is_disjoint_with_map(cla,map)){
          cla = cla->next;
	  continue;
      }
      flag = true;
      opt_pk_t * oak = poly_a[k];
      intersect_map[k] = 1;
      num_vertex_a[k] = opt_generator_rearrange(oak->F,oak->satF);
      if(!num_vertex){
	 num_vertex = num_vertex_a[k];
      }
      else{
	 num_vertex = num_vertex*num_vertex_a[k];
      }
      nbgen = nbgen + oak->F->nbrows - num_vertex_a[k];
      cla = cla->next;
  }
  if(!flag){
	elina_interval_set_top(interval);
	free(map);
	free(intersect_map);
	free(num_vertex_a);
	return interval;
  }
  // compute the factor by fusing all intersecting blocks
  comp_list_t * cl = acl->head;
  unsigned short int num_comp = acl->size;
  unsigned short int comp_size=0;
  unsigned short int * ca = NULL;
  elina_linexpr0_t * new_expr=NULL; 
  for(k=0; k < num_comp; k++){
        if(!is_disjoint(cl,clb,maxcols)){
		comp_size = cl->size;
		ca = to_sorted_array(cl,maxcols);
		new_expr = copy_linexpr0_with_comp_list(opk,expr,ca, comp_size);
		break;
	}       
        cl = cl->next;
  }

  opt_matrix_t *F = opt_matrix_alloc(nbgen + num_vertex, comp_size + 2, false);
  fuse_generators_intersecting_blocks(F,poly_a,acla,ca,num_vertex_a,intersect_map,maxcols);
  free(ca);
  free(map);
  free(intersect_map);

  /* we fill the vector with the expression, taking lower bound of the interval
     constant */
  elina_rat_t *inf = (elina_rat_t*)malloc(sizeof(elina_rat_t));
  elina_rat_t *sup = (elina_rat_t*)malloc(sizeof(elina_rat_t));
  opt_generator_bound_elina_linexpr0(opk,inf,sup,new_expr,F);
  elina_scalar_set_elina_rat(interval->inf,inf);
  elina_scalar_set_elina_rat(interval->sup,sup);
  man->result.flag_exact = man->result.flag_best = 
    ( (opk->funopt->flag_exact_wanted || opk->funopt->flag_best_wanted) &&
      elina_linexpr0_is_real(new_expr,comp_size) ) ? 
    exact :
    false;
  elina_linexpr0_clear(new_expr);
  free(new_expr);
  free_array_comp_list(aclb);
  free_array_comp_list(acl);
  free(inf);
  free(sup);
  free(num_vertex_a);
  opt_matrix_free(F);
  return interval;
}



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
  elina_interval_t** tinterval;
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
	tinterval = opt_generator_to_box(opk,ok->F);	
	for (i=0; i<comp_size; i++){
	     unsigned short int var = ca[i] - opk->dec;
	     elina_interval_set(interval[var],tinterval[i]);
    	}
	elina_interval_array_free(tinterval,comp_size);
	free(ca);
	cl = cl->next;
  }
  man->result.flag_exact = man->result.flag_best = true;

  #if defined(TIMING)
       record_timing(poly_to_box_time);
  #endif
  return interval;
}



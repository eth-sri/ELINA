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

#include "opt_pk_config.h"
#include "opt_pk.h"
#include "opt_pk_matrix.h"
#include "opt_pk_test.h"
#include "opt_pk_vector.h"
#include "opt_pk_representation.h"
#include "opt_pk_user.h"
#include "opt_pk_constructor.h"
#include "opt_pk_widening.h"
#include "opt_pk_meetjoin.h"


/* ====================================================================== */
/* Emptiness test */
/* ====================================================================== */

bool opt_pk_is_bottom(elina_manager_t* man, opt_pk_array_t* op)
{
  start_timing();
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_IS_BOTTOM);
  array_comp_list_t * acl = op->acl;
  if(op->is_bottom || !acl){
	man->result.flag_exact = man->result.flag_best = true;
    	//opt_poly_set_bottom(opk,op);
	op->is_bottom = true;
	//printf("is bottom end true 1\n");
  	//fflush(stdout);
	record_timing(is_bottom_time);
	return true;
  }
 
  opt_pk_t ** poly = op->poly;
  unsigned short int num_comp = acl->size;
  unsigned short int k;
  for(k=0; k < num_comp; k++){
	opt_pk_t * op_k = poly[k];
	if(!op_k->C && !op_k->F){
		record_timing(is_bottom_time);
		op->is_bottom = true;
		return true;
	}
	if(op_k->F){
		continue;
	}
	else{
	     if (opk->funopt->algorithm>=0){
     			opt_poly_chernikova(man,op_k,NULL);
			if (opk->exn){
				man->result.flag_exact = man->result.flag_best = false;
				opk->exn = ELINA_EXC_NONE;
	 			record_timing(is_bottom_time);
				return false;
      			}
			man->result.flag_exact = man->result.flag_best = op_k->intdim>0 && op_k->F ? false : true;
			if(!op_k->F){
	 			 record_timing(is_bottom_time);
				 op->is_bottom = true;
				 return true;
			}
    	     }
	     else {
		if(!op_k->C){
			man->result.flag_exact = man->result.flag_best = false;
	       		record_timing(is_bottom_time);
			op->is_bottom = true;
	      		return true;
		}
	    }
       }
  }
  record_timing(is_bottom_time);
  return false;
}




void swelina_index(unsigned short int *idx1, unsigned short int *idx2){
  unsigned short int *tmp = idx1;
  idx1 = idx2;
  idx2 = tmp;
}

unsigned short int *sorted_index (opt_numint_t *arr, unsigned short int n)
{
  unsigned short int *idx, i, j;
 
  idx = malloc (sizeof (unsigned short int) * n);
 
  for (i=0; i<n; i++)
  {
    idx[i] = i;
  }
 
  for (i=0; i<n; i++)
  {
    for (j=i+1; j<n; j++)
    {
      if (arr[idx[i]] > arr[idx[j]])
      {
        swelina_index(&idx[i], &idx[j]);
      }
    }
  }
 
  return idx;
}


bool opt_pk_is_top(elina_manager_t* man, opt_pk_array_t* o){
	//printf("is top %p\n",o);
	//elina_lincons0_array_t arr1 = opt_pk_to_lincons_array(man,o);
  	//elina_lincons0_array_fprint(stdout,&arr1,NULL);
	//fflush(stdout);
	if(o==NULL){
		return false;
	}
	#if defined(TIMING)
 	 	start_timing();
   	#endif
	opt_pk_internal_t *opk = opt_pk_init_from_manager(man,ELINA_FUNID_IS_TOP);
	man->result.flag_exact = man->result.flag_best = true;
	array_comp_list_t * acl = o->acl;
	if(o->is_bottom || !acl){
		#if defined(TIMING)
 	 		record_timing(is_top_time);
   		#endif
		return false;
	}
	if(acl->size==0){
		#if defined(TIMING)
 	 		record_timing(is_top_time);
   		#endif
		return true;
	}
        unsigned short int num_comp = acl->size;
	unsigned short int k;
	opt_pk_t ** poly = o->poly;
	//sort components in increasing order
	opt_numint_t * comps = (opt_numint_t *)malloc(num_comp*sizeof(opt_numint_t));
	for(k=0; k < num_comp; k++){
		 if(poly[k]->C && poly[k]->F){
			comps[k] = 0;
		 }
		 else if(poly[k]->C){
			comps[k] = poly[k]->C->nbrows*poly[k]->C->nbcolumns;
		 }
		 else if(poly[k]->F){
			comps[k] = poly[k]->F->nbrows*poly[k]->F->nbcolumns;
		 }
		 else{
			free(comps);
			return false;
		 }
	}

	unsigned short int *comp_ind = sorted_index(comps,num_comp);
	free(comps);
        for(k=0; k < num_comp; k++){
		unsigned short int k1 = comp_ind[k];
		opt_pk_t *ok = poly[k1];
		if (opk->funopt->algorithm>=0){
    		    opt_poly_chernikova(man,ok,NULL);
		}
		if (!ok->C && !ok->F){
		     free(comp_ind);
		     #if defined(TIMING)
			 record_timing(is_top_time);
		     #endif    
		     return false;
		}
		else if (ok->C && ok->F){
		      if(ok->C->nbrows != (opk->dec - 1)){
			 free(comp_ind);
			 #if defined(TIMING)
			     record_timing(is_top_time);
			 #endif
			 return false;
		      }
		}
		else {	    
		      man->result.flag_exact = man->result.flag_best = false;
		      free(comp_ind);
		      #if defined(TIMING)
			  record_timing(is_top_time);
		      #endif
		      return false;
	        }
	}
	free(comp_ind);
	#if defined(TIMING)
 	 	record_timing(is_top_time);
   	#endif
	
	return true;
}


unsigned short int * insertion_sort(unsigned short int *arr, unsigned short int size){
	unsigned short int * res = (unsigned short int *)malloc(size*sizeof(unsigned short int));
	unsigned short int i,j;
	for(i = 0; i < size; i++){
		res[i] = i;
	}
	for (i=0; i<size; i++){
    		for (j=i+1; j<size; j++){
      			if (arr[res[i]] > arr[res[j]]){
        			unsigned short int tmp  = res[i];
				res[i] = res[j];
				res[j] = tmp;
      			}
    		}
  	}
	return res;
}



bool opt_generators_sat_vector(opt_pk_internal_t* opk, opt_matrix_t* F,  
			       opt_numint_t* tab, bool is_strict)
{
  size_t i;
  if (opt_numint_sgn(tab[0])==0){
    /* 1. constraint is an equality */
    for (i=0; i<F->nbrows; i++){	
      //opt_numint_t * fpi = opt_melina_vector(F->p[i],mapF,comp_size,F->nbcolumns);
      
      opt_numint_t prod = opt_vector_product_strict(opk, F->p[i], tab,F->nbcolumns);
      if(opk->exn){
	opk->exn = ELINA_EXC_NONE;
	return false;
      }	
      if (opt_numint_sgn(prod)) {
		//free(fpi);
		return false;
      }
	//free(fpi);
    }
    
    return true;
  }
  else {
    /* 2. constraint is an inequality */
    int sign;      /* sign of the scalar product */

    for (i=0; i<F->nbrows; i++){
	
      //opt_numint_t * fpi = opt_melina_vector(F->p[i],mapF,comp_size,F->nbcolumns);
      opt_numint_t prod;
      prod = opt_vector_product_strict(opk, F->p[i], tab, F->nbcolumns);
      if(opk->exn){
	opk->exn = ELINA_EXC_NONE;
	return false;
      }	
      sign = opt_numint_sgn(prod);

      if (sign<0){
        //free(fpi);
	return false;
      }
      else {
	if (opt_numint_sgn(F->p[i][0])==0){
	  /* line */
	  if (sign!=0) {
		//free(fpi);
		return false;
	  }
	}
	else {
	  /* ray or vertex */
	  if (is_strict && sign==0 &&
	      (opk->strict ? opt_numint_sgn(F->p[i][opt_polka_eps])>0 : true)){
		    //free(fpi);
		    return false;
	    }
	}
      }
      //free(fpi);
    }
    
    return true;
  }
}

bool opt_poly_leq(opt_pk_internal_t * opk, opt_matrix_t * C, opt_matrix_t * F){
	size_t nbrows = C->nbrows;
        unsigned short int comp_size = C->nbcolumns - 2;
        size_t i;
	for(i=0; i < nbrows; i++){
		if(!opt_generators_sat_vector(opk,F,C->p[i],  opk->strict && opt_numint_sgn(C->p[i][opt_polka_eps])<0)){
			return false;
		}
	}
	return true;
}


bool opt_pk_is_leq_gen(elina_manager_t * man, opt_pk_array_t *oa, opt_pk_array_t *ob){
 
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_IS_LEQ);
  
  man->result.flag_exact = man->result.flag_best = false;
  unsigned short int k, ka,kb;
  array_comp_list_t *acla = oa->acl;
  unsigned short int maxcols = oa->maxcols;
  unsigned short int num_compa = acla->size;
  opt_pk_t ** poly_a = oa->poly;
  for(ka=0; ka < num_compa; ka++){
	  if (opk->funopt->algorithm>=0)
	    opt_poly_chernikova(man,poly_a[ka],"is leq first argument");
	  else
	    opt_poly_obtain_F(man,poly_a[ka],"is leq first argument");

	  if (opk->exn){
	    opk->exn = ELINA_EXC_NONE;
	    return false;
	  }
	  if (!poly_a[ka]->F){ /* pa is empty */
	    man->result.flag_exact = man->result.flag_best = true;
	    return true;
	  }
  }
  
  array_comp_list_t *aclb = ob->acl;
  unsigned short int num_compb = aclb->size;
  opt_pk_t ** poly_b = ob->poly; 
  for(kb=0; kb < num_compb; kb++){
	  if (opk->funopt->algorithm>=0)
	    opt_poly_chernikova(man,poly_b[kb],"is leq second argument");
	  else
	    opt_poly_obtain_C(man,poly_b[kb],"is leq second argument");

	  if (opk->exn){
	    opk->exn = ELINA_EXC_NONE;
	    return false;
	  }
	  if (!poly_b[kb]->C){/* pb is empty */
	    return false;
	  }
  }
  
  bool sat;
  
  /*******************************
	Compute union of independent components
  *******************************/
  
  array_comp_list_t * acl = union_array_comp_list(acla,aclb,maxcols);
  
  unsigned short int num_comp = acl->size;
  
  unsigned short int ** ca_arr = (unsigned short int **)malloc(num_comp*sizeof(unsigned short int *));
  unsigned short int * comp_size_arr = (unsigned short int *)malloc(num_comp*sizeof(unsigned short int));
  opt_pk_t ** tmp = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
  size_t * nbgena = (size_t *)calloc(num_comp,sizeof(size_t));
  comp_list_t *cl = acl->head;
  size_t * num_vertex1 = (size_t *)calloc(num_comp,sizeof(size_t));  
  char ** var_map = (char **)malloc(num_comp*sizeof(char *)); 

  for(k=0; k < num_comp; k++){
	unsigned short int comp_size = cl->size;
	ca_arr[k] = to_sorted_array(cl,maxcols);
	comp_size_arr[k] = comp_size;
	num_vertex1[k] = 1;
	var_map[k] = (char *)calloc(comp_size,sizeof(char));
	cl = cl->next;
  }
 
  /********************************
	Data structures for A
  *********************************/
  unsigned short int * rmapa = (unsigned short int *)calloc(num_compa, sizeof(unsigned short int));	
  unsigned short int ** ind_melina_a = (unsigned short int **)malloc(num_compa*sizeof(unsigned short int *));
  comp_list_t * cla = acla->head; 
  size_t * num_vertex_a = (size_t *)calloc(num_compa,sizeof(size_t));
 
  for(ka=0; ka < num_compa; ka++){
	short int inda = is_comp_list_included(acl,cla,maxcols);
	rmapa[ka] = inda;
	unsigned short int * ca_a = to_sorted_array(cla,maxcols);
	unsigned short int *ca = ca_arr[inda];
	ind_melina_a[ka] = map_index(ca_a,ca, cla->size);
	opt_pk_t *oak = poly_a[ka];
	nbgena[inda] = nbgena[inda] + oak->F->nbrows;
	opt_poly_obtain_satF(oak);
	num_vertex_a[ka] = opt_generator_rearrange(oak->F,oak->satF);
	if(num_vertex_a[ka]){
		num_vertex1[inda] = num_vertex1[inda] * num_vertex_a[ka];
	}
	unsigned short int k1;
	for(k1=0; k1 < cla->size; k1++){
		unsigned short int num = ind_melina_a[ka][k1];
		var_map[inda][num]++;
	}
	free(ca_a);
	cla = cla->next;
  }

   /********************************
	Data structures for B
  *********************************/
 
  unsigned short int * rmapb = (unsigned short int *)calloc(num_compb, sizeof(unsigned short int));	
  unsigned short int ** ind_melina_b = (unsigned short int **)malloc(num_compb*sizeof(unsigned short int *));
  comp_list_t * clb = aclb->head;
  for(kb=0; kb < num_compb; kb++){
	short int indb = is_comp_list_included(acl,clb,maxcols);
	rmapb[kb] = indb;
	unsigned short int * ca_b = to_sorted_array(clb,maxcols);
	unsigned short int *ca = ca_arr[indb];
	ind_melina_b[kb] = map_index(ca_b,ca, clb->size);
	free(ca_b);
        opt_pk_t *obk = poly_b[kb];
        // if (opk->funopt->algorithm>0){
    	//	opt_poly_chernikova(man,obk,"of the second argument");
	 //}
  	//else{
    	//	opt_poly_obtain_C(man,obk,"of the second argument");
        //}
	
        //if (opk->exn){
    	//	opk->exn = ELINA_EXC_NONE;
    	//	sat = false;
	//	goto opt_pk_is_leq_gen_exit;
  	//}
  	//if (!obk->C){/* ob is empty */
    	//	man->result.flag_exact = man->result.flag_best = (obk->intdim==0);
    	//	sat = false;
	//	goto opt_pk_is_leq_gen_exit;
  	//}
	clb = clb->next;
  }   
	 
  	cl = acl->head;
	 
	for(k=0; k < num_comp; k++){
		unsigned short int comp_size = cl->size; 
		tmp[k] = opt_poly_alloc(comp_size,0);
		unsigned short int k1, nblines = 0;
		for(k1=0; k1< comp_size; k1++){
			if(!var_map[k][k1]){
				nblines++;
			}
		}
		tmp[k]->F = opt_matrix_alloc(nbgena[k] + 2*num_vertex1[k]+nblines,comp_size+2, false);
		num_vertex1[k] = 0;
		
		cl = cl->next;
	}

        cla = acla->head;
        size_t *counter = (size_t *)calloc(num_comp,sizeof(size_t));
	/******************
		Cartesian product of vertices of A
	*******************/
	cartesian_product_vertices(oa,tmp,rmapa,ca_arr,num_vertex_a, num_vertex1,counter);

	/******************
		Meet rays of A
	*******************/
	meet_rays(oa,tmp,rmapa,ca_arr,num_vertex_a,counter);
	cl = acl->head;
	for(k=0; k < num_comp; k++){
		size_t count = counter[k];
		unsigned short int k1;
		unsigned short int comp_size = cl->size;
		opt_matrix_t * F = tmp[k]->F;
		for(k1=0; k1< comp_size; k1++){
			if(!var_map[k][k1]){
				F->p[count][k1+2] = 1;
				count++;
			}
		}
		F->nbrows = count;
		cl = cl->next;
	}
	
  //man->result.flag_exact = man->result.flag_best =
    //oa->intdim==0;
  /* if both are mininmal, check the dimensions */
  /*if (pa->C && pa->F && pb->C && pb->F
      && (pa->nbeq < pb->nbeq || pa->nbline > pb->nbline))
    {
     man->result.flag_exact = man->result.flag_best = true;
     return false;
    }
  else {*/
   
    /* does the frames of pa satisfy constraints of pb ? */
    size_t i;
    sat = true;
    
    //cla = acla->head;
    for(kb=0; kb < num_compb; kb++){
	//opt_pk_t * oak = poly_a[ka];
	//unsigned short int inda = rmapa[ka];
	//man->result.flag_exact = man->result.flag_best = (oak->intdim==0);
        //clb = aclb->head;
        unsigned short int indb = rmapb[kb];
	//for(k=0; k < num_comp; k++){	
		//if(k==indb){
			opt_pk_t * obk = poly_b[kb];
			size_t nbconsb = obk->C->nbrows;
			for(i=0; i<nbconsb;i++){
				
				opt_numint_t * cpi = opt_map_vector(obk->C->p[i],ind_melina_b[kb], comp_size_arr[indb],obk->C->nbcolumns);
				
				sat = opt_generators_sat_vector(opk,
					  tmp[indb]->F, cpi,
					  opk->strict &&
					  opt_numint_sgn(obk->C->p[i][opt_polka_eps])<0);

                                // free(cpi);
                                if(!sat){ 
	 				goto opt_pk_is_leq_gen_exit;
      				}
			}
			
		//}
	   //clb = clb->next;
	//}
        //cla = cla->next;
    }
  
  opt_pk_is_leq_gen_exit:
	
	free(rmapa);
	free(rmapb);
	free(comp_size_arr);
	free(num_vertex_a);
	free(num_vertex1);
	for(k=0; k < num_comp; k++){
		free(ca_arr[k]);
                free(tmp[k]->F);
                free(tmp[k]);
		free(var_map[k]);
	}
	free(tmp);
	free(ca_arr);
	free(counter);
	free(var_map);
	for(k=0; k < num_compa; k++){
		free(ind_melina_a[k]);
	}
	free(ind_melina_a);
	free(nbgena);
	for(k=0; k < num_compb; k++){
		free(ind_melina_b[k]);
	}
	free(ind_melina_b);

        return sat;
}


/*****************************************************
		Inequality test
******************************************************/
bool opt_pk_is_leq(elina_manager_t *man, opt_pk_array_t *oa, opt_pk_array_t *ob){
	#if defined(TIMING)
		start_timing();
     	 #endif
                opt_pk_internal_t *opk =
                    opt_pk_init_from_manager(man, ELINA_FUNID_IS_LEQ);
                array_comp_list_t *acla = oa->acl;
                if (oa->is_bottom || !acla) {
		#if defined(TIMING)
			record_timing(is_lequal_time);
		#endif
		return true;
	}
	
	array_comp_list_t * aclb = ob->acl;
	if(ob->is_bottom || !aclb){
		#if defined(TIMING)
			record_timing(is_lequal_time);
		#endif
		return false;
	}
	bool res = opt_pk_is_leq_gen(man,oa,ob);
	#if defined(TIMING)
		record_timing(is_lequal_time);
	#endif
	return res;
}



/*****************************************************
		Equality test
******************************************************/
bool opt_pk_is_eq(elina_manager_t* man, opt_pk_array_t* oa, opt_pk_array_t* ob)
{
  opt_pk_internal_t *opk = opt_pk_init_from_manager(man, ELINA_FUNID_IS_EQ);
  array_comp_list_t * acla = oa->acl;
  array_comp_list_t * aclb = ob->acl;
  if(oa->is_bottom || !acla){
	if(ob->is_bottom || !aclb){
		return true;
	}
	else{
		return false;
	}
  }
  else{
	if(ob->is_bottom || !aclb){
		return false;
	}
	else{
		/**************************
			Both A and B are not bottom, compare the number of equalities
		**************************/
		unsigned short int ka, kb;
		opt_pk_t ** poly_a = oa->poly;
		size_t nbeqa = 0;
		unsigned short int num_compa = acla->size;
		bool flaga = true;
		for(ka=0; ka < num_compa; ka++){
			opt_pk_t * oak = poly_a[ka];
			if(oak->is_minimized){
				nbeqa = nbeqa + oak->nbeq;
			}
			else{
				flaga = false;
			}
		}
		
		opt_pk_t ** poly_b = ob->poly;
		size_t nbeqb = 0;
		bool flagb = true;
		unsigned short int num_compb = aclb->size;
		for(kb = 0; kb < num_compb; kb++){
			opt_pk_t * obk = poly_b[kb];
			if(obk->is_minimized){
				nbeqb = nbeqb + obk->nbeq;
			}
			else{
				flagb = false;
			}
		} 
		if(flaga && flagb && (nbeqa != nbeqb)){
			return false;
		}
		 man->result.flag_exact = man->result.flag_best = true;
    		bool res1 = opt_pk_is_leq(man,oa,ob);
		if(!res1){
			return false;
		}
    		bool res2 = opt_pk_is_leq(man,ob,oa);
    		bool res = res1 && res2;
		return res;
	}
  }
  
  
 
   
}


/********************************************************
		Is dimension unconstrained
********************************************************/
bool opt_pk_is_dimension_unconstrained(elina_manager_t* man, opt_pk_array_t* oa,
					elina_dim_t dim)
{
  size_t i,j;
  bool res;
  opt_matrix_t* F;
  opt_matrix_t* C;
  elina_dim_t var = dim+2;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED);
  if(oa->is_bottom || !oa->acl){
	man->result.flag_exact = man->result.flag_best = true;
    	return false;
  }

  comp_list_t * cl = find(oa->acl,var);
  if(cl==NULL){
	man->result.flag_exact = man->result.flag_best = true;
    	return true;
  }
  short int k = find_index(oa->acl,cl);
  opt_pk_t ** poly_a = oa->poly;
  opt_pk_t * oak = poly_a[k];
  opt_poly_chernikova(man,oak,NULL);
  if (opk->exn){
    opk->exn = ELINA_EXC_NONE;
    return false;
  }
  if (!oak->F){ /* po is empty */
    man->result.flag_exact = man->result.flag_best = true;
    return false;
  }
  unsigned short int * ca = to_sorted_array(cl,oa->maxcols);
  unsigned short int l, new_dim=0;
  unsigned short int comp_size = cl->size;
  for(l=0; l < comp_size; l++){
	if(ca[l]==var){
		new_dim = l;
	}
  }
  free(ca);
  /* We test if there exists the line of direction dim */
  F = oak->F;
  res = false;
  for (i=0; i<oak->nbline; i++){
    if (opt_numint_sgn(F->p[i][opk->dec+new_dim])){
      res = true;
      for(j=opk->dec; j<F->nbcolumns; j++){
	if (j!=opk->dec+new_dim && opt_numint_sgn(F->p[i][j])){
	  res = false;
	  break;
	}
      }
      break;
    }
  }
  man->result.flag_exact = man->result.flag_best = true;
  return res;
}


/* ====================================================================== */
/* Satisfiability of a linear constraint */
/* ====================================================================== */
void fuse_generators_intersecting_blocks(opt_matrix_t *F, opt_pk_t ** poly_a, array_comp_list_t *acla, unsigned short int *ca, 
					 size_t * num_vertex_a, char *intersect_map, unsigned short int maxcols){
	unsigned short int num_compa = acla->size;
	unsigned short int j,k;
	size_t i;
	size_t num_vertex = 0;
	size_t counter = 0;
	comp_list_t * cla = acla->head;
	unsigned short int ** ca_a_arr = (unsigned short int **)malloc(num_compa*sizeof(unsigned short int *));
	for(k=0; k < num_compa; k++){
		if(!intersect_map[k]){
			cla = cla->next;
			continue;
		}
		opt_pk_t * src = poly_a[k];
		ca_a_arr[k] = to_sorted_array(cla,maxcols);
		unsigned short int * ca_a = ca_a_arr[k];
		unsigned short int comp_size = cla->size;
		opt_matrix_t * src_mat = src->F;
                size_t nbconsa = src_mat->nbrows;
                opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = F->p;
		if(!num_vertex){
			for(i=0; i < num_vertex_a[k]; i++){
				opt_numint_t * src_pi = src_p[i];
				opt_numint_t * dst_pi = dst_p[counter];
				dst_pi[0] = src_pi[0];
				dst_pi[1] = src_pi[1];
				unsigned short int l = 0;
				for(j = 0; j < comp_size; j++){
					while(ca[l] != ca_a[j]){
						l++;
					}
					dst_pi[l+2] = src_pi[j+2];
					l++;
				}
				counter++;
				num_vertex++;
			}
			
			F->nbrows = counter;
		}
		else{
			size_t start = counter - num_vertex;
			size_t i2, i3 = start;
			size_t end = counter;			
			for(i=0; i < num_vertex_a[k]; i++){
				opt_numint_t * src_pi = src_p[i];
				for(i2=start; i2 < end; i2++){
					opt_numint_t * vec = dst_p[i2];
					opt_numint_t * dst_pi = dst_p[i3];
					opt_numint_t lcm = opt_numint_lcm(src_pi[1],vec[1]);
					opt_numint_t lcm1 = lcm/vec[1];
					opt_numint_t lcm2 = lcm/src_pi[1];
					if(i3 >= counter){
						num_vertex++;
						opt_vector_copy(dst_pi,vec,F->nbcolumns);
					}
					dst_pi[1] = lcm;
					unsigned short int l = 0;
					
					
					if(lcm1!=1){
						
						for(j=2; j < F->nbcolumns; j++){
							dst_pi[j] = lcm1*dst_pi[j];
						}
					}
					
					
					if(lcm2==1){
						for(j = 0; j < comp_size; j++){
							while(ca[l] != ca_a[j]){
								l++;
							}
							dst_pi[l+2] = src_pi[j+2];
							l++;
						}
					}
					else{
						for(j = 0; j < comp_size; j++){
							while(ca[l] != ca_a[j]){
								l++;
							}
							dst_pi[l+2] = lcm2*src_pi[j+2];
							l++;
						}
					}
					i3++;
							
				}
			}
			counter = i3;
		}
		//free(ca_a);
		cla = cla->next;
	}
	
	cla = acla->head;
        for(k=0; k < num_compa;k++){
		if(!intersect_map[k]){
			cla = cla->next;
			continue;
		}
		opt_pk_t * src = poly_a[k];
		opt_matrix_t * src_mat = src->F;
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = F->p;
		size_t nbconsa = src_mat->nbrows;
		unsigned short int comp_size = cla->size;
		unsigned short int * ca_a = ca_a_arr[k];
		for(i=num_vertex_a[k]; i < nbconsa; i++){
			opt_numint_t * src_pi = src_p[i];
			opt_numint_t * dst_pi = dst_p[counter];
			dst_pi[0] = src_pi[0];
			dst_pi[1] = src_pi[1];
			unsigned short int l=0;
			for(j=0; j < comp_size; j++){
				while(ca[l]!=ca_a[j]){
					l++;
				}
				dst_pi[l+2] = src_pi[j+2];
				l++;
			}
			counter++;
		}
		cla = cla->next;
		free(ca_a_arr[k]);
	}
	F->nbrows = counter;
	free(ca_a_arr);
}


bool opt_pk_sat_lincons(elina_manager_t* man, opt_pk_array_t* oa, elina_lincons0_t* lincons0)
{
  bool exact=true;
  bool sat;
  size_t dim;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_SAT_LINCONS);

  if(oa->is_bottom || !oa->acl){
	man->result.flag_exact = man->result.flag_best = true;
	return true;
  }
 
  array_comp_list_t * acla = oa->acl;
  opt_pk_t ** poly_a = oa->poly;
  unsigned short int num_compa = acla->size;
  unsigned short int maxcols = oa->maxcols;
  unsigned short int k;
  for(k=0; k < num_compa; k++){
      opt_pk_t * oak = poly_a[k];
      if (opk->funopt->algorithm>0){
	  opt_poly_chernikova(man,oak,NULL);
      }
      else{
	  opt_poly_obtain_F(man,oak,NULL);
      }
      if (opk->exn){
    	  opk->exn = ELINA_EXC_NONE;
    	  return false;
      }
      if (!oak->C && !oak->F){ /* one of the factor is empty */
           man->result.flag_exact = man->result.flag_best = true;
           return true;
      }
  }
  switch (lincons0->constyp){
          case ELINA_CONS_EQ:
          case ELINA_CONS_SUPEQ:
          case ELINA_CONS_SUP:
               break;
          default:
               man->result.flag_exact = man->result.flag_best = false;
               return false;
  }  

  elina_linexpr0_t * expr = lincons0->linexpr0;
  if(is_linexpr_zero(expr)){
	int sgn = elina_coeff_sgn(&expr->cst);
	if(lincons0->constyp==ELINA_CONS_EQ){
		if(sgn){
			return false;
		}
		else{
			return true;
		}
	}
	else if(lincons0->constyp==ELINA_CONS_SUPEQ){
		if(sgn< 0){
			return false;
		}
		else{
			return true;
		}
	}
	else if(lincons0->constyp==ELINA_CONS_SUP){
		if(sgn != 1){
			return false;
		}
		else{
			return true;
		}
	}
  }

  comp_list_t * clb = lincons0_to_comp_list(opk,lincons0);
  array_comp_list_t *aclb = create_array_comp_list();
  insert_comp_list(aclb,clb);
  array_comp_list_t * acl = union_array_comp_list(acla,aclb,maxcols);
  //short int res = is_comp_list_included(acl,clb,maxcols);
  char * map = create_map(clb,maxcols);
  comp_list_t * cla = acla->head;
  char * intersect_map = (char *)calloc(num_compa,sizeof(char));
  size_t * num_vertex_a = (size_t *)calloc(num_compa,sizeof(size_t));
  size_t num_vertex=0;
  size_t nbgen=0;
  for(k=0; k < num_compa; k++){
      if(is_disjoint_with_map(cla,map)){
          cla = cla->next;
	  continue;
      }
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
  
  comp_list_t * cl = acl->head;
  unsigned short int num_comp = acl->size;
  unsigned short int comp_size = 0;
  unsigned short int * ca = NULL;
  elina_lincons0_t * new_lincons0 = (elina_lincons0_t *)malloc(sizeof(elina_lincons0_t));
  new_lincons0->scalar = NULL;
  for(k=0; k < num_comp; k++){
        if(!is_disjoint(cl,clb,maxcols)){
		comp_size = cl->size;
		ca = to_sorted_array(cl,maxcols);
		copy_lincons0_with_comp_list(opk,new_lincons0,lincons0,ca, comp_size);
		break;
	}       
        cl = cl->next;
  }
  opt_matrix_t *F = opt_matrix_alloc(nbgen + num_vertex, comp_size + 2, false);
  fuse_generators_intersecting_blocks(F,poly_a,acla,ca,num_vertex_a,intersect_map,maxcols);
  free(ca);
  free(map);
  free(intersect_map);
  //dim = po->intdim + po->realdim;
  dim = comp_size;
  if (!elina_linexpr0_is_quasilinear(lincons0->linexpr0)){
    elina_interval_t** env = opt_matrix_to_box(opk,F);
    exact = quasilinearize_elina_lincons0(new_lincons0, env, false,ELINA_SCALAR_MPQ) && exact;
    elina_interval_array_free(env,dim);
  }

  sat = opt_vector_set_elina_lincons0_sat(opk,
				   opk->poly_numintp,
				   new_lincons0,
				   comp_size, 0, true);
  if (sat){
    if(F->nbrows ==0){
	sat = false;
    }
    else{
	    sat = opt_generators_sat_vector(opk,F,
					   opk->poly_numintp,
					   new_lincons0->constyp==ELINA_CONS_SUP);
    }
  }
  man->result.flag_exact = man->result.flag_best =
    sat ?
    true :
    (
     ( (opk->funopt->flag_exact_wanted || opk->funopt->flag_best_wanted) &&
       exact && elina_linexpr0_is_real(lincons0->linexpr0,dim) ) ?
     true :
     false );
  elina_lincons0_clear(new_lincons0);
  free(new_lincons0);
  free_array_comp_list(aclb);
  free_array_comp_list(acl);
  free(num_vertex_a);
  opt_matrix_free(F);
  return sat;
}


bool opt_pk_sat_tcons(elina_manager_t* man, opt_pk_array_t* oa, elina_tcons0_t* cons)
{
  size_t dim;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_SAT_LINCONS);
  array_comp_list_t *acla = oa->acl;
  if (oa->is_bottom || !acla){ /* oa is empty */
    man->result.flag_exact = man->result.flag_best = true;
    return true;
  }
  opt_pk_t ** poly_a = oa->poly;
  unsigned short int num_compa = acla->size;
  unsigned short int maxcols = oa->maxcols; 
  unsigned short int k;
  for(k=0; k < num_compa; k++){
	opt_pk_t * oak = poly_a[k];
	if (opk->funopt->algorithm>0){
	    opt_poly_chernikova(man,oak,"sat tcons input");
	}
	else{
	    opt_poly_obtain_F(man,oak,"sat tcons input");
	}
	if (opk->exn){
	    opk->exn = ELINA_EXC_NONE;
	    return false;
	}
	if (!oak->F){ /* the factor is empty */
	    man->result.flag_exact = man->result.flag_best = true;
	    return true;
	}
  }
  switch (cons->constyp){
  case ELINA_CONS_EQ:
  case ELINA_CONS_SUPEQ:
  case ELINA_CONS_SUP:
    break;
  default:
    man->result.flag_exact = man->result.flag_best = false;
    return false;
  }
  dim = maxcols - opk->dec;

  elina_interval_t** env = opt_pk_to_box(man,oa);
  /*******
	Linearize the tree expression
  ********/
  elina_lincons0_t *lincons0 = (elina_lincons0_t *)malloc(sizeof(elina_lincons0_t));
  elina_linexpr0_t *linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,0);
  elina_scalar_t *scalar = elina_scalar_alloc();
  lincons0->constyp = cons->constyp;
  lincons0->linexpr0 = linexpr0;
  lincons0->scalar = scalar;
  elina_intlinearize_elina_tcons0(lincons0,
			     cons,env,dim,ELINA_SCALAR_MPQ);
  elina_lincons0_fprint(stdout, lincons0, NULL);
  quasilinearize_elina_lincons0(lincons0,env,false,ELINA_SCALAR_MPQ);
  elina_interval_array_free(env,dim);
  elina_linexpr0_t * expr = lincons0->linexpr0;
  if(is_linexpr_zero(expr)){
	int sgn = elina_coeff_sgn(&expr->cst);
	if(lincons0->constyp==ELINA_CONS_EQ){
		elina_lincons0_clear(lincons0);
		free(lincons0);
		if(sgn){
			return false;
		}
		else{
			return true;
		}
	}
	else if(lincons0->constyp==ELINA_CONS_SUPEQ){
		elina_lincons0_clear(lincons0);
		free(lincons0);
		if(sgn< 0){
			return false;
		}
		else{
			return true;
		}
	}
	else if(lincons0->constyp==ELINA_CONS_SUP){
		elina_lincons0_clear(lincons0);
		free(lincons0);
		if(sgn != 1){
			return false;
		}
		else{
			return true;
		}
	}
  }

  /**************************
		compute partition for the linearized constraint
  ***************************/
  comp_list_t * clb = lincons0_to_comp_list(opk,lincons0);
  array_comp_list_t *aclb = create_array_comp_list();
  insert_comp_list(aclb,clb);
  array_comp_list_t * acl = union_array_comp_list(acla,aclb,maxcols);
  char * map = create_map(clb,maxcols);
  comp_list_t * cla = acla->head;
  char * intersect_map = (char *)calloc(num_compa,sizeof(char));
  size_t * num_vertex_a = (size_t *)calloc(num_compa,sizeof(size_t));
  size_t num_vertex=0;
  size_t nbgen=0;
  for(k=0; k < num_compa; k++){
      if(is_disjoint_with_map(cla,map)){
          cla = cla->next;
	  continue;
      }
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
  
  comp_list_t * cl = acl->head;
  unsigned short int num_comp = acl->size;
  unsigned short int comp_size=0;
  unsigned short int * ca = NULL;
  elina_lincons0_t * new_lincons0 = (elina_lincons0_t *)malloc(sizeof(elina_lincons0_t));
  new_lincons0->scalar = NULL;
  for(k=0; k < num_comp; k++){
        if(!is_disjoint(cl,clb,maxcols)){
		comp_size = cl->size;
		ca = to_sorted_array(cl,maxcols);
		copy_lincons0_with_comp_list(opk,new_lincons0,lincons0,ca, comp_size);
		//free(ca);
		break;
	}       
        cl = cl->next;
  }
  opt_matrix_t *F = opt_matrix_alloc(nbgen + num_vertex, comp_size + 2, false);

  fuse_generators_intersecting_blocks(F,poly_a,acla,ca,num_vertex_a,intersect_map,maxcols);
  free(map);
  free(ca);
  free(intersect_map);
  bool sat = opt_vector_set_elina_lincons0_sat(opk,
					opk->poly_numintp,
					new_lincons0,
					comp_size, 0, true);
  
  if (sat){
    if(F->nbrows==0){
	sat = false;
    }
    else{
	    sat = opt_generators_sat_vector(opk,F,
					   opk->poly_numintp,
					   cons->constyp==ELINA_CONS_SUP);
   }
  }
  man->result.flag_exact = man->result.flag_best = sat;
  elina_lincons0_clear(lincons0);
  elina_lincons0_clear(new_lincons0);
  free(lincons0);
  free(new_lincons0);
  free_array_comp_list(aclb);
  free_array_comp_list(acl);
  free(num_vertex_a);
  opt_matrix_free(F);
  return sat;
}

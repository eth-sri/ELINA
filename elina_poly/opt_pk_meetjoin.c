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
/* opt_pk_meetjoin.c: Meet and join operations */
/* ********************************************************************** */


#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"
#include "opt_pk_representation.h"
#include "opt_pk_user.h"
#include "opt_pk_constructor.h"
#include "opt_pk_test.h"
#include "opt_pk_meetjoin.h"
#include "opt_pk_project.h"
#include "opt_pk_cherni.h"
#include "rdtsc.h"

/* ====================================================================== */
/* II.2 Meet with (array of) linear constraint(s) */
/* ====================================================================== */

/* ---------------------------------------------------------------------- */
/* Factorized version */



bool opt_poly_meet_matrix(bool meet,
		      bool lazy,
		      elina_manager_t* man,
		      opt_pk_t* op,
		      opt_pk_t* oa, opt_matrix_t* oc)
{
  
  opt_pk_internal_t* opk = (opt_pk_internal_t*)man->internal;
  man->result.flag_best = (oa->intdim==0);
  man->result.flag_exact = meet;
  size_t start = oa->C->nbrows;
  // assert(pa->satC);
  if (lazy) {

    opt_poly_obtain_sorted_C(opk, oa);
    if (op != oa) {
      op->C = opt_matrix_merge_sort(opk, oa->C, oc);
      // op->nbeq = oa->nbeq;
    } else {
      opt_matrix_merge_sort_with(opk, oa->C, oc);
      if (oa->F) {
        opt_matrix_free(oa->F);
        oa->F = NULL;
      }
      if (oa->satC) {
        opt_satmat_free(oa->satC);
        oa->satC = NULL;
      }
      if (oa->satF) {
        opt_satmat_free(oa->satF);
        oa->satF = NULL;
      }
      oa->nbline = 0;
      op->nbeq = oa->nbeq;
    }
    oa->status = 0;
    return false;
     }
     else{
	size_t start = oa->C->nbrows;
    	assert(oa->satC);
		
	if (op != oa){
	      	op->C = opt_matrix_append(oa->C,oc);
	      	op->satC = opt_satmat_copy_resize_cols(oa->satC,
					opt_bitindex_size(op->C->nbrows));
	      	op->nbline = oa->nbline;
	      	op->nbeq = oa->nbeq;
	}
    	else {
      		opt_matrix_append_with(oa->C,oc);
      		opt_satmat_resize_cols(oa->satC,
			 opt_bitindex_size(oa->C->nbrows));
    	}
	opt_matrix_sort_rows(opk,op->F);
    	combine_satmat(opk,op,op->C->nbcolumns -2,start,true);
    	opt_cherni_add_and_minimize(opk,meet,op,start);
	if(opk->exn){
          // opk->exn = ELINA_EXC_NONE;

          return false;
	}
	if(!op->F){
		return true;
	}
	return false;
     }
    
  op->is_minimized = false;
  return false;
 
}




bool opt_poly_meet_elina_lincons_array(bool lazy,
				 elina_manager_t* man,
				 opt_pk_t* op, opt_pk_t* oa,
				 elina_lincons0_array_t* array)
{
  opt_matrix_t* mat;
  bool quasilinear;
  opt_pk_internal_t* opk = (opt_pk_internal_t*)man->internal;
  unsigned short int num_var = oa->intdim + oa->realdim;
  unsigned short int k; 
  quasilinear = elina_lincons0_array_is_quasilinear(array);
  /* quasilinearize if needed */
  if (!quasilinear){
    	elina_interval_t ** env = opt_generator_to_box(opk,oa->F);
    	quasilinearize_elina_lincons0_array(array,env,true,ELINA_SCALAR_MPQ);
	
	for(k=0; k < num_var; k++){
		elina_interval_free(env[k]);
	}
	free(env);
  }
  linearize_elina_lincons0_array(array,true,ELINA_SCALAR_MPQ);
  elina_lincons0_array_reduce_integer(array,op->intdim,ELINA_SCALAR_MPQ);
  bool exact = opt_matrix_set_elina_lincons0_array(opk,&mat,array,op->intdim,op->realdim,true);
  opt_matrix_sort_rows(opk,mat);
  size_t i;
  for(i=0; i < mat->nbrows; i++){
	if(!mat->p[i][0]){
		op->nbeq++;
	}
  }
  
  bool is_bottom = opt_poly_meet_matrix(true,lazy,man,op,oa,mat);
  
  opt_matrix_free(mat);
  if (opk->exn){
    // opk->exn = ELINA_EXC_NONE;
    man->result.flag_exact = man->result.flag_best = false;
  }
  else {
    	man->result.flag_best = man->result.flag_exact = exact ? true : false;
  }
  return is_bottom;
}

comp_list_t * linexpr0_to_comp_list(opt_pk_internal_t * opk, elina_linexpr0_t * expr){
	comp_list_t * cl = create_comp_list();
	size_t size = expr->size;
	size_t j;
	elina_linexpr_discr_t discr = expr->discr;
	if(discr==ELINA_LINEXPR_DENSE){
		elina_coeff_t* coeff = expr->p.coeff;
		for(j=0; j < size; j++){
			if(!elina_coeff_zero(&coeff[j])){
				insert_comp(cl,j + opk->dec);
			}
		}
	}
	else{
		elina_linterm_t* linterm = expr->p.linterm;
		for(j = 0; j < size; j++){
			elina_dim_t dim = linterm[j].dim;
			insert_comp(cl,dim + opk->dec);	
		}
	}
	return cl;
}

comp_list_t * lincons0_to_comp_list(opt_pk_internal_t * opk, elina_lincons0_t * cons){
	elina_linexpr0_t * expr = cons->linexpr0;
	return linexpr0_to_comp_list(opk,expr);
}

array_comp_list_t * lincons0_array_to_array_comp_list(opt_pk_internal_t *opk, elina_lincons0_array_t* array, unsigned short int n, char * is_trivial){
	size_t size = array->size;
	size_t i;
	array_comp_list_t * acl = create_array_comp_list();
	for(i = 0; i < size; i++){
		if(is_trivial[i]){
			continue;
		}
		elina_lincons0_t * cons = array->p + i;
		comp_list_t * cl = lincons0_to_comp_list(opk,cons);
		if(cl->size!=0){
			insert_comp_list_with_union(acl,cl,n);
		}			
	}
	return acl;
}

elina_linexpr0_t * copy_linexpr0_with_comp_list(opt_pk_internal_t *opk, elina_linexpr0_t * src, unsigned short int * ca, unsigned short int comp_size){
	elina_linexpr0_t * dst;
	elina_linexpr_discr_t discr = src->discr; 
	size_t size = src->size;
	unsigned short int j, k;
	if(discr==ELINA_LINEXPR_DENSE){
		dst = elina_linexpr0_alloc(discr,comp_size);
		elina_coeff_set(&dst->cst,&src->cst);
		elina_coeff_t* src_coeff = src->p.coeff;
		elina_coeff_t* dst_coeff = dst->p.coeff;
		for(k=0; k < comp_size; k++){
			unsigned short int num = ca[k] - opk->dec;
			if(!elina_coeff_zero(&src_coeff[num])){
				elina_coeff_set(&dst_coeff[k],&src_coeff[num]);
			}
		}
		
	}
	else{
		dst = elina_linexpr0_alloc(discr,size);
		elina_coeff_set(&dst->cst,&src->cst);
		elina_linterm_t* src_linterm = src->p.linterm;
		elina_linterm_t * dst_linterm = dst->p.linterm;
		
		for(j = 0; j < size; j++){
			elina_dim_t src_dim = src_linterm[j].dim + opk->dec;
			elina_coeff_t src_coeff = src_linterm[j].coeff;
			/************************
				Linexpr dimensions may not be sorted as observed in SeaHorn
			***********************/
			k = 0;
			while(ca[k] != src_dim){
				k++;
			}
			dst_linterm[j].dim = k;
			elina_coeff_set(&dst_linterm[j].coeff, &src_coeff);
                        k++;
                }
	}
	return dst;
}

void copy_lincons0_with_comp_list(opt_pk_internal_t *opk, elina_lincons0_t * dst, elina_lincons0_t * src, unsigned short int * ca, unsigned short int comp_size){
	dst->constyp = src->constyp;
        if(src->scalar!=NULL){
		if(dst->scalar==NULL){
			dst->scalar = elina_scalar_alloc();
		}
		elina_scalar_set(dst->scalar,src->scalar);
        }
	elina_linexpr0_t * src_expr = src->linexpr0;
	dst->linexpr0 = copy_linexpr0_with_comp_list(opk,src_expr,ca,comp_size);
}

bool is_linexpr_zero(elina_linexpr0_t * expr){
	elina_linexpr_discr_t discr = expr->discr;
	unsigned short int size = expr->size;
	unsigned short int j;	
	if(discr==ELINA_LINEXPR_DENSE){
		elina_coeff_t * coeff = expr->p.coeff;
		for(j=0; j < size; j++){
			if(!elina_coeff_zero(&coeff[j])){
				return false;
			}			
		}
	}
	else{
		elina_linterm_t * linterm = expr->p.linterm;
		for(j=0; j < size; j++){
			elina_coeff_t coeff = linterm[j].coeff;
			if(!elina_coeff_zero(&coeff)){
				return false;
			}
		}
	}
	return true;
}

int elina_coeff_sgn(elina_coeff_t * coeff){
    if(coeff->discr==ELINA_COEFF_SCALAR){
	return elina_scalar_sgn(coeff->val.scalar);
    }
    else{
	int si = elina_scalar_sgn(coeff->val.interval->inf);
	int ss = elina_scalar_sgn(coeff->val.interval->sup);
	if(!si){
		if(!ss){
			return 0;
		}
		else{
			return 2;
		}
	}
	else if(si > 0){
		return 1;
	}
	else{
		if(!ss){
			return -2;
		}
		else if(ss > 0){
			return -3;
		}
		else{
			return -1;
		}
	}
    }
	
}

opt_pk_array_t* opt_pk_meet_lincons_array_cons(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, elina_lincons0_array_t* array)
{
  printf(".");
  //printf("meet lincons input\n");
  //elina_lincons0_array_t arr2 = opt_pk_to_lincons_array(man,oa);
  //elina_lincons0_array_fprint(stdout,&arr2,NULL);
  //elina_lincons0_array_clear(&arr2);
  //elina_lincons0_array_fprint(stdout,array,NULL);
  fflush(stdout);
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_MEET_LINCONS_ARRAY);
  size_t i;
  size_t size = array->size; 
  opt_pk_array_t * op;
  if(size==0){
	
	if(destructive){
		op = oa;
		return op;
	}
	else{
		op = opt_pk_copy(man,oa);
		return op;
	}
  }

  op = destructive ? oa :  opt_pk_array_alloc(NULL,NULL,oa->maxcols);
  array_comp_list_t * acla = oa->acl;
  
  /* if oa is bottom, return bottom */
  if(oa->is_bottom || !acla){
	 man->result.flag_best = man->result.flag_exact = true;
    	opt_poly_set_bottom(opk,op);
    	return op;
  }
  unsigned short int k;
  unsigned short int num_compa = acla->size;
  opt_pk_t ** poly_a = oa->poly;
  for(k=0; k < num_compa; k++){
      opt_pk_t * oak = poly_a[k];
      if(opk->funopt->algorithm>=0){
	 opt_poly_chernikova(man,oak,"meet lincons input");
      }
      else{
	 opt_poly_obtain_C(man,oak,"meet lincons input");
      }
      if(opk->exn){
	 opk->exn = ELINA_EXC_NONE;
	 if(!oak->C){
	    man->result.flag_best = man->result.flag_exact = false;
	    opt_poly_set_top(opk,op);
	    return op;
	 }
      }
      if(!oak->C && !oak->F){
	 man->result.flag_best = man->result.flag_exact = true;
   	 opt_poly_set_bottom(opk,op);
	 return op;
      }
  }
  bool is_bottom = false;
  char * is_trivial = (char *)calloc(size, sizeof(char));
  /*****************************************
	Handle Trivial constraints
  *****************************************/
  for(i=0; i < size; i++){
	elina_lincons0_t * cons = array->p + i;
	elina_linexpr0_t * expr = cons->linexpr0;
	elina_constyp_t constyp = cons->constyp;
	if(constyp==ELINA_CONS_DISEQ){
		is_trivial[i] = 1;
	}
	if(is_linexpr_zero(expr)){
		is_trivial[i] = 1;
                int sgn = elina_coeff_sgn(&expr->cst);
                if(constyp==ELINA_CONS_EQ){
                  if (sgn) {
                    is_bottom = true;
                    break;
                  }
                }
		else if(constyp==ELINA_CONS_SUPEQ){
                  if (sgn < 0) {
                    is_bottom = true;
                    break;
                  }
                }
		else if(constyp==ELINA_CONS_SUP){
                  if (sgn != 1) {
                    is_bottom = true;
                    break;
                  }
                }
	}
  }

  if(is_bottom){
	man->result.flag_best = man->result.flag_exact = true;
    	opt_poly_set_bottom(opk,op);
    	return op;
  }


  array_comp_list_t * aclb = lincons0_array_to_array_comp_list(opk,array,oa->maxcols, is_trivial);
  
  /**********************************
	Compute Union of partitions
  **********************************/
  unsigned short int maxcols = oa->maxcols;
  array_comp_list_t * acl = union_array_comp_list(acla,aclb,maxcols);
  unsigned short int num_comp = acl->size;
 
  unsigned short int ** ca_arr = (unsigned short int **)malloc(num_comp*sizeof(unsigned short int *));
  unsigned short int * comp_size_map = (unsigned short int *)calloc(num_comp,sizeof(unsigned short int));
  comp_list_t * cl = acl->head;
  for(k=0; k < num_comp; k++){
	unsigned short int comp_size = cl->size; 
	ca_arr[k] = to_sorted_array(cl,maxcols);
	comp_size_map[k] = comp_size;
	cl = cl->next;
  }

  /**********************************
	Factor linear constraints according to union
  ***********************************/
  elina_lincons0_array_t * arr = (elina_lincons0_array_t *)malloc(num_comp*sizeof(elina_lincons0_array_t ));
  unsigned short int * rmapb = (unsigned short int *)calloc(size, sizeof(unsigned short int));
  size_t * nbmapb = (size_t *)calloc(num_comp,sizeof(size_t));
  for(i=0; i < size; i++){
	if(is_trivial[i]){
		continue;
	}
	elina_lincons0_t * cons = array->p + i;
	comp_list_t * clb = lincons0_to_comp_list(opk,cons);
	short int res = is_comp_list_included(acl,clb,maxcols);
	rmapb[i] = res;
        nbmapb[res]++;
	free_comp_list(clb);
  }
  
  for(k=0; k < num_comp; k++){
	arr[k] = elina_lincons0_array_make(nbmapb[k]);
  }
 
  size_t * counter = (size_t *)calloc(num_comp,sizeof(size_t));
  for(i =0; i < size; i++){
	if(is_trivial[i]){
		continue;
	}
	short int k1 = rmapb[i];
	unsigned short int * ca = ca_arr[k1];
        unsigned short int i1 = counter[k1];
	counter[k1]++;
	elina_lincons0_t * src = array->p + i;
	elina_lincons0_t * dst = (arr+k1)->p + i1;
	copy_lincons0_with_comp_list(opk,dst,src,ca, comp_size_map[k1]);
  }

  /*********************************
	Factor A according to union
  **********************************/
  
  unsigned short int * rmapa = (unsigned short int *)calloc(num_compa, sizeof(unsigned short int));
  size_t * nbmapa = (size_t *)calloc(num_comp,sizeof(size_t));
  size_t * nbeqmapa = (size_t *)calloc(num_comp,sizeof(size_t));
  char * disjoint_map = (char *)calloc(num_compa,sizeof(char));
  size_t * nbgenmapa = (size_t *)calloc(num_comp,sizeof(size_t));
  size_t * nblinemapa = (size_t *)calloc(num_comp,sizeof(size_t));
  size_t * num_vertex_a = (size_t *)calloc(num_compa,sizeof(size_t));
  size_t * num_vertex = (size_t *)calloc(num_comp,sizeof(size_t));
  unsigned short int * array_map_a = create_array_map(acla,maxcols);
  comp_list_t * cla = acla->head;
  //unsigned short int * array_map_b = create_array_map(aclb,maxcols);
  for(k=0; k < num_compa; k++){
	opt_pk_t * oak = poly_a[k];
	short int res = is_comp_list_included(acl,cla,maxcols);
	rmapa[k] = res;
	nbmapa[res] = nbmapa[res] + oak->C->nbrows;
	nbeqmapa[res] = nbeqmapa[res] + oak->nbeq;
	opt_poly_obtain_satF(oak);
	num_vertex_a[k] = opt_generator_rearrange(oak->F,oak->satF);
	if(num_vertex_a[k]){
		if(!num_vertex[res]){
			num_vertex[res] = num_vertex_a[k];
		}
		else{
			num_vertex[res] = num_vertex[res] * num_vertex_a[k];
		}
	}
	nbgenmapa[res] = nbgenmapa[res] + oak->F->nbrows;
	nblinemapa[res] = nblinemapa[res] + oak->nbline;
	cla = cla->next;
  }
  
  opt_pk_t ** poly = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
  size_t * counterC = (size_t *)calloc(num_comp,sizeof(size_t));
  char * pos_con_map = (char *)calloc(num_comp,sizeof(char));
  size_t * counterF = (size_t *)calloc(num_comp,sizeof(size_t));
  cl = acl->head;
  for(k=0; k < num_comp; k++){
	unsigned short int comp_size = cl->size;
	poly[k] = opt_poly_alloc(comp_size,0);
	pos_con_map[k] = 1;
	cl = cl->next;
  }

  cl = acl->head;  	
  for(k=0; k < num_comp; k++){
	unsigned short int l,j;
	if(!nbmapb[k]){
		for(l=0; l < num_compa;l++){
			if(rmapa[l]==k){
				break;
			}
		}
		disjoint_map[l] = 1;
		if(destructive){
			poly[k]->C = poly_a[l]->C;
			poly[k]->nbeq = poly_a[l]->nbeq;
			poly[k]->F = poly_a[l]->F;
			poly[k]->satF = poly_a[l]->satF;
			poly[k]->satC = poly_a[l]->satC;
			poly[k]->nbline = poly_a[l]->nbline;
		}
		else{
			
			poly[k]->C = opt_matrix_copy(poly_a[l]->C);
			poly[k]->nbeq = poly_a[l]->nbeq;
			poly[k]->F = poly_a[l]->F ? opt_matrix_copy(poly_a[l]->F) : NULL; 
			poly[k]->satF = poly_a[l]->satF ? opt_satmat_copy(poly_a[l]->satF) : NULL; 
			poly[k]->satC = poly_a[l]->satC ? opt_satmat_copy(poly_a[l]->satC) : NULL; 
			poly[k]->nbline = poly_a[l]->nbline;
		}
	}
	else{
		unsigned short int comp_size = comp_size_map[k];
		poly[k]->C = opt_matrix_alloc(nbmapa[k]+1, comp_size+2,false);
		poly[k]->nbeq = nbeqmapa[k];		
		unsigned short int k1;
		unsigned short int nblines = 0;
		unsigned short int * ca = ca_arr[k];
		for(k1=0; k1 < comp_size; k1++){
			unsigned short int var = ca[k1];
			if(!array_map_a[var]){
				nblines++;
			}
		}
		poly[k]->F = opt_matrix_alloc(nbgenmapa[k]+2*num_vertex[k]+nblines+1, comp_size+2,false);
		num_vertex[k] = 0;
		nblines = 0;
		for(k1=0; k1 < comp_size; k1++){
			unsigned short int var = ca[k1];
			if(!array_map_a[var]){
				poly[k]->F->p[nblines][k1+2] = 1;
				nblines++;
			}
		}
		poly[k]->nbline = nblinemapa[k] + nblines;
		if(nblines==comp_size){
			poly[k]->F->p[nblines][0] = 1;
			poly[k]->F->p[nblines][1] = 1;
			nblines++;
		}
		counterF[k] = nblines;
		poly[k]->F->nbrows = nblines;
	}
	cl = cl->next;
  }
  
  free(array_map_a);
  // meet the constraints
  meet_cons_with_map(opk,oa,poly,rmapa, ca_arr, counterC, pos_con_map, disjoint_map);
  /*************************
		Add positivity constraint of A
  **************************/      
  for(k=0; k < num_comp; k++){
	size_t count = counterC[k];
        if (pos_con_map[k] && nbmapb[k]) {
          poly[k]->C->p[count][0] = 1;
          poly[k]->C->p[count][1] = 1;
          counterC[k]++;
          poly[k]->C->nbrows = counterC[k];
        }
  }	
	
  // cartesian product of vertices
  cartesian_product_vertices_with_map(oa, poly, rmapa, ca_arr, num_vertex_a, num_vertex, counterF, disjoint_map);

  // meet of rays
  meet_rays_with_map(oa, poly, rmapa, ca_arr, num_vertex_a, counterF, disjoint_map);
  char * exc_map = (char *)calloc(num_comp,sizeof(char));
  is_bottom = false;	
  for(k=0; k < num_comp; k++){
	if(nbmapb[k] && !is_bottom){
		poly[k]->satC = opt_satmat_alloc(poly[k]->F->nbrows,opt_bitindex_size(poly[k]->C->nbrows));
		is_bottom = opt_poly_meet_elina_lincons_array(opk->funopt->algorithm<0,
				      man,poly[k],poly[k],arr+k);
		if(opk->exn){

                        opk->exn = ELINA_EXC_NONE;
			exc_map[k]=1;
		}
	  	
	}
	free(ca_arr[k]);
	elina_lincons0_array_clear(arr+k);
  }

  if(destructive){
	for(k=0; k < num_compa; k++){
		//unsigned short int ind = rmapa[k];
		if(!disjoint_map[k]){
			opt_poly_clear(poly_a[k]);
		}
		free(poly_a[k]);
	}
	free(poly_a);
  }
  
  if(!is_bottom){
	array_comp_list_t * tmpa = oa->acl;
	
	if(destructive){
		free_array_comp_list(tmpa);
	}
        unsigned short int k1 = 0;
        unsigned short int bound = num_comp;
        cl = acl->head;
        for (k = 0; k < num_comp; k++) {
          opt_pk_t *oak = poly[k1];
          if (exc_map[k]) {
            comp_list_t *tmp = cl;
            cl = cl->next;
            remove_comp_list(acl, tmp);
            unsigned short int k2;
            for (k2 = k1; k2 < bound - 1; k2++) {
              poly[k2] = poly[k2 + 1];
            }
            opt_poly_clear(oak);
            bound--;
          } else {
            k1++;
            cl = cl->next;
          }
        }
        op->acl = acl;
	op->poly = poly;
	op->is_bottom = false;
  }
  else{
	for(k = 0; k < num_comp; k++){
		opt_pk_t * op_k = poly[k];
		if(op_k){
			opt_poly_clear(op_k);
			free(op_k);
		}
	}
	free(poly);
	free_array_comp_list(acl);
	op->is_bottom = true;
	
  }
  free(comp_size_map);
  free(rmapa);
  free(rmapb);
  free(nbmapa);
  free(nbmapb);
  free(nbeqmapa);
  free(counter);
  free(counterC);
  free(counterF);
  free(ca_arr);
  free(arr);
  free(pos_con_map);
  free(is_trivial);
  free(disjoint_map);
  free(num_vertex);
  free(num_vertex_a);
  free(nbgenmapa);
  free(nblinemapa);
  free(exc_map);
  free_array_comp_list(aclb);
  //for(k=0; k<num_comp; k++){
//	opt_matrix_fprint(stdout,poly[k]->F);
 // }
  //printf("meet lincons output\n");
  //elina_lincons0_array_t arr1 = opt_pk_to_lincons_array(man,op);
  //elina_lincons0_array_fprint(stdout,&arr1,NULL);
  //elina_lincons0_array_clear(&arr1);
  //	fflush(stdout);
  return op;
}

opt_pk_array_t* opt_pk_meet_lincons_array(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, elina_lincons0_array_t* array){
	#if defined(TIMING)
		start_timing();
	#endif
	opt_pk_array_t * op = opt_pk_meet_lincons_array_cons(man,destructive,oa,array);
	#if defined(TIMING)
		record_timing(meet_lincons_time);
	#endif
	return op;
}
/************************************************
	Meet with tcons array
*************************************************/
opt_pk_array_t* opt_pk_meet_tcons_array(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, elina_tcons0_array_t* array)
{
  opt_pk_array_t *op;
  op= elina_generic_meet_intlinearize_tcons_array(man,destructive,oa,array,
						  ELINA_SCALAR_MPQ, ELINA_LINEXPR_LINEAR,
						  &opt_pk_meet_lincons_array);
  return op;
}



/* ********************************************************************** */
/* III. Join */
/* ********************************************************************** */

/* ====================================================================== */
/* III.1 Join of two or more polyhedra, functional and side-effect versions */
/* ====================================================================== */


bool are_vertex_disjoint(opt_pk_internal_t *opk, opt_numint_t * v1, opt_numint_t * v2, unsigned short int size){
	unsigned short int j;
	char * map1 = (char *)calloc(size - 2, sizeof(char));
	char * map2 = (char *)calloc(size - 2, sizeof(char));
	for(j = opk->dec; j < size; j++){
		map1[j-2] = (v1[j]> 0);
		map2[j-2] = (v2[j]> 0);
	}
	for(j=0; j < size-2; j++){
		if(map1[j]&&map2[j]){
			free(map1);
			free(map2);
			return false;
		}
	}
	free(map1);
	free(map2);
	return true;
}

bool is_vertex_included(opt_pk_internal_t *opk, opt_numint_t * v1, opt_numint_t * v2, unsigned short int size){
	unsigned short int j;
	char * map1 = (char *)calloc(size - 2, sizeof(char));
	char * map2 = (char *)calloc(size - 2, sizeof(char));
	for(j = opk->dec; j < size; j++){
		map1[j-2] = (v1[j]> 0);
		map2[j-2] = (v2[j]> 0);
	}
	for(j=0; j < size-2; j++){
		if(map1[j]&&map2[j]){
			if(v1[j+2]!=v2[j+2]){
				free(map1);
				free(map2);
				return false;
			}
		}
		if(map1[j]&&!map2[j]){
			free(map1);
			free(map2);
			return false;
		}
	}
	free(map1);
	free(map2);
	return true;
}


void combine_vertex(opt_pk_internal_t *opk, opt_numint_t *v1, opt_numint_t *v2, unsigned short int size){
	unsigned short int j;
	for(j=opk->dec; j < size; j++){
		v1[j] = v1[j] + v2[j];
	}
}

/*****************************
		Join based on generator representation
*******************************/
opt_pk_array_t * opt_poly_join_gen(elina_manager_t *man, opt_pk_array_t *oa, opt_pk_array_t *ob, bool destructive){

        join_count++;
	opt_pk_internal_t *opk = opt_pk_init_from_manager(man, ELINA_FUNID_JOIN);
	size_t i;
	unsigned short int j,k;
	unsigned short int maxcols = oa->maxcols;
	array_comp_list_t * acla = oa->acl;
	array_comp_list_t * aclb = ob->acl;
	opt_pk_array_t *op = destructive ? oa : opt_pk_array_alloc(NULL,NULL,maxcols);
	if(oa->is_bottom || !acla){
		if(destructive){
			opt_poly_array_clear(opk,op);
			free(op);
		}
		else{
			free(op);
		}
		op = opt_pk_copy(man,ob);
		
		return op;
	}
	
	if(ob->is_bottom || !aclb){
		if(destructive){
			op = oa;
		}
		else{
			free(op);
			op = opt_pk_copy(man,oa);
		}
		
		return op;
	}
	
	/***********************
			Minimize A and compute components that could be connected by sigmas
	************************/
        opt_pk_t **poly_a = oa->poly;
        unsigned short int num_compa = acla->size;
        comp_list_t *cla = acla->head;
        size_t * num_vertex_a = (size_t *)calloc(num_compa,sizeof(size_t));
	for(k=0; k < num_compa; k++){
		opt_pk_t * oak = poly_a[k];
		
		if(opk->funopt->algorithm>=0){
			opt_poly_chernikova(man,oak,"cons to gen");
		}
		else{	
			opt_poly_obtain_F(man,oak,"cons to gen");
		}
		if(opk->exn){
			opt_poly_set_top(opk,op);
			
			return op;
		}
		if(!oak->F){
			if(destructive){
				opt_poly_array_clear(opk,op);	
  				free(op);
			}
			else{
				free(op);
			}
			op = opt_pk_copy(man,ob);
			return op;
		}
		/*if(destructive){
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
		}*/
		//printf("A \n");
		//opt_matrix_fprint(stdout,oak->F);
		//opt_matrix_fprint(stdout,oak->C);
		//opt_satmat_fprint(stdout,oak->satF);
		//fflush(stdout);
		opt_poly_obtain_satF(oak);
		num_vertex_a[k] = opt_generator_rearrange(oak->F,oak->satF);
		cla = cla->next;
	}	
	
	
	/***********************
			Minimize B and compute components that could be connected by sigmas
	************************/
        opt_pk_t **poly_b = ob->poly;
        unsigned short int num_compb = aclb->size;
        comp_list_t *clb = aclb->head;
        size_t * num_vertex_b = (size_t *)calloc(num_compb,sizeof(size_t));
	for(k=0; k < num_compb; k++){
		opt_pk_t * obk = poly_b[k];
		if(opk->funopt->algorithm >=0){
			opt_poly_chernikova(man,obk,"cons to gen");
		}
		else{
			opt_poly_obtain_F(man,obk,"cons to gen");
		}
		if(opk->exn){
			opt_poly_set_top(opk,op);
			
			return op;
		}
		if(!obk->F){
			if(destructive){
				return oa;
			}
			else{
				free(op);
				op = opt_pk_copy(man,oa);
				return op;
			}
		}
		opt_poly_obtain_satF(obk);
		num_vertex_b[k] = opt_generator_rearrange(obk->F,obk->satF);
		clb = clb->next;
	}
	
	/**************************
			Compute union of independent components
	***************************/
	array_comp_list_t * acl = union_array_comp_list(acla,aclb,maxcols);
	unsigned short int num_comp = acl->size;
	char ** var_map_a = (char **)malloc(num_comp*sizeof(char*));
	char ** var_map_b = (char **)malloc(num_comp*sizeof(char*));
	opt_pk_t ** poly1 = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
 	opt_pk_t ** poly2 = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
	size_t * num_vertex1 = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * num_vertex2 = (size_t *)calloc(num_comp,sizeof(size_t));	
	size_t * nbmapCa = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * nblinemapa = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * nbeqmapa = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * nbmapCb = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * nblinemapb = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * nbeqmapb = (size_t *)calloc(num_comp,sizeof(size_t));		
	unsigned short int ** ca_arr = (unsigned short int **)malloc(num_comp*sizeof(unsigned short int *));


	comp_list_t * cl = acl->head;
        for(k=0; k < num_comp; k++){
		var_map_a[k] = (char *)calloc(cl->size, sizeof(char)); 
		var_map_b[k] = (char *)calloc(cl->size, sizeof(char)); 
		ca_arr[k] = to_sorted_array(cl,maxcols);
		cl = cl->next;
	}
	/**************************
		Partition A according to union 
	***************************/
	unsigned short int * rmapa = (unsigned short int *)calloc(num_compa,sizeof(unsigned short int));
	size_t * nbmapFa = (size_t *)calloc(num_comp,sizeof(size_t));
	
	cla = acla->head;
	for(k=0; k < num_compa; k++){
		opt_pk_t *oak = poly_a[k];
		unsigned short int inda = is_comp_list_included(acl,cla,maxcols);
		rmapa[k] = inda;	
		if(num_vertex_a[k]){
			if(!num_vertex1[inda]){
				num_vertex1[inda] = num_vertex_a[k];
			}
			else{
				num_vertex1[inda] = num_vertex1[inda] * num_vertex_a[k];
			}
		}
		
		nbmapFa[inda] = nbmapFa[inda] + oak->F->nbrows; 
		nbmapCa[inda] = nbmapCa[inda] + oak->C->nbrows; 
		nblinemapa[inda] = nblinemapa[inda] + oak->nbline;
		nbeqmapa[inda] = nbeqmapa[inda] + oak->nbeq;
		unsigned short int *ca = to_sorted_array(cla,maxcols);
		unsigned short int * ind_map = map_index(ca,ca_arr[inda],cla->size);
		unsigned short int k1;
		for(k1=0; k1 < cla->size; k1++){
			unsigned short int num = ind_map[k1];
			var_map_a[inda][num]++;
		}
		free(ca);
		free(ind_map);
		cla = cla->next;
	}
 
	/**************************
		Partition B according to union
	***************************/
	unsigned short int * rmapb = (unsigned short int *)calloc(num_compb,sizeof(unsigned short int));
	
	size_t * nbmapFb = (size_t *)calloc(num_comp,sizeof(size_t));
	clb = aclb->head;
	for(k=0; k < num_compb; k++){
		opt_pk_t * obk = poly_b[k];
		unsigned short int indb = is_comp_list_included(acl,clb,maxcols);
		rmapb[k] = indb;
		
		
		if(num_vertex_b[k]){
			if(!num_vertex2[indb]){
				num_vertex2[indb] = num_vertex_b[k];
			}
			else{
				num_vertex2[indb] = num_vertex2[indb] * num_vertex_b[k];
			}
		}
		nbmapCb[indb] = nbmapCb[indb] + obk->C->nbrows;
		nbeqmapb[indb] = nbeqmapb[indb] + obk->nbeq;
		nbmapFb[indb] = nbmapFb[indb] + obk->F->nbrows;
		nblinemapb[indb] = nblinemapb[indb] + obk->nbline;
		//printf("B %d\n",num_vertex_b[k]);
		//opt_matrix_fprint(stdout,obk->F);
		//opt_matrix_fprint(stdout,obk->C);
		//opt_satmat_fprint(stdout,obk->satF);
		//fflush(stdout);
		unsigned short int *ca = to_sorted_array(clb,maxcols);
		unsigned short int * ind_map = map_index(ca,ca_arr[indb],clb->size);
		unsigned short int k1;
		for(k1=0; k1 < clb->size; k1++){
			unsigned short int num = ind_map[k1];
			var_map_b[indb][num]++;
		}
		free(ca);
		free(ind_map);
		clb = clb->next;
	}
	
	cl = acl->head;
	size_t * counterFa = (size_t *)calloc(num_comp, sizeof(size_t));
	size_t * counterFb = (size_t *)calloc(num_comp, sizeof(size_t));
	size_t * counterCa = (size_t *)calloc(num_comp, sizeof(size_t));
	size_t * counterCb = (size_t *)calloc(num_comp, sizeof(size_t));
	//char ** vertex_map = (char **)malloc(num_comp*sizeof(char *));
	for(k = 0; k < num_comp; k++){
		
		unsigned short int comp_size = cl->size;
		size_t gen_size, nblines=0, nbeq=0;
		unsigned short int k1;
		
		// A
		poly1[k] = opt_poly_alloc(comp_size,0);
		if(!nbmapFa[k]){
		    gen_size = 1;
		}
		else{
			nblines = 0;			
			for(k1=0; k1 < cl->size; k1++){
				if(!var_map_a[k][k1]){
					nblines++;
					nbeq++;
				}
			}
			gen_size = nbmapFa[k] + 2*num_vertex1[k] + nblines;
		}
		poly1[k]->F = opt_matrix_alloc(gen_size,comp_size+2,false);
		poly1[k]->C = opt_matrix_alloc(nbmapCa[k]+1+nbeq,comp_size+2,false);
		poly1[k]->nbline = nblinemapa[k]+nblines;
		poly1[k]->nbeq = nbeqmapa[k]+nbeq;
		
		// B
		poly2[k] = opt_poly_alloc(comp_size,0);
		if(!nbmapFb[k]){
		    gen_size =  1;
		}
		else{
		    nblines = 0;
		    for(k1=0; k1 < cl->size; k1++){
			if(!var_map_b[k][k1]){
				nblines++;
			}
		    }
		    gen_size = nbmapFb[k] + 2*num_vertex2[k] + nblines;
		}
		
		poly2[k]->F = opt_matrix_alloc(gen_size,comp_size+2,false);
		poly2[k]->C = opt_matrix_alloc(nbmapCb[k]+1,comp_size+2,false);
		num_vertex1[k] = 0;
		num_vertex2[k] = 0;
		poly2[k]->nbline = nblinemapb[k]+nblines;
		poly2[k]->nbeq = nbeqmapb[k];
		cl = cl->next;
		
	}
	
	
	/*************************
		Cartesian Product of Vertices from  A
	************************/

        cartesian_product_vertices(oa, poly1, rmapa, ca_arr, num_vertex_a,
                                   num_vertex1, counterFa);

        /************************
		Now Consider rays of A
	************************/
        meet_rays(oa, poly1, rmapa, ca_arr, num_vertex_a, counterFa);

        //for(k=0; k < num_comp; k++){
	//	size_t count = counterF[k];
	//	start_counter[k] = count;
	//}

	/***********************
		Meet constraints of A
	************************/
	char * pos_con_map = (char *)calloc(num_comp, sizeof(char)); 
	for(k=0; k < num_comp; k++){
		pos_con_map[k] = 1;
	}
        meet_cons(opk, oa, poly1, rmapa, ca_arr, counterCa, pos_con_map);

        /*************************
		Add positivity constraint of A
	**************************/
        cl = acl->head;
	for(k=0; k < num_comp; k++){
		size_t countC = counterCa[k];
		size_t countF = counterFa[k];
		nblinemapa[k] = 0;
		if(pos_con_map[k]){
			poly1[k]->C->p[countC][0] = 1;
			poly1[k]->C->p[countC][1] = 1;
			countC++;
		}
		if(nbmapFa[k]){
			unsigned short int k1;
			for(k1=0; k1 < cl->size; k1++){
				if(!var_map_a[k][k1]){
					poly1[k]->F->p[countF][k1+2]=1;
					poly1[k]->C->p[countC][k1+2]=1;
					nblinemapa[k]++;
					countC++;
					countF++;
				}
			}
		}
		counterCa[k] = countC;
		counterFa[k] = countF;
		poly1[k]->C->nbrows = counterCa[k];
		poly1[k]->F->nbrows = counterFa[k];
		cl = cl->next;
	}	
	
	
	
	/*************************
		Cartesian Product of Vertices from B
	************************/
        cartesian_product_vertices(ob, poly2, rmapb, ca_arr, num_vertex_b,
                                   num_vertex2, counterFb);

        /************************
		 Consider rays of B
	************************/
        meet_rays(ob, poly2, rmapb, ca_arr, num_vertex_b, counterFb);

        /***********************
		Meet constraints of B
	************************/ 
	for(k=0; k < num_comp; k++){
		pos_con_map[k] = 1;
	}
        meet_cons(opk, ob, poly2, rmapb, ca_arr, counterCb, pos_con_map);

        /*************************
		Add positivity constraint of B
	**************************/
        cl = acl->head;
	for(k=0; k < num_comp; k++){
		size_t countC = counterCb[k];
		size_t countF = counterFb[k];
		if(pos_con_map[k]){
			poly2[k]->C->p[countC][0] = 1;
			poly2[k]->C->p[countC][1] = 1;
			countC++;
		}
		if(nbmapFb[k]){
			unsigned short int k1;
			for(k1=0; k1 < cl->size; k1++){
				if(!var_map_b[k][k1]){
					poly2[k]->F->p[countF][k1+2]=1;
					countF++;
				}
			}
		}
		counterCb[k] = countC;
		counterFb[k] = countF;
		poly2[k]->C->nbrows = counterCb[k];
		poly2[k]->F->nbrows = counterFb[k];
		cl = cl->next;
	}	
	

	unsigned short int num_comp_res = 0;
	bool flag = false;
	cl = acl->head;
	for(k=0; k < num_comp; k++){
		pos_con_map[k] = 0;
		opt_matrix_t * Ca = poly1[k]->C;
		opt_matrix_t * Fa = poly1[k]->F;
		
		opt_matrix_t * Cb = poly2[k]->C;
		opt_matrix_t * Fb = poly2[k]->F;
		if(!nbmapFa[k] || !nbmapFb[k]){
			pos_con_map[k] = 3;
		}
		else if(opt_poly_leq(opk,Ca,Fb)&&opt_poly_leq(opk,Cb,Fa)){
			num_comp_res++;
			// take B
			pos_con_map[k] = 2;
		}
		//else if(opt_poly_leq(opk,Cb,Fa)){
		//	num_comp_res++;
			//take A
			
			//pos_con_map[k] = 1;
		//}
		else{
			flag = true;
		}
		cl =cl->next;
	}

	if(flag){
		num_comp_res++;
	}
	opt_pk_t ** poly = (opt_pk_t **)malloc(num_comp_res*sizeof(opt_pk_t *));
	unsigned short int k1 = 0;
	comp_list_t * clp = create_comp_list();
	array_comp_list_t * res = create_array_comp_list();
	char * clp_map = (char *)calloc(maxcols, sizeof(char));
	size_t num_vertex1a = 0, num_vertex2b=0;
	size_t nbF = 0, nbC = 0, nbeq = 0, nbline = 0, nblinea = 0;
	cl = acl->head;

        for(k=0; k < num_comp; k++){
		unsigned short int k2 = num_comp_res - k1 - 1;
		unsigned short int comp_size = cl->size;
		if(pos_con_map[k]==1){	
			poly[k2] = opt_poly_alloc(comp_size,0);		
			poly[k2]->C = poly1[k]->C;
			poly[k2]->F = poly1[k]->F;
			poly[k2]->nbeq = poly1[k]->nbeq;
			poly[k2]->nbline = poly1[k]->nbline;
			poly[k2]->satF = opt_satmat_alloc(poly[k2]->C->nbrows,opt_bitindex_size(poly[k2]->F->nbrows));
			combine_satmat(opk,poly[k2],comp_size,poly[k2]->F->nbrows,false);
			insert_comp_list(res,copy_comp_list(cl));
			
			k1++;
			
		}
		else if(pos_con_map[k]==2){
			poly[k2] = opt_poly_alloc(comp_size,0);
			poly[k2]->C = poly2[k]->C;
			poly[k2]->F = poly2[k]->F;
			opt_matrix_sort_rows(opk,poly[k2]->F);
			poly[k2]->nbeq = poly2[k]->nbeq;
			poly[k2]->nbline = poly2[k]->nbline;
			//poly[num_comp_res - k1 - 1]->satC = poly2[k]->satC;
			poly[k2]->satF = opt_satmat_alloc(poly[k2]->C->nbrows,opt_bitindex_size(poly[k2]->F->nbrows));
			combine_satmat(opk,poly[k2],comp_size,poly[k2]->F->nbrows,false);
			insert_comp_list(res,copy_comp_list(cl));
			
			k1++;
		}
		else if(!pos_con_map[k]){
			
			union_comp_list(clp,cl,clp_map);
			if(!num_vertex1a){
				num_vertex1a = num_vertex1[k];
			}
			else{
				num_vertex1a *= num_vertex1[k];
			}
			if(!num_vertex2b){
				num_vertex2b = num_vertex2[k];
			}
			else{
				num_vertex2b *= num_vertex2[k];
			}
			nbF = nbF + poly1[k]->F->nbrows + poly2[k]->F->nbrows;
			nbC = nbC + poly1[k]->C->nbrows;
			nbeq = nbeq + poly1[k]->nbeq;
			nbline = nbline + poly1[k]->nbline + poly2[k]->nbline;
			nblinea += nblinemapa[k];  
			//comp_list_t * tcl = cl;
			
			//free_comp_list(tcl);
		}
		cl = cl->next;
	}
	
	free(clp_map);

        if(flag){
          // printf("start %d %d\n",clp->size,oa->maxcols-2);
          // fflush(stdout);
          unsigned short int comp_size = clp->size;
          poly[0] = opt_poly_alloc(comp_size, 0);
          poly[0]->F = opt_matrix_alloc(nbF + 2 * (num_vertex1a + num_vertex2b),
                                        comp_size + 2, false);
          poly[0]->C = opt_matrix_alloc(nbC, comp_size + 2, false);
          poly[0]->nbeq = nbeq;
          poly[0]->nbline = nbline;
          unsigned short int *ca = to_sorted_array(clp, maxcols);
          size_t count = 0;
          // combine all overlapping vertices of A into one
          count = cartesian_product_vertices_one_comp(
              poly1, acl, ca_arr, num_vertex1, poly[0], count, ca, pos_con_map);
          // combine all overlapping rays of A into one

          meet_rays_one_comp(poly1, acl, ca_arr, nblinemapa, poly[0],
                             num_vertex1, ca, count, pos_con_map);
          size_t begin = poly[0]->F->nbrows - nblinea;
          // combine all overlapping constraints of A into one
          bool is_pos_con = meet_cons_one_comp(opk, poly1, acl, ca_arr, poly[0],
                                               ca, pos_con_map);
          if (is_pos_con) {
            size_t count = poly[0]->C->nbrows;
            poly[0]->C->p[count][0] = 1;
            poly[0]->C->p[count][1] = 1;
            count++;
            poly[0]->C->nbrows = count;
                }
		
		//combine all overlapping vertices of B into one
		count = cartesian_product_vertices_one_comp(poly2,acl,ca_arr,num_vertex2,poly[0],poly[0]->F->nbrows,ca,pos_con_map);
		//combine all overlapping rays of B into one

                meet_rays_one_comp(poly2,acl,ca_arr,NULL,poly[0],num_vertex2,ca,count,pos_con_map);
                free(ca);
                opt_matrix_t * F = poly[0]->F;
		opt_matrix_t * C = poly[0]->C;
		//remove_common_gen(opk,F,begin);
		opt_matrix_sort_rows(opk,C);
		/************************
			Combine satmat of A
		*************************/
		poly[0]->satF = opt_satmat_alloc(C->nbrows,opt_bitindex_size(F->nbrows));
		combine_satmat(opk,poly[0],comp_size,begin,false);			
		size_t num = F->nbrows - begin;
		opt_matrix_sort_rows_from(opk,F,begin,num);
		opt_poly_dual(poly[0]);
		opt_cherni_add_and_minimize(opk,false,poly[0],begin);
		opt_poly_dual(poly[0]);
		if(opk->exn){
                  opk->exn = ELINA_EXC_NONE;
                  opt_pk_t *tmp = poly[0];
                  unsigned short int k1;
                  for (k1 = 0; k1 < num_comp_res - 1; k1++) {
                    poly[k1] = poly[k1 + 1];
			}
			opt_poly_clear(tmp);
		}
		else{
			insert_comp_list(res,clp);
		}
                // printf("finish %lld\n",poly[0]->F->nbrows);
                // fflush(stdout);
        }

        if(destructive){
		for(k=0; k < num_compa; k++){
			opt_poly_clear(poly_a[k]);
		}
		array_comp_list_t * tmp = acla;
		free(poly_a);
		op->poly = poly;
		op->acl = res;
		free(tmp);
	}
	else{
		op->poly = poly;
		op->acl = res;
	}

        //cl = acl->head;
	//unsigned short int k2;
	//unsigned short int nc = num_comp;
	for(k=0; k< num_comp; k++){
		if(pos_con_map[k]==1){
			opt_poly_clear(poly2[k]);
		}
		else if(pos_con_map[k]==2){
			opt_poly_clear(poly1[k]);
		}
		else{
			opt_poly_clear(poly1[k]);
			opt_poly_clear(poly2[k]);
			
		}
		free(poly1[k]);
		free(poly2[k]);
		free(ca_arr[k]);
		free(var_map_a[k]);
		free(var_map_b[k]);
	}
	free(rmapa);
	free(rmapb);
	free(counterFa);
	free(counterFb);
	free(counterCa);
	free(counterCb);
	free(nbmapCa);
	free(nbmapCb);
	free(nbmapFa);
	free(nbmapFb);
	free(nbeqmapa);
	free(nbeqmapb);
	free(nblinemapa);
	free(nblinemapb);
	free(var_map_a);
	free(var_map_b);
	free(ca_arr);
	free(poly1);
	free(poly2);
	free(num_vertex1);
	free(num_vertex2);
	free(num_vertex_a);
	free(num_vertex_b);
	free(pos_con_map);
	free_array_comp_list(acl);

        return op;
}


opt_pk_array_t* opt_pk_join(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, opt_pk_array_t* ob)
{
  opt_pk_array_t * res;
  #if defined (TIMING)
	start_timing();
  #endif
  res = opt_poly_join_gen(man,oa,ob,destructive);
  #if defined (TIMING)
	record_timing(join_time);
  #endif
  return res;
}


void opt_poly_meet(bool meet,
	       bool lazy,
	       elina_manager_t* man,
	       opt_pk_array_t* op, opt_pk_array_t* oa, opt_pk_array_t* ob){
	
	opt_pk_internal_t *opk = opt_pk_init_from_manager(man, ELINA_FUNID_MEET);
	array_comp_list_t *acla = oa->acl;
	unsigned short int num_compa = acla->size;
	array_comp_list_t *aclb = ob->acl;
	unsigned short int num_compb = aclb->size;
	unsigned short int maxcols = oa->maxcols;
	/*************************
		Compute union of independent components
	*************************/
	array_comp_list_t *acl = union_array_comp_list(acla, aclb, maxcols);
	unsigned short int num_comp = acl->size;
	opt_pk_t **poly = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
	size_t * nbeqmap = (size_t *)calloc(num_comp,sizeof(size_t));
	/*****************************
		Factor A according to union
	*****************************/	
	opt_pk_t ** poly_a = oa->poly;	
	unsigned short int * rmapa = (unsigned short int *)calloc(num_compa, sizeof(unsigned short int));
	size_t * nbmapa = (size_t *)calloc(num_comp,sizeof(size_t));
	comp_list_t * cla = acla->head;
	size_t i;
	unsigned short int k,j;
	for(k = 0; k < num_compa; k++){
		opt_pk_t * oak = poly_a[k];	
		opt_matrix_t * ocak = oak->C;
		short int ind = is_comp_list_included(acl,cla,maxcols);
		rmapa[k] = ind;
		nbmapa[ind] = nbmapa[ind] + ocak->nbrows;
		nbeqmap[ind] = nbeqmap[ind] + oak->nbeq;
		cla = cla->next;
	}

	/*****************************
		Factor B according to union
	*****************************/		
	opt_pk_t ** poly_b = ob->poly;
	unsigned short int * rmapb = (unsigned short int *)calloc(num_compb, sizeof(unsigned short int));
	size_t * nbmapb = (size_t *)calloc(num_comp,sizeof(size_t));
	comp_list_t * clb = aclb->head;
	for(k = 0; k < num_compb; k++){
		opt_pk_t * obk = poly_b[k];
		opt_matrix_t * ocbk = obk->C;
		short int ind = is_comp_list_included(acl,clb,maxcols);
		rmapb[k] = ind;
		nbmapb[ind] = nbmapb[ind] + ocbk->nbrows;
		nbeqmap[ind] = nbeqmap[ind] + obk->nbeq;
		//opt_matrix_sort_rows(opk,ocbk);
		clb = clb->next;
	}

	comp_list_t * cl = acl->head;
	size_t * counterC = (size_t *)calloc(num_comp, sizeof(size_t));
	unsigned short int ** ca_arr = (unsigned short int **)malloc(num_comp*sizeof(unsigned short int *));
	for(k = 0; k < num_comp; k++){
		unsigned short int comp_size = cl->size;
		poly[k] = opt_poly_alloc(comp_size,0);
		poly[k]->C = opt_matrix_alloc(nbmapa[k] + nbmapb[k]+1,comp_size+2,false);
		poly[k]->C->p[0][0] = 1;
		poly[k]->C->p[0][1] = 1;
		poly[k]->nbeq = nbeqmap[k];
		counterC[k] = 1;
		ca_arr[k] = to_sorted_array(cl,maxcols);
		cl = cl->next;
	}

	cla = acla->head;
	for(k = 0; k < num_compa; k++){
		opt_pk_t * src = poly_a[k];
		unsigned short int ind = rmapa[k];
		opt_pk_t * dst = poly[ind];
		unsigned short int * ca_a = to_sorted_array(cla,maxcols);
		unsigned short int * ca = ca_arr[ind];
		opt_matrix_t * src_mat;
		opt_matrix_t * dst_mat;
		size_t * counter;
		src_mat = src->C;
		dst_mat = dst->C;
		counter = counterC;
		unsigned short int comp_size = cla->size;
		size_t nbconsa = src_mat->nbrows;
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = dst_mat->p;
		size_t i1 = counter[ind];
		counter[ind] = counter[ind] + nbconsa;
		for(i = 0; i < nbconsa; i++){
		    opt_numint_t * src_pi = src_p[i];
		    opt_numint_t * dst_pi = dst_p[i1];
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
		   i1++;
	        }
		free(ca_a);
		cla = cla->next;
	}
	
	clb = aclb->head;
	for(k = 0; k < num_compb; k++){
		opt_pk_t * src = poly_b[k];
		unsigned short int ind = rmapb[k];
		opt_pk_t * dst = poly[ind];
		unsigned short int * ca_b = to_sorted_array(clb, maxcols);
		unsigned short int * ca = ca_arr[ind];
		opt_matrix_t * src_mat = src->C;
		opt_matrix_t * dst_mat = dst->C;
		unsigned short int comp_size = clb->size;
		size_t nbconsb = src_mat->nbrows;
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = dst_mat->p;
		size_t *counter = counterC;
		size_t i1 = counter[ind];
		counter[ind] = counter[ind] + nbconsb;
		for(i = 0; i < nbconsb; i++){
			opt_numint_t * src_pi = src_p[i];
			opt_numint_t * dst_pi = dst_p[i1];
			dst_pi[0] = src_pi[0];
			dst_pi[1] = src_pi[1];
			unsigned short int l = 0;
			for(j = 0; j < comp_size; j++){
				while(ca[l] != ca_b[j]){
					l++;
				}
				dst_pi[l+2] = src_pi[j+2];
				l++;
			}
			i1++;
		}
		free(ca_b);
		clb = clb->next;
	}
	
	/*****************************
		Sort rows of each block of op
	*****************************/
	bool is_bottom = false;
	for(k=0; k < num_comp; k++){
		opt_pk_t * src = poly[k];
		opt_poly_chernikova(man,src,"meet abstract");
		if(src->F==NULL){
		   opt_poly_set_bottom(opk,op);
		   is_bottom = true;
		   break;
		}
	}
	if(!is_bottom){
		array_comp_list_t * tmp = oa->acl;
		if(op==oa){
			free_array_comp_list(tmp);
		}
		op->acl = acl;
		op->poly = poly;
	}
	else{
		for(k = 0; k < num_comp; k++){
			opt_pk_t * op_k = poly[k];
			if(op_k){
				opt_poly_clear(op_k);
			}
		}
		free(poly);
		free_array_comp_list(acl);
	}

	for(k=0; k < num_comp; k++){
		free(ca_arr[k]);
	}
	
	free(rmapa);
	free(rmapb);
	free(nbmapa);
	free(nbmapb);
	free(nbeqmap);
	free(ca_arr);
	free(counterC);
	//#if defined (CONVERT)
	//	free(nbgenmap);
	//	free(nblinemap);
	//	free(counterF);
	//	free(counterS);
	//	free(colmapS);
	//#endif
}


/* ********************************************************************** */
/* II. Meet */
/* ********************************************************************** */

/* ********************************************************************** */
/* II.1 Meet of two or more polyhedra */
/* ********************************************************************** */

opt_pk_array_t* opt_pk_meet_cons(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, opt_pk_array_t* ob){
	opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_MEET);
	opt_pk_array_t* op = destructive ? oa :  opt_pk_array_alloc(NULL,NULL,oa->maxcols);
	opt_pk_t ** poly_a = oa->poly;
	array_comp_list_t * acla = oa->acl;
	array_comp_list_t * aclb = ob->acl;

	if(oa->is_bottom || !acla){
		if(destructive){
			return oa;
		}
		else{
			free(op);
			return opt_pk_bottom(man,oa->maxcols - 2, 0);
		}
        }
	if(ob->is_bottom || !aclb){
		if(destructive){
			opt_poly_set_bottom(opk,oa);
			return oa;
		}
		else{
			free(op);
			return opt_pk_bottom(man,oa->maxcols - 2, 0);
		}
	}
	unsigned short int num_compa = acla->size;
	unsigned short int k;
	for(k = 0; k < num_compa; k++){
		opt_pk_t * oak = poly_a[k];
		if(opk->funopt->algorithm>=0){
		   opt_poly_chernikova(man,oak,"meet abstract");
		}
		else{	
		    opt_poly_obtain_C(man,oak,"meet abstract");
		}
		if(opk->exn){
		   opk->exn = ELINA_EXC_NONE;
		   if(destructive){
			opt_poly_array_clear(opk,op);
			free(op);
		   }
		   else{
		   	free(op);
		   }
		   op = opt_pk_copy(man,ob);
		   return op;
		}
		if(!oak->C){
		   opt_poly_set_bottom(opk,op);
		   return op;
		}
	}
	
	opt_pk_t ** poly_b = ob->poly;
	
	unsigned short int num_compb = aclb->size;
	for(k = 0; k < num_compb; k++){
	    opt_pk_t * obk = poly_b[k];
	    if(opk->funopt->algorithm>=0){
	       opt_poly_chernikova(man,obk,"meet abstract");
	    }
	    else{	
		opt_poly_obtain_C(man,obk,"meet abstract");
	    }
	    if(opk->exn){
	       opk->exn = ELINA_EXC_NONE;
	       if(destructive){
		  return oa;
	       }
	       else{
		  free(op);
		  op = opt_pk_copy(man,oa);
		  return op;
	       }
	   }
	   if(!obk->C){
	      opt_poly_set_bottom(opk,op);
	      return op;
	   }
	}
	
	opt_poly_meet(true, opk->funopt->algorithm < 0,
		  man, op,oa,ob);
	
	return op;
}


opt_pk_array_t* opt_pk_meet(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, opt_pk_array_t* ob){
	#if defined (TIMING)
		start_timing();
	#endif
	opt_pk_array_t *op = opt_pk_meet_cons(man,destructive,oa,ob);
	#if defined (TIMING)
		record_timing(meet_time);
	#endif
	return op;
}




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

#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"
#include "opt_pk_test.h"
#include "opt_pk_representation.h"


opt_pk_array_t * opt_pk_array_alloc(opt_pk_t ** poly, array_comp_list_t *acl, unsigned short int maxcols){
	opt_pk_array_t * op = (opt_pk_array_t *)malloc(sizeof(opt_pk_array_t));
	op->poly = poly;
	op->acl = acl;
	op->maxcols = maxcols;
	op->is_bottom = false;
	return op;
}

opt_pk_t* opt_poly_alloc(size_t intdim, size_t realdim)
{
  opt_pk_t* op = (opt_pk_t*)malloc(sizeof(opt_pk_t));
  op->C = NULL;
  op->F = NULL;
  op->satC = op->satF = NULL;
  op->intdim = intdim;
  op->realdim = realdim;
  op->nbeq = 0;
  op->nbline = 0;
  op->status = 0;
  op->is_minimized = false;
  return op;
}

/* Clearing a polyhedron */
void opt_poly_clear(opt_pk_t* opo)
{
  if(opo->C){
  	opt_matrix_free(opo->C);
  }

  if(opo->F){
	opt_matrix_free(opo->F);
  }
  
  if(opo->satC){
	opt_satmat_free(opo->satC);
  }
  if(opo->satF){
	opt_satmat_free(opo->satF);
  }
  opo->status = 0;
  opo->nbeq = 0;
  opo->nbline = 0;
}

void opt_poly_array_clear(opt_pk_internal_t *opk, opt_pk_array_t * op){
  array_comp_list_t * acl = op->acl;
  if(acl){
	opt_pk_t ** poly = op->poly;
	unsigned short int num_comp = acl->size;
	unsigned short int k;
	for(k=0; k < num_comp; k++){
		opt_pk_t * opp = poly[k];
		if(opp){
			opt_poly_clear(opp);
			free(opp);
			opp=NULL;
		}
	}
	free(poly);
	poly=NULL;
	free_array_comp_list(acl);
  }
}


/* Finalization function for polyhedra, which frees
   recursively the members of the structure. */
void opt_pk_free(elina_manager_t* man, opt_pk_array_t* op)
{
  #if defined(TIMING)
 	 start_timing();
   #endif 
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_FREE);
  opt_poly_array_clear(opk,op);
  free(op);
  op=NULL;
  #if defined(TIMING)
 	 record_timing(free_time);
   #endif 
}

void opt_poly_copy(opt_pk_t *dst, opt_pk_t *src){
	dst->nbeq = src->nbeq;
	dst->C = src->C ? opt_matrix_copy(src->C) : NULL;
	dst->status = src->status;
	dst->F = src->F ? opt_matrix_copy(src->F) : NULL;
	dst->satC = src->satC ? opt_satmat_copy(src->satC) : NULL;
	dst->satF = src->satF ? opt_satmat_copy(src->satF) : NULL;
	dst->nbline = src->nbline;
}

/* Duplicate (recursively) a polyhedron. */
opt_pk_array_t* opt_pk_copy(elina_manager_t* man, opt_pk_array_t* src)
{
  #if defined(TIMING)
 	 start_timing();
  #endif
  if(src->is_bottom || !src->acl){
	#if defined(TIMING)
		record_timing(copy_time);
 	#endif
	return opt_pk_bottom(man,src->maxcols - 2,0);
  }
  array_comp_list_t *acl = copy_array_comp_list(src->acl);
  unsigned short int maxcols = src->maxcols;
  unsigned short int num_comp = acl->size;
  opt_pk_t ** poly = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
  unsigned short int k;
  opt_pk_t ** spoly = src->poly;
  for(k = 0; k < num_comp; k++){ 
	poly[k] = opt_poly_alloc(spoly[num_comp-k-1]->intdim,spoly[num_comp-k-1]->realdim);
	poly[k]->C = spoly[num_comp-k-1]->C ? opt_matrix_copy(spoly[num_comp-k-1]->C) : NULL;	
        poly[k]->F = spoly[num_comp-k-1]->F ? opt_matrix_copy(spoly[num_comp-k-1]->F) : NULL;
	poly[k]->satC = spoly[num_comp-k-1]->satC ? opt_satmat_copy(spoly[num_comp-k-1]->satC) : NULL;
	poly[k]->satF = spoly[num_comp-k-1]->satF ? opt_satmat_copy(spoly[num_comp-k-1]->satF) : NULL;
	poly[k]->nbeq = spoly[num_comp-k-1]->nbeq;
        poly[k]->nbline = spoly[num_comp-k-1]->nbline;
	poly[k]->status = spoly[num_comp-k-1]->status;
	poly[k]->is_minimized = spoly[num_comp-k-1]->is_minimized;
  }
  
  opt_pk_array_t * dst = opt_pk_array_alloc(poly,acl,maxcols);
  dst->is_bottom = src->is_bottom;
  #if defined(TIMING)
		record_timing(copy_time);
  #endif
  return dst;
}



void opt_poly_set(opt_pk_t* oa, opt_pk_t* ob)
{
  if (oa!=ob){
    opt_poly_clear(oa);
    oa->intdim = ob->intdim;
    oa->realdim = ob->realdim;
    oa->C = ob->C ? opt_matrix_copy(ob->C) : NULL;
    oa->status = ob->status;
    oa->nbeq = ob->nbeq;
  }
}


double abs_diff(double a, opt_numint_t b){
	double tmp = (double)b;
	double res = tmp - a;
	if(res < 0){
		res = -res;
	}
	return res;
}


void opt_poly_chernikova(elina_manager_t* man,
		     opt_pk_t* poly,
		     char* msg)
{
  opt_pk_internal_t* opk = (opt_pk_internal_t*)man->internal;
  if ((poly->C && poly->F) || (!poly->C && !poly->F)){
    
    return;
  }
  else {
    if (poly->C){
      //if (!poly_is_conseps(pk,po) ){
	opt_matrix_normalize_constraint(opk,poly->C,poly->intdim,poly->realdim);
      //}
     
      opt_matrix_sort_rows(opk,poly->C);
      opt_cherni_minimize(opk,true,poly);
      if (opk->exn) goto poly_chernikova_exit0;
      poly->status = opt_pk_status_consgauss;
    }
    else {
      poly->C = poly->F; poly->F = NULL;
      opt_matrix_sort_rows(opk,poly->C);
      opt_cherni_minimize(opk,false,poly);
      opt_poly_dual(poly);
      if (opk->exn) goto poly_chernikova_exit0;
      poly->status = opt_pk_status_gengauss;
    }
  }
  return;
 poly_chernikova_exit0:
  poly->status = 0;
  {
    char str[160];
    sprintf(str,"conversion from %s %s\n",
	    poly->C ? "constraints to generators" : "generators to constraints",
	    msg);
    elina_manager_raise_exception(man,opk->exn,opk->funid,str);
  }
  return;
}


void opt_pk_convert(elina_manager_t *man, opt_pk_array_t * op, char *msg){
	array_comp_list_t *acl = op->acl;
	if(op->is_bottom || !acl){
		return;
	}
	unsigned short int num_comp = acl->size;
	if(num_comp==0){
		return;
	} 
	unsigned short int k;
	opt_pk_t ** poly = op->poly;
	for(k=0; k < num_comp; k++){
		opt_pk_t * ok = poly[k];
		opt_poly_chernikova(man,ok,msg);
	}
}

//void poly_swap(opt_pk_t *poly1, opt_pk_t *poly2){
//	opt_pk_t *tmp = poly1;
//	poly1 = poly2;
//	poly2 = tmp;
//}

void replace_ineq_with_eq(opt_pk_internal_t *opk, opt_pk_t * op){
	opt_matrix_t * oc = op->C;
	size_t nbcons = oc->nbrows;
	size_t nbeq = op->nbeq;
	unsigned short int nbcolumns = oc->nbcolumns;
	size_t i,k;
	unsigned short int j;
	opt_numint_t ** p = oc->p;
	for(i=0; i < nbcons; i++){
		opt_numint_t * pi = p[i];
		if(!pi[0]){
			continue;
		}
		k = i+1;
		while(k < nbcons){
			opt_numint_t * pk = p[k];
			if(!pk[0]){
				k++;
				continue;
			}
			bool flag = true;
			for(j=1; j < nbcolumns; j++){
				if(pi[j] != (-pk[j])){
					flag = false;
					break;
				}
			}
			if(flag){
				pi[0] = 0;
				nbcons--;
				nbeq++;
				opt_matrix_exch_rows(oc,k,nbcons);
			}
			else{
				k++;
			}
		}
	}
	oc->nbrows = nbcons;
	op->nbeq = nbeq;
}

/************************************
	Remove syntactically equivalent inequalities
***********************************/

void quasi_removal(opt_pk_internal_t *opk, opt_pk_t * o){
	opt_matrix_t *oc = o->C;
	size_t nbcons, nbcolumns;
	nbcons = oc->nbrows;
	nbcolumns = oc->nbcolumns;
	size_t nb = nbcons;
	size_t nbeq = o->nbeq;
        size_t *rmap = (size_t *)calloc(nbcons,sizeof(size_t));
	int i,j;
	//opt_matrix_fprint(stdout,oc);
	for(i = 0; i < nbcons; i++){
		opt_numint_t *pi = oc->p[i];
		for(j = i+1; j < nbcons; j++){
			opt_numint_t* pj = oc->p[j];
			if((pi[0]==pj[0]) && (!opt_vector_compare_coeff(opk,pi,pj,nbcolumns))){
				if(!pi[0]){
					nbeq--;
				}
				//nbcons--;
				if(pi[1] > pj[1]){
                    			rmap[i] = 1;
					//opt_matrix_exch_rows(oc,i,nbcons);
				}
				else{
                   			 rmap[j] = 1;
					//opt_matrix_exch_rows(oc,j,nbcons);
				}	
				//opt_vector_print(pi,nbcolumns);
				//opt_vector_print(pj,nbcolumns);
			}
			else if(!pi[0]){
				opt_numint_t *npj = opt_vector_neg(opk,pj,nbcolumns);
				if(!opt_vector_compare(opk,pi,npj,nbcolumns)){
					if(pi[1]==npj[1]){
						rmap[j] = 1;
						nbeq--;
					}
				}
				free(npj);
			}
		}
	}
    j = nbcons - 1;
    for(i=0; i < nb; i++){
	//printf("i: %d rmap[i]: %d\n",i,rmap[i]);
        if(rmap[i]){
            nbcons--;
            while((j>i) && rmap[j]){
                j--;
            }
            if(j>i){
                opt_matrix_exch_rows(oc,i,j);
	    }
	    j--;
        }
    }
    free(rmap);
   //opt_matrix_fprint(stdout,oc);
   oc->nbrows = nbcons;
   o->nbeq = nbeq;	

}


/* Return the abstract size of a polyhedron, which is the number of
   coefficients of its current representation, possibly redundant. */
size_t opt_pk_size(elina_manager_t* man, opt_pk_array_t* op)
{
  size_t res;
  array_comp_list_t * acl = op->acl;
  if(op->is_bottom || !acl){
	res = 0;
  }
  else{
	unsigned short int num_compa = acl->size;
	opt_pk_t ** poly = op->poly; 
	unsigned short int k;
	for(k=0; k < num_compa; k++){
		opt_pk_t * src = poly[k];
		size_t s1,s2;
		s1 = src->C ? src->C->nbrows : 0;
		s2 = src->F ? src->F->nbrows : 0;
		res =  res + (s1+s2)*(src->intdim + src->realdim);
	} 
  }
  return res;
}


/*******************
	Combine inequalities into one component
********************/
char meet_cons_one_comp(opt_pk_internal_t *opk, opt_pk_t **poly_a, array_comp_list_t *acla,
			unsigned short int **ca_arr, opt_pk_t *poly, unsigned short int * ca,  char * map){
	unsigned short int num_compa = acla->size;
	unsigned short int j,k;
	size_t i;
	comp_list_t * cla = acla->head;
	char is_pos_con = 1;
	size_t  counter = 0;
	for(k = 0; k < num_compa; k++){
		if(map[k]){
			cla = cla->next;
			continue;
		}
		opt_pk_t * src = poly_a[k];
		opt_pk_t * dst = poly;
		unsigned short int * ca_a = ca_arr[k];
		unsigned short int comp_size = cla->size;
		opt_matrix_t * src_mat = src->C;
		opt_matrix_t * dst_mat = dst->C;
		size_t nbconsa = src_mat->nbrows;
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = dst_mat->p;
		size_t i1 = counter;
		bool flag = false;		
		for(i = 0; i < nbconsa; i++){
			opt_numint_t * src_pi = src_p[i];
			opt_numint_t * dst_pi = dst_p[i1];
			if(opt_vector_is_positivity_constraint(opk,src_pi,comp_size+2)){
				//map[ind] = 1;
				flag = true;
				continue;
			}
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
		if(!flag){
			is_pos_con = 0;
		}
		counter = i1;
		dst_mat->nbrows = counter;
		cla = cla->next;
	}
	return is_pos_con;
}



/*******************
	Combine inequalities
********************/
void meet_cons_with_map(opt_pk_internal_t *opk, opt_pk_array_t *oa, opt_pk_t **poly, unsigned short int *rmapa, 
		unsigned short int **ca_arr, size_t *counterC, char *map, char * exclusion_map){
	array_comp_list_t *acla = oa->acl;
	unsigned short int num_compa = acla->size;
	unsigned short int maxcols = oa->maxcols;
	unsigned short int j,k;
	size_t i;
	opt_pk_t ** poly_a = oa->poly;
	comp_list_t * cla = acla->head;
	for(k = 0; k < num_compa; k++){
		if(exclusion_map[k]){
			cla = cla->next;
			continue;
		}
		opt_pk_t * src = poly_a[k];
		unsigned short int ind = rmapa[k];
		opt_pk_t * dst = poly[ind];
		unsigned short int * ca_a = to_sorted_array(cla,maxcols);
		unsigned short int * ca = ca_arr[ind];
		unsigned short int comp_size = cla->size;
		opt_matrix_t * src_mat = src->C;
		opt_matrix_t * dst_mat = dst->C;
		size_t nbconsa = src_mat->nbrows;
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = dst_mat->p;
		size_t * counter = counterC;
		size_t i1 = counter[ind];
		bool flag = false;		
		for(i = 0; i < nbconsa; i++){
			opt_numint_t * src_pi = src_p[i];
			opt_numint_t * dst_pi = dst_p[i1];
			if(opt_vector_is_positivity_constraint(opk,src_pi,comp_size+2)){
				//map[ind] = 1;
				flag = true;
				continue;
			}
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
		if(!flag){
			map[ind] = 0;
		}
		counter[ind] = i1;
		dst_mat->nbrows = counter[ind];
		free(ca_a);
		cla = cla->next;
	}
	
}



void meet_cons(opt_pk_internal_t *opk, opt_pk_array_t *oa, opt_pk_t **poly, unsigned short int *rmapa, 
		unsigned short int **ca_arr, size_t *counterC, char *map){
	array_comp_list_t *acla = oa->acl;
	unsigned short int num_compa = acla->size;
	unsigned short int maxcols = oa->maxcols;
	unsigned short int j,k;
	size_t i;
	opt_pk_t ** poly_a = oa->poly;
	comp_list_t * cla = acla->head;
	for(k = 0; k < num_compa; k++){
		opt_pk_t * src = poly_a[k];
		unsigned short int ind = rmapa[k];
		opt_pk_t * dst = poly[ind];
		unsigned short int * ca_a = to_sorted_array(cla,maxcols);
		unsigned short int * ca = ca_arr[ind];
		unsigned short int comp_size = cla->size;
		opt_matrix_t * src_mat = src->C;
		opt_matrix_t * dst_mat = dst->C;
		size_t nbconsa = src_mat->nbrows;
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = dst_mat->p;
		size_t * counter = counterC;
		size_t i1 = counter[ind];
		bool flag = false;		
		for(i = 0; i < nbconsa; i++){
			opt_numint_t * src_pi = src_p[i];
			opt_numint_t * dst_pi = dst_p[i1];
			if(opt_vector_is_positivity_constraint(opk,src_pi,comp_size+2)){
				//map[ind] = 1;
				flag = true;
				continue;
			}
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
		if(!flag){
			map[ind] = 0;
		}
		counter[ind] = i1;
		dst_mat->nbrows = counter[ind];
		free(ca_a);
		cla = cla->next;
	}
	
}


//combine all vertices into one component
size_t cartesian_product_vertices_one_comp(opt_pk_t **poly_a, array_comp_list_t *acla,
					 unsigned short int ** ca_arr, size_t * num_vertex_a,
					 opt_pk_t * poly, size_t counter, unsigned short int *ca, char * map){

	unsigned short int j,k;
	size_t num_vertex = 0, i;
	unsigned short int num_compa = acla->size;
	comp_list_t *cla = acla->head;
	opt_pk_t * dst = poly;
	for(k=0; k < num_compa; k++){
		if(map[k]){
			cla = cla->next;
			continue;
		}
		opt_pk_t * src = poly_a[k];
		
		
		opt_matrix_t * src_mat = src->F;
		opt_matrix_t * dst_mat = dst->F;
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = dst_mat->p;
		unsigned short int comp_size = cla->size;
		unsigned short int * ca_a = ca_arr[k];
		size_t i1 = counter;
		if(!num_vertex){
			
			for(i=0; i < num_vertex_a[k]; i++){
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
				num_vertex++;
			}
			counter = i1;
			dst_mat->nbrows = counter;
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
					if(i3 >= end){
						num_vertex++;
						opt_vector_copy(dst_pi,vec,dst_mat->nbcolumns);
					}
					dst_pi[1] = lcm;
					unsigned short int l = 0;
					
					
					if(lcm1!=1){
						
						for(j=2; j < dst_mat->nbcolumns; j++){
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
			dst_mat->nbrows = counter;
		}
		cla = cla->next;
		
	}
	
	return counter;
}



void cartesian_product_vertices_with_map(opt_pk_array_t *oa, opt_pk_t ** poly, 
				  unsigned short int *rmapa, unsigned short int ** ca_arr, 
				  size_t * num_vertex_a, size_t *num_vertex, size_t *counterF, char * exclusion_map){
	/*************************
		Cartesian Product of Vertices from A
	************************/
	array_comp_list_t * acla = oa->acl;
	unsigned short int num_compa = acla->size;
	unsigned short int maxcols = oa->maxcols;
	unsigned short int j,k;
	size_t i;
	opt_pk_t **poly_a = oa->poly;
	comp_list_t *cla = acla->head;
	for(k=0; k < num_compa; k++){
		if(exclusion_map[k]){
			cla = cla->next;
			continue;
		}
		if(!num_vertex_a[k]){
			cla = cla->next;
			continue;
		}
		opt_pk_t * src = poly_a[k];
		unsigned short int ind = rmapa[k];
		opt_pk_t * dst = poly[ind];
		unsigned short int * ca_a = to_sorted_array(cla,maxcols);
		unsigned short int * ca = ca_arr[ind];
		unsigned short int comp_size = cla->size;
		opt_matrix_t * src_mat = src->F;
		opt_matrix_t * dst_mat = dst->F;
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = dst_mat->p;
		size_t i1 = counterF[ind];
		if(!num_vertex[ind]){
			for(i=0; i < num_vertex_a[k]; i++){
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
				num_vertex[ind]++;
			}
			counterF[ind] = i1;
			dst_mat->nbrows = i1;
		}
		else{
			size_t start = counterF[ind] - num_vertex[ind];
			size_t i2, i3 = start;
			size_t end = counterF[ind];			
			for(i=0; i < num_vertex_a[k]; i++){
				opt_numint_t * src_pi = src_p[i];
				for(i2=start; i2 < end; i2++){
					opt_numint_t * vec = dst_p[i2];
					opt_numint_t * dst_pi = dst_p[i3];
					opt_numint_t lcm = opt_numint_lcm(src_pi[1],vec[1]);
					opt_numint_t lcm1 = lcm/vec[1];
					opt_numint_t lcm2 = lcm/src_pi[1];
					if(i3 >= counterF[ind]){
						num_vertex[ind]++;
						opt_vector_copy(dst_pi,vec,dst_mat->nbcolumns);
					}
					dst_pi[1] = lcm;
					unsigned short int l = 0;
					
					
					if(lcm1!=1){
						
						for(j=2; j < dst_mat->nbcolumns; j++){
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
			counterF[ind] = i3;
			dst_mat->nbrows = i3;
		}
		free(ca_a);
		cla = cla->next;		
	}
	
}





void cartesian_product_vertices(opt_pk_array_t *oa, opt_pk_t ** poly, 
				  unsigned short int *rmapa, unsigned short int ** ca_arr, 
				  size_t * num_vertex_a, size_t *num_vertex, size_t *counterF){
	/*************************
		Cartesian Product of Vertices from A
	************************/
	array_comp_list_t * acla = oa->acl;
	unsigned short int num_compa = acla->size;
	unsigned short int maxcols = oa->maxcols;
	unsigned short int j,k;
	size_t i;
	opt_pk_t **poly_a = oa->poly;
	comp_list_t *cla = acla->head;
	for(k=0; k < num_compa; k++){
		if(!num_vertex_a[k]){
			cla = cla->next;
			continue;
		}
		opt_pk_t * src = poly_a[k];
		unsigned short int ind = rmapa[k];
		opt_pk_t * dst = poly[ind];
		unsigned short int * ca_a = to_sorted_array(cla,maxcols);
		unsigned short int * ca = ca_arr[ind];
		unsigned short int comp_size = cla->size;
		opt_matrix_t * src_mat = src->F;
		opt_matrix_t * dst_mat = dst->F;
		//size_t nbconsa = src_mat->nbrows;
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = dst_mat->p;
		size_t i1 = counterF[ind];
		
		if(!num_vertex[ind]){
			for(i=0; i < num_vertex_a[k]; i++){
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
				num_vertex[ind]++;
			}
			counterF[ind] = i1;
			dst_mat->nbrows = i1;
		}
		else{
			size_t start = counterF[ind] - num_vertex[ind];
			size_t i2, i3 = start;
			size_t end = counterF[ind];			
			for(i=0; i < num_vertex_a[k]; i++){
				opt_numint_t * src_pi = src_p[i];
				for(i2=start; i2 < end; i2++){
					opt_numint_t * vec = dst_p[i2];
					opt_numint_t * dst_pi = dst_p[i3];
					opt_numint_t lcm = opt_numint_lcm(src_pi[1],vec[1]);
					opt_numint_t lcm1 = lcm/vec[1];
					opt_numint_t lcm2 = lcm/src_pi[1];
					if(i3 >= counterF[ind]){
						num_vertex[ind]++;
						opt_vector_copy(dst_pi,vec,dst_mat->nbcolumns);
					}
					dst_pi[1] = lcm;
					unsigned short int l = 0;
					
					
					if(lcm1!=1){
						
						for(j=2; j < dst_mat->nbcolumns; j++){
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
			counterF[ind] = i3;
			dst_mat->nbrows = i3;
		}
		free(ca_a);
		
		cla = cla->next;		
	}
	
}

//meet rays in one component
void meet_rays_one_comp(opt_pk_t **poly_a, array_comp_list_t *acla, unsigned short int **ca_arr,
	       		 size_t * nblinemap, opt_pk_t *poly, size_t * num_vertex_a,
	        	unsigned short int  * ca, size_t start, char * map){
	/************************
		Consider rays of A
	************************/
	unsigned short int num_compa = acla->size;
	unsigned short int j,k;
	size_t i;
	comp_list_t * cla = acla->head;
	size_t i1 = start;
	opt_pk_t * dst = poly;
	opt_matrix_t * dst_mat = dst->F;
	for(k = 0; k < num_compa; k++){
		if(map[k]){
			cla = cla->next;
			continue;
		}
		opt_pk_t * src = poly_a[k];
		
		unsigned short int * ca_a = ca_arr[k];
		unsigned short int comp_size = cla->size;
		opt_matrix_t * src_mat = src->F;
		size_t nbconsa = nblinemap==NULL? src_mat->nbrows : src_mat->nbrows - nblinemap[k];
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = dst_mat->p;
		
		for(i = num_vertex_a[k]; i < nbconsa; i++){
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
		
		cla = cla->next;
	}
	if(nblinemap){
		cla = acla->head;
		for(k = 0; k < num_compa; k++){
			if(map[k]){
				cla = cla->next;
				continue;
			}
			opt_pk_t * src = poly_a[k];
		
			unsigned short int * ca_a = ca_arr[k];
			unsigned short int comp_size = cla->size;
			opt_matrix_t * src_mat = src->F;
			size_t start = src_mat->nbrows - nblinemap[k];
			size_t nbconsa = src_mat->nbrows;
			opt_numint_t ** src_p = src_mat->p;
			opt_numint_t ** dst_p = dst_mat->p;
		
			for(i = start; i < nbconsa; i++){
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
		
			cla = cla->next;
		}
	}
	dst_mat->nbrows = i1;
}





void meet_rays_with_map(opt_pk_array_t *oa, opt_pk_t **poly, unsigned short int *rmapa, 
	       unsigned short int **ca_arr, size_t * num_vertex_a, size_t * counterF, char * exclusion_map){
	/************************
		Consider rays of A
	************************/
	opt_pk_t **poly_a = oa->poly;
	array_comp_list_t *acla = oa->acl;
	unsigned short int num_compa = acla->size;
	unsigned short int maxcols = oa->maxcols;
	unsigned short int j,k;
	size_t i;
	comp_list_t * cla = acla->head;
	for(k = 0; k < num_compa; k++){
		if(exclusion_map[k]){
			cla = cla->next;
			continue;
		}
		opt_pk_t * src = poly_a[k];
		unsigned short int ind = rmapa[k];
		opt_pk_t * dst = poly[ind];
		unsigned short int * ca_a = to_sorted_array(cla,maxcols);
		unsigned short int * ca = ca_arr[ind];
		unsigned short int comp_size = cla->size;
		opt_matrix_t * src_mat = src->F;
		opt_matrix_t * dst_mat = dst->F;
		size_t nbconsa = src_mat->nbrows;
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = dst_mat->p;
		//size_t * counter = loop==0 ? counterF : counterC;
		size_t * counter = counterF;
		size_t i1 = counter[ind];
				
		for(i = num_vertex_a[k]; i < nbconsa; i++){
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
		counter[ind] = i1;
		dst_mat->nbrows = counter[ind];
		free(ca_a);
		cla = cla->next;
	}

}


void meet_rays(opt_pk_array_t *oa, opt_pk_t **poly, unsigned short int *rmapa, 
	       unsigned short int **ca_arr, size_t * num_vertex_a, size_t * counterF){
	/************************
		Consider rays of A
	************************/
	opt_pk_t **poly_a = oa->poly;
	array_comp_list_t *acla = oa->acl;
	unsigned short int num_compa = acla->size;
	unsigned short int maxcols = oa->maxcols;
	unsigned short int j,k;
	size_t i;
	comp_list_t * cla = acla->head;
	for(k = 0; k < num_compa; k++){
		opt_pk_t * src = poly_a[k];
		unsigned short int ind = rmapa[k];
		opt_pk_t * dst = poly[ind];
		unsigned short int * ca_a = to_sorted_array(cla,maxcols);
		unsigned short int * ca = ca_arr[ind];
		unsigned short int comp_size = cla->size;
		opt_matrix_t * src_mat = src->F;
		opt_matrix_t * dst_mat = dst->F;
		size_t nbconsa = src_mat->nbrows;
		opt_numint_t ** src_p = src_mat->p;
		opt_numint_t ** dst_p = dst_mat->p;
		//size_t * counter = loop==0 ? counterF : counterC;
		size_t * counter = counterF;
		size_t i1 = counter[ind];
				
		for(i = num_vertex_a[k]; i < nbconsa; i++){
			opt_numint_t * src_pi = src_p[i];
			opt_numint_t * dst_pi = dst_p[i1];
			//if(is_vertex(src_pi)){
			//	vertex_map[ind][i1] = 1;
			//}
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
		//}
		counter[ind] = i1;
		dst_mat->nbrows = counter[ind];
		// handle the satF matrix
		//size_t i1 = counterS[ind];
		//size_t j1 = colmapS[ind];
		//size_t nbconsa = poly_a[k]->C->nbrows;
		
		//for(i=0; i < nbconsa; i++){
		//	opt_bitstring_copy(poly[ind]->satF->p[i1]+j1,poly_a[k]->satF->p[i],poly_a[k]->F->nbrows);
		//	i1++;
		//}
		//counterS[ind]=counterS[ind] + nbconsa;
		//colmapS[ind] = colmapS[ind] + poly_a[k]->F->nbrows;
		free(ca_a);
		cla = cla->next;
	}

}

void combine_satmat(opt_pk_internal_t *opk, opt_pk_t *poly, unsigned short int comp_size, size_t end, bool con_to_ray){
	size_t i,j;
	opt_matrix_t *rmat = con_to_ray ? poly->F : poly->C;
	opt_matrix_t *cmat = con_to_ray ? poly->C : poly->F;
	size_t nbrows = rmat->nbrows;
	//opt_satmat_t *src = con_to_ray ? poly_a[k]->satC : poly_a[k]->satF;
	opt_satmat_t *sat = con_to_ray ? poly->satC : poly->satF;
	for(i=0; i < nbrows; i++){
		opt_bitstring_t * dst = sat->p[i];
		opt_numint_t *ri = rmat->p[i];
		opt_bitindex_t ind = opt_bitindex_init(0);
		for(j=0; j < end; j++){
			opt_numint_t *ci = cmat->p[j];
			if(opt_vector_product(opk,ri,ci,comp_size+2)){
				opt_bitstring_set(dst,ind);
			}
			else if(opk->exn){
				opk->exn = ELINA_EXC_NONE;
				opt_bitstring_set(dst,ind);
			}
			opt_bitindex_inc(&ind);
		}
	}
}

void opt_poly_minimize(elina_manager_t* man, opt_pk_t* op){
	if(op->is_minimized){
		return;
	}
        //printf("Minimize IN:\n");
	//fflush(stdout);
	opt_pk_internal_t *opk = (opt_pk_internal_t *)man->internal;
	opt_matrix_t * oc = op->C;
	if(oc){
		quasi_removal(opk,op);
		size_t nbeq = op->nbeq;
		size_t nbcons = oc->nbrows;
		if(nbcons==1){
			return;
		}
		size_t i,j;
		//printf("nbeq: %d\n",nbeq);
		
		opt_matrix_rearrange(oc,nbeq);
		//
		size_t rank = opt_matrix_gauss_elimination(opk,oc,nbeq);
		if(rank < nbeq){
			i = rank;
			//j = nbco;
			while(i<nbeq){
				nbcons --;
				opt_matrix_exch_rows(oc,i,nbcons);
				i++;
				//j++;
			}
			//fprintf(stdout,"nbrows: %d rank: %d nbeq: %d\n", nbcons, rank, nbeq);
			//nbcons = nbcons + rank - nbeq;
			op->nbeq = rank;
			oc->nbrows = nbcons;
		}
		gauss_backsubstitute(opk,oc,rank);
		//remove_redundancy(opk,op);
		op->is_minimized = true;
		
	}
	//printf("Minimize OUT:\n");
        //fflush(stdout);
}

/* Not implemented */
void opt_pk_minimize(elina_manager_t* man, opt_pk_array_t* op)
{
    
}

void opt_poly_obtain_sorted_C(opt_pk_internal_t* opk, opt_pk_t* op)
{
  assert (op->C);

  if (!opt_matrix_is_sorted(op->C)){
    if (op->satF){
      if (op->satC){ opt_satmat_free(op->satC); op->satC = NULL; }
      opt_matrix_sort_rows_with_sat(opk,op->C,op->satF);
    }
    else if (op->satC){
      op->satF = opt_satmat_transpose(op->satC,op->C->nbrows);
      opt_satmat_free(op->satC); op->satC = NULL;
      opt_matrix_sort_rows_with_sat(opk,op->C,op->satF);
    }
    else {
      opt_matrix_sort_rows(opk,op->C);
    }
  }
}



void opt_pk_array_canonicalize(elina_manager_t *man, opt_pk_array_t * o){
	//printf("canonicalize\n");
	//fflush(stdout);
	array_comp_list_t *acl = o->acl;
	if(o->is_bottom || !acl){
		return;
	}
	if(acl->size==0){
		return;
	}
	unsigned short int num_comp = acl->size;
	opt_pk_t ** poly = o->poly;
	unsigned short int k;
	for(k=0; k < num_comp; k++){
	    opt_pk_t * ok = poly[k];
	    opt_poly_chernikova(man,ok,"canonicalize");
	}
}

/* ********************************************************************** */
/* III Printing */
/* ********************************************************************** */
void opt_pk_fprint(FILE* stream, elina_manager_t *man, opt_pk_t* op,
	       char** name_of_dim)
{ 
  opt_pk_internal_t *opk = (opt_pk_internal_t*)man->internal;
  opt_poly_minimize(man,op);
  
  if (!op->C){
    assert(opk->exn == ELINA_EXC_NONE);
    fprintf(stream,"empty polyhedron of dim (%lu,%lu)\n",
	    (unsigned long)op->intdim,(unsigned long)op->realdim);
  }
  else {
    fprintf(stream,"polyhedron of dim (%lu,%lu)\n",
	    (unsigned long)op->intdim,(unsigned long)op->realdim);
    if (opk->exn){
      opk->exn = ELINA_EXC_NONE;
      fprintf(stream,"cannot display due to exception\n");
    }
    else {
      //elina_lincons0_array_t cons = opt_pk_to_lincons_array(man,op);
      //elina_lincons0_array_fprint(stream,&cons,name_of_dim);
      //elina_lincons0_array_clear(&cons);
	opt_matrix_fprint(stdout,op->C);
    }
  }
}


void opt_pk_array_fprint(FILE* stream, elina_manager_t * man, opt_pk_array_t * oa, char ** name_of_dim){
	/*opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_FPRINT);
	array_comp_list_t * acl = oa->acl;
	unsigned short int maxcols = oa->maxcols;
	if(oa->is_bottom || !acl){
		fprintf(stream,"empty polyhedron of dim (%lu)\n",maxcols);
		return ;
	}
	unsigned short int num_comp = acl->size;
	
	comp_list_t * cl = acl->head;
	unsigned short int k;
	opt_pk_t ** poly = oa->poly;
	for(k=0; k < num_comp; k++){
		fprint_comp_list(stream,cl,maxcols);
		opt_pk_t *oak = poly[k];
		opt_pk_fprint(stream,man,oak,name_of_dim); 
		cl = cl->next;
	}*/
        #if defined(TIMING)
		fprintf(stdout,"Times are in CPU Cycles\n");
		fprintf(stdout,"Top: %g\n",top_time);
		fprintf(stdout,"Bottom: %g\n",bottom_time);
		fprintf(stdout,"Free: %g\n",free_time);
		fprintf(stdout,"Copy: %g\n",copy_time);
		fprintf(stdout,"Is_Lequal: %g\n",is_lequal_time);
		fprintf(stdout,"Meet_Abstract: %g\n",meet_time);
		fprintf(stdout,"Join: %g\n",join_time);
		fprintf(stdout,"Widening: %g\n",widening_time);
		fprintf(stdout,"Add_dimension: %g\n",add_dimension_time);
		fprintf(stdout,"remove: %g\n",remove_dimension_time);
		fprintf(stdout,"Permute_dimension: %g\n",permute_dimension_time);
		fprintf(stdout,"Meet_Lincons_Array: %g\n",meet_lincons_time);
		fprintf(stdout,"Forget_Array %g\n",forget_array_time);
		fprintf(stdout,"Poly_to_Box: %g\n",poly_to_box_time);
		fprintf(stdout,"Is_Top: %g\n",is_top_time);
		fprintf(stdout,"Is_Bottom: %g\n",is_bottom_time);
		fprintf(stdout,"Expand: %g\n",expand_time);
		fprintf(stdout,"Fold: %g\n",fold_time);
		fprintf(stdout,"Sat_Lincons: %g\n",sat_lincons_time);
		fprintf(stdout,"Assign Linexpr: %g\n",assign_linexpr_time);
		fprintf(stdout,"Substitute Linexpr: %g\n",substitute_linexpr_time);
		fprintf(stdout,"Bound Dimension: %g\n",bound_dimension_time);
		//fprintf(stdout,"Conversion time: %g\n",opt_conversion_time);
		fprintf(stdout,"Poly is unconstrained: %g\n",poly_is_unconstrained_time);
		//fprintf(stdout,"Join Count: %lld\n",join_count);
		double total_time = top_time + free_time + copy_time + bottom_time + remove_dimension_time + is_lequal_time + meet_time + join_time + widening_time + add_dimension_time 
					+ permute_dimension_time + meet_lincons_time + bound_dimension_time + 
					forget_array_time + poly_to_box_time  + is_top_time + is_bottom_time + expand_time + fold_time + sat_lincons_time + assign_linexpr_time + substitute_linexpr_time+
					poly_is_unconstrained_time;
		fprintf(stdout,"Total OptPoly Analysis: %g\n",total_time);
		fflush(stdout);
	#endif
}



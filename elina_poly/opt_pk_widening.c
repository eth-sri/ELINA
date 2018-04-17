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

#include "opt_pk.h"
#include "opt_pk_config.h"
#include "opt_pk_internal.h"
#include "opt_pk_representation.h"
#include "opt_pk_vector.h"
#include "opt_pk_constructor.h"
#include "opt_pk_test.h"
#include "opt_pk_cherni.h"
#include "opt_mf_qsort.h"

typedef struct opt_satmat_row_t {
  opt_bitstring_t* p; 
  int index;
} opt_satmat_row_t;

static opt_satmat_row_t* opt_esatmat_of_satmat(opt_satmat_t* sat)
{
  size_t i;
  opt_satmat_row_t* tab;

  tab = malloc(sat->nbrows * sizeof(opt_satmat_row_t));
  for (i=0; i<sat->nbrows; i++){
    tab[i].p = sat->p[i];
    tab[i].index = i;
  }
  return tab;
}

/* Row sorting.

We use here the insertion sort. 
The array tab is supposed to be of size sat->nbrows.
*/
static void opt_esatmat_isort_rows(opt_satmat_row_t* tab, opt_satmat_t* sat)
{
  size_t i,j;

  for (i=1; i<sat->nbrows; i++){
    opt_satmat_row_t row = tab[i];
    j = i;
    while (j > 0 && opt_bitstring_cmp(tab[j-1].p, row.p, sat->nbcolumns) > 0){
      tab[j] = tab[j-1];
      j--;
    }
    tab[j] = row;
  }
}

/* Row sorting.

We use here the quick sort. */
typedef struct opt_qsort_man_t {
  size_t size;
} opt_qsort_man_t;

static int opt_qsort_rows_compar(void* qsort_man, const void* p1, const void* p2)
{
  opt_qsort_man_t* qm = (opt_qsort_man_t*)qsort_man;
  return (opt_bitstring_cmp( ((opt_satmat_row_t*)p1)->p,
			 ((opt_satmat_row_t*)p2)->p,
			 qm->size));
}

static void opt_esatmat_sort_rows(opt_satmat_row_t* tab, opt_satmat_t* sat)
{
  if (sat->nbrows>=6){
    opt_qsort_man_t qsort_man;
    qsort_man.size = sat->nbcolumns;
    opt_qsort2(tab,
	  (size_t)sat->nbrows, sizeof(opt_satmat_row_t),
	   opt_qsort_rows_compar, &qsort_man);
  }
  else {
    opt_esatmat_isort_rows(tab,sat);
  }
}

/* Membership test.

The following function tests if the given row belongs to the sorted saturation
matrix. If it is the case, it returns its rank in the saturation
matrix. Otherwise, it returns -1 */

typedef struct opt_bsearch_man_t {
  opt_bitstring_t* satline;
  opt_satmat_row_t* tab;
  size_t size;
} opt_bsearch_man_t;

static bool opt_bsearch2(opt_bsearch_man_t* man, size_t low, size_t high)
{
  if (high - low <= 4){
    size_t i;
    int res=-1;
    for (i=low; i<high; i++){
      int cmp = opt_bitstring_cmp(man->tab[i].p, man->satline, man->size);
      if (cmp==0){
	res=i; break;
      }
      else if (cmp>0) break;
    }
    return res;
  }
  else {
    size_t mid = low+(high-low)/2;
    int cmp = opt_bitstring_cmp(man->tab[mid].p,man->satline,man->size);
    if (cmp<0)
      return (opt_bsearch2(man,mid+1,high));
    else if (cmp>0)
      return (opt_bsearch2(man,low,mid));
    else
      return mid;
  }
}

static
int opt_esatmat_index_in_sorted_rows(opt_bitstring_t* satline, 
				 opt_satmat_row_t* tab, 
				 opt_satmat_t* sat)
{
  opt_bsearch_man_t man;
  man.satline = satline;
  man.tab = tab;
  man.size = sat->nbcolumns;
  return opt_bsearch2(&man,0,sat->nbrows);
}


bool is_vectors_coeff_equal_comp_list(opt_pk_internal_t *opk, opt_numint_t *v1,
                                      opt_numint_t *v2, unsigned short int *ca1,
                                      unsigned short int *ca2, unsigned short int comp_size1,
                                      unsigned short int comp_size2){
    bool res;
    unsigned short int i,j;
    unsigned short int s = opk->dec;
    for(i=0; i < comp_size1; i++){
        unsigned short int num1 = ca1[i];
        j = 0;
        while((j < comp_size2) && (ca2[j]!=num1)){
            j++;
        }
        /******************
         there is a match
         *******************/
        if(j < comp_size2){
            if(v1[i+s] != v2[j+s] ){
                //printf("false\n");
                return false;
            }
        }
        else{
            if(v1[i+s]!=0){
                //printf("false 2 ind: %d num1: %d\n",i+s,num1);
                return false;
            }
        }
    }
    
    return true;
    
}

bool is_vectors_equal_comp_list(opt_pk_internal_t *opk, opt_numint_t * v1, 
				opt_numint_t * v2, unsigned short int * ca1, 
				unsigned short int * ca2, unsigned short int comp_size1, 
				unsigned short int comp_size2){
		
	if(v1[0] != v2[0]){
		return false;
	}
	if(v1[1] != v2[1]){
		return false;
	}
    return is_vectors_coeff_equal_comp_list(opk,v1,v2,ca1,ca2,comp_size1, comp_size2);
}


void vector_copy_comp_list(opt_pk_internal_t *opk, opt_numint_t * dst, opt_numint_t * src, 
			   unsigned short int * ind_map, unsigned short int comp_size){
		unsigned short int i;
		dst[0] = src[0];
		dst[1] = src[1];
		unsigned short int s = opk->dec;
		for(i = 0; i < comp_size; i++){
			unsigned short int j = ind_map[i];
			dst[j+s] = src[i+s];
			j++;
		}
}


/*********************************
		Required for widening
**********************************/
void opt_build_satline(opt_pk_internal_t *opk, opt_matrix_t *mat, 
		       opt_numint_t *cons, opt_bitstring_t * satline){
	opt_bitindex_t ind = opt_bitindex_init(0);
	unsigned short int nbcolumns = mat->nbcolumns;
	while(ind.index < mat->nbrows){
		opt_numint_t prod = opt_vector_product(opk,mat->p[ind.index],cons,nbcolumns);
		if(prod){
			opt_bitstring_set(satline,ind);
		}
		opt_bitindex_inc(&ind);
	}
}


comp_list_t * vector_to_comp_list(opt_pk_internal_t *opk, opt_numint_t *v, unsigned short int * ca, unsigned short int comp_size){
	unsigned short int j;
	comp_list_t *cl = create_comp_list();
	for(j=0; j < comp_size; j++){
		if(v[j+2]){
			insert_comp(cl,ca[j]);
		}
	}
	return cl;
}

/*************************
		Standard Polyhedra widening with generators
**************************/

void opt_poly_widening_gen(elina_manager_t *man, opt_pk_array_t *op, opt_pk_array_t *oa, opt_pk_array_t *ob){
	opt_pk_internal_t * opk = opt_pk_init_from_manager(man,ELINA_FUNID_WIDENING);
	unsigned short int maxcols = oa->maxcols;
	array_comp_list_t * acla = oa->acl;
	if(oa->is_bottom || !acla){
		op = opt_pk_copy(man,ob);
		return;
	}
	unsigned short int num_compa = acla->size;
	if(num_compa==0){
		opt_poly_set_top(opk,op);
		return;
	}
	
	array_comp_list_t * aclb = ob->acl;
	if(ob->is_bottom || !aclb){
		op = opt_pk_copy(man,oa);
		return;
	}
	unsigned short int num_compb = aclb->size;
	if(num_compb==0){
		opt_poly_set_top(opk,op);
		return;
	}
	unsigned short int k, ka, kb;
	size_t i;
	/**************************
		Minimize A
	***************************/
	opt_pk_t ** poly_a = oa->poly;
	for(ka=0; ka < num_compa; ka++){
		opt_pk_t *oak = poly_a[ka];
		opt_poly_chernikova(man,oak,"cons to gen");
		if(opk->exn){
			opk->exn = ELINA_EXC_NONE;
			man->result.flag_best = man->result.flag_exact = false;
			opt_poly_set_top(opk,op);
			return;
		}
		if(!oak->C && !oak->F){
			op = opt_pk_copy(man,ob);
			return;
		}
	}

	/**************************
		Minimize B
	***************************/
	opt_pk_t ** poly_b = ob->poly;
	for(kb=0; kb < num_compb; kb++){
		opt_pk_t * obk = poly_b[kb];
		opt_poly_chernikova(man,obk,"cons to gen");
		if(opk->exn){
			opk->exn = ELINA_EXC_NONE;
			man->result.flag_best = man->result.flag_exact = false;
			opt_poly_set_top(opk,op);
			return;
		}
	}

	/******************************
			Compute union of independent components
	*******************************/
	
	array_comp_list_t * acl = union_array_comp_list(acla,aclb,maxcols);
	unsigned short int num_comp = acl->size;
	
	
	/*******************************
			Factor A according to union
	********************************/
	unsigned short int * rmapa = (unsigned short int *)calloc(num_compa,sizeof(unsigned short int));
	size_t * nbgena = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * nbmapa = (size_t *)calloc(num_comp,sizeof(size_t));
	size_t * num_vertex_a = (size_t *)calloc(num_compa,sizeof(size_t));
	size_t * num_vertex = (size_t *)calloc(num_comp,sizeof(size_t));
	comp_list_t *cla = acla->head;
	for(ka=0; ka < num_compa; ka++){
		opt_pk_t * oak = poly_a[ka];
		
		unsigned short int inda = is_comp_list_included(acl,cla,maxcols);
		rmapa[ka] = inda;
		nbgena[inda] = nbgena[inda] + oak->F->nbrows;
		nbmapa[inda] = nbmapa[inda] + oak->C->nbrows;
		opt_poly_obtain_satF(oak);
		num_vertex_a[ka] = opt_generator_rearrange(oak->F,oak->satF);
		if(num_vertex_a[ka]){
			if(!num_vertex[inda]){
				num_vertex[inda] = num_vertex_a[ka];
			}
			else{
				num_vertex[inda] = num_vertex[inda] * num_vertex_a[ka];
			}
		}
		
		cla = cla->next;
	}	

	/********************************
			Factor B according to union
	*********************************/
	unsigned short int * rmapb = (unsigned short int *)calloc(num_compb, sizeof(unsigned short int));
	size_t * nbmapb = (size_t *)calloc(num_comp,sizeof(size_t));
	comp_list_t *clb = aclb->head;
	for(kb=0; kb < num_compb; kb++){
		opt_pk_t * obk = poly_b[kb];
		
		unsigned short int indb = is_comp_list_included(acl,clb,maxcols);
		rmapb[kb] = indb;
		nbmapb[indb] = nbmapb[indb] + obk->C->nbrows;
		//printf("B\n");
		//opt_matrix_fprint(stdout,obk->C);
		//fflush(stdout);
		clb = clb->next;
	}
	
	/*************************************
			Initialize data structures for the result
	**************************************/
	unsigned short int ** ca_arr = (unsigned short int **)malloc(num_comp*sizeof(unsigned short int *));
	unsigned short int * comp_size_arr = (unsigned short int *)calloc(num_comp, sizeof(unsigned short int));
	opt_pk_t ** poly = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));	
	size_t * counter = (size_t *)calloc(num_comp, sizeof(size_t));
	size_t * counterC = (size_t *)calloc(num_comp, sizeof(size_t));
	size_t * counterF = (size_t *)calloc(num_comp, sizeof(size_t));
	//size_t *rowS = (size_t *)calloc(num_comp,sizeof(size_t));
	//opt_satmat_t ** satF_arr = (opt_satmat_t **)malloc(num_comp*sizeof(opt_satmat_t*));
	comp_list_t *cl = acl->head;
	opt_pk_t ** tmp = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t*));
	
	for(k=0; k < num_comp; k++){
		unsigned short int comp_size = cl->size;
		ca_arr[k] = to_sorted_array(cl,maxcols);
		poly[k] = opt_poly_alloc(comp_size,0);
		tmp[k] = opt_poly_alloc(comp_size,0);
		tmp[k]->C = opt_matrix_alloc(nbmapa[k]+1,comp_size+2,false);
		tmp[k]->F = opt_matrix_alloc(nbgena[k]+2*num_vertex[k],comp_size+2,false);
		poly[k]->C = opt_matrix_alloc(opk->dec - 1 + nbmapb[k],comp_size+2,false);
		opt_matrix_fill_constraint_top(opk,poly[k]->C,0);
		//satF_arr[k] = opt_satmat_alloc(nbmapa[k],opt_bitindex_size(nbgena[k]));
		comp_size_arr[k] = comp_size;
		counter[k] = opk->dec - 1;
		num_vertex[k] = 0;
		cl = cl->next;
	}	
	
	
	
	//unsigned short int ** ind_map_a_arr = (unsigned short int **)malloc(num_compa*sizeof(unsigned short int *));
	//cla = acla->head;
	
	//for(ka=0; ka < num_compa; ka++){
		//unsigned short int inda = rmapa[ka];
		//unsigned short int * ca = to_sorted_array(cla,maxcols);
		//ind_map_a_arr[ka] = map_index(ca,ca_arr[inda],cla->size);
		//free(ca);
		//size_t nka;
		//opt_satmat_t * satF = satF_arr[inda];	
		//opt_matrix_t * C = poly_a[ka]->C;
		//size_t counter = rowS[inda];
		//for(i=0; i < C->nbrows; i++){
		//	size_t colS = 0;
		//	opt_numint_t * pci = C->p[i];
		//	if(opt_vector_is_positivity_constraint(opk, pci,poly_a[ka]->C->nbcolumns)){
		//		continue;
		//	}
		//	for(nka=0; nka < num_compa; nka++){
		//		unsigned short int n_inda = rmapa[nka];
		//		if(n_inda==inda){
		//			opt_matrix_t * F = poly_a[nka]->F;
					
		//			if(ka==nka){
		//				opt_bitstring_move(satF->p[counter],poly_a[ka]->satF->p[i],F->nbrows,colS);
		//			}
		//			colS =  colS + F->nbrows;
		//		}
		//	}
		//	counter++;
		//}
		//rowS[inda] = counter;
		//cla = cla->next;
	//}
	/************************
		Now Consider rays of A
	************************/
	meet_rays(oa,tmp,rmapa,ca_arr,num_vertex_a,counterF);

	/*************************
		Cartesian Product of Vertices from  A
	************************/
	
	cartesian_product_vertices(oa,tmp,rmapa,ca_arr,num_vertex_a,num_vertex,counterF);
	
	/***********************
		Meet constraints of A
	************************/
	char * pos_con_map = (char *)calloc(num_comp, sizeof(char)); 
	//for(k=0; k < num_comp; k++){
	//	pos_con_map[k] = 1;
	//}
	meet_cons(opk,oa,tmp,rmapa,ca_arr,counterC,pos_con_map);
		
	/*************************
		Add positivity constraint of A
	**************************/
        
	//for(k=0; k < num_comp; k++){
		//size_t count = counterC[k];
		//if(pos_con_map[k]){
		//	poly[k]->C->p[count][0] = 1;
		//	poly[k]->C->p[count][1] = 1;
		//	counterC[k]++;
		//}
		//poly[k]->C->nbrows = counterC[k];
	//}	

	char * exc_map = (char *)calloc(num_comp,sizeof(char));
	for(k=0; k < num_comp; k++){
		opt_matrix_t * C = tmp[k]->C;
		opt_matrix_t * F = tmp[k]->F;
		tmp[k]->satF = opt_satmat_alloc(C->nbrows,opt_bitindex_size(F->nbrows));
		combine_satmat(opk,tmp[k],comp_size_arr[k],F->nbrows,false);
	}

	clb = aclb->head;
	for(kb=0; kb < num_compb; kb++){
		opt_pk_t * obk = poly_b[kb];
		opt_matrix_t * C = obk->C;
		
		size_t nbcons = C->nbrows;
		opt_numint_t ** pc = C->p;
		unsigned short int indb = rmapb[kb];
		
		unsigned short int * ca = to_sorted_array(clb,maxcols);
		unsigned short int * ind_map_b = map_index(ca,ca_arr[indb], clb->size); 
		opt_satmat_t * satF = tmp[indb]->satF;
		/********************
				Sort rows of saturation matrix
		*********************/
		opt_satmat_row_t * satrow = opt_esatmat_of_satmat(satF);
		opt_esatmat_sort_rows(satrow,satF);		
		unsigned short int nbcols_sat = satF->nbcolumns;
		opt_bitstring_t * satline = opt_bitstring_alloc(nbcols_sat);
		opt_matrix_t * F = tmp[indb]->F;
		for(i=0; i < nbcons; i++){
			opt_numint_t * pci = pc[i];
			//Fix Me
			//comp_list_t * res = vector_to_comp_list(opk,pci,ca,clb->size);
			opt_numint_t * npci = opt_map_vector(pci,ind_map_b, comp_size_arr[indb],obk->C->nbcolumns);
			//size_t start = 0;
			opt_bitstring_clear(satline,nbcols_sat);
			int index;
			//cla = acla->head;
			//for(ka=0; ka < num_compa; ka++){
				//unsigned short int inda = rmapa[ka];
				//if(inda==indb){
					//unsigned short int * ind_map_a = ind_map_a_arr[ka];
					
					//if(!is_disjoint(res,cla,maxcols)){
						opt_build_satline(opk,F,npci,satline);
						if(opk->exn){
							opk->exn = ELINA_EXC_NONE;
							exc_map[indb] = 1;
						}
					//}
					//start = start + F->nbrows;
					
				//}
				//cla = cla->next;
			//}
			
			index = opt_esatmat_index_in_sorted_rows(satline,satrow,satF);
			if(index>=0){
				index = satrow[index].index;
				//if(opk->funopt->algorithm > 0 || !opt_vector_is_positivity_constraint(opk, oak->C->p[index],oak->C->nbcolumns)){
					size_t i1 = counter[indb];
					opt_vector_copy(poly[indb]->C->p[i1],npci,comp_size_arr[indb]+2);
					counter[indb]++;
				//}
			}
			free(npci);
			//free(res);
		}
		
		free(ca);
		free(ind_map_b);
		clb = clb->next;
	}


	cl = acl->head;
	
	for(k=0; k < num_comp; k++){
		free(ca_arr[k]);
		opt_poly_clear(tmp[k]);
	}
	unsigned short int k1=0;
	unsigned short int bound = num_comp;
	for(k=0; k < num_comp; k++){
		opt_pk_t *oak = poly[k1];
		if(exc_map[k]){
			comp_list_t * tmp = cl;
			cl = cl->next;
			remove_comp_list(acl,tmp);
			unsigned short int k2;
			for(k2=k1; k2 < bound - 1; k2++){
				poly[k2] = poly[k2+1];
			}
			opt_poly_clear(oak);
			bound--;
		}
		else{
			poly[k1]->C->nbrows = counter[k1];
			opt_poly_chernikova(man,poly[k1],"widening result");
			k1++;
			cl=cl->next;
		}
		
	}
	
        op->acl = acl;
	op->poly = poly;
	free(rmapa);
	free(rmapb);
	free(exc_map);
	free(nbmapb);
	free(nbgena);
	free(counter);
	free(counterC);
	free(counterF);
	free(tmp);
	free(nbmapa);
	free(comp_size_arr);
	//free(ind_map_a_arr);
	free(ca_arr);
	free(pos_con_map);
}

opt_pk_array_t* opt_pk_widening(elina_manager_t* man, opt_pk_array_t* oa, opt_pk_array_t* ob){
	#if defined(TIMING)
 	    start_timing();
   	#endif 
	//printf("Widening INPUT\n");
	//elina_lincons0_array_t arr1 = opt_pk_to_lincons_array(man,oa);
	//elina_lincons0_array_fprint(stdout,&arr1,NULL);
	//elina_lincons0_array_t arr2 = opt_pk_to_lincons_array(man,ob);
	//elina_lincons0_array_fprint(stdout,&arr2,NULL);
	//elina_lincons0_array_clear(&arr1);
	//elina_lincons0_array_clear(&arr2);
	//fflush(stdout);
	opt_pk_array_t *op;
	op = opt_pk_array_alloc(NULL,NULL,oa->maxcols);
	opt_poly_widening_gen(man,op,oa,ob);
	#if defined(TIMING)
 	    record_timing(widening_time);
   	#endif 
	//printf("Widening OUTPUT\n");
	//elina_lincons0_array_t arr3 = opt_pk_to_lincons_array(man,op);
	//elina_lincons0_array_fprint(stdout,&arr3,NULL);
	//elina_lincons0_array_clear(&arr3);
	//fflush(stdout);
	return op;
}

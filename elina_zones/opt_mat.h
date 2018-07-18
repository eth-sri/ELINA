/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
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

#ifndef __OPT_MAT_H__
#define __OPT_MAT_H__

#ifdef __cplusplus
extern "C" {
#endif

#if defined(TIMING)
	#define start_timing()				\
		tsc_counter start, end;				\
		double cycles;					\
		CPUID();					\
		RDTSC(start)						

  	#define record_timing(counter)		\
		RDTSC(end);				\
  		CPUID();				\
  		cycles = (double)(COUNTER_DIFF(end, start));	\
  		counter += cycles	
	extern double zones_closure_time;
	extern double zones_copy_time;
	extern double zones_is_equal_time;
	extern double zones_is_lequal_time;
	extern double zones_permute_dimension_time;
	extern double zones_top_time;
	extern double zones_meet_time;
	extern double zones_join_time;
	extern double zones_add_dimension_time;
	extern double zones_widening_time;
	extern double zones_free_time;
	extern double zones_forget_array_time;
	extern double zones_meet_lincons_time;
	extern double zones_to_box_time;
	extern double zones_alloc_time;
	extern double zones_is_top_time;
	extern double zones_expand_time;
	extern double zones_fold_time;
	extern double zones_sat_lincons_time;
	extern double zones_assign_linexpr_time;
        extern double zones_narrowing_time;
	extern double zones_is_unconstrained_time;
#endif


#include "opt_zones_internal.h"
#include "opt_zones_closure.h"

/****************** Basic Opeartors **********************/
opt_zones_mat_t* opt_zones_mat_alloc_top(unsigned short int dim);
opt_zones_mat_t* opt_zones_mat_alloc(unsigned short int dim);
opt_zones_mat_t* opt_zones_mat_copy(opt_zones_mat_t * oz, unsigned short int dim);
bool opt_zones_mat_closure(opt_zones_mat_t * oz, unsigned short int dim);
void opt_zones_mat_free(opt_zones_mat_t * oz);

/********************* Predciate operators *****************/
bool is_top_zones_mat(opt_zones_mat_t *oz, unsigned short int dim);
bool is_lequal_zones_mat(opt_zones_mat_t *oz1, opt_zones_mat_t *oz2, unsigned short int dim);
bool is_equal_zones_mat(opt_zones_mat_t *oz1, opt_zones_mat_t *oz2, unsigned short int dim);


/********************** Binary Operators   *******************/
void meet_zones_mat(opt_zones_mat_t *oz, opt_zones_mat_t *oz1, opt_zones_mat_t *oz2, unsigned short int dim, bool destructive);
void sparse_join_zones_mat(opt_zones_mat_t *oo, opt_zones_mat_t *oo1, opt_zones_mat_t *oo2, unsigned short int dim, bool destructive);
void join_zones_mat(opt_zones_mat_t *oz, opt_zones_mat_t *oz1, opt_zones_mat_t *oz2, unsigned short int dim, bool destructive);
void sparse_widening_zones_mat(opt_zones_mat_t *oz, opt_zones_mat_t *oz1, opt_zones_mat_t *oz2, unsigned short int dim);
void widening_zones_mat(opt_zones_mat_t *oz, opt_zones_mat_t *oz1, opt_zones_mat_t *oz2, unsigned short int dim);

/********************** Resize operators  *********************/
void forget_array_zones_mat(opt_zones_mat_t *oz, elina_dim_t * tdim, unsigned short int dim, unsigned short int size, bool project);
void opt_zones_mat_addrem_dimensions(opt_zones_mat_t * dz, opt_zones_mat_t *sz, elina_dim_t *tdim, 
				     unsigned short int size, unsigned short int mult, unsigned short int dim, bool add);
void opt_zones_mat_permute(opt_zones_mat_t * dst, opt_zones_mat_t *src, unsigned short int dst_dim, unsigned short int src_dim, elina_dim_t* permutation);

/********************** Transfer function ********************/
bool opt_zones_mat_add_lincons(opt_zones_internal_t * pr,opt_zones_mat_t *oz, unsigned short int intdim, 
			       unsigned short int dim, elina_lincons0_array_t *array, bool* exact, bool* respect_closure);

static inline bool check_trivial_relation_zones(double *m, unsigned short int i, unsigned short int j, unsigned short int dim){
	unsigned short int n = dim + 1;
	int ind1 = n*i+j;
	int ind2 = n*j+i;
	if(i==j){
		if(m[ind1] != 0){
			return false;
		}
	}
	else{
		if(m[ind1] != INFINITY){
			return false;
		}
		if(m[ind2] != INFINITY){
			return false;
		}
	}
	return true;
}

static inline void ini_unary_relation_zones(double *m, unsigned short int i, unsigned short int dim){
	if(i>=dim){
		return;
	}
	unsigned short int n = dim + 1;
	unsigned short int i1 = i;
	m[n*i1+i1] = 0;
	m[n*i1] = INFINITY;
	m[i1] = INFINITY;
}

//assumes i!=j
static inline void ini_binary_relation_zones(double *m, unsigned short int i, unsigned short int j, unsigned short int dim){
	if((i>dim) || (j > dim)){
		return;
	}
	
	unsigned short int n = dim + 1;
	unsigned short int i1 = i;
	unsigned short int j1 = j;
	int ind1 = n*i1+j1;
	int ind2 = n*j1+i1;
	m[ind1] = INFINITY;
	m[ind2] = INFINITY;
}

static inline void ini_relation_zones(double *m, unsigned short int i, unsigned short int j, unsigned short int dim){
	
	if((i>dim) || (j > dim)){
		return;
	}
	if(i==j){
		ini_unary_relation_zones(m,i,dim);
	}
	else{
		ini_binary_relation_zones(m,i,j,dim);
	}
}

static inline void ini_comp_relations_zones(double * result, comp_list_t * cl1, comp_list_t *cl2,unsigned short int dim){
	comp_t *c1 = cl1->head;
	while(c1!=NULL){
		comp_t *c2 = cl2->head;
		unsigned short int i = c1->num;
		while(c2!=NULL){
			unsigned short int j = c2->num;
			if(i!=j){
				
				ini_binary_relation_zones(result,i+1,j+1,dim);
			}
			c2 = c2->next;
		}
		c1 = c1->next;
	}
}

static inline void ini_comp_elem_relation_zones(double * m, comp_list_t * cl1, unsigned short int j, unsigned short int dim){
	comp_t *c1 = cl1->head;
	while(c1!=NULL){
		int i = c1->num;
		if(i!=j){
			ini_binary_relation_zones(m,i+1,j+1,dim);
		}
		c1 = c1->next;
	}
}

static inline void print_mat(opt_zones_mat_t *oz, unsigned short int dim){
	if(!oz){	
		printf("null matrix\n");
		return;
	}
	unsigned short int n = dim+1;
	unsigned short int i,j;
	double *m = oz->mat;
	if(!oz->is_dense){
		array_comp_list_t *acl = oz->acl;
		comp_list_t *cl = acl->head;
		while(cl!=NULL){
			print_comp_list(cl,n);
			unsigned short int comp_size = cl->size;	
			unsigned short int *ca = to_sorted_array(cl,n);
			for(i=0; i < comp_size; i++){
				unsigned short int i1 = ca[i]+1;
				for(j=0; j < comp_size; j++){
					unsigned short int j1 = ca[j]+1;
					printf("%g ",m[n*i1+j1]);
				}
				printf("\n");
			}
			for(i=0; i < comp_size; i++){
				unsigned short int i1 = ca[i]+1;
				printf("[%g, %g]\n",m[n*i1],m[i1]);
			}
			free(ca);
			cl = cl->next;
		}
		
	}
	else{
		printf("dim: %d\n",dim);
		for(i=0;i<n; i++){
			for(j=0; j < n; j++){
				printf("%g ",m[n*i+j]);
			}
			printf("\n");
		}
	}
}

void texpr0_to_comp_list_zones(comp_list_t *res,elina_texpr0_t * expr);
void texpr0_node_to_comp_zones(comp_list_t *res, elina_texpr0_node_t* node);

#ifdef __cplusplus
}
#endif

#endif

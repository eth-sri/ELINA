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

#include "opt_mat.h"
#include "opt_zones_incr_closure.h"

#if defined(TIMING)
	double zones_closure_time = 0;
	double zones_copy_time = 0;
	double zones_is_equal_time = 0;
	double zones_is_lequal_time = 0;
	double zones_permute_dimension_time = 0;
	double zones_top_time = 0;
	double zones_meet_time = 0;
	double zones_join_time = 0;
	double zones_add_dimension_time = 0;
	double zones_widening_time = 0;
	double zones_free_time = 0;
	double zones_forget_array_time = 0;
	double zones_meet_lincons_time = 0;
	double zones_to_box_time = 0;
	double zones_alloc_time = 0;
	double zones_is_top_time = 0;
	double zones_expand_time = 0;
	double zones_fold_time = 0;
	double zones_sat_lincons_time = 0;
	double zones_assign_linexpr_time = 0;
        double zones_narrowing_time = 0;
	double zones_is_unconstrained_time = 0;
#endif

void texpr0_node_to_comp_zones(comp_list_t *res, elina_texpr0_node_t* node){
	texpr0_to_comp_list_zones(res,node->exprA);
	texpr0_to_comp_list_zones(res,node->exprB);
}

void texpr0_to_comp_list_zones(comp_list_t *res,elina_texpr0_t * expr){
  switch (expr->discr) {
  case ELINA_TEXPR_CST:
    break;
  case ELINA_TEXPR_DIM:
    insert_comp(res,expr->val.dim);
    break;
  case ELINA_TEXPR_NODE:
    texpr0_node_to_comp_zones(res,expr->val.node);
    break;
  default:
    assert(false);
  }
}

/********* extract independent components ***************/
array_comp_list_t * extract_comps_zones(double *m, unsigned short int dim){
	comp_list_t **comp = (comp_list_t **)calloc(dim,sizeof(comp_list_t *));
	char *map = (char *)calloc(dim*dim,sizeof(char));
	unsigned short int * owner = (unsigned short int *)calloc(dim,sizeof(unsigned short int));
	for(unsigned short int i = 0; i < dim; i++){
		comp[i]= create_comp_list();
		owner[i] = i;
	}
	unsigned short int n = dim + 1;
	for(unsigned short int i =1; i < dim; i++){
		for(unsigned short int j = 1; j <=i; j++){
			bool flag = false;
			int ind1 = n*i+j;
			int ind2 = n*j+i;
			if(i != j){
				if((m[ind1]!=INFINITY) || (m[ind2]!=INFINITY)){
					flag = true;
				}
			}
			else{
				if((m[i]!=INFINITY) || (m[n*i]!=INFINITY)){
					flag = true;
				}
			}
			if(flag){
				unsigned short int oi,oj;
			
				////printf("%d\t%d\n",i,j);
				oi = owner[i-1];
				//int s1 = comp[oi].size;
				if(!map[dim*oi + i-1]){
					insert_comp(comp[oi],i-1);
					map[dim*oi + i-1] = 1;
				}				
				
				oj  = owner[j-1];
				if(!map[dim*oj + j-1]){
					insert_comp(comp[oj],j-1);
					map[dim*oj + j-1] = 1;
				}
				unite_comp_lists(comp[oi],comp[oj],map,oi,oj,dim);
				comp_t *c = comp[oi]->head;
				for(int k = 0; k < comp[oi]->size; k++){
					unsigned short int num = c->num;
					if(num!=oi ){
						owner[num] = owner[oi];
						//free((blk + k)->ind);
					}
					c = c->next;
				}
			}
		}
	}
	free(map);
	int size = 0;
	//double * map = (double *)calloc(n,sizeof(double));
	array_comp_list_t * res = create_array_comp_list();
	for(unsigned short int i = 0; i < dim; i++){
		if((owner[i]==i) && (comp[i]->size > 0)){
			comp_list_t * dst = copy_comp_list(comp[i]);
			insert_comp_list(res,dst);
		}
	}
	for(unsigned short int i = 0; i < dim; i++){
		free_comp_list(comp[i]);
	}
	free(comp);
	free(owner);
	//printf("Extracted:\n");
	return res;
}

void set_mat_top(double *m, unsigned short int n){
	int size = n*n;
	int i;
	#if defined(VECTOR)
		v_double_type inf = v_set1_double(INFINITY);
		for(i=0;i < size/v_length; i++){
			v_store_double(m+i*4,inf);
		}
		for(i=(size/v_length)*v_length; i < size; i++ ){
			m[i] = INFINITY;
		}
	#else
		for(i=0; i < size; i++){
			m[i] = INFINITY;
		}
	#endif	
	for(i=0; i < n; i++){
		m[n*i+i] = 0;
	}
}

opt_zones_mat_t* opt_zones_mat_alloc_top(unsigned short int dim){
	#if defined(TIMING)
		start_timing();
    	#endif
	double * m;
	unsigned short int n = dim+1;
	int size = n*n;
	m = (double *)malloc(size*sizeof(double));
	//set_mat_top(m,n);
	m[0] = 0;
	array_comp_list_t * acl = create_array_comp_list();
	opt_zones_mat_t * res = (opt_zones_mat_t *)malloc(sizeof(opt_zones_mat_t));
	
	res->mat = m;	
	res->acl = acl;
	res->nni = 2*n;
	res->is_top = true;
	res->is_dense = false;
	res->ti = false;
	res->ind = -1;
	
        #if defined(TIMING)
		record_timing(zones_top_time);
	#endif
	return res;
}

opt_zones_mat_t* opt_zones_mat_alloc(unsigned short int dim){
	#if defined(TIMING)
		start_timing();
	#endif
	
	double *m;
	unsigned short int n = dim + 1;
	int size = n*n;
	m = (double *)malloc(size*sizeof(double));
	opt_zones_mat_t * res = (opt_zones_mat_t *)malloc(sizeof(opt_zones_mat_t));
	res->acl = create_array_comp_list();
	res->mat = m;
	res->ti = false;
	res->is_top = false;
	res->ind = -1;
	
	#if defined(TIMING)
		record_timing(zones_alloc_time);
	#endif
	return res;
}


void opt_zones_comp_copy(double * dest, double * src, comp_list_t * cl, unsigned short int dim){
	unsigned short int n = dim+1;
	unsigned short int * ca = to_sorted_array(cl,n);
	unsigned short int i,j;
	unsigned short int comp_size = cl->size;
	for(i=0; i < comp_size; i++){
		unsigned short int i1 = ca[i]+1;
		for(j=0; j < i;j++){
			unsigned short int j1 = ca[j]+1;
			dest[n*i1+j1] = src[n*i1+j1];
			dest[n*j1+i1] = src[n*j1+i1];
		}
		dest[n*i1] = src[n*i1];
		dest[i1] = src[i1];
		dest[n*i1+i1] = 0;
	}
			
	free(ca);
	
}

void opt_zones_sparse_copy(opt_zones_mat_t * dst_mat, opt_zones_mat_t *src_mat, unsigned short int dim){
	dst_mat->ti = false;
	comp_list_t * cl = src_mat->acl->head;
	int size = (dim+1)*(dim+1);
	double * dest = (double *)malloc(size*sizeof(double));
	double * src = src_mat->mat;
	dest[0] = src[0];
	while(cl!=NULL){
		opt_zones_comp_copy(dest,src,cl,dim);
		cl = cl->next;
	}
	dst_mat->acl = copy_array_comp_list(src_mat->acl);
	dst_mat->nni = src_mat->nni;
  	dst_mat->is_dense = false;
	dst_mat->is_top = src_mat->is_top;
	dst_mat->ti = false;
	dst_mat->mat = dest;
}

opt_zones_mat_t* opt_zones_mat_copy(opt_zones_mat_t * oz, unsigned short int dim){
	if(!oz){
		return NULL;
	}
	
	#if defined(TIMING)
		start_timing();
	#endif
	
	unsigned short int n = dim + 1;
	int size = n*n;
	opt_zones_mat_t * res = (opt_zones_mat_t *)malloc(sizeof(opt_zones_mat_t));
	double * m = (double *)malloc(size*sizeof(double));
	//set_mat_top(m,n);
	if(!oz->is_dense){
		//printf("sparse\n");
		//fflush(stdout);
		array_comp_list_t * acl = oz->acl;
		comp_list_t *cl = acl->head;
		m[0] = 0;
		while(cl!=NULL){
			opt_zones_comp_copy(m,oz->mat,cl,dim);
			cl = cl->next;
		}
		res->ti = false;
		res->acl = copy_array_comp_list(acl);
	}
	else{
		//printf("dense\n");
		//fflush(stdout);
		int i;
		#if defined (VECTOR)
			for(i=0; i < size/v_length; i++){
				v_double_type tmp = v_load_double(oz->mat + i*v_length);
				v_store_double(m+i*v_length,tmp);
			}
			for(i=(size/v_length)*v_length; i < size; i++){
				m[i] = oz->mat[i];
			}
		#else
			for(i=0; i < size; i++){
				m[i] = oz->mat[i];
			}
		#endif
		res->ti = true;
		res->acl = NULL;
	}
	res->mat = m;
	res->nni = oz->nni;
	
	res->is_dense = oz->is_dense;
	res->is_top = oz->is_top;
	res->ind = oz->ind;
	#if defined(TIMING)
		record_timing(zones_copy_time);
	#endif
	
	return res;
}

double recalculate_zones_sparsity(opt_zones_mat_t *oz, unsigned short int dim){
	unsigned short int n = dim + 1;
	int size = n*n;
	int count = 0;
	double *m = oz->mat;
	if(!oz->is_dense){
		array_comp_list_t * acl = oz->acl;	
		comp_list_t *cl = acl->head;
		while(cl!=NULL){
			unsigned short int * ca = to_sorted_array(cl,n);
			unsigned short int i,j;
			for(i=0; i < dim; i++){
				unsigned short int i1 = ca[i]+1;
				for(j=0; j < dim; j++){
					unsigned short int j1 = ca[j]+1;
					if(oz->mat[n*i1+j1]!=INFINITY){
						count++;
					}
				}
				if(oz->mat[i1]!=INFINITY){
					count++;
				}
				if(oz->mat[n*i1]!=INFINITY){
					count++;
				}
			}
			free(ca);
			cl = cl->next;
		}
	}
	else{
		int i;
		for(i=0; i < size; i++){
			if(m[i]!=INFINITY){
				count++;
			}
		}
	}
	oz->nni = count;
	return 1- ((double)(oz->nni/(double)size));
}

/****
	Fully initialize the decomposed type DBM
*****/
void convert_to_dense_zones(opt_zones_mat_t * oz, unsigned short int dim, bool flag){
	double *src = oz->mat;
	src[0] = 0;
	
	for(unsigned short int i = 0; i < dim; i++){
		comp_list_t * li = find(oz->acl,i);
		if(li==NULL){
			ini_unary_relation_zones(src,i+1,dim);
		}
	}
	for(unsigned short int i = 0; i < dim; i++){
		comp_list_t * li = find(oz->acl,i);
		for(unsigned short int j = 0; j <i; j++){
			comp_list_t * lj = find(oz->acl,j);
			if((li!=lj) || ((li==NULL) && (lj==NULL))){
				ini_relation_zones(src,i+1,j+1,dim);
			}
			
		}
	}
}

bool opt_zones_mat_closure(opt_zones_mat_t * oz, unsigned short int dim){
        #if defined(TIMING)
		start_timing();
	#endif
	unsigned short int n = dim + 1;
	int size = n*n;
	//printf("closure required\n");
	double sparsity = 1- ((double)(oz->nni/size));
	bool res;
	if(sparsity >= zone_sparse_threshold){
		if(oz->is_dense){
			/*****
				If matrix is dense, convert it into decomposed type,
				extract the independent components
			*****/
			oz->is_dense = false;
			oz->acl = extract_comps_zones(oz->mat,dim);
		}
		closure_comp_sparse(oz,dim);
		res = strengthening_zones(oz,dim);
	}
	else{
		sparsity = recalculate_zones_sparsity(oz,dim);
		if(sparsity >= zone_sparse_threshold){
			if(oz->is_dense){
				/*****
					If matrix is dense, convert it into decomposed type,
					extract the independent components
				*****/
				oz->is_dense = false;
				oz->acl = extract_comps_zones(oz->mat,dim);
			}
			closure_comp_sparse(oz,dim);
			res = strengthening_zones(oz,dim);	
		}
		else{
			if(!oz->is_dense){
				oz->is_dense = true;
				if(!oz->ti){
					oz->ti = true;
					convert_to_dense_zones(oz,dim,false);	
				}
				free_array_comp_list(oz->acl);
			}
			#if defined (VECTOR)
				res = closure_dense(oz,dim);		
			#else
				res = closure_dense_scalar(oz,dim);
			#endif
		}
	}
	#if defined(TIMING)
		record_timing(zones_closure_time);
	#endif
	return res;
}

void opt_zones_mat_free(opt_zones_mat_t * oz){
	#if defined(TIMING)
		start_timing();
	#endif
	free(oz->mat);
	if(!oz->is_dense){
		free_array_comp_list(oz->acl);
	}
	free(oz);
	#if defined(TIMING)
		record_timing(zones_free_time);
	#endif
}

bool is_top_zones_mat(opt_zones_mat_t *oz, unsigned short int dim){
	double *m = oz->mat;
	unsigned short int n = dim + 1;
	int size = n*n;
	if(!oz->is_dense){
		array_comp_list_t * acl = oz->acl;
		if(acl->size==0){
			return true;
		}
		comp_list_t * cl = acl->head;
		while(cl!=NULL){
			unsigned short int * ca = to_sorted_array(cl,n);
			unsigned short int comp_size = cl->size;
			unsigned short int i,j;
			for(i=0; i < comp_size; i++){
				unsigned short int i1 = ca[i]+1;
				for(j=0; j < comp_size; j++){
					unsigned short int j1 = ca[j]+1;
					if(i1!=j1){
						if(m[n*i1+j1]!=INFINITY){
							free(ca);
							return false;
						}
					}
				}
				if(m[n*i1]!=INFINITY){
					free(ca);
					return false;
				}
				if(m[i1]!=INFINITY){
					free(ca);
					return false;
				}
			}
			free(ca);
			cl = cl->next;
		}
		clear_array_comp_list(acl);
		return true;
	}
	else{
		bool flag = true;
		int i;
		/*small trick just replace 0 at diagonal with infinity temporarily*/
		for(i=0; i < n; i++){
			m[n*i+i] = INFINITY;
		}
		#if defined (VECTOR)
			v_double_type infty = v_set1_double(INFINITY);
			v_int_type one = v_set1_int(1);
			for(int i = 0; i < size/v_length; i++){
				v_double_type t1 = v_load_double(m + i*v_length);
				v_double_type res = v_cmp_double(t1,infty, _CMP_EQ_OQ);
				v_int_type op = v_double_to_int(res);
				if(!v_test_int(op,one)){
					flag = false;
					break;
				}
		
			}
		#else
			for(int i = 0; i < (size/v_length)*v_length; i++){
				if(m[i]!=INFINITY){
					flag = false;
					break;
				}
			}
		#endif
		for(int i = (size/v_length)*v_length; (i <size) && flag; i++){
			if(m[i] != INFINITY){
				flag = false;
				break;
			}
		}
		
		/* now make diagonal elements OPT_ZERO again*/	
		for(int i = 0; i < n; i++){
			int ind = n*i+i;
			m[ind] = 0;
		}
		return flag;
	}
}

bool check_is_leq_zones(double * m1 , double *m2, unsigned short int i, unsigned short int j, unsigned short int dim){
	unsigned short int n = dim + 1;
	int ind = n*i+j;
	int ind1 = n*i;
	int ind2 = j;
	double val = min(m1[ind],m1[ind1] + m1[ind2]);	
	if(val > m2[ind]){
		return false;
	}

	ind = n*j+i;
	ind1 = n*j;
	ind2 = i;
	val = min(m1[ind],m1[ind1] + m1[ind2]);
	if(val > m2[ind]){
		return false;
	}

	return true;
}


bool is_lequal_zones_mat(opt_zones_mat_t *oz1, opt_zones_mat_t *oz2, unsigned short int dim){
	double *m1 = oz1->mat;
	double *m2 = oz2->mat;
	unsigned short int n = dim + 1;
	int size = n*n;
	//print_mat(oz2,dim);
	if(!oz1->is_dense && !oz2->is_dense){
		/****
			If both oo1 and oo2 are decomposed type, apply decomposed type operator
			Compute intersection of components and check if oo1 is <= oo2 in lattice as follows:
			1. Elements corresponding to components in oo2 not in intersection should all be oo
			2. Elements cooresponding to components in intersection should be equal  
		****/
		unsigned short int * arr_map1 = create_array_map(oz1->acl,dim);
		
		comp_list_t * cl2 = oz2->acl->head;
		while(cl2!=NULL){
			unsigned short int comp_size = cl2->size;
			unsigned short int * ca = to_sorted_array(cl2,dim);
			for(int i = 0; i < comp_size; i++){
				unsigned short int i1 = ca[i]+1;
				unsigned short int ci = arr_map1[i1-1];
				for(int j = 0; j <i; j++){
					unsigned short int j1 = ca[j]+1;
					unsigned short int cj = arr_map1[j1-1];
					if(!ci || !cj || ci!=cj){
						if(!check_trivial_relation_zones(m2,i1,j1,dim)){
                            				free(ca);
							//free_array_comp_list(acl);
							free(arr_map1);
							return false;
						}
					}
					else{
						if(!check_is_leq_zones(m1,m2,i1,j1,dim)){
							free(ca);
							//free_array_comp_list(acl);
							free(arr_map1);
							return false;
						}
					}
				}
				int ind = n*i1;
				if(!ci){
					if((m2[i1]!=INFINITY) || (m2[ind]!=INFINITY) ){
						free(arr_map1);
						return false;
					}
				}
				else{
					if((m1[i1] > m2[i1]) || (m1[ind] > m2[ind]) ){
						free(arr_map1);
						return false;
					}
				}
				
			}
			free(ca);
			cl2 = cl2->next;
		}	
		free(arr_map1);
	}
	else{
		/***
			Apply dense operator.
			If the operand is sparse, we fully initialize it,
			but no need to change type.
		****/
		if(!oz1->is_dense){
			if(!oz1->ti){
				oz1->ti = true;
				convert_to_dense_zones(oz1,dim,false);
			}
			
		}
		if(!oz2->is_dense){
			if(!oz2->ti){
				oz2->ti = true;
				convert_to_dense_zones(oz2,dim,false);
			}
			
		}
		if(oz1->ind!=-1){
			unsigned short int ind = oz1->ind;
			if(m1[ind] > m2[ind]){
				return false;
			}
		}
		#if defined(VECTOR)
			v_int_type one = v_set1_int(1);
			
			for(int i = 0; i < size/v_length; i++){
				v_double_type t1 = v_load_double(m1 + i*v_length);
				v_double_type t2 = v_load_double(m2 + i*v_length);
				v_double_type res = v_cmp_double(t1,t2, _CMP_LE_OQ);
				v_int_type op = v_double_to_int(res);
				if(!v_test_int(op,one)){
					oz1->ind = i*v_length;
					return false;
				}
			}
		#else
			for(int i = 0; i < (size/v_length)*v_length;i++){
				if(m1[i] > m2[i]){
					oz1->ind = i;
					return false;
				}
			}
		#endif
		for(int i = (size/v_length)*v_length; i < size; i++){
			if(m1[i] > m2[i]){
				return false;
			}
		}
	}
	
	return true;
}


bool check_is_eq_zones(double * m1 , double *m2, unsigned short int i, unsigned short int j, unsigned short int dim){
	unsigned short int n = dim + 1;
	int ind = n*i+j;
	int ind1 = n*i;
	int ind2 = j;
	double val1 = min(m1[ind],m1[ind1] + m1[ind2]);
	double val2 = min(m2[ind],m2[ind1] + m2[ind2]);	
	if(val1 != val2){
		return false;
	}

	ind = n*j+i;
	ind1 = n*j;
	ind2 = i;
	val1 = min(m1[ind],m1[ind1] + m1[ind2]);
	val2 = min(m2[ind],m2[ind1] + m2[ind2]);	
	if(val1!= val2){
		return false;
	}

	return true;
}



bool is_equal_zones_mat(opt_zones_mat_t *oz1, opt_zones_mat_t *oz2, unsigned short int dim){
	double *m1= oz1->mat;
	double *m2 = oz2->mat;
	unsigned short int n = dim+1;
	int size = n*n;
	if(!oz1->is_dense && !oz2->is_dense){
		/*****
			If both oo1 and oo2 are decomposed type, we use decomposed operator
			
		******/
		if(!is_equal_array_comp_list(oz1->acl,oz2->acl,dim)){
			/***
				oo1 and oo2 do not contain same set of independent components. 
				Compute intersection of components and check the elements for equality as follows:
				1. Elements corresponding to components in oo1 not in intersection should all be oo
				2. Elements corresponding to components in oo2 not in intersection should all be oo
				3. Elements cooresponding to components in intersection should be equal
			****/
			array_comp_list_t * acl = intersection_array_comp_list(oz1->acl,oz2->acl,dim);
			/****
				Check 1
			****/
			comp_list_t *cl1 = oz1->acl->head;
			unsigned short int * arr_map2 = create_array_map(oz2->acl,dim);
			while(cl1!=NULL){
				unsigned short int comp_size = cl1->size;
				unsigned short int *ca = to_sorted_array(cl1,dim);
				for(int i = 0; i < comp_size; i++){
					unsigned short int i1 = ca[i]+1;
					unsigned short int ci = arr_map2[i1-1];
					for(int j = 0; j < i; j++){
						unsigned short int j1 = ca[j]+1;
						unsigned short int cj = arr_map2[j1-1];
						if(!ci || !cj){
							if(!check_trivial_relation_zones(m1,i1,j1,dim)){
                                				free(ca);
								free(arr_map2);
								free_array_comp_list(acl);
								return false;
							}
						}
					}
					if(!ci){
						if(m1[i1]!=INFINITY){
							free(ca);
							free(arr_map2);
							free_array_comp_list(acl);
							return false;
						}
						if(m1[n*i1]!=INFINITY){
							free(ca);
							free(arr_map2);
							free_array_comp_list(acl);
							return false;
						}
					}
				}
				free(ca);	
				cl1 = cl1->next;
			}
			free(arr_map2);
			unsigned short int * arr_map1 = create_array_map(oz1->acl,dim);
			/*****
				Check 2
			******/
			comp_list_t *cl2 = oz2->acl->head;
			while(cl2!=NULL){
				unsigned short int comp_size = cl2->size;
				unsigned short int * ca = to_sorted_array(cl2,dim);
				for(int i = 0; i < comp_size; i++){
					unsigned short int i1 = ca[i]+1;
					unsigned short int ci = arr_map1[i1-1];
					for(int j = 0; j < i; j++){
						unsigned short int j1 = ca[j]+1;
						unsigned short int cj = arr_map1[j1-1];	
						if(!ci || !cj){
							if(!check_trivial_relation_zones(m2,i1,j1,dim)){
								free_array_comp_list(acl);
                                				free(ca);
								free(arr_map1);
								return false;
							}
						}  
					}
					if(!ci){
						if(m2[i1]!=INFINITY){
							free(ca);
                                                        free(arr_map2);
                                                        free_array_comp_list(acl);
							return false;
						}
						if(m2[n*i1]!=INFINITY){
							free(ca);
                                                        free(arr_map2);
                                                        free_array_comp_list(acl);
							return false;
						}
					}
				}
				free(ca);
				cl2 = cl2->next;
			}
			free(arr_map1);
			
			/******
				Check 3
			*******/	
			comp_list_t * cl = acl->head;
			while(cl!=NULL){
				unsigned short int comp_size = cl->size;
				unsigned short int * ca = to_sorted_array(cl,dim);
				for(int i = 0; i < comp_size; i++){
					unsigned short int i1 = ca[i] + 1;
					for(int j = 0; j < i; j++){
						unsigned short int j1 = ca[j] + 1;
						if(!check_is_eq_zones(m1,m2,i1,j1,dim)){
                            				free(ca);
							free_array_comp_list(acl);
							return false;
						}
					}
					if((m1[i1]!=m2[i1]) || (m1[n*i1]!=m2[n*i1]) ){
						free(ca);
						free_array_comp_list(acl);
						return false;
					}
				}
				free(ca);
				cl = cl->next;
			}
			
						
						
			free_array_comp_list(acl);
		}
		else{
			/****
				If we come here we know that both oo1 and oo2 contain the same set of Independent components
			*****/
		
			comp_list_t *cl = oz1->acl->head;
			while(cl!=NULL){
				unsigned short int comp_size = cl->size;
				unsigned short int * ca = to_sorted_array(cl,dim);
				for(int i = 0; i < comp_size; i++){
					unsigned short int i1 = ca[i] + 1;
					for(int j = 0; j < i; j++){
						unsigned short int j1 = ca[j]+1;
						if(!check_is_eq_zones(m1,m2,i1,j1,dim)){
							free(ca);
							return false;
						}
					}
					if((m1[i1]!=m2[i1]) || (m1[n*i1]!=m2[n*i1])){
						free(ca);
						return false;
					}
				}
				free(ca);
				cl= cl->next;
			}
		}
  		
  		
		
	}
	else{
		/***
			If the operand is sparse, we fully initialize it,
			but no need to change type.
		****/
		if(!oz1->is_dense){
			if(!oz1->ti){
				oz1->ti = true;
				convert_to_dense_zones(oz1,dim,false);
			}
			
		}
		if(!oz2->is_dense){
			if(!oz2->ti){
				oz2->ti = true;
				convert_to_dense_zones(oz2,dim,false);
			}
			
		}
		
		#if defined(VECTOR) 
			v_int_type one = v_set1_int(1);
			for(int i = 0; i < size/v_length; i++){
				v_double_type t1 = v_load_double(m1 + i*v_length);
				v_double_type t2 = v_load_double(m2 + i*v_length);
				v_double_type res = v_cmp_double(t1,t2, _CMP_EQ_OQ);
				v_int_type op = v_double_to_int(res);
				if(!v_test_int(op,one)){
					return false;
				}
			}
	
		#else
			for(int i = 0; i < (size/v_length)*v_length;i++){
				if(m1[i] != m2[i]){
					return false;
				}
			}
		#endif
		for(int i = (size/v_length)*v_length; i < size; i++){
			if(m1[i] != m2[i]){
				return false;
			}
		}
		
	}
	return true;
}


void meet_zones_mat(opt_zones_mat_t *oz, opt_zones_mat_t *oz1, 
				 opt_zones_mat_t *oz2, unsigned short int dim, 
				 bool destructive){
	double *m = oz->mat;
	double *m1 = oz1->mat;
	double *m2 = oz2->mat;
	unsigned short int n = dim + 1;
        int size = n*n;
	if(!oz1->is_dense && !oz2->is_dense){
		/*****
			If both oo1 and oo2 are decomposed type, then apply the decomposed type operator,
			1. Compute union of corresponding sets of independent components.
			2. Initialize oo1 elements corresponding to components in union not in oo1 
			3. Initialize oo2 elements corresponding to components in union not in oo2
			4. Compute meet of oo1 and oo2 by operating on elements corresponding to union
		*****/
		array_comp_list_t * temp = oz1->acl;
		/*****
			Step 1
		******/
		array_comp_list_t * acl = union_array_comp_list(oz1->acl,oz2->acl,dim);		
		oz->is_dense = false;
		if(!destructive){
			oz->ti = false;
			m[0] = 0;
			//printf("m: %d dim: %d\n",m==NULL,dim);
			//fflush(stdout);
			//set_mat_top(m,dim+1);
		}
		comp_list_t * cl = acl->head;
		unsigned short int * arr_map1 = create_array_map(oz1->acl,dim);
		unsigned short int * arr_map2 = create_array_map(oz2->acl,dim);
		while(cl!=NULL){
			unsigned short int comp_size = cl->size;
			unsigned short int * ca = to_sorted_array(cl,dim);
			/*****
				Step 2
			******/
			if(!oz1->ti){
				
				for(int i = 0; i < comp_size; i++){
					unsigned short int i1 = ca[i]+1;
					unsigned short int ci = arr_map1[i1-1];
					for(int j = 0; j <i; j++){
						unsigned short int j1 = ca[j]+1;
						unsigned short int cj = arr_map1[j1-1];
						if(!ci || !cj || ci!=cj){
							ini_relation_zones(m1,i1,j1,dim);
							//handle_binary_relation(m1,oo1->acl,i1,j1,dim);
						}
					}
					if(!ci){
						m1[n*i1] = INFINITY;
						m1[i1] = INFINITY;
						m1[n*i1+i1] = 0;
					}
				}
			}
			/*****
				Step 3
			******/
			if(!oz2->ti){
				for(int i = 0; i < comp_size; i++){
					unsigned short int i1 = ca[i]+1;
					unsigned short int ci = arr_map2[i1-1];
					for(int j = 0; j <i; j++){
						unsigned short int j1 = ca[j]+1;
						unsigned short int cj = arr_map2[j1-1];
						if(!ci || !cj || ci!=cj){
							ini_relation_zones(m2,i1,j1,dim);
							//handle_binary_relation(m2,oo2->acl,i1,j1,dim);
						}
					}
					if(!ci){
						m2[n*i1] = INFINITY;
						m2[i1] = INFINITY;
						m2[n*i1+i1] = 0;
					}
				}
			}
			/*****
				Step 4
			******/
			for(unsigned short int i = 0; i < comp_size; i++){
				unsigned short int i1 = ca[i]+1;
				for(unsigned short int j = 0; j < comp_size; j++){
					unsigned short int j1 = ca[j]+1;
					int ind = n*i1+j1;	
					m[ind] = min(m1[ind],m2[ind]);
				}
				m[n*i1] = min(m1[n*i1],m2[n*i1]);
				m[i1] = min(m1[i1],m2[i1]); 
				m[n*i1+i1] = 0;
			}
			free(ca);
			cl = cl->next;
		}
		free(arr_map1);
		free(arr_map2);
		if(!destructive && oz->acl){
			free_array_comp_list(oz->acl);
		}
		oz->acl = acl;
		if(destructive){
			free_array_comp_list(temp);
		}
	}
	else{
		/*****
			Apply the dense operator,
			If the operand is decomposed type, we fully initialize it,
			no need to change type.
		******/
		
		if(!oz2->is_dense){
			if(!oz2->ti){
				oz2->ti = true;
				convert_to_dense_zones(oz2,dim,false);
			}
			
		}
		if(!oz1->is_dense){
			if(!oz1->ti){
				oz1->ti = true;
				convert_to_dense_zones(oz1,dim,false);
			}
			
		}
		
		if(!destructive){
			oz->ti = true;
			free_array_comp_list(oz->acl);
		}
		if(destructive && !(oz1->is_dense)){
			free_array_comp_list(oz1->acl);
		}
		oz->is_dense = true;
		
		#if defined(VECTOR)
			for(int i = 0; i < size/v_length; i++){
				v_double_type t1 = v_load_double(m1 + i*v_length);
				v_double_type t2 = v_load_double(m2 + i*v_length);
				v_double_type t3 = v_min_double(t1,t2);
				v_store_double(m + i*v_length,t3);	
				//count = count+4;
			}
			
		#else
			for(int i = 0; i < (size/v_length)*v_length;i++){
				m[i] = min(m1[i],m2[i]);
			}
		#endif
			for(int i = (size/v_length)*v_length; i < size; i++){
				m[i] = min(m1[i],m2[i]);
			}
	}
	
	
	if(oz->is_dense){
		oz->nni = size;
	}
	else{
		comp_list_t *cl = oz->acl->head;
		int count = 0;
		while(cl!=NULL){
			count = count + (cl->size)*(cl->size);
			cl = cl->next;
		}
		oz->nni = count;
	}	
}



bool is_equal_comp_zones(opt_zones_mat_t *oz1, opt_zones_mat_t *oz2, unsigned short int *arr_map1, unsigned short int *arr_map2, comp_list_t *cl, unsigned short int dim){
	double * m1 = oz1->mat;
	double * m2 = oz2->mat;
	unsigned short int comp_size = cl->size;
	unsigned short int * ca = to_sorted_array(cl,dim);
	unsigned short int i,j;
	unsigned short int n = dim+1;
	bool res = true;
	/************
		Initialize oo1 according to cl
	*************/
	if(!oz1->ti){	
		/*for(unsigned short int i = 0; i < comp_size; i++){
			unsigned short int i1 = ca[i];
			if(!arr_map1[i1]){
				ini_self_relation(m1,i1,dim);
			}
		}*/	
		for(unsigned short int i = 0; i < comp_size; i++){
			unsigned short int i1 = ca[i]+1;
			unsigned short int ci = arr_map1[i1-1];
			for(unsigned short int j = 0; j <i; j++){
				unsigned short int j1 = ca[j]+1;
				unsigned short int cj = arr_map1[j1-1];
				
				if(!ci||!cj|| ci!=cj){
					ini_relation_zones(m1,i1,j1,dim);
				}
			}
			if(!ci){
				m1[n*i1] = INFINITY;
				m1[i1] = INFINITY;
				m1[n*i1+i1] = 0;
			}
		}
	}

	/************
		Initialize oo2 according to cl
	*************/
	if(!oz2->ti){
		/*for(unsigned short int i = 0; i < comp_size; i++){
			unsigned short int i1 = ca[i];
			if(!arr_map2[i1]){
				ini_self_relation(m2,i1,dim);
			}
		}*/	
		for(unsigned short int i = 0; i < comp_size; i++){
			unsigned short int i1 = ca[i]+1;
			unsigned short int ci = arr_map2[i1-1];
			for(int j = 0; j <i; j++){
				unsigned short int j1 = ca[j]+1;
				unsigned short int cj = arr_map2[j1-1];
				if(!ci || !cj || ci!=cj){
					ini_relation_zones(m2,i1,j1,dim);
				}
			}
			if(!ci){
				m2[n*i1] = INFINITY;
				m2[i1] = INFINITY;
				m2[n*i1+i1] = 0;
			}
		}
	}

	
	
		
	
	/**************
		Perform equality test
	***************/
	for(i = 0; i < comp_size; i++){
		unsigned short int i1 = ca[i] + 1;
		for(j = 0; j < i; j++){
			unsigned short int j1 = ca[j]+1;
			int ind = n*i1 + j1;	
			int ind1 = n*i1;
			int ind2 = j1;
		
			int ind3 = n*j1 + i1;
			int ind4 = n*j1;
			int ind5 = i1;

			if((arr_map1[i1-1] * arr_map1[j1-1]) && arr_map1[i1-1]!=arr_map1[j1-1]){
				m1[ind] = min(m1[ind],m1[ind1]+m1[ind2]);
				m1[ind3] = min(m1[ind3],m1[ind4]+m1[ind5]);
			}
			if((arr_map2[i1-1] * arr_map2[j1-1]) && arr_map2[i1-1]!=arr_map2[j1-1]){
				m2[ind] = min(m2[ind],m2[ind1]+m2[ind2]);
				m2[ind3] = min(m2[ind3],m2[ind4]+m2[ind5]);
			}
			if(m1[ind]!=m2[ind]){
				res = false;
			}
			if(m1[ind3]!=m2[ind3]){
				res = false;
			}
			
		}
		if(m1[n*i1]!=m2[n*i1]){
			res = false;
		}
		if(m1[i1]!=m2[i1]){
			res = false;
		}
	}
	free(ca);
	
	return res;
}




void sparse_join_zones_mat(opt_zones_mat_t *oz, opt_zones_mat_t *oz1, opt_zones_mat_t *oz2, unsigned short int dim, bool destructive){
	
	array_comp_list_t *acl1 = oz1->acl;
	array_comp_list_t *acl2 = oz2->acl;
	
	array_comp_list_t * acl = union_array_comp_list(acl1,acl2,dim);
	
        comp_list_t * cl = acl->head;
	unsigned short int n = dim+1;
	unsigned short int * arr_map1 = create_array_map(acl1,dim);
	unsigned short int * arr_map2 = create_array_map(acl2,dim);
	//cl = acl->head;
	while(cl!=NULL){
		comp_t * c = cl->head;
		while(c!=NULL){
			unsigned short int num = c->num;
			if(!arr_map1[num] || !arr_map2[num]){
				c = c->next; 
				remove_comp(cl,num);
			}
			else{
				c = c->next;
			}
		}
		if(cl->size==0){
			comp_list_t * tmp = cl;
			cl = cl->next;
			remove_comp_list(acl,tmp);
		}
		else{
			cl = cl->next;
		}
	}
	
	cl = acl->head;
	unsigned short int num_comp = acl->size;
	char * map = (char *)calloc(num_comp,sizeof(char));
	unsigned short int k = 0, count = 0;
	oz->mat[0] = 0;
	while(cl!=NULL){
		if(is_equal_comp_zones(oz1,oz2,arr_map1,arr_map2,cl,dim)){
		   if(!destructive){
		   	opt_zones_comp_copy(oz->mat,oz1->mat,cl,dim);
		   }
		}
		else{
			 map[k] = 1;
			 count++;
		}
		k++;
		cl = cl->next;
	} 
	
	if(count>=1){
		
		comp_list_t * res=NULL;
		if(count>1){
			strengthening_incr_init_zones(oz1,acl,map,dim);
			strengthening_incr_init_zones(oz2,acl,map,dim);
			//unsigned short int * index = (unsigned short int *)calloc(2*(2*dim + 1),sizeof(unsigned short int));
			//double * temp = (double *)malloc(2*dim*sizeof(double));
			strengthening_inter_comp_zones(oz1,acl,map,dim);
			strengthening_inter_comp_zones(oz2,acl,map,dim);
			
			//	Update components
			char *res_map = (char *)calloc(dim,sizeof(char)); 
			res = create_comp_list();
			comp_list_t * cl = acl->head;
			for(k=0; k < num_comp; k++){
				if(map[k]){
					/*comp_t *c = cl->head;
					while(c!=NULL){
						unsigned short int num = c->num;
						if(arr_map1[num] && arr_map2[num]){
							insert_comp(res,num);
						}
						c = c->next;
					}*/
					union_comp_list(res,cl,res_map);
					comp_list_t *temp = cl;
					cl = cl->next;
					remove_comp_list(acl,temp);
				}
				else{
					cl = cl->next;
				}
			}
			if(res->size){
				insert_comp_list(acl,res);
			}
			free(res_map);
		}
		else{
			comp_list_t * cl = acl->head;
			for(k=0; k < num_comp; k++){
				if(map[k]){
					res = cl;
					break;
				}
				cl = cl->next;
			}
		}
		
		//	Perform the actual join
		
		
		double *m = oz->mat;
		double *m1 = oz1->mat;
		double *m2 = oz2->mat;
		
		unsigned short int comp_size = res->size;
		unsigned short int * ca = to_sorted_array(res,dim);
		for(unsigned short int i = 0; i < comp_size; i++){
			unsigned short int i1 = ca[i] + 1;
			for(unsigned short int j = 0; j < i; j++){
				unsigned short int j1 = ca[j]+1;
				int ind = n*i1 + j1;	
				m[ind] = max(m1[ind],m2[ind]);
				ind = n*j1+i1;
				m[ind] = max(m1[ind],m2[ind]);
			}
			m[n*i1] = max(m1[n*i1],m2[n*i1]);
			m[i1] = max(m1[i1],m2[i1]);
			m[n*i1+i1] = 0;
		}
		
		free(ca);
	}
	
	oz->is_dense = false;
	oz->ti = false;
	oz->nni = min(oz1->nni,oz2->nni);
	free(map);
	array_comp_list_t *tacl = oz->acl;
	
	if(destructive){
		free_array_comp_list(tacl);	
	}
	oz->acl = acl;
	free(arr_map1);
	free(arr_map2);
 	
	return;
}


void join_zones_mat(opt_zones_mat_t *oz, opt_zones_mat_t *oz1, 
				 opt_zones_mat_t *oz2, unsigned short int dim, 
				 bool destructive){
	double *m = oz->mat;
	double *m1 = oz1->mat;
	double *m2 = oz2->mat;
	//int count = 0;
	unsigned short int n = dim + 1;
	int size = n*n;
	if((!oz1->is_dense) || (!oz2->is_dense)){
		/******
			If either oo1 or oo2 is decomposed type, apply the decomposed type operator
			1. Compute the intersection of corresponding sets of independent components.
			2. Compute join of oo1 and oo2 by operating on elements corresponding to intersection.
		******/
		array_comp_list_t * temp = oz1->acl;
		/*****
			Step 1
		******/
		if(oz1->is_dense){
			if(!destructive && oz->acl){
				free_array_comp_list(oz->acl);
			}
			oz->acl = copy_array_comp_list(oz2->acl);
		}
		else if(oz2->is_dense){
			if(!destructive){
				if(oz->acl){
					free_array_comp_list(oz->acl);
				}
				oz->acl = copy_array_comp_list(oz1->acl);
			}
		}
		else{
			if(!destructive){
				if(oz->acl){
					free_array_comp_list(oz->acl);
				}
			}
			oz->acl = intersection_array_comp_list(oz1->acl,oz2->acl,dim);
			// Can only destroy after computing intersection
			if(destructive){
				free_array_comp_list(temp);
			}
		}
		oz->is_dense = false;
		
		if(!destructive){
			oz->ti = false;
			m[0] = 0;
			//set_mat_top(m,dim+1);
		}
		/******
			Step 2
		******/
		comp_list_t * cl = oz->acl->head;
		m[0] = 0;
		while(cl!=NULL){
			unsigned short int comp_size = cl->size;
			unsigned short int * ca = to_sorted_array(cl,dim);
			for(int i = 0; i < comp_size; i++){
				int i1 = ca[i]+1;
				for(int j = 0; j < comp_size; j++){
					int j1 = ca[j]+1;
					int ind = n*i1 + j1;	
					m[ind] = max(m1[ind],m2[ind]);
				}
				m[n*i1] = max(m1[n*i1],m2[n*i1]);
				m[i1] = max(m1[i1],m2[i1]);
				m[n*i1+i1] = 0;
			}
			free(ca);
			cl = cl->next;
		}
		
	}
	else{
		/******
			Apply dense operator if both operands are dense
		*******/
		if(destructive && !(oz1->is_dense)){
			free_array_comp_list(oz1->acl);
		}
		if(!destructive){
			oz->ti = true;
			free_array_comp_list(oz->acl);
		}
		oz->is_dense = true;
		
		#if defined(VECTOR)
			for(int i = 0; i < size/v_length; i++){
				v_double_type t1 = v_load_double(m1 + i*v_length);
				v_double_type t2 = v_load_double(m2 + i*v_length);
				v_double_type t3 = v_max_double(t1,t2);
				v_store_double(m + i*v_length,t3);
				
			}
			
		#else
			for(int i = 0; i < (size/v_length)*v_length;i++){
				m[i] = max(m1[i],m2[i]);
			}
		#endif
			for(int i = (size/v_length)*v_length; i <size; i++){
				m[i] = max(m1[i],m2[i]);
			}
	}
	oz->nni = min(oz1->nni,oz2->nni);
}


void set_infinity_widening_zones(double *m, unsigned short int i, unsigned short int j, unsigned short int dim){
	unsigned short int n = dim + 1;
	int ind = n*i+j;
	m[ind] = INFINITY;

	ind = n*j+i;
	m[ind] = INFINITY;
}

void set_widening_zones(double *m, double *m1, double *m2, unsigned short int i, unsigned short int j, unsigned short int dim, int * count){
	unsigned short int n = dim+1;
	int ind = n*i+j;
	int ind1 = n*i;
	int ind2 = j;
	double val1 = m1[ind];
       // if(val1!=INFINITY){
		
	//	val1 = min(val1,(m1[ind1] + m1[ind2])/2);
	//}
	double val2 = min(m2[ind],m2[ind1] + m2[ind2]);	
	if(val1 >= val2){
		m[ind] = val1;
		if(m[ind]!=INFINITY){
			*count++;
		}
	}
	else{
		m[ind] = INFINITY;
	}

	ind = n*j+i;
	ind1 = n*j;
	ind2 = i;
	val1 = m1[ind];
	//if(val1!=INFINITY){
		
	//	val1 = min(val1,(m1[ind1] + m1[ind2])/2);
	//}
	val2 = min(m2[ind],m2[ind1] + m2[ind2]);	
	if(val1 >= val2){
		
		m[ind] = val1;
		if(m[ind]!=INFINITY){
			*count++;
		}
	}
	else{
		m[ind] = INFINITY;
	}
}

void sparse_widening_zones_mat(opt_zones_mat_t *oz, opt_zones_mat_t *oz1, opt_zones_mat_t *oz2, unsigned short int dim){
	double *m = oz->mat;
	double *m1 = oz1->mat;
	double *m2 = oz2->mat;
	int count = 0;
	unsigned short int n = dim+1;
	
	/******
		If either oo1 or oo2 is decomposed type, apply the decomposed type operator
		1. Compute the intersection of corresponding sets of independent components.
		2. Operate on elements corresponding to intersection.
	******/
	oz->is_dense = false;
	oz->ti = false;
	
	oz->acl = copy_array_comp_list(oz1->acl);
			
	/*******
		Step 2
	*******/
	comp_list_t *cl = oz->acl->head;
	unsigned short int * arr_map2 = create_array_map(oz2->acl,dim);
	m[0] = 0;	
	while(cl!=NULL){
		unsigned short int comp_size = cl->size;
		unsigned short int * ca = to_sorted_array(cl,dim);
		for(unsigned short int i =0; i < comp_size; i++){
			unsigned short int i1 = ca[i]+1;
			for(unsigned short int j = 0; j < i; j++){
				unsigned short int j1 = ca[j]+1;
				if(!arr_map2[i1-1] || !arr_map2[j1-1]){
					set_infinity_widening_zones(m,i1,j1,dim);
				}
				else{
					set_widening_zones(m,m1,m2,i1,j1,dim,&count);
				}
			}
			int ind = n*i1;
			if(!arr_map2[i1-1]){
				m[ind] = INFINITY;
				m[i1] = INFINITY;
			}
			else{
				double val1 = m1[ind];
				double val2 = m2[ind];
				if(val1 >= val2){
					m[ind] = val1;
					if(m[ind]!=INFINITY){
						count++;
					}
				}
				else{
					m[ind] = INFINITY;
				}
				val1 = m1[i1];
				val2 = m2[i1];
				if(val1 >= val2){
					m[i1] = val1;
					if(m[i1]!=INFINITY){
						count++;
					}
				}
				else{
					m[i1] = INFINITY;
				}
			}
			m[n*i1+i1] = 0;
		}
		
		free(ca);
		cl = cl->next;
	}
	free(arr_map2);
	oz->nni = count;
	
	
}


void widening_zones_mat(opt_zones_mat_t *oz, opt_zones_mat_t *oz1, 
				     opt_zones_mat_t *oz2, unsigned short int dim){
	double *m = oz->mat;
	double *m1 = oz1->mat;
	double *m2 = oz2->mat;
	unsigned short int n = dim + 1;
	int count = 0;
	int size = n*n;
	if(!oz1->is_dense || !oz2->is_dense){
		/******
			If either oo1 or oo2 is decomposed type, apply the decomposed type operator
			1. Compute the intersection of corresponding sets of independent components.
			2. Operate on elements corresponding to intersection.
		******/
		oz->is_dense = false;
		oz->ti = false;
		//set_mat_top(m,dim+1);
		/*******
			Step 1
		*******/
		if(oz->acl!=NULL){
			free_array_comp_list(oz->acl);
		}
		if(oz1->is_dense){
			oz->acl = copy_array_comp_list(oz2->acl);
		}
		else if(oz2->is_dense){
			oz->acl = copy_array_comp_list(oz1->acl);
		}
		else{
			
			oz->acl = intersection_array_comp_list(oz1->acl,oz2->acl,dim);
			
		}
		/*******
			Step 2
		*******/
		comp_list_t *cl = oz->acl->head;
		m[0] = 0;
		while(cl!=NULL){
			unsigned short int comp_size = cl->size;
			unsigned short int * ca = to_sorted_array(cl,dim);
			for(int i = 0; i < comp_size; i++){
				int i1 = ca[i]+1;
				for(int j = 0; j < comp_size; j++){
					int j1 = ca[j]+1;
					int ind = n*i1 + j1;
					if(m1[ind] >=m2[ind]){
						m[ind] = m1[ind];
					}
					else{
						m[ind] = INFINITY;
					}
					if(m[ind]!=INFINITY){
						count++;
					}
				}
				m[n*i1] = (m1[n*i1] >= m2[n*i1]) ? m1[n*i1] : INFINITY;
				m[i1] =   (m1[i1] >= m2[i1]) ? m1[i1] : INFINITY;
				m[n*i1+i1] = 0;
			}
			free(ca);
			cl = cl->next;
		}
	}
	else{
		/*****
			Apply the dense operator in case both operands are dense.
		******/
		oz->is_dense = true;
		oz->ti = true;
		free_array_comp_list(oz->acl);
		for(int i = 0; i <size; i++){
			if(m1[i] >= m2[i]){
				m[i] = m1[i];
			}
			else{
				m[i] = INFINITY;
			}
			if(m[i] != INFINITY){
				count++;
			}
		}
	}
	oz->nni = count;
	
}	


void forget_array_zones_mat(opt_zones_mat_t *oz, elina_dim_t * tdim, 
			    unsigned short int dim, unsigned short int size,
			    bool project){
	double *m = oz->mat;
	array_comp_list_t * acl = oz->acl;
	unsigned short int n = dim + 1;
	for(unsigned short int i = 0; i < size; i++){
		elina_dim_t d = tdim[i] + 1;
		/*****
			If the matrix is decomposed type,remove
			the component part of arr if it exists
		******/
		if(!oz->is_dense){
			comp_list_t *cl = find(acl,tdim[i]);
			if(cl!=NULL){
				remove_comp(cl,tdim[i]);
			}
		}
	
		double *pd = m + n*d;
		unsigned short int i, j;	
		#if defined(VECTOR)
			v_double_type infty = v_set1_double(INFINITY);
			for(j = 0; j < n/v_length; j++){
				v_store_double(pd + j*v_length,infty);
			}
		#else
			for(j = 0; j < (n/v_length)*v_length;j++){
				pd[j] = INFINITY;
				
			}
		#endif
		for(j = (n/v_length)*v_length; j < n; j++){
			pd[j] = INFINITY;
		}
	
		for(i = 0; i < n; i++){
			m[n*i + d] = INFINITY;
			
		}

		if(project){
			m[d] = 0;
			m[n*d] = 0;
			/*****
				Handle Independent Components in case of Project
			******/
			if(!oz->is_dense){
				comp_list_t * cj = create_comp_list();
				insert_comp(cj,tdim[i]);
				insert_comp_list(acl,cj);
			}
			
			oz->nni = oz->nni + 2;
		}
		else{
			m[d] = INFINITY;
			m[n*d] = INFINITY;			
		}
		m[n*d + d] = 0;
	}
}


void opt_zones_mat_addrem_dimensions(opt_zones_mat_t * dz, opt_zones_mat_t *sz,
				     elina_dim_t *pos, unsigned short int nb_pos, 
				     unsigned short int mult, unsigned short int dim,
				     bool add){
  
  unsigned short int i,j;
  double * dst = dz->mat;
  double * src = sz->mat;
  unsigned short int * map = (unsigned short int *)calloc(dim+1,sizeof(unsigned short int));
  elina_dim_t * add_pos = (elina_dim_t *)calloc(nb_pos,sizeof(elina_dim_t));
  unsigned short int new_dim=0;
  
  if(nb_pos){
	  
	  if(add){
		/****
			Exact number of non infinities for add 
		****/
		new_dim = dim + (nb_pos > mult ? nb_pos : mult);
		
		int max_nni = (new_dim+1)*(new_dim +1);
		dz->nni = min(max_nni,sz->nni + nb_pos);
		/******
			Because of add, the independent components are shifted in destination matrix
			Build Map of Independent Components for add in destination matrix
			Map contains new positions of components in the destination matrix
		******/
	
		int l = 0,k= 0,p =0,ac=0;
		
		for(int i = 0; i <=dim; i++){
			//
			while((l < nb_pos) &&(i==pos[l])){
				/***
					For expand mult can be greater than 1 in which case nb_pos will be 1
				****/
				add_pos[ac] = p;
				p = p + mult;
				ac++;
				l++;
			}
		
			map[k] = p;
			k++;
			p++;
		}
	
	  }
	  else{
		/**** 
			Approximation of number of non infinities for remove 
		****/
	
		new_dim = dim - nb_pos;
		int new_size = (new_dim+1)*(new_dim+1);
		dz->nni = min(sz->nni - nb_pos,new_size);
		dz->nni = max(new_dim+1,dz->nni);
		/******
			Because of remove, the independent components are shifted in destination matrix
			Build Map of Independent Components for add in destination matrix
			Map contains new positions of components in the destination matrix
		******/
		unsigned short int l = 0;
		for(unsigned short int i = 0; i < dim; i++){
			if((l < nb_pos) && (i==pos[l])){
				map[i] = new_dim;
				l++;
			}
			else{
				map[i] = i - l;
			}
		}
	  }
		
		/*****
			Handle Independent Components for decomposed type
		******/
		
		if(!sz->is_dense){
			array_comp_list_t * acl1 = sz->acl;	
			array_comp_list_t * acl2 = dz->acl;
			comp_list_t *cl1 = acl1->head;
			while(cl1!=NULL){
				comp_list_t *cl2 = create_comp_list();
				comp_t * c1 = cl1->head;
				while(c1!=NULL){
					unsigned short int num = c1->num;
					if(map[num]!=new_dim){
						insert_comp(cl2,map[num]);
					}
					c1 = c1->next;
				}
				if(cl2->size > 0){
					insert_comp_list(acl2,cl2);
				}
				else{
					free_comp_list(cl2);
				}
				cl1 = cl1->next;
			}
		}
		
     }
     dz->is_dense = sz->is_dense;
     unsigned short int dn = new_dim + 1;
     unsigned short int sn = dim + 1;
     //set_mat_top(dst,new_dim+1);
     if(!sz->is_dense){
		/*****
			If the source matrix is decomposed type,
			apply the decomposed type operator.
		*****/
		if(nb_pos){
			dst[0] = 0;
			array_comp_list_t * acl = sz->acl;
			comp_list_t *cl = acl->head;
			while(cl!=NULL){
				unsigned short int * ca = to_sorted_array(cl,dim);
				unsigned short int comp_size = cl->size;
				for(i = 0; i < comp_size; i++){
					unsigned short int i1 = ca[i];
					unsigned short int ni = map[i1];
					if(ni==new_dim){
						continue;
					}
				
					for(j = 0; j < comp_size; j++){
						unsigned short int j1 = ca[j];
						unsigned short int nj = map[j1];
						if(nj==new_dim){
							continue;
						}
						dst[dn*(ni+1) + nj + 1] = src[sn*(i1+1) + j1 + 1];
					}
					dst[dn*(ni+1)] = src[sn*(i1+1)];
					dst[ni+1] = src[i1+1];
					dst[dn*(ni+1) + ni+1] = 0;
				}
                		free(ca);
				cl = cl->next;
			}
			
		}
		dz->ti = false;
     }
     else{
	  /*****
		If the source matrix is dense type,
		apply the dense type operator.
	 *****/
	 if(!add){
		dst[0] = 0;
	 }
	 for(i=0; i < dim; i++){
		unsigned short int ni = map[i];
		if(ni==new_dim){
			continue;
		}
		for(j=0; j < dim; j++){
			unsigned short int nj = map[j];
			if(nj==new_dim){
				continue;
			}
			dst[dn*(ni+1) + nj + 1] = src[sn*(i+1) + j + 1];
		}
		dst[dn*(ni+1)] = src[sn*(i+1)];
		dst[ni+1] = src[i+1];
		if(!add){
			dst[dn*(ni+1) + ni+1] = 0;
		}
	 }
	 
	 if(add){
		forget_array_zones_mat(dz,add_pos,new_dim,nb_pos,false);		
	 }	
	dz->ti = true;
	
      }
	
	if(sz->ind==-1){
		dz->ind = -1;
	}
	else{
		unsigned short int i1 = sz->ind/sn;
		unsigned short int j1 = sz->ind % sn;
		unsigned short int ni, nj; 
		if(i1==0){
			nj = map[j1-1];
			dz->ind = nj+1; 
		}
		else if(j1==0){
			ni = map[i1-1];
			dz->ind = dn*(ni+1);
		}
		else{
			ni = map[i1-1];
			nj = map[j1-1];
			dz->ind = dn*(ni+1) + nj+1;
		}
	}
  	free(map);
        free(add_pos);
	
}

void opt_zones_mat_permute(opt_zones_mat_t * dz, opt_zones_mat_t *sz, 
			   unsigned short int dst_dim, unsigned short int src_dim,
		  	   elina_dim_t* permutation){
  double *dst = dz->mat;
  double *src = sz->mat; 
  unsigned short int sn = src_dim + 1;
  unsigned short int dn = dst_dim + 1;
  if(!sz->is_dense){
	  /******
		If source matrix is decomposed type, then apply decomposed operator.
	  ******/
	  dz->ti = false;
	  //set_mat_top(dst,dst_dim+1);
	  array_comp_list_t * acl1 = sz->acl;
  	  comp_list_t * cl1 = acl1->head;
	  dst[0] = 0;
	  while(cl1 != NULL){
		unsigned short int comp_size = cl1->size;
		unsigned short int * ca1 = to_sorted_array(cl1,src_dim);
		for(unsigned short int i = 0; i < comp_size; i++){
			unsigned short int i1 = ca1[i];
			unsigned short int new_ii = permutation[i1];
			
			if(new_ii >= dst_dim){
				continue;
			}
			for(unsigned short int j = 0; j < i; j++){
				unsigned short int j1 = ca1[j];
				//if(j1 > (i1|1)){
					//break;
				//}
				unsigned short int new_jj = permutation[j1];
				if(new_jj >= dst_dim){
					continue;
				}
				
				int d_ind = dn*(new_ii+1) + new_jj + 1;
				int s_ind = sn*(i1+1) + (j1+1);
				dst[d_ind] = src[s_ind];	
				d_ind = dn*(new_jj+1) + new_ii + 1;
				s_ind = sn*(j1+1) + i1 + 1;
				dst[d_ind] = src[s_ind];
			}
			dst[dn*(new_ii+1)] = src[sn*(i1+1)];
			dst[new_ii+1] = src[i1+1];
			dst[dn*(new_ii+1) + (new_ii+1)] = 0;
		}
		free(ca1);
		cl1 = cl1->next;
	  }
	/*****
		Handle Independent Components.
		The set of independent components gets permute 
		according to permutation specified by permutation array 
  	******/
 	array_comp_list_t * acl2 = dz->acl;
  	cl1 = acl1->head;
  	
  	while(cl1 != NULL){
		comp_list_t * cl2 = create_comp_list();
		comp_t * c1 = cl1->head;
		while(c1 != NULL){
			unsigned short int num = c1->num;
			insert_comp(cl2,permutation[num]);
			c1 = c1->next;
		}
		
		insert_comp_list(acl2,cl2);
		cl1 = cl1->next;
 	 } 
  }
  else{
	  /*****
		If source matrix is dense, apply the dense operator.
	  *****/
	  dz->ti = true;
          unsigned short int i,j;
	  dst[0] = 0;
	  for (i=0;i<src_dim;i++) {
	    unsigned short int new_ii = permutation[i];
	    if (new_ii >= dst_dim)  { continue; }
	    for (j=0;j<i;j++) {
		      unsigned short int new_jj = permutation[j];
		      if (new_jj >= dst_dim) continue;
		      int d_ind = dn*(new_ii+1) + new_jj + 1;
		      int s_ind = sn*(i+1) + (j+1);
		      dst[d_ind] = src[s_ind];
		      d_ind = dn*(new_jj+1) + new_ii + 1;
		      s_ind = sn*(j+1) + (i+1);
		      dst[d_ind] = src[s_ind];
	    }
	    dst[dn*(new_ii+1)] = src[sn*(i+1)];
	    dst[new_ii+1] = src[i+1];
	    dst[dn*(new_ii+1) + (new_ii+1)] = 0;
	 }
	 
  }
  if(sz->ind==-1){
		dz->ind = -1;
  }
  else{
	unsigned short int i1 = sz->ind/sn;
	unsigned short int j1 = sz->ind % sn;
	unsigned short int ni, nj; 
	if(i1==0){
		nj = permutation[j1-1];
		dz->ind = nj+1; 
	}
	else if(j1==0){
		ni = permutation[i1-1];
		dz->ind = dn*(ni+1);
	}
	else{
		ni = permutation[i1-1];
		nj = permutation[j1-1];
		dz->ind = dn*(ni+1) + nj+1;
	}
  }
  dz->nni = sz->nni;
  dz->is_dense = sz->is_dense;
  dz->is_top = sz->is_top;
}

/***********************
	Meet with Linear Constraint
************************/
zone_expr zone_expr_of_linexpr(opt_zones_internal_t* pr, double* dst,
			   elina_linexpr0_t* e, unsigned short int intdim, unsigned short int dim)
{
#define CLASS_COEFF(idx,coef)						\
  if ((dst[2*idx+2] == -coef) &&				\
      (dst[2*idx+3] == coef)) {				\
    if (z.type==OPT_ZERO) { z.type = OPT_UNARY;  z.i = idx; z.coef_i = coef; }	\
    else if(z.type==OPT_UNARY){\
      if(z.coef_i==-coef) { z.type = OPT_BINARY; z.j = idx; z.coef_j = coef; }	\
      else {z.type= OPT_OTHER;}}\
    continue;								\
  }
  
#define CLASS_VAR(idx)							\
  if (idx>=intdim) z.is_int = false;                                        \
  if (z.type==OPT_EMPTY) continue;						\
  if ((dst[2*idx+2] == 0) && (dst[2*idx+3] == 0)) continue;	\
  if (z.type>=OPT_BINARY) { z.type = OPT_OTHER; continue; } 			\
  CLASS_COEFF(idx,1);							\
  CLASS_COEFF(idx,-1);							\
  z.type = OPT_OTHER;
  
#define COEFF(c,i)                                                      \
  elina_coeff_reduce(&c);                                                   \
  if (opt_bounds_of_coeff(pr,dst + i,dst + i + 1,c)) z.type = OPT_EMPTY;    \
  if (c.discr!=ELINA_COEFF_SCALAR || !is_integer(dst[i])) z.is_int = false;
  zone_expr z = { OPT_ZERO, 0, 0, 0, 0, true};
  unsigned short int i;
  COEFF(e->cst,0);
  switch (e->discr) {
  case ELINA_LINEXPR_DENSE:
    if(e->size > dim)return z;
    for (i=0;i<e->size;i++) {
      
      COEFF(e->p.coeff[i],2*i+2);
      CLASS_VAR(i);
    }
    for (;i<dim;i++) {
      dst[2*i+2] = 0;
      dst[2*i+3] = 0;
    }
    break;
  case ELINA_LINEXPR_SPARSE:
	
    for (i=0;i<dim;i++) {
      dst[2*i+2] = 0;
      dst[2*i+3] = 0;
    }
	
    for (i=0;i<e->size;i++) {
      elina_dim_t d = e->p.linterm[i].dim;
      if (d==ELINA_DIM_MAX) continue;
      if(d>=dim)return z;
      COEFF(e->p.linterm[i].coeff,2*d+2);
      CLASS_VAR(d);
    }
	
    break;
  default: break;/********TODO: handle arg_assert arg_assert(0,return u;);*****/
  }
  return z;
}

comp_list_t * linexpr0_to_comp_list_zones(elina_linexpr0_t * expr){
	comp_list_t * cl = create_comp_list();
	size_t size = expr->size;
	size_t j;
	elina_linexpr_discr_t discr = expr->discr;
	if(discr==ELINA_LINEXPR_DENSE){
		elina_coeff_t* coeff = expr->p.coeff;
		for(j=0; j < size; j++){
			if(!elina_coeff_zero(&coeff[j])){
				insert_comp(cl,j);
			}
		}
	}
	else{
		elina_linterm_t* linterm = expr->p.linterm;
		for(j = 0; j < size; j++){
			elina_dim_t dim = linterm[j].dim;
			insert_comp(cl,dim);	
		}
	}
	return cl;
}

bool opt_zones_mat_add_lincons(opt_zones_internal_t * pr,opt_zones_mat_t *oz, unsigned short int intdim, 
			       unsigned short int dim, elina_lincons0_array_t *array, bool* exact, bool* respect_closure){
  double *m = oz->mat;
  unsigned short int i, j, k;
  int zi,zj;
  unsigned short int var_pending = 0; /* delay incremental closure as long as possible */
  unsigned short int closure_pending = 0;
  *exact = 1;
  unsigned short int n = dim + 1;
  int max_nni = n*n;
  int *ind1, *ind2;
  //posix_memalign((void **)&temp1, 32, 2*dim*sizeof(double));
  //posix_memalign((void **)&temp2, 32, 2*dim*sizeof(double));
  bool (*incr_closure)(opt_zones_mat_t * ,...);
  double size = n*n;
  //double sparsity = 1- ((double)(oz->nni)/size);
  
  /******
	Measure sparsity to decide on whether to use dense or decomposed type incremental closure.
	We do not recalculate sparsity here (if estimate of nni is too imprecise) as it increases overhead.
  ******/
  //if(sparsity >=zone_sparse_threshold){
	if(oz->is_dense){
		
		oz->is_dense = false;
		//printf("convert to sparse %g\n",sparsity);
		
		//fflush(stdout);
		oz->acl = extract_comps_zones(oz->mat,dim);
	}
	incr_closure = &incr_closure_comp_sparse;
  //}
  //else{ 
	//if(!oz->is_dense){
		//if(!oz->ti){
	//		oz->ti = true;
	//		convert_to_dense_zones(oz,dim,false);
	//	}
	//	oz->is_dense = true;
	//	free_array_comp_list(oz->acl);
	//}
  	//#if defined(VECTOR)
	//	incr_closure = &incr_closure_dense;
  //}
  	//#else
	//	incr_closure = &incr_closure_dense_scalar;
  //	#endif
 // }
  
  for (i=0;i<array->size;i++) {
   elina_constyp_t c = array->p[i].constyp;
    zone_expr z;

    switch (c) {

      /* skipped */
    case ELINA_CONS_EQMOD:
    case ELINA_CONS_DISEQ:
      *exact = 0;
      continue;

      /* handled */
    case ELINA_CONS_EQ:
    case ELINA_CONS_SUPEQ:
    case ELINA_CONS_SUP:
      break;

      /* error */
    default:
      assert(0);
    }
    
	
    /* now handle ==, >=, > */
     if(!oz->is_dense){
	//print_opt_oct_mat(oo,dim);
    	comp_list_t * clb = linexpr0_to_comp_list_zones(array->p[i].linexpr0);
	if(clb->size){
		if(!oz->ti){
			m[0] = 0;
			unsigned short int * arr_map = create_array_map(oz->acl,dim);
			unsigned short int *ca_b = to_sorted_array(clb,dim);
			unsigned short int comp_size_b = clb->size;
			
			for(unsigned short int i = 0; i < comp_size_b; i++){
				unsigned short int i1 = ca_b[i]+1;
				unsigned short int ci = arr_map[i1-1];
				for(unsigned short int j=0; j < i; j++){
					unsigned short int j1 = ca_b[j]+1;
					unsigned short int cj = arr_map[j1-1];
					if(!ci || !cj || ci!=cj){
						ini_relation_zones(m,i1,j1,dim);
					}
				}
				
				if(!ci){
					m[n*i1] = INFINITY;
					m[i1] = INFINITY;
					m[n*i1+i1] = 0;	
				}
			}
			comp_list_t * cl = oz->acl->head;
			unsigned short int k = 0;
			unsigned short int num_comp = oz->acl->size;
			char * map = (char *)calloc(num_comp,sizeof(char));
			while(cl!=NULL){
				if(!is_disjoint(cl,clb,dim)){
					map[k] = 1;
					
					unsigned short int * ca = to_sorted_array(cl,dim);
					unsigned short int comp_size = cl->size;
					for(unsigned short int i =0; i < comp_size_b; i++){
						unsigned short int i1 = ca_b[i]+1;
						unsigned short int ci = arr_map[i1-1];
						if(!ci){
							for(unsigned short int j = 0; j < comp_size; j++){
								unsigned short int j1 = ca[j]+1;
								ini_relation_zones(m,i1,j1,dim);
							}	
						}
					}
					free(ca);
				}
				k++;
				cl = cl->next;
			}
			strengthening_incr_init_zones(oz,oz->acl,map,dim);
			free(map);
			free(arr_map);
			free(ca_b);
		}
		//comp_list_t * cli = copy_comp_list(clb);
		//unsigned short int num = clb->head->num;
		//char * map = (char *)calloc(dim,sizeof(char));
		//comp_t * c = clb->head;
		//while(c!=NULL){
		//	unsigned short int num = c->num;
		//	if(find(oo->acl)){
		//	} 
		//	c=c->next;
		//}
		insert_comp_list_with_union(oz->acl,clb,dim);
		
	}
    }
	
    z = zone_expr_of_linexpr(pr,pr->tmp,array->p[i].linexpr0,intdim,dim);
     	
    /* transform e+[-a,b] > 0 into >= e+[-(a+1),b-1] >= 0 on integer constraints */
    if (z.is_int && c==ELINA_CONS_SUP) {
      c = ELINA_CONS_SUPEQ;
      pr->tmp[0] = pr->tmp[0] + 1;
      pr->tmp[1] = pr->tmp[1] - 1;
    }
    
    int count;
   
    switch (z.type) {
	
    case OPT_EMPTY:
      /* OPT_EMPTY constraint: no added information */
      break;

    case OPT_ZERO:
      if ((c==ELINA_CONS_SUPEQ && (pr->tmp[1]>=0)) ||
	  /* [-a,b] >= 0 <=> b >= 0 */
	  (c==ELINA_CONS_SUP && (pr->tmp[1]>0)) ||
	  /* [-a,b] > 0 <=> b > 0 */
	  (c==ELINA_CONS_EQ && (pr->tmp[0]>=0) && (pr->tmp[1]>=0))
	  /* [-a,b] = 0 <=> a >= 0 && b >= 0 */
	  )
	; /* trivial */
	
      else{ return true;} /* unsatisfiable */
      break;
    array_comp_list_t *acl;
    case OPT_UNARY:
	
      /* can we delay incremental closure further? */
      if (*respect_closure && closure_pending && var_pending!=z.i) {
          if (incr_closure(oz,dim,var_pending)){
                return true;
          }
      }
     
      count = oz->nni;
      closure_pending = 1;
      var_pending = z.i;
      /******
		Handle Indepenedent Components for Unary Case.
		If z.i if not part of any component then create 
		a new component containing z.i and also possibly initialize
		relation between z.i with itself.
      *******/
     
      if (z.coef_i==1){
		 zi = n*(z.i+1);
		 zj = z.i+1; 
      }else{
		 zi = z.i+1;
		 zj = n*(z.i+1);
      }
      if(m[zi]==INFINITY){
      	m[zi] = pr->tmp[1];
	count++;
      }
      else{
	m[zi] = min(m[zi], pr->tmp[1]);
      }
      /*  c_i X_i + [-a,b] >= 0 <=> -c_i X_i <= b */
      if (c==ELINA_CONS_EQ) {
	if(m[zj]==INFINITY){
		m[zj] = pr->tmp[0];
		count++;
	}
        else{
		m[zj] = min(m[zj], pr->tmp[0]);
	}
      }
      /*  c_i X_i + [-a,b] <= 0 <=>  c_i X_i <= a */
      if (c==ELINA_CONS_SUP) *exact = 0; /* not exact for strict constraints */
      oz->nni = min(max_nni,count);	
      break;

    case OPT_BINARY:
      
      /* can we delay incremental closure further? */
      if (*respect_closure && closure_pending &&
	  var_pending!=z.i && var_pending!=z.j) {
	
          if (incr_closure(oz,dim,var_pending)) {
              return true;
          }
      }
	
      closure_pending = 1;
      var_pending = (var_pending==z.j) ? z.j : z.i;
      count = oz->nni;
      /******
	Handle the Independent components for Binary Case for decomposed type.
	There are four cases to handle:
	1. If u.i and u.j are not present in any component. In this case, create
	   a component conatining u.i and u.j
	2. If u.j is present in a componenet and u.i is not, then add u.i to 
	   component containing u.j, (similar for the other case)
	3. If u.i and u.j are present in different components, then take the union 
	   of containing components.
        4. If u.i and u.j are already contained in the same component, 
	   then do nothing  
      ****/
 
      
      if ( z.coef_i==1){
		zi = n*(z.i+1) + z.j+1;
		zj = n*(z.j+1) + z.i+1;
      } 
      else{
 	 	zi = n*(z.j+1) + z.i+1;
		zj = n*(z.i+1) + z.j+1;
      }
      
      if(m[zi]==INFINITY){
		m[zi]=pr->tmp[1];
		count++;
      }
      else{
      	m[zi] = min(m[zi], pr->tmp[1]);
      }
	
      /*  c_i X_i + c_j X_j + [-a,b] >= 0 <=> -c_i X_i - c_j X_j <= b */
      if (c==ELINA_CONS_EQ){
	if(m[zj]== INFINITY){
		m[zj] = pr->tmp[0];
		count++;
	}
        else{
		m[zj] = min(m[zj], pr->tmp[0]); 
	}
      }
     
      /*  c_i X_i + c_j X_j + [-a,b] <= 0 <=>  c_i X_i + c_j X_j <= a */
      if (c==ELINA_CONS_SUP) *exact = 0; /* not exact for strict constraints */
      oz->nni = min(max_nni,count);	
      break;

    case OPT_OTHER:
      {
	/* general, approximated case */
	
	double tmpa = 0, tmpb = 0, Cb = 0, cb = 0;
	int Cinf = 0;            /* number of infinite upper bounds */
	int Cj1 = 0, Cj2 = 0; /* variable index with infinite bound */
	int cinf = 0;            /* number of infinite lower bounds */
	int cj1 = 0, cj2 = 0; /* variable index with infinite bound */

	*respect_closure = false; /* do not respect closure */

	count = oz->nni;
	/* compute 2 * upper bound, ignoring components leading to +oo */
	cb = pr->tmp[0];
	Cb = pr->tmp[1];
	for (j=0;j<dim;j++) {
	  double tmp[8];
	  
	  double b_inf = pr->tmp[2*j + 2];
	  double b_sup = pr->tmp[2*j + 3];
	  if((b_inf==0) && (b_sup==0)){
		continue;
	  }
	  
	  double a_inf = m[n*(j+1)];
	  double a_sup = m[j+1];
	  if((a_sup == 0) || (b_sup == 0)){
		tmp[0] = 0;
		tmp[4] = 0;
	  }
	  else{
		  tmp[0] = a_sup * b_sup;
	  	  tmp[4] = -a_sup; 
		  tmp[4] = tmp[4] * b_sup;
	  }
	
	  if((a_inf == 0) || (b_inf == 0)){
		  tmp[1] = 0;
		  tmp[5] = 0;		
	  }
	  else{
	  	  tmp[1] = a_inf * b_inf;
	  	  tmp[5] = - a_inf;  
		  tmp[5] = tmp[5] * b_inf;
	 }

	 if((a_sup== 0) || (b_inf == 0)){
		  tmp[2] = 0;
		  tmp[6] = 0;
	 }
	 else{
	  	  tmp[6] = a_sup * b_inf;
	  	  tmp[2] = -a_sup;  
		  tmp[2] = tmp[2] * b_inf;
	}
	
	if((a_inf == 0) || (b_sup == 0)){
		  tmp[3] = 0;
		  tmp[7] = 0;
	}
	else{
	  	  tmp[7] = a_inf * b_sup;
	  	  tmp[3] = -a_inf;  
		  tmp[3] = tmp[3] * b_sup;
	}

  	  tmpb = max(tmp[0],tmp[1]);
  	  tmpb = max(tmpb,tmp[2]);
  	  tmpb = max(tmpb,tmp[3]);

  	  tmpa = max(tmp[4],tmp[5]);
  	  tmpa = max(tmpa,tmp[6]);
  	  tmpa = max(tmpa,tmp[7]);
	  
	  if (tmpa == INFINITY) { 
		cinf++; 
		cj2 = cj1; 
		cj1 = j; 
	  }
	  else {
		cb = cb + tmpa;
	  }
	  if (tmpb == INFINITY) { 
		Cinf++; 
		Cj2 = Cj1; 
		Cj1 = j; 
	  }
	  else {
		Cb = Cb + tmpb;
	  }
	  
	}
	
	/* upper bound */
	if (Cb == INFINITY) ;
	
	else if (!Cinf) {
	  /* no infinite bound: derive quadratic number of bounds */
	 /******
		Handle Independent Components in case Quadratic Number of Bounds are Created.
		Find the component containing each j and merge it with other components using 
		strengthening, this creates over approximation.
	  *****/
          
	  //acl = oo->acl;
	  for (j=0;j<dim;j++) {
	    //contribution by vj is positive
	    if ((pr->tmp[2*j+2] <= -1) &&
		(m[j+1] != INFINITY)) {
	      /* -vj <= expr-vj <= max(expr) - max(vj) */
		
	      tmpa = Cb - m[j+1];
	      zj = -1;
	    }
	    //contribution by vj is negative
	    else if ((pr->tmp[2*j+3]<=-1) &&
		     (m[n*(j+1)] != INFINITY)) {
	      /* vj <= expr+vj <= max(expr) - max(-vj) */
		
	      tmpa = Cb - m[n*(j+1)];
	      zj = 1;
	    }
	    else continue;
	    
	    for (k=j+1;k<dim;k++) {
	      
	      if ((pr->tmp[2*k+2]<=-1) &&
		  (m[k+1] != INFINITY) && (zj==1)) {
		/* vj - vk <= max(expr) - max(-vj) - max (vk) */
		tmpb = tmpa - m[k+1];
		if(m[n*(k+1)+j+1] ==INFINITY){
			m[n*(k+1)+j+1] = tmpb;
			count++;
		}
		else{
			m[n*(k+1) + j+1] = min(m[n*(k+1) + j+1], tmpb);
		}
	      }
	      else if ((pr->tmp[2*k+3] <=-1) &&
		       (m[n*(k+1)] !=  INFINITY) &&(zj==-1)) {
		/* vk - vj <= max(expr) - max(vj) - max (-vk) */
		tmpb = tmpa - m[n*(k+1)];
		if(m[n*(j+1) + k+1]==INFINITY){
			m[n*(j+1)+k+1] = tmpb;
			count++;
		}
		else{
			m[n*(j+1)+k+1] = min(m[n*(j+1)+k+1],tmpb);
		}
	      }
	       /******
			Add Components to the new list
	       *******/
		
	    }
		
	  }
	  oz->nni = min(max_nni,count);
	}

	else if (Cinf==1) {
	  /* one infinite bound: derive linear number of bounds */
	  
	  if ((pr->tmp[2*Cj1+3] == -1) &&
	      (pr->tmp[2*Cj1+2] == 1)) zj = -1;
	  else if ((pr->tmp[2*Cj1+3] == 1) &&
		   (pr->tmp[2*Cj1+2] == -1)) zj = 1;
	  else goto Cbrk;
	  /****
		Handle the Independent Components in case Linear Number of Bounds are created. 
		Find the component containing each Cj1 and merge it with other components using 
		strengthening, this creates over approximation for the set of independent components.
	  *****/
	  //acl = oo->acl;
	 
	  for (k=0;k<dim;k++) {
	    if (k==Cj1) continue;

            if ((pr->tmp[2 * k + 2] <= -1) && (m[k + 1] != INFINITY) &&
                (zj == 1)) {
              /* vj - vk <= max(expr) - max(-vj) - max(vk) */
	      tmpb = Cb - m[k+1];
              if (m[n * (j + 1) + k + 1] == INFINITY) {
                m[n * (j + 1) + k + 1] = tmpb;
                count++;
              } else {
                m[n * (j + 1) + k + 1] = min(m[n * (j + 1) + k + 1], tmpb);
              }

            } else if ((pr->tmp[2 * k + 3] <= -1) &&
                       (m[n * (k + 1)] != INFINITY) && (zj == -1)) {
              /* vk-vj <= max(expr) - max(-vk) - max (vj) */
	      tmpb = Cb - m[n*(k+1)];
              tmpb = tmpb / 2;
              if (m[n * (k + 1) + j + 1] == INFINITY) {
                m[n * (k + 1) + j + 1] = tmpb;
                count++;
              } else {
                m[n * (k + 1) + j + 1] = min(m[n * (k + 1) + j + 1], tmpb);
              }
            }
            /*****
		Add components to new List
	    *****/
	    
	  }
	   
	   
          oz->nni = min(max_nni,count);
	}

	else if (Cinf==2) {
	  /* two infinite bounds: derive just one bound */
	  if ((pr->tmp[2*Cj1+3]==-1) &&
	      (pr->tmp[2*Cj1+2]==1)) zi = -1;
	  else if ((pr->tmp[2*Cj1+3] == 1) &&
		   (pr->tmp[2*Cj1+2] ==-1)) zi = 1;
	  else goto Cbrk;
	  if ((pr->tmp[2*Cj2+3] == -1) &&
	      (pr->tmp[2*Cj2+2] == 1)) zj = -1;
	  else if ((pr->tmp[2*Cj2+3] == 1) &&
		   (pr->tmp[2*Cj2+2] == -1)) zj = 1;
	  else goto Cbrk;
	  /*****
		Handle Independent Components in case only one bound is created.
		Handling is similar for the case of binary type.
	  ******/
		
	       
          int ind;
	  if((zi==-1) && (zj==1)){
            ind = n * (Cj1 + 1) + Cj2 + 1;
            if (m[ind] == INFINITY) {
              m[ind] = tmpa;
              count++;
	  	}
		else{
		  	m[ind] = min(m[ind],tmpa);
		}
	  }
	  else if((zi==1) && (zj==-1)){
            ind = n * (Cj2 + 1) + Cj1 + 1;
            if (m[ind] == INFINITY) {
              m[ind] = tmpa;
              count++;
	  	}
		else{
		  	m[ind] = min(m[ind],tmpa);
		}
	  }
	  oz->nni = min(max_nni,count);
	  
	}
	
	/* if more than two infinite bounds: do nothing */

      Cbrk:
	/* lower bound */
	count = oz->nni;
	if (c==ELINA_CONS_EQ) {
		
	if (cb == INFINITY) ;
	else if (!cinf) {
	   /******
		Handle Independent Components in case Quadratic Number of Bounds are Created.
		Handling is similar to quadratic case above.
	  *****/
         
	 // acl = oo->acl;
	  for (j=0;j<dim;j++) {
	    if ((pr->tmp[2*j+3] <= -1) &&
		(m[j+1] != INFINITY)) {
	      tmpa = cb - m[j+1];
	      zj = -1;
	    }
	    else if ((pr->tmp[2*j+2] <=-1) &&
		     (m[n*(j+1)] != INFINITY)) {
	      tmpa = cb - m[n*(j+1)];
	      zj = 1;
	    }
	    else continue;
            
	   
	    for (k=j+1;k<dim;k++) {
	      
	      if ((pr->tmp[2*k+3] <=-1) &&
		  (m[k+1] != INFINITY)&&(zj==1)) {
		tmpb = tmpa - m[k+1];
		if(m[n*(k+1) + j+1]==INFINITY){
			m[n*(k+1) + j+1] = tmpb;
			count++;
		}
		else{
			m[n*(k+1) + j+1] = min(m[n*(k+1) + j+1], tmpb);
		}
	      }
	      else if ((pr->tmp[2*k+2] <=-1) &&
		       (m[n*(k+1)] != INFINITY)&&(zj==-1)) {
		tmpb = tmpa - m[n*(j+1) + k+1];
		if(m[n*(j+1) + k+1]==INFINITY){
			m[n*(j+1) + k+1] = tmpb;
			count++;
		}
		else{
			m[n*(j+1) + k+1] = min(m[n*(j+1) + k+1], tmpb);
		}
	      }
		
		
	    }
		
	  }
	  oz->nni = min(max_nni,count);
	}
	else if (cinf==1) {
	  
	  if ((pr->tmp[2*cj1+2] ==-1) &&
	      (pr->tmp[2*cj1+3] == 1)) zj = 1;
	  else if ((pr->tmp[2*cj1+2] == 1) &&
		   (pr->tmp[2*cj1+3] ==-1)) zj = -1;
	  else goto cbrk;
	   /****
		Handle the Independent Components in case Linear Number of Bounds are Created.
		Handling is simialr to Linear case above. 
	  *****/
	  //acl = oo->acl;	
	   
	   
	   for (k=0;k<dim;k++) {
	     if (k==cj1) continue;
	   
	     if ((pr->tmp[2*k+3] <= -1) &&
		(m[k+1] != INFINITY)&&(zj==1)) {
		//vk-vj
	      tmpb = cb - m[k+1];

              if (m[n * (j + 1) + k + 1] == INFINITY) {
                m[n * (j + 1) + k + 1] = tmpb;
                count++;
              } else {
                m[n * (j + 1) + k + 1] = min(m[n * (j + 1) + k + 1], tmpb);
              }

            }
	    else if ((pr->tmp[2*k+2] <= -1) &&
		     (m[n*(k+1)] != INFINITY)&&(zj==-1)) {
	      tmpb = cb - m[n*(k+1)];
	      //Incremental initialization

              if (m[n * (k + 1) + j + 1] == INFINITY) {
                m[n * (k + 1) + j + 1] = tmpb;
                count++;
              } else {
                m[n * (k + 1) + j + 1] = min(m[n * (k + 1) + j + 1], tmpb);
              }
            }
	       
		
	  }
	   oz->nni = min(max_nni,count);
	}
	else if (cinf==2) {
	  if ((pr->tmp[2*cj1+2]==-1) &&
	      (pr->tmp[2*cj1+3] == 1)) zi = 1;
	  else if ((pr->tmp[2*cj1+2] == 1) &&
		   (pr->tmp[2*cj1+3] ==-1)) zi = -1;
	  else goto cbrk;
	  if ((pr->tmp[2*cj2+2] == -1) &&
	      (pr->tmp[2*cj2+3] == 1)) zj = 1;
	  else if ((pr->tmp[2*cj2+2] == 1) &&
		   (pr->tmp[2*cj2+3] == -1)) zj = -1;
	  else goto cbrk;
	 
	  
	  
	  int ind;
	  if((zi==1) && (zj==-1)){
            int ind = n * (cj1 + 1) + cj2 + 1;
            if (m[ind] == INFINITY) {
              m[ind] = tmpa;
              count++;
		}
		else{
		  	m[ind] = min(m[ind], tmpa);
		}
	  }
	  else if((zi==-1)&&(zj==1)){
            int ind = n * (cj2 + 1) + cj1 + 1;
            if (m[ind] == INFINITY) {
              m[ind] = tmpa;
              count++;
		}
		else{
		  	m[ind] = min(m[ind], tmpa);
		}
	  }
          
	}
	 oz->nni = min(max_nni,count);
	 
	}

      cbrk:
	*exact = 0;
      }
      break;

    default: assert(0);
    }
  }
    
  /* apply pending incremental closure now */
   
  if (*respect_closure && closure_pending)
      if (incr_closure(oz,dim,var_pending)) {
          return true;
      }
  
  return false;
}

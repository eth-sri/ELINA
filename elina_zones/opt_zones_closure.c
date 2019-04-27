/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
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


#include "opt_zones_closure.h"
#include "opt_mat.h"

/******
	Perform strengthening on independent components.
	Joins set cd with other sets which have at least
	one variable with a finite unary inequality.
******/
void strengthening_comp_zones(opt_zones_mat_t *oz,comp_list_t *cd, unsigned short int dim,int sgn){
	double *m = oz->mat;
	array_comp_list_t *acl = oz->acl;
	//char *cm = (char *)calloc(dim,sizeof(char));
	comp_list_t * cl = oz->acl->head;
	unsigned short int num_comp = acl->size;
	unsigned short int n = dim + 1;
	char * jm = (char *)calloc(num_comp,sizeof(char));
	int jc = 0;
	int l = 0;
	double * temp = (double *)calloc(dim,sizeof(double));
	/*****
		Find the sets to be merged with cd. put all
		such sets including cd into join set.
	******/
	while(cl!=NULL){
		if(cl==cd){
			cl = cl->next;
			jm[l] = 1;
			jc++;
			l++;
			continue;
		}
		int comp_size = cl->size;
		comp_t *c = cl->head;
		for(unsigned i = 0; i < comp_size; i++){
			int i1 = c->num;
			int ind = sgn==1 ? n*(i1+1) : i1+1;
			//temp[i] = m[n*(i^1) + i];
			temp[i1] = m[ind];
			if(temp[i1] != INFINITY){
				jm[l] = 1;	
				jc++;
				break;
			}
		}
		l++;
		cl=cl->next;
	}
	//free(cm);
	/*****
		Code for Incremental Initialization
	******/
	if(jc > 1 && !(oz->ti)){
		comp_list_t * ci = acl->head;
		for(unsigned short int i = 0; i < num_comp; i++){
			if(jm[i]){
				comp_list_t * cj = acl->head; 
				for(unsigned short int j = 0; j < i; j++){
					if(!jm[j]){
						cj = cj->next;
						continue;
					}
					else{
						ini_comp_relations_zones(m,ci,cj,dim);
						cj = cj->next;
					}
				}
			}
			ci = ci->next;
		}
	}
	if(jc > 1){
		/******
			If the size of join set is more than 1, then
			merge all the sets.
		******/
		cl = acl->head;
		unsigned short int nr = 0;
		comp_list_t *jl = create_comp_list();
		char *map = (char *)calloc(dim,sizeof(char)); 
		for(unsigned short int l = 0; l < num_comp; l++){
			if(jm[l]){
				union_comp_list(jl,cl,map);
				comp_list_t *temp = cl;
				cl = cl->next;
				remove_comp_list(acl,temp);
				nr++;
			}
			else{
				cl = cl->next;
			}
		}
		insert_comp_list(acl,jl);
		free(map);
	}
	free(temp);
	free(jm);		
}


// check for negative cycle

bool check_negative_cycle(opt_zones_mat_t * oz, unsigned short int dim){
	unsigned short int n = dim + 1;
	double * m = oz->mat;
	unsigned short int i;
	for(i=0; i < n; i++){
		if(m[i*n+i] < 0){
			return true;
		}
		else{
			m[i*n+i] = 0;
		}
	}
	return false;
}

/*********************
	vectorized dense closure
***********************/
void floyd_warshall_vector(double *m, unsigned short int n){
	#if defined (VECTOR)
	unsigned short int i,j,k;
	for(k=0; k < n; k++){
		// load the k-th row
		double * p2 = m + k*n;
		for(i=0; i < n; i++){
			//load the i-th row
			double *p1 = m + i*n;
			// load the first operand (m[i,k])
			double m_ik = p1[k];
			v_double_type op1 = v_set1_double(m_ik);
			//vectorized j-loop 
			for(j=0; j < n/v_length; j++){
				//load the second operand (m[k,j])
				v_double_type op2 = v_load_double(p2 + j*v_length);
				//compute the sum (m[i,k] + m[k,j])
				v_double_type op3 = v_add_double(op1,op2);
				//load (m[i,j])
				v_double_type op4 = v_load_double(p1 + j*v_length);
				//compute min (m[i,j], m[i,k] + m[k,j])
				v_double_type res = v_min_double(op3,op4);
				//store the result
				v_store_double(p1 + j*v_length,res); 
			}
			// non-vectorized clean up code
			for(j=n/v_length; j < n; j++){
				p1[j] = min(p1[j], m_ik + p2[j]);
			}
		}
	}

	#endif
}

bool closure_dense(opt_zones_mat_t * oz, unsigned short int dim){
	double * m = oz->mat;
	unsigned short int n = dim + 1;
	floyd_warshall_vector(m,n);
	oz->nni = n*n;
	return check_negative_cycle(oz,dim);
}


/***********************
	scalar dense closure
************************/
void floyd_warshall_scalar(double *m, unsigned short int n){
	unsigned short int i,j,k;
	for(k=0; k < n; k++){
		//load the k-th row
		double * p2 = m + k*n;
		for(i=0; i < n; i++){
			//load the i-th row
			double *p1 = m + i*n;
			//load the first operand (m[i,k])
			double m_ik = p1[k];
			for(j=0; j < n; j++){
				// load the second opearnd (m[k,j])
				double m_kj = p2[j];
				//compute min(m[i,j],m[i,k] + m[k,j])
				p1[j] = min(p1[j],m_ik + m_kj);
			}
		}
	}
	
}

bool closure_dense_scalar(opt_zones_mat_t *oz, unsigned short int dim){
	double *m = oz->mat;
	unsigned short int n = dim + 1;
	floyd_warshall_scalar(m,n);
	oz->nni = n*n;
	return check_negative_cycle(oz,dim);
}

/**************
	 Decomposition based Closure
***************/
int update_bounds_comp(double *m, unsigned short int * ca, 
			unsigned short int comp_size, unsigned short int n){
	unsigned short int i,j;
	int count = 0;
	for(i=0; i < comp_size; i++){
		unsigned short int i1 = ca[i]+1;
		bool flag1 = (m[n*i1]==INFINITY);
		bool flag2 = (m[i1]==INFINITY);
		for(j=0; j < comp_size; j++){
			unsigned short int j1 = ca[j]+1;
			m[n*i1] = min(m[n*i1],m[n*i1 + j1] + m[n*j1]);
			m[i1] = min(m[i1],m[j1] + m[n*j1+i1]);
		}
		if(flag1 && (m[n*i1]!=INFINITY)){
			count++;
		}
		if(flag2 && (m[i1]!=INFINITY)){
			count++;
		}
	}
	return count;
}


bool strengthening_single_comp_zones(opt_zones_mat_t * oz, comp_list_t * cl, unsigned short int dim){
	double * m = oz->mat;
	int count = oz->nni;
	unsigned short int *ca = to_sorted_array(cl,dim);
	unsigned short int i,j;
	unsigned short int comp_size = cl->size;
	unsigned short int n = dim+1;
	for(i=0; i < comp_size; i++){
		unsigned short int i1 = ca[i] + 1;
		for(j=0; j < comp_size; j++){
			unsigned short int j1 = ca[j] + 1;
			if(m[n*i1+j1]==INFINITY){
				m[n*i1+j1] = m[n*i1] + m[j1];
				count++;
			}
			else{
				m[n*i1+j1] = min(m[n*i1+j1],m[n*i1] + m[j1]);
			}
		}
	}
	
	m[0] = 0;
	for(i = 0; i < comp_size; i++){
		int i1 = ca[i]+1;
		int ind = i1*n + i1;
		if(m[ind] < 0){
			return true;
		}
		else{
			m[ind] = 0;
		}
	}
	free(ca);
	oz->nni = min(n*n,count);
	return false;
}

void strengthening_incr_init_zones(opt_zones_mat_t *oz, array_comp_list_t *acl, char *map , unsigned short int dim){
	unsigned short int num_comp = acl->size;
	unsigned short int i, j;
	double *m = oz->mat;
	/***********************
		Incremental initilaization
	************************/
	comp_list_t * ci = acl->head;
	for(i = 0; i < num_comp; i++){
		if(map[i]){
			comp_list_t * cj = acl->head; 
			for( j = 0; j < i; j++){
				if(!map[j]){
					cj = cj->next;
					continue;
				}
				else{
					ini_comp_relations_zones(m,ci,cj,dim);
					cj = cj->next;
				}
			}
		}
		ci = ci->next;
	}

}

void strengthening_inter_comp_zones(opt_zones_mat_t *oz, array_comp_list_t *acl, char *map , unsigned short int dim){
	unsigned short int num_comp = acl->size;
	unsigned short int i,j;
	unsigned short int i1,j1;
	double *m = oz->mat;
	comp_list_t * cli = acl->head; 
	unsigned short int n = dim + 1;
	unsigned short int ** ca_arr= (unsigned short int **)malloc(num_comp*sizeof(unsigned short int *));
	unsigned short int k = 0;
	while(cli!=NULL){
		ca_arr[k] = to_sorted_array(cli,dim);
		cli = cli->next;
		k++;
	}
	cli = acl->head;
	for(i=0; i < num_comp; i++){
		if(map[i]){
			unsigned short int comp_size_i = cli->size;
			unsigned short int * ca_i = ca_arr[i];
			comp_list_t *clj = acl->head;
			for(j=0; j < i; j++){
				if(!map[j]){
					clj = clj->next;
					continue;
				}
				else{
					unsigned short int comp_size_j = clj->size;
					unsigned short int * ca_j = ca_arr[j];
					for(i1=0; i1 < comp_size_i; i1++){
						unsigned short int i2 = ca_i[i1]+1;
						int ind1 = n*i2;
						for(j1=0; j1 < comp_size_j; j1++){
							unsigned short int j2 = ca_j[j1]+1;
							int ind2 = j2;
							int ind =  n*i2+j2;
							m[ind] = min(m[ind],m[ind1]+m[ind2]);
						}
					}
					clj = clj->next;
				}
			}
		}
		cli = cli->next;
	}
	
	for(k=0; k < num_comp; k++){
		free(ca_arr[k]);
	}
	free(ca_arr);
	return;
}


bool strengthening_intra_comp_zones(opt_zones_mat_t * oz, unsigned short int dim){
	array_comp_list_t *acl = oz->acl;
	comp_list_t * cl = acl->head;
	double *m = oz->mat;
	if(m[0] < 0){
		return true;
	}
	else{
		m[0] = 0;
	}
	while(cl!=NULL){
		if(strengthening_single_comp_zones(oz,cl,dim)){
			return true;
		}
		cl = cl->next;
	}
	return false;
}


bool strengthening_zones(opt_zones_mat_t *oz, unsigned short int dim){
	//printf("strengthening start\n");
	//print_mat(oz,dim);
	//fflush(stdout);
	double *m = oz->mat;
	array_comp_list_t *acl = oz->acl;
	int s = 0;
	unsigned short int n = dim+1;
	//int count = oz->nni;
	
	comp_list_t * cl = acl->head;
	unsigned short int num_comp = acl->size;
	char * jm = (char *)calloc(num_comp,sizeof(char));
	
	/****
		Corresponding to each component, store the index of set containing it. 
	*****/
	unsigned short int ** ca_arr = (unsigned short int **)malloc(num_comp*sizeof(unsigned short int *));
	unsigned short int k,k1;
	bool jc = false;
	
	/***************
		Strengthening for entries in same component
	**************/	
	for(k=0; k < num_comp; k++){
		ca_arr[k] = to_sorted_array(cl,n);
		/*unsigned short int * ca = ca_arr[k];
		unsigned short int i,j;
		unsigned short int comp_size = cl->size;
		for(i=0; i < comp_size; i++){
			unsigned short int i1 = ca[i]+1;
			for(j=0; j < comp_size; j++){
				unsigned short int j1 = ca[j]+1;
				if(m[n*i1+j1]==INFINITY){
					m[n*i1+j1] = m[n*i1] + m[j1];
					count++;
				}
				else{
					m[n*i1+j1] = min(m[n*i1+j1],m[n*i1] + m[j1]);
				}
			}
		}*/
		
		cl = cl->next;	
	}
	
	/***************
		Check which components will get connected
	**************/
	cl = acl->head;
        for(k=0; k < num_comp; k++){
		unsigned short int comp_size = cl->size;
		comp_list_t * cl1 = acl->head;
		unsigned short int *ca = ca_arr[k];
		for(k1=0; k1 < num_comp; k1++){
			if(k1<=k){
				cl1= cl1->next;
				continue;
			}
			unsigned short int *ca1 = ca_arr[k1];
			unsigned short int comp_size1 = cl1->size;
			unsigned short int i,j;
			for(i=0; i < comp_size; i++){
				unsigned short int i1 = ca[i]+1;
				double op1 = m[i1*n];
				double op3 = m[i1];	
				for(j=0; j < comp_size1; j++){
					unsigned short int j1 = ca1[j]+1;
					//ini_binary_relation_zones(m,i1,j1,dim);
					double op2 = m[j1];
					double op4 = m[n*j1];
					if((op1!=INFINITY) && (op2!=INFINITY)){
						jm[k] = 1;
						jm[k1] = 1;
						jc = true;
						goto diff_comp_loop;
						//m[i1*n + j1] = op1 + op2;
						//count++;
					}
					if((op3!=INFINITY)&&(op4!=INFINITY)){
						jm[k] = 1;
						jm[k1] = 1;
						jc = true;
						goto diff_comp_loop;
						//m[n*j1+i1] = op3 + op4;
						//count++;
					}
				}
			}
			diff_comp_loop:
				cl1 = cl1->next;
		}
		cl = cl->next;
	}
	
	
	
        /*****
		Code for Incremental Initialization
	******/
	if(jc){
		
		comp_list_t * ci = acl->head;
		for(int i = 0; i < num_comp; i++){
			if(jm[i]){
				comp_list_t * cj = acl->head; 
				for(int j = 0; j < i; j++){
					if(!jm[j]){
						cj = cj->next;
						continue;
					}
					else{
						ini_comp_relations_zones(m,ci,cj,dim);
						cj = cj->next;
					}
				}
			}
			ci = ci->next;
		}
	}
	
	
	for(k=0; k < num_comp; k++){
		free(ca_arr[k]);
	}
	
	
	free(ca_arr);
	if(jc){
		/******
			If the size of join set is more than 1, then
			merge all the sets.
		******/
		cl = acl->head;
		int nr = 0;
		comp_list_t *jl = create_comp_list();
		char *map = (char *)calloc(n,sizeof(char)); 
		for(unsigned short int l = 0; l < num_comp; l++){
			if(jm[l]){
				union_comp_list(jl,cl,map);
				comp_list_t *temp = cl;
				cl = cl->next;
				remove_comp_list(acl,temp);
				nr++;
			}
			else{
				cl = cl->next;
			}
		}
		insert_comp_list(acl,jl);
		free(map);
	}
	
	free(jm);
	if(m[0] < 0){
		return true;
	}
	else{
		m[0] = 0;
	}
	cl = acl->head;
	while(cl!=NULL){
		if(strengthening_single_comp_zones(oz,cl,dim)){
			return true;
		}
		cl = cl->next;
	}
	//printf("strengthening end\n");
	//print_mat(oz,dim);
	//fflush(stdout);
	return false;
}

int calculate_zones_comp_sparsity(opt_zones_mat_t *oz, comp_list_t * cl, unsigned short int dim){
	unsigned short int comp_size = cl->size;
	int count = 0;
	unsigned short int n = dim + 1;
	unsigned short int i,j;
	double *m = oz->mat;
	unsigned short int *ca = to_sorted_array(cl,n);
	for(i=0; i < comp_size; i++){
		unsigned short int i1 = ca[i]+1;
		double * pi = m + i1*n;
		for(j=0; j < comp_size; j++){
			unsigned short int j1 = ca[j]+1;
			if(i1!=j1){
				if(pi[j1]!=INFINITY){
					count++;
				}
			}
		}
		if(m[n*i1]!=INFINITY){
			count++;
		}
		if(m[i1]!=INFINITY){
			count++;
		}
	} 
	free(ca);
	return count;
}

void floyd_warshall_comp_zones(opt_zones_mat_t *oz,comp_list_t *cl, unsigned short int dim){
	
	unsigned short int comp_size = cl->size;
	int size = comp_size*comp_size;
	unsigned short int i,j;
	unsigned short int n = dim + 1;
	unsigned short int *ca = to_sorted_array(cl,n);
	double * m = oz->mat;
	int ind = 0;
	double * tmp = (double*)malloc(size*sizeof(double));
	for(i=0; i < comp_size; i++){
		unsigned short int i1 = ca[i]+1;
		for(j=0; j < comp_size;j++){
			unsigned short int j1 = ca[j]+1;
			tmp[ind] = m[n*i1+j1];
			ind++;
		}
	}
	
	#if defined (VECTOR)
		floyd_warshall_vector(tmp,comp_size);
	#else
		floyd_warshall_scalar(tmp,comp_size);
	#endif
	
	ind = 0;
	for(i=0; i < comp_size; i++){
		unsigned short int i1 = ca[i]+1;
		for(j=0; j < comp_size;j++){
			unsigned short int j1 = ca[j]+1;
			m[n*i1+j1] = tmp[ind];
			ind++;
		}
	}
	
	free(ca);
	free(tmp);
}

void compute_sparse_index(double *m, unsigned short int *ca, 
			  unsigned short int comp_size, unsigned short int *index1,
			  unsigned short int * index2, unsigned short int k, unsigned short int dim){
	
	unsigned short int s1=0, s2=0;
	unsigned short int i,j;
	unsigned short int n = dim + 1;
	unsigned short int k1 = k+1;
	//load the k-th row
	double * pk = m + n*k1;
	for(i=0; i < comp_size; i++){
		unsigned short int i1 = ca[i]+1;
		double * pi = m + i1*n;
		if(pi[k1]!=INFINITY){
			index1[s1+1] = i1;
			s1++;
		}
		if(pk[i1]!=INFINITY){
			index2[s2+1] = i1;
			s2++;
		}
	}
	index1[0] = s1;
	index2[0] = s2;
}

void closure_comp_sparse(opt_zones_mat_t *oz, unsigned short int dim){
	//printf("closure input\n");
	//print_mat(oz,dim);
	//fflush(stdout);
	unsigned short int n = dim+1;
	double *m = oz->mat;
	array_comp_list_t * acl = oz->acl;
	unsigned short int num_comp = acl->size;
	unsigned short int i,j,k,l;
	comp_list_t * cl = acl->head;
	unsigned short int *index1 = (unsigned short int *)calloc(dim+1,sizeof(unsigned short int));
	unsigned short int *index2 = (unsigned short int *)calloc(dim+1,sizeof(unsigned short int)); 
	int count = oz->nni;
	for(l=0; l < num_comp; l++){
		unsigned short int comp_size = cl->size;
        if(comp_size==0){
            cl = cl->next;
            continue;
        }
		/******
			Calculate precise sparsity of each component set.
			If it is less than threshold use dense Floyd Warshall,
			otherwise use decomposition based Floyd Warshall.
	    	*******/
		int nni = calculate_zones_comp_sparsity(oz,cl,dim);
		int comp_mat_size = comp_size*comp_size; 
		count = count + comp_mat_size - nni;
	    	double sparsity = 1- ((double)(nni/comp_mat_size));
	    	if(sparsity < zone_sparse_threshold){
				floyd_warshall_comp_zones(oz,cl,dim);
				
				cl = cl->next;
				
				continue;
	    	}
		
		unsigned short int * ca = to_sorted_array(cl,n);
		for(k=0; k < comp_size; k++){
			unsigned short int k1 = ca[k]+1;
			//load the k-th row
			double * p2 = m + k1*n;
			compute_sparse_index(m,ca,comp_size,index1,index2,ca[k],dim);
			unsigned short int n_ik = index1[0];
			unsigned short int n_kj = index2[0];
			for(i=0; i < n_ik; i++){
				unsigned short int i1 = index1[i+1];
				//load the i-th row
				double *p1 = m + i1*n;
				//load the first operand (m[i,k])
				double m_ik = p1[k1];
				for(j=0; j < n_kj; j++){
					unsigned short int j1 = index2[j+1];
					// load the second opearnd (m[k,j])
					double m_kj = p2[j1];
					if(p1[j1]==INFINITY){
						p1[j1] = m_ik + m_kj;
						count++;
					}
					else{
					//compute min(m[i,j],m[i,k] + m[k,j])
						p1[j1] = min(p1[j1],m_ik + m_kj);
					}
				}
				
			}
		}
		/***************
			Update the variable bounds
		**************/	
		nni = update_bounds_comp(m,ca,comp_size,n);
		count = count + nni;
		free(ca);
		
		cl = cl->next;
	}

	
	free(index1);
	free(index2);
	oz->nni = count;
    	return;
	// similar to strengthening for octagons 
	//bool res = strengthening_zones(oz,dim);
        //printf("closure output\n");
	//print_mat(oz,dim);
	//fflush(stdout);
        //return res;
}

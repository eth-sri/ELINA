/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
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


#include "opt_oct_closure_comp_sparse.h"


/******
	Perform strengthening on independent components.
	Joins set cd with other sets which have at least
	one variable with a finite unary inequality.
******/
void strengthening_comp_list(opt_oct_mat_t *oo,comp_list_t *cd, unsigned short int dim){
	double *m = oo->mat;
	array_comp_list_t *acl = oo->acl;
	//char *cm = (char *)calloc(dim,sizeof(char));
	comp_list_t * cl = acl->head;
	unsigned short int num_comp = acl->size;
	char * jm = (char *)calloc(num_comp,sizeof(char));
	int jc = 0;
        cl = oo->acl->head;
	int l = 0;
	double * temp = (double *)calloc(2*dim,sizeof(double));
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
		for(unsigned i = 0; i < 2*comp_size; i++){
			int i1 = (i%2==0)? 2*c->num : 2*c->num+1;
			int ind = i1 + ((((i1^1) + 1)*((i1^1) + 1))/2);
			//temp[i] = m[n*(i^1) + i];
			temp[i1] = m[ind];
			if(temp[i1] != INFINITY){
				jm[l] = 1;
				
				jc++;
				break;
			}
			if(i%2==1){
				c = c->next;
			}
		}
		l++;
		cl=cl->next;
	}
	//free(cm);
	/*****
		Code for Incremental Initialization
	******/
	if(jc > 1 && !(oo->ti)){
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
						ini_comp_relations(m,ci,cj,dim);
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
		int nr = 0;
		comp_list_t *jl = create_comp_list();
		char *map = (char *)calloc(dim,sizeof(char)); 
		for(int l = 0; l < num_comp; l++){
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


bool strengthning_int_comp_sparse(opt_oct_mat_t * oo,  unsigned short int * ind1, double *temp, int n){
	double *m = oo->mat;
	array_comp_list_t *acl = oo->acl;
	int count = oo->nni;
	int s = 0;
	char *cm = (char *)calloc(n/2, sizeof(char));
	comp_list_t *cl = acl->head;
	unsigned short int num_comp = acl->size;
	char * jm = (char *)calloc(num_comp,sizeof(char));
	/****
		Corresponding to each component, store the index of set containing it. 
	*****/
	for(int l = 0; l < num_comp; l++){
		comp_t *c = cl->head;
		while(c!=NULL){
			unsigned short int num = c->num;
			cm[num] = l;
			c = c->next;
		}
		cl = cl->next;
	}
	int jc = 0;
	for(unsigned i = 0; i < n; i++){
		int ind = i + ((((i^1) + 1)*((i^1) + 1))/2);
		temp[i] = ceil(m[ind]/2);
		if(temp[i] != INFINITY){
			ind1[s+1] = i;
			/***
				look for the index of set that contains the component i 
				and then add that set to the join set
			***/
			int cn = i/2;
			int cln = cm[cn];
			if(!jm[cln]){
				jm[cln] = 1;
				jc++;
			}
			s++;
		}
	}
	free(cm);
	/*****
		Code for Incremental Initialization
	******/
	if(jc > 1){
		
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
						ini_comp_relations(m,ci,cj,n/2);
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
		int nr = 0;
		comp_list_t *jl = create_comp_list();
		char *map = (char *)calloc(n/2,sizeof(char)); 
		for(int l = 0; l < num_comp; l++){
			if(jm[l]){
				union_comp_list(jl,cl,map);
				comp_list_t * temp = cl;
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
	ind1[0] = s;
	/******
		Perform sparse strengthening
	******/
	for(unsigned i = 0; i < ind1[0]; i++){
		unsigned i1 = ind1[i + 1];
		double t1 = temp[i1];
		for(unsigned j = 0; j < ind1[0];j++){
			unsigned j1 = ind1[j + 1];
			double t2 = temp[j1];
			int ind = j1 + ((((i1^1) + 1)*((i1^1) + 1))/2);
			//m[n*(i1^1) + j1] = min(m[n*(i1^1) + j1], t1 + t2);
			if(m[ind]!=INFINITY){
				m[ind] = min(m[ind], t1 + t2);
			}
			else{
				m[ind] = t1 + t2;
				count++;
			}
		}	
	}
	oo->nni = count;
	/******
		Check for negative cycle
	*******/
	cl = oo->acl->head;
	
	while(cl != NULL){
		comp_t * c = cl->head;
		unsigned short int comp_size = cl->size;
		for(unsigned short int i = 0; i < 2*comp_size; i++){
			int i1 = (i%2==0)? 2*c->num : 2*c->num+1;
			int ind = i1 + (((i1+1)*(i1+1))/2);
			if(m[ind] < 0){
				
				return true;
			}
			else{
				m[ind] = 0;
			}
			if(i1%2==1){
				c = c->next;
			}
		}
		cl = cl->next;
	}
	return false;
}


bool strengthning_comp_sparse(opt_oct_mat_t *oo, unsigned short int * ind1, double *temp, int n){
	double *m = oo->mat;
	array_comp_list_t *acl = oo->acl;
	int s = 0;
	int count = oo->nni;
	char *cm = (char *)calloc(n/2,sizeof(char));
	
	comp_list_t * cl = acl->head;
	unsigned short int num_comp = acl->size;
	char * jm = (char *)calloc(num_comp,sizeof(char));
	/****
		Corresponding to each component, store the index of set containing it. 
	*****/
	for(int l = 0; l < num_comp; l++){
		comp_t *c = cl->head;
		while(c!=NULL){
			unsigned short int num = c->num;
			cm[num] = l;
			c = c->next;
		}
		cl = cl->next;
	}
	int jc = 0;
        cl = oo->acl->head;
        while(cl!=NULL){
		int comp_size = cl->size;
		comp_t *c = cl->head;
		for(unsigned i = 0; i < 2*comp_size; i++){
			int i1 = (i%2==0)? 2*c->num : 2*c->num+1;
			int ind = i1 + ((((i1^1) + 1)*((i1^1) + 1))/2);
			//temp[i] = m[n*(i^1) + i];
			temp[i1] = m[ind];
			if(temp[i1] != INFINITY){
				ind1[s+1] = i1;
				/***
					look for the index of set that contains the component i 
					and then add that set to the join set
				***/
				int cn = i1/2;
				int cln = cm[cn];
				if(!jm[cln]){
					jm[cln] = 1;
					jc++;
				}
				s++;
			}
			if(i%2==1){
				c = c->next;
			}
		}
		cl = cl->next;
	}
	free(cm);
        /*****
		Code for Incremental Initialization
	******/
	if(jc > 1){
		
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
						ini_comp_relations(m,ci,cj,n/2);
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
		int nr = 0;
		comp_list_t *jl = create_comp_list();
		char *map = (char *)calloc(n/2,sizeof(char)); 
		for(int l = 0; l < num_comp; l++){
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
	
	ind1[0] = s;
	/******
		Perform sparse strengthening
	******/
	for(unsigned i = 0; i < ind1[0]; i++){
		unsigned i1 = ind1[i + 1];
		double t1 = temp[i1];
		for(unsigned j = 0; j < ind1[0];j++){
			unsigned j1 = ind1[j + 1];
			if(j1 > (i1|1)){
				continue;
			}
			double t2 = temp[j1];
			int ind = j1 + ((((i1^1) + 1)*((i1^1) + 1))/2);
			//m[n*(i1^1) + j1] = min(m[n*(i1^1) + j1], (t1 + t2)/2);
			if(m[ind]!=INFINITY){
				m[ind] = min(m[ind], (t1 + t2)/2);
			}
			else{
				m[ind] = (t1+t2)/2;
				count++;
			}
		}	
	}
	oo->nni = count;
	/******
		Check for negative cycle
	*******/
	cl = oo->acl->head;
	
	while(cl != NULL){
		comp_t * c = cl->head;
		unsigned short int comp_size = cl->size;
		for(unsigned short int i = 0; i < 2*comp_size; i++){
			int i1 = (i%2==0)? 2*c->num : 2*c->num+1;
			int ind = i1 + (((i1+1)*(i1+1))/2);
			if(m[ind] < 0){
				
				return true;
			}
			else{
				m[ind] = 0;
			}
			if(i1%2==1){
				c = c->next;
			}
		}
		cl = cl->next;
	}
	
	return false;
}

/******
	Compute the index containing location of finite entries
	in 2k and 2k+1-th row and column for the given component set/
******/
void compute_index_comp_sparse(double *m, unsigned short int *ca, unsigned short int comp_size, unsigned short int *index1, unsigned short int *index2, unsigned short int k, int dim){
    int n = 2*dim;
    index1[0] = 0;
    index1[n + 1] = 0;
    index2[0] = 0;
    index2[n + 1] = 0;
    int s1 = 0, s2 = 0;
    int ind;
    //comp_t *ci = cl->head;
    //unsigned short int comp_size = cl->size;
   // unsigned short int * ca = to_sorted_array(cl,dim);
    for(int i = 0; i < 2*comp_size; i++){
	int i1 = (i%2==0) ? 2*ca[i/2] : 2*ca[i/2]+1;
	if(i1>=(2*k+2)){
		ind = (2*k) + (((i1 + 1)*(i1 + 1))/2); 
		//if(m[n*i + 2*k] != INFINITY){
		if(m[ind] != INFINITY){
		    index2[s1 + 1] = i1;
		    s1++;
		}
		ind = ((2*k)^1) + (((i1 + 1)*(i1 + 1))/2);  
		//if(m[n*i + ((2*k)^1)] != INFINITY){
		if(m[ind] != INFINITY){
		    index2[n + s2 + 2] = i1;
		    s2++;
		}
	}
    }
    index2[0] = s1;
    index2[n + 1] = s2;

    s1 = 0, s2 = 0;
    for(int j  = 0; j < 2*comp_size; j++){
	int j1 = (j%2==0) ? 2*ca[j/2] : 2*ca[j/2]+1;
	if(j1 < (2*k)){
		ind = j1 + ((((2*k) + 1)*((2*k) + 1))/2);
		//if(m[n*2*k + j] != INFINITY){
		if(m[ind] != INFINITY){
		    index1[s1 + 1] = j1;
		    s1++;
		}
		ind = j1 + (((((2*k)^1) + 1)*(((2*k)^1) + 1))/2);
		//if(m[n*((2*k)^1) + j] != INFINITY){
		if(m[ind] != INFINITY){
		    index1[n + s2 + 2] = j1;
		    s2++;
		}
	}
    }
    index1[0] = s1;
    index1[n + 1] = s2;
    //free(ca);
}



void print_index(unsigned short int *ind,int dim){
	int n = 2*dim;
	int s = ind[0];
	fprintf(stdout,"Size is %d\n",s);
	for(int i = 0; i < s; i++){
		fprintf(stdout,"%d ",ind[i+1]);
	}
	fprintf(stdout,"\n");
	s = ind[n+1];
	fprintf(stdout,"Size is %d\n",s);
	for(int i = 0; i < s; i++){
		fprintf(stdout,"%d ",ind[n + i +2]);
	}
	fprintf(stdout,"\n");
	fflush(stdout);
}

/*******
	Calculate sparsity of component set.
*******/
double calculate_comp_sparsity(opt_oct_mat_t *oo, comp_list_t *cl, int dim){
	unsigned short int *ca = to_sorted_array(cl,dim);
	unsigned short int comp_size = cl->size;
	int count = 0;
	double *m = oo->mat;
	int size = 2*comp_size*(comp_size+1);
	for(int i = 0; i < 2*comp_size; i++){
		int i1 = (i%2==0) ? 2*ca[i/2] : 2*ca[i/2]+1;
		for(int j = 0; j < 2*comp_size;j++){
			int j1 = (j%2==0) ? 2*ca[j/2] : 2*ca[j/2] + 1;
			if(j1 > (i1 | 1)){
				break;
			}
			int ind = j1 + (((i1+1)*(i1+1))/2);
			if(m[ind]!=INFINITY){
				count++;
			}
		}
	}
	free(ca);
	return 1- ((double)(count/(double)size));
}

bool floyd_warshall_comp_dense(opt_oct_mat_t * oo, comp_list_t * cl, int dim){
	unsigned short int comp_size = cl->size;
	double *m = oo->mat;
	int size = 2*comp_size*(comp_size+1);
	opt_oct_mat_t * ot = opt_hmat_alloc(size);
	double * temp = ot->mat;
	ot->is_dense = true;
	/******
		Copy the component set to temporary dense matrix.
	******/
	unsigned short int *ca = to_sorted_array(cl,dim);
	int ind = 0;
	for(int i = 0; i < 2*comp_size; i++){
		int i1 = (i%2==0) ? 2*ca[i/2] : 2*ca[i/2]+1;
		for(int j = 0; j < 2*comp_size;j++){
			int j1 = (j%2==0) ? 2*ca[j/2] : 2*ca[j/2] + 1;
			if(j1 > (i1 | 1)){
				break;
			}
			int ind1 = j1 + (((i1+1)*(i1+1))/2);
			temp[ind] = m[ind1];
			ind++;
		}
	}
	bool flag = is_int_flag ? true : false;
	double * temp1 = (double *)malloc(2*comp_size*sizeof(double));
	double * temp2 = (double *)malloc(2*comp_size*sizeof(double));
	/******
		Apply Floyd-Warshall on temporary matrix.
	*******/
	#if defined(VECTOR)
			floyd_warshall_dense(ot,temp1,temp2,comp_size, flag);
	#else
			floyd_warshall_dense_scalar(ot,temp1,temp2,comp_size, flag);
	#endif
	free(temp1);
	free(temp2);
	ind = 0;
	/*******
		Copy back the result of applying strong closure
		on temporary matrix to component set.
	*******/
	for(int i = 0; i < 2*comp_size; i++){
		int i1 = (i%2==0) ? 2*ca[i/2] : 2*ca[i/2]+1;
		for(int j = 0; j < 2*comp_size;j++){
			int j1 = (j%2==0) ? 2*ca[j/2] : 2*ca[j/2] + 1;
			if(j1 > (i1 | 1)){
				break;
			}
			int ind1 = j1 + (((i1+1)*(i1+1))/2);
			m[ind1] = temp[ind];
			ind++;
		}
	}
	free(ca);
	free_array_comp_list(ot->acl);
	opt_hmat_free(ot);
}

bool strong_closure_comp_sparse(opt_oct_mat_t *oo, double *temp1, double *temp2, unsigned short int *index1, unsigned short int *index2, int dim, bool is_int){
    double *m = oo->mat;
    array_comp_list_t *acl = oo->acl;
    int count = oo->nni;
    int n = 2*dim; 
    for(int i = 0; i < n; i++){
        temp1[i] = 0;
        temp2[i] = 0;
    }
    unsigned short int num_comp = acl->size;
    comp_list_t * cl = acl->head;
    
    for(int l = 0; l < num_comp; l++){
	    /******
			Calculate precise sparsity of each component set.
			If it is less than threshold use dense Floyd Warshall,
			otherwise use decomposition based Floyd Warshall.
	    *******/
	    double sparsity = calculate_comp_sparsity(oo,cl,dim);
	    if(sparsity < sparse_threshold){
			floyd_warshall_comp_dense(oo,cl,dim);
			cl = cl->next;
			continue;
	    }
	    unsigned short int comp_size = cl->size;
	    unsigned short int * ca = to_sorted_array(cl,dim);
	    //comp_t * ck = cl->head;
	    /******
			Floyd-Warshall step for each set independently
	    ******/
	    for(int k = 0; k < comp_size; k++){
		//Compute index at start of iteration
		unsigned short int k1 = ca[k];
		//ck = ck->next;
		/******
			Compute index for k-th iteration
		*******/
		compute_index_comp_sparse(m,ca, comp_size, index1, index2, k1, dim);
		
		int s1 = index1[0],s2 = index1[n + 1];
		int s3 = index2[0], s4 = index2[n + 1];
		int pos1 = ((2*k1)^1) + ((((2*k1) + 1)*((2*k1) + 1))/2);
		//int pos2 = matpos2((2*k)^1, 2*k);
		int pos2 = (2*k1) + (((((2*k1)^1) + 1)*(((2*k1)^1) + 1))/2);
		//Compute (2k) and (2k+1)-th row and column, update the index
		if(m[pos1]!= INFINITY){
			  for(int i = 0; i < s3;i++){
				int i2 = index2[i + 1];
				int i1 = i2;
				//int ind2 = n*i2 + ((2*k)^1);
				//int ind1 = n*i1 + (2*k);
				int ind1 = (2*k1) + (((i1 + 1)*(i1 + 1))/2);
				int ind2 = ((2*k1)^1) + (((i2 + 1)*(i2 + 1))/2);

				if(m[ind2]!= INFINITY){
					m[ind2] = min(m[ind2],  m[pos1] + m[ind1]);
				}
				else{
					m[ind2] =  m[pos1] + m[ind1];
					//index1[m*i2 + index1[m*i2] + 1] = ((2*k)^1);
					index2[n + s4 + 2] = i2;
					//index1[m*i2]++;
					s4++;
					count++;
				}
				//temp2[i2] = m[ind2];
			}
			index2[n + 1] = s4;
		}
		
		for(int i = 2*k1+2; i < n; i++){
			int ind = ((2*k1)^1) + (((i+1)*(i+1))/2);
			//temp2[i] = m[n*i + ((2*k)^1)];
			temp2[i] = m[ind];
		}
	

		if(m[pos2] != INFINITY){
		
			for(int i = 0; i < s4; i++){
				int i2 = index2[n + i + 2];
				int i1 = i2;
				//int ind2 = n*i2 + ((2*k)^1);
				//int ind1 = n*i1 + (2*k);
				int ind1 = (2*k1) + (((i1 + 1)*(i1 + 1))/2);
				int ind2 = ((2*k1)^1) + (((i2 + 1)*(i2 + 1))/2);
				if(m[ind1] != INFINITY){
					m[ind1] = min(m[ind1], m[pos2] + m[ind2]);
				}
				else{
					m[ind1] = m[pos2] + m[ind2];
					//index1[m*i1 + index1[m*i1] + 1] = 2*k;
					index2[s3 + 1] = i1;
					//index1[m*i1]++;
					s3++;
					count++;
				}
				//temp1[i1] = m[ind1];
			}

			index2[0] = s3;
		}
		
		for(int i = 2*k1+2; i < n; i++){
			int ind = (2*k1) + (((i+1)*(i+1))/2);
			//temp1[i] = m[n*i + (2*k)];
			temp1[i] = m[ind];
		}

		if(m[pos2] != INFINITY){

			for(int j = 0; j < s1; j++){
				//int ind4 = get_index(n, 2*k,j);
				//int j1 = index1[m*(2*k) + j + 1];
				int j1 = index1[j + 1];
				int ind1 = j1 + ((((2*k1) + 1)*((2*k1) + 1))/2);
				int ind2 = j1 + (((((2*k1)^1) + 1)*(((2*k1)^1) + 1))/2);
				if(m[ind2] != INFINITY ){
					m[ind2] = min(m[ind2], m[pos2] + m[ind1]);
				}
				else{
					m[ind2] =  m[pos2] + m[ind1];
					index1[n + s2 + 2] = j1;
					//index2[m*j1 + index2[m*j1] + 1] = ((2*k)^1);
					s2++;
					count++;
					//index2[m*j1]++;
				} 
			}
			index1[n + 1] = s2;
		}

		if(m[pos1] != INFINITY){

			for(int j = 0; j < s2; j++){
				//int ind4 = get_index(n, 2*k,j);
				//int j1 = index1[m*((2*k)^1) + j + 1];
				int j1 = index1[n + j + 2];
				int ind1 = j1 + ((((2*k1) + 1)*((2*k1) + 1))/2);
				int ind2 = j1 + (((((2*k1)^1) + 1)*(((2*k1)^1) + 1))/2);
				//if(m[ind2] != std::numeric_limits<double>::infinity(
				if(m[ind1] != INFINITY ){
					m[ind1] = min(m[ind1], m[pos1] + m[ind2]);
				}
				else{
					m[ind1] =  m[pos1] + m[ind2];
					index1[s1 + 1] = j1;
					//index2[m*j1 + index2[m*j1] + 1] = (2*k);
					s1++;
					count++;
					//index2[m*j1]++;
				}
			}
			index1[0] = s1;
		}
	        
		//This is the Main Loop Divided into four parts
		//First Part through 2k+1
		int ind1_k = index1[0];
		int ind2_k = index1[n + 1];
		int ind3_k = index2[0];
		int ind4_k = index2[n + 1];
		for(int i = 0; i < ind1_k; i++){
			int i1 = index1[i + 1];
			int i2 = (i1%2==0) ? (i1 + 1): i1;
			int br = i2 < 2*k1 ? i2 : 2*k1 - 1;
			int ind1 = i1 + ((((2*k1) + 1)*((2*k1) + 1))/2);
			//double t1 = m[n*(2*k) + i1];
			double t1 = m[ind1];
			//double t2 = m[n*((2*k)^1) + (i^1)];
			//int j2 = (j/2)*2;
			for(int j = 0;j < ind2_k ; j++){
				//int ind2 = get_index(k,j);
		        	//int j1 = index1[m*((2*k)^1) + j + 1];
				int j1 = index1[n + j + 2];
				if(j1 > br){
					break;
					//continue;
				}
				int ind2 = j1 + (((((2*k1)^1) + 1)*(((2*k1)^1) + 1))/2);
				//double op1 = t1 + m[n*((2*k)^1) + j1];RDTSC(end);
        
				double op1 = t1 + m[ind2];
				//double op2 = t2 + m[n*(2*k) + j];
				//double op3 = min(op1, op2);
				int ind3 = j1 + ((((i1^1) + 1)*((i1^1) + 1))/2);
				//m[n*(i1^1) + j1] = min(m[n*(i1^1) + j1],op1 );
				if(m[ind3]!=INFINITY){
					m[ind3] = min(m[ind3],op1 );
				}
				else{
					m[ind3] = op1;
					count++;
				}
			}
			//for(int j = 0; j < index2[m*2*k]; j++){
			for(int j = 0;j < ind3_k; j++){
				int j1 = index2[j + 1];
				if(j1>i2){
					 break;
					//continue;
			    	}
				double op1 = t1 + temp1[j1];
				int ind3 = (j1^1) + ((((i1^1) + 1)*((i1^1) + 1))/2);
				//m[n*(i1^1) + (j1^1)] = min(m[n*(i1^1) + (j1^1)],op1 );
				if(m[ind3]!=INFINITY){
					m[ind3] = min(m[ind3],op1 );
				}
				else{
					m[ind3] = op1;
					count++;
				}
			}
			//}
		}
		
		//Second Part through 2k
		for(int i = 0; i < ind2_k; i++){
		    int i1 = index1[n + i + 2];
		    int i2 = (i1%2==0) ? (i1 + 1): i1;
		    int br = i2 < 2*k1 ? i2 : 2*k1 - 1;
		    //double t1 = m[n*(2*k) + i1];
		    int ind1 = i1 + (((((2*k1)^1) + 1)*(((2*k1)^1) + 1))/2);
		    //double t2 = m[n*((2*k)^1) + i1];
		    double t2 = m[ind1];
		    //int j2 = (j/2)*2;
		    for(int j = 0; j < ind1_k; j++){
			    int j1 = index1[j + 1];
			    if(j1 > br){
				break;
				//continue;
			    }
			    int ind2 = j1 + ((((2*k1) + 1)*((2*k1) + 1))/2);
		            //double op2 = t2 + m[n*(2*k) + j1];
			    double op2 = t2 + m[ind2];
			    int ind3 = j1 + ((((i1^1) + 1)*((i1^1) + 1))/2);
		            //m[n*(i1^1) + j1] = min(m[n*(i1^1) + j1],op2 );
			    if(m[ind3] !=INFINITY){
			    	m[ind3] = min(m[ind3],op2 );
			    }
			    else{
				m[ind3] = op2;
				count++;
			    }
		        }
		        //for(int j = 0; j < index2[m*((2*k)^1)]; j++){
			  for(int j = 0; j < ind4_k; j++){
			     int j1 = index2[n + j + 2];
			     if(j1>i2){
				break;
				//continue;
			    }
		            double op2 = t2 + temp2[j1];
			    int ind3 = (j1^1) + ((((i1^1) + 1)*((i1^1) + 1))/2);
		            //m[n*(i1^1) + (j1^1)] = min(m[n*(i1^1) + (j1^1)],op2 );
			    if(m[ind3]!=INFINITY){
			    	m[ind3] = min(m[ind3],op2 );
			    }
			    else{
				m[ind3] = op2;
				count++;
			    }
		        }
		    //}
		}
	    
		//Third Part i >= (2*k+1)
		for(int i = 0; i < ind4_k; i++){
		    int i1 = index2[n + i + 2];
		    int i2 = (i1%2==0) ? (i1 + 1): i1;
		    int br = i2 < 2*k1 ? i2 : 2*k1 - 1;
		    int ind1 = ((2*k1)^1) + (((i1 + 1)*(i1 + 1))/2);
		    //double t1 = m[n*i1 + ((2*k)^1)];
		    double t1 = m[ind1];
		   
		    for(int j = 0; j < ind2_k; j++){
		   	    //int j1 = index1[m*((2*k)^1) + j + 1];
			    int j1 = index1[n + j + 2];
			    if(j1 > br){
				break;
				//continue;
			    }
			    int ind2 = j1 + (((((2*k1)^1) + 1)*(((2*k1)^1) + 1))/2);
		            //double op1 = t1 + m[n*((2*k)^1) + j1];
			    double op1 = t1 + m[ind2];
			    int ind3 = j1 + (((i1 + 1)*(i1 + 1))/2);
		            //m[n*i1 + j1] = min(m[n*i1 + j1],op1 );
			    if(m[ind3]!=INFINITY){
			    	m[ind3] = min(m[ind3],op1 );
			    }
			    else{
				m[ind3] = op1;
				count++;
			    }
		        }
		        
		 	for(int j = 0; j < ind3_k ; j++){
		            
			    int j1 = index2[j + 1];
			    if(j1>i2){
				break;
				//continue;
			     }
		            double op1 = t1 + temp1[j1];
		            //double op1 = t1 + m[n*(j1) + 2*k];
		     	    int ind3 = (j1^1) + (((i1 + 1)*(i1 + 1))/2);
			    if(m[ind3]!=INFINITY){
			    	m[ind3] = min(m[ind3],op1 );
			    }
			    else{
				m[ind3] = op1;
				count++;
			    }
		        }
		}
		
		//Fourth Part i >= 2*k
		for(int i = 0; i < ind3_k; i++){
		    //int i1 = index2[m*(2*k) + i + 1];
		    int i1 = index2[i + 1];
		    int i2 = (i1%2==0) ? (i1 + 1): i1;
		    int br = i2 < 2*k1 ? i2 : 2*k1 - 1;
		    //double t2 = m[n*i1 + (2*k)];
		    int ind1 = (2*k1) + (((i1 + 1)*(i1 + 1))/2);
		    double t2 = m[ind1];
		    for(int j = 0; j < ind1_k; j++){
			    int j1 = index1[j + 1];
			    if(j1 > br){
				break;
				//continue;
			    }
			    
		            //double op2 = t2 + m[n*(2*k) + j1];
			    int ind2 = j1 + ((((2*k1) + 1)*((2*k1) + 1))/2);
			    double op2 = t2 + m[ind2];
			    int ind3 = j1 + (((i1 + 1)*(i1 + 1))/2);
		            //m[n*i1 + j1] = min(m[n*i1 + j1],op2 );
			    if(m[ind3]!=INFINITY){
			    	m[ind3] = min(m[ind3],op2 );
			    }
			    else{
				m[ind3] = op2;
				count++;
			    }
		        }
		       
		       // j >= 2*k
		       for(int j = 0; j < ind4_k ; j++){
			    int j1 = index2[n + j + 2];
			    if(j1>i2){
				break;
				//continue;
			    }
			    
		            double op2 = t2 + temp2[j1];
			    int ind3 = (j1^1) + (((i1 + 1)*(i1 + 1))/2);
		            //m[n*i1 + (j1^1)] = min(m[n*i1 + (j1^1)],op2 );
			    if(m[ind3]!=INFINITY){
			    	m[ind3] = min(m[ind3],op2 );
			    }
			    else{
				m[ind3] = op2;
				count++;
			    }
		        }
		}
		
	    }
	free(ca);
	cl = cl->next;
    }
    oo->nni = count;
    
    if(is_int){
		if(strengthning_int_comp_sparse(oo,index1,temp1,n)){
			return 1;
		}
    	}
    	else{
    		if(strengthning_comp_sparse(oo,index1,temp1,n)){
			return 1;
		}
    }
    //return strengthning_dense_scalar(m,temp1,n);
}


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


#include "opt_zones_incr_closure.h"

/******************************
	Vector Incremental Closure
*******************************/

void floyd_warshall_incr_vector(double *m, unsigned short int n, unsigned short int v){
	#if defined(VECTOR)
	unsigned short int i,j,k;
	unsigned short int v1 = v +1 ;
	//load the v-th row
	double *pv = m + n*v1;

	double * tmpv = (double *)calloc(n,sizeof(double));
	// store v-th column in array to improve cache performance and enable vectorization
	for(i=0; i < n; i++){
		tmpv[i] = m[n*i+v1];
	}
	
	double * tmpk = (double *)calloc(n,sizeof(double));

	for(k=0; k < n; k++){
		v_double_type op1 = v_set1_double(pv[k]);
		//load the k-th row
		double * pk = m + n*k;
		//update the v-th row
		for(j=0; j < n/v_length; j++){
			v_double_type op2 = v_load_double(pk+j*v_length);
			v_double_type op3 = v_add_double(op1,op2);
			v_double_type op4 = v_load_double(pv+j*v_length);
			v_double_type res = v_min_double(op3,op4);
			v_store_double(pv+j*v_length,res);
		}
		
		for(j=n/v_length; j < n; j++){
			pv[j] = min(pv[j], pv[k] + pk[j]);
		}
	
		//store the k-th column in array to allow vectorization
		for(i=0; i < n; i++){
			tmpk[i] = m[n*i+k];
		}
		//update the v-th column
		op1 = v_set1_double(pk[v1]);
		for(i=0; i < n/v_length; i++){
			v_double_type op2 = v_load_double(tmpk + i*v_length);
			v_double_type op3 = v_add_double(op1,op2);
			v_double_type op4 = v_load_double(tmpv+i*v_length);
			v_double_type res = v_min_double(op3,op4);
			v_store_double(tmpv+i*v_length,res);
			
		}
		
		for(i=n/v_length; i < n; i++){
			tmpv[i] = min(tmpv[i], tmpk[i] + pk[v1]);
		}		
	}
	
	for(i=0; i < n; i++){
		m[n*i+v1] = tmpv[i];
	}

	//update remaining elements
	for(i=0; i < n; i++){
		double *pi = m + n*i;
		v_double_type op1 = v_set1_double(pi[v1]);
		for(j=0; j < n/v_length; j++){
			v_double_type op2 = v_load_double(pv+j*v_length);
			v_double_type op3 = v_add_double(op1,op2);
			v_double_type op4 = v_load_double(pi + j*v_length);
			v_double_type res = v_min_double(op3,op4);
			v_store_double(pi+j*v_length,res);	
		}
		
		for(j=n/v_length; j < n; j++){
			pi[j] = min(pi[j],pi[v1] + pv[j]);		
		}
	}
	free(tmpv);
	free(tmpk);
	#endif
}

bool incr_closure_dense(opt_zones_mat_t *oz, unsigned short int dim, unsigned short int v){
	double *m = oz->mat;
	unsigned short int n = dim + 1;
	floyd_warshall_incr_vector(m,n,v);
	oz->nni = n*n;
	return check_negative_cycle(oz,dim);
}

/************************
	Scalar Incremental Closure
************************/

void floyd_warshall_incr_scalar(double *m, unsigned short int n, unsigned short int v){
	unsigned short int i,j,k;
	//load the v-th row
	unsigned short int v1 = v + 1;
	double * pv = m + n*v1;
	// store v-th column in array to improve cache performance
	double * tmpv = (double *)calloc(n,sizeof(double));
	for(i=0; i < n; i++){
		tmpv[i] = m[n*i+v1];
	}

	for(k=0; k < n; k++){
		//load the k-th row
		double * pk = m + n*k;
		//update the v-th row
		for(j=0; j < n; j++){
			pv[j] = min(pv[j], pv[k] + pk[j]);
		}
		//update the v-th column
		for(i=0; i < n; i++){
			tmpv[i] = min(tmpv[i], m[n*i+k] + pk[v1]);
		}		
	}
	
	for(i=0; i < n; i++){
		m[n*i+v1] = tmpv[i];
	}

	//update remaining elements
	for(i=0; i < n; i++){
		double *pi = m + n*i;
		for(j=0; j < n; j++){
			pi[j] = min(pi[j],pi[v1] + pv[j]);			
		}
	}
	free(tmpv);
	
}

bool incr_closure_dense_scalar(opt_zones_mat_t *oz, unsigned short int dim, unsigned short int v){
	double *m = oz->mat;
	unsigned short int n = dim + 1;
	floyd_warshall_incr_scalar(m,n,v);
	oz->nni = n*n;
	return check_negative_cycle(oz,dim);
}

/*********************************
	Decomposition based incremental closure
*********************************/

bool incr_closure_comp_sparse(opt_zones_mat_t *oz, unsigned short int dim, unsigned short int v){
	
	double *m = oz->mat;
	array_comp_list_t * acl = oz->acl;
	comp_list_t *cl = find(acl,v);
	unsigned short int v1 = v+1;
	unsigned short int n = dim + 1;
	int count = oz->nni;
	if(cl!=NULL){
		unsigned short int comp_size = cl->size;
		unsigned short int * ca = to_sorted_array(cl,n);
		unsigned short int i,j,k; 
		//load the v-th row
		double * pv = m + n*v1;
		// store v-th column in array to improve cache performance
		double * tmpv = (double *)calloc(comp_size,sizeof(double));
		for(i=0; i < comp_size; i++){
			unsigned short int i1 = ca[i]+1;
			tmpv[i] = m[n*i1 + v1];			
		}
		
		for(k=0; k < comp_size; k++){
			unsigned short int k1 = ca[k]+1;
			//load the k-th row
			double * pk = m + n*k1;

			//update the v-th row
			for(j=0; j < comp_size; j++){
				unsigned short int j1 = ca[j]+1;
				//if(pv[j1]==INFINITY){
				//	count++;
				//	pv[j1] = pv[k1] + pk[j1];
				//}
				//else{
					
					pv[j1] = min(pv[j1], pv[k1] + pk[j1]);
				//}
			}

			//update the v-th column
			for(i=0; i < comp_size; i++){
				unsigned short int i1 = ca[i]+1;
				//if(tmpv[i]==INFINITY){
				//	count++;
				//	tmpv[i] = m[n*i1+k1] + pk[v1];
				//}
				//else{
					tmpv[i] = min(tmpv[i], m[n*i1+k1] + pk[v1]);
				//}
			}
			
		}
		count = count + 2*n;
		for(i=0; i < comp_size; i++){
			unsigned short int i1 = ca[i]+1;
			m[n*i1+v1] = tmpv[i];
		}
		
		unsigned short int * index1 = (unsigned short int *)calloc(comp_size+1,sizeof(unsigned short int));
		unsigned short int * index2 = (unsigned short int *)calloc(comp_size+1,sizeof(unsigned short int));
		compute_sparse_index(m,ca,comp_size,index1,index2,v,dim);
		unsigned short int n_ik = index1[0];
		unsigned short int n_kj = index2[0];
	
		//update remaining elements
		for(i=0; i < n_ik; i++){
			unsigned short int i1 = index1[i+1];
			double *pi = m + n*i1;
			for(j=0; j < n_kj; j++){
				unsigned short int j1 = index2[j+1];
				if(pi[j1]==INFINITY){
					pi[j1] = pi[v1] + pv[j1];
					count++;
				}
				else{
					pi[j1] = min(pi[j1],pi[v1] + pv[j1]);
				}			
			}
		}
		
		
		
		// update the 0-th column
		for(i=0; i < comp_size; i++){
			unsigned short int i1 = ca[i]+1;
			if(m[n*i1]==INFINITY){
				m[n*i1] = m[n*i1+v1] + pv[0];
				count++;
			}
			else{
				m[n*i1] = min(m[n*i1], m[n*i1 + v1] + pv[0]);
			}
		}
		
		
		// update the 0-th row
		for(j=0; j < comp_size; j++){
			unsigned short int j1 = ca[j]+1; 
			
			if(m[j1]==INFINITY){
				m[j1] = m[v1] + pv[j1];
				count++;
			}
			else{
				m[j1] = min(m[j1],m[v1] + pv[j1]);
			}
			
		}
		free(ca);
		free(tmpv);
		free(index1);
		free(index2);
		oz->nni = count;
		
		return strengthening_single_comp_zones(oz,cl,dim);
	}
	else{
		return check_negative_cycle(oz,dim);
	}
}


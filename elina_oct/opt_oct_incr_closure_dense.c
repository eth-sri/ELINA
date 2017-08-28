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

#include <stdio.h>
#include <immintrin.h>
#include <math.h>
#include "opt_oct_incr_closure_dense.h"
#include "opt_oct_closure_dense.h"


double incremental_closure_calc_perf_dense(double cycles, int dim){
  int n = 2*dim;
  return  (7*n*n)/cycles;
}


bool incremental_closure_opt_dense(opt_oct_mat_t *oo, int dim, int v, bool is_int){
	double *m = oo->mat;
	int n = 2*dim;
	int ii = 2*v + 1;
	int j;
	double *temp1, *temp2;
	temp1 = (double *)calloc(2*dim, sizeof(double));
  	temp2 = (double *)calloc(2*dim, sizeof(double));
	/******
		incremental Floyd-Warshall : v in end-point position 
	******/
	for(unsigned k = 0; k < n; k=k + 2){
		
		int v1 = 2*v;
		int v2 = 2*v + 1;
		int v1v2 = v2 + (((v1 + 1)*(v1 + 1))/2);
		int v2v1 = v1 + (((v2 + 1)*(v2 + 1))/2);
		int kk = (k^1);
		int br1 = k < v1 ? k : v1;
		int br2 = kk < v1 ? kk : v1;
		/*int v1v2 = v2 + (((v1 + 1)*(v1 + 1))/2);
		int v2v1 = v1 + (((v2 + 1)*(v2 + 1))/2);
		int v1k = k + (((v1 + 1)*(v1 + 1))/2);
		int v2k = k +  (((v2 + 1)*(v2 + 1))/2);
		int kv1 = v1 + (((k + 1)*(k + 1))/2);
		int kv2 = v2 +  (((k + 1)*(k + 1))/2);*/
		//m[n*2*v + (2*v + 1)] = min(m[n*2*v + (2*v + 1)], m[n*2*v + k] + m[n*k + 2*v + 1]);
		//m[v1v2] = min(m[v1v2], m[v1k] + m[kv2]);
		//m[n*(2*v + 1) + 2*v] = min(m[n*(2*v + 1) + 2*v], m[n*(2*v + 1) + k] + m[n*k + 2*v]);
		//m[v2v1] = min(m[v2v1], m[v2k] + m[kv1]);
		for(unsigned i = 2*v; i < 2*v + 2; i++){
			//double ik = m[n*i + k];
			int ind_ik, ind_ikk;
			if(k <=i){
				ind_ik = k + (((i + 1)*(i + 1))/2);
			}
			else{
				ind_ik = (i^1) + ((((k^1) + 1)*((k^1) + 1))/2);
			}
			
			if(kk <=i){
				ind_ikk = kk + (((i + 1)*(i + 1))/2);
			}
			else{
				ind_ikk = (i^1) + ((((kk^1) + 1)*((kk^1) + 1))/2);
			}
			double ik = m[ind_ik];
			double ikk = m[ind_ikk];
			v_double_type vik = v_set1_double(m[ind_ik]);
			v_double_type vikk = v_set1_double(m[ind_ikk]);
			//double ki = m[n*k + i];
			int ind_ki, ind_kki;
			if ( k <= i){
				ind_ki = (k^1) + ((((i^1) + 1)*((i^1) + 1))/2);
			}
			else{
				ind_ki = i + (((k + 1)*(k + 1))/2);
			}

			if ( kk <= i){
				ind_kki = (kk^1) + ((((i^1) + 1)*((i^1) + 1))/2);
			}
			else{
				ind_kki = i + (((kk + 1)*(kk + 1))/2);
			}
			
			//int ind_ki = i + (((k + 1)*(k + 1))/2);
			double ki = m[ind_ki];	
			double kki = m[ind_kki];
			/******* 
				v in first end-point position.
				This part is vectorized. 
			********/		
			for(j = 0; j <br1/v_length; j++){
				//double kj = m[n*k + j];
				int ind_kj = j*v_length + (((k + 1)*(k + 1))/2);
				v_double_type vkj = v_load_double(m+ind_kj);
				int ind_kkj = j*v_length + (((kk + 1)*(kk + 1))/2);
				v_double_type vkkj = v_load_double(m+ind_kkj);
				//double jk = m[n*j + k];
				//int ind_jk = k + (((j + 1)*(j + 1))/2);
				//double jk = m[ind_jk];
				//m[n*i + j] = min(m[n*i + j], ik + kj);
				int ind_ij = j*v_length + (((i + 1)*(i + 1))/2);
				v_double_type vij = v_load_double(m+ind_ij);
				v_double_type op1 = v_add_double(vik,vkj);
				v_double_type op2 = v_add_double(vikk,vkkj);
				v_double_type op3 = v_min_double(op1, op2);
				v_double_type res = v_min_double(vij,op3);
				v_store_double(m+ind_ij,res);
				//m[n*j + i] = min(m[n*j + i], jk + ki);
			}
			for(j = (br1/v_length)*v_length; j < br1; j++){
				int ind_kj = j + (((k + 1)*(k + 1))/2);
				double kj = m[ind_kj];
				int ind_kkj = j + (((kk + 1)*(kk + 1))/2);
				double kkj = m[ind_kkj];
				int ind_ij = j + (((i + 1)*(i + 1))/2);
				m[ind_ij] = min(m[ind_ij], ik + kj);
				m[ind_ij] = min(m[ind_ij], ikk + kkj);
			}
			for(; j < v1; j++){
				//double kj = m[n*k + j];
				int ind_kj = (k^1) + ((((j^1) + 1)*((j^1) + 1))/2);
				double kj = m[ind_kj];
				int ind_kkj = (kk^1) + ((((j^1) + 1)*((j^1) + 1))/2);
				double kkj = m[ind_kkj];
				//double jk = m[n*j + k];
				int ind_ij = j + (((i + 1)*(i + 1))/2);
				//m[n*i + j] = min(m[n*i + j], ik + kj);
				m[ind_ij] = min(m[ind_ij], ik + kj);
				m[ind_ij] = min(m[ind_ij], ikk + kkj);
				//m[n*j + i] = min(m[n*j + i], jk + ki);
			}
			/****** 
				v in second end-point position.
				This part is scalar. 
			*****/
			for(j= v1+2; j < k; j++ ){
				int ind_jk = (j^1) + ((((k^1) + 1)*((k^1) + 1))/2);
				double jk = m[ind_jk];
				int ind_jkk = (j^1) + ((((kk^1) + 1)*((kk^1) + 1))/2);
				double jkk = m[ind_jkk];
				int ind_ji = i + (((j + 1)*(j + 1))/2);
				m[ind_ji] = min(m[ind_ji], jk + ki);
				m[ind_ji] = min(m[ind_ji], jkk + kki);
			}
			
			for(; j < 2*dim; j++){
				int ind_jk = k + (((j + 1)*(j + 1))/2);
				double jk = m[ind_jk];
				int ind_jkk = kk + (((j + 1)*(j + 1))/2);
				double jkk = m[ind_jkk];
				int ind_ji = i + (((j + 1)*(j + 1))/2);
				m[ind_ji] = min(m[ind_ji], jk + ki);
				m[ind_ji] = min(m[ind_ji], jkk + kki);
			}
			
		}
		m[v1v2] = min(m[v1v2],m[opt_matpos2(v1,k)] + m[opt_matpos2(k,v2)]);
		m[v1v2] = min(m[v1v2],m[opt_matpos2(v1,kk)] + m[opt_matpos2(kk,v2)]);
		m[v2v1] = min(m[v2v1],m[opt_matpos2(v2,k)] + m[opt_matpos2(k,v1)]);
		m[v2v1] = min(m[v2v1],m[opt_matpos2(v2,kk)] + m[opt_matpos2(kk,v1)]);
		
	}
	
	int v1 = (2*v);
	int v2 = (2*v)^1;
	int vi = (((v1 + 1)*(v1 + 1))/2);
	int vvi = (((v2 + 1)*(v2 + 1))/2);
	int pos1 = v2 + vi;
	//int pos2 = opt_matpos2((2*k)^1, 2*k);
	int pos2 = v1 + vvi;
	int l = (v1 + 2);
	int mod = l%v_length;
	if(mod){
		l = l + (v_length - mod);
	}
	//variable v in pivot position
	for(int i = v1 + 2; i < n;i++){
		int ind1 = v2 + (((i+1)*(i+1))/2);
		int ind2 = v1 + (((i+1)*(i+1))/2);
		m[ind1] = min(m[ind1], m[ind2] + m[pos1] );
		temp2[i^1] = m[ind1];
	}
	for(int i = v1 + 2; i < n; i++){
		int ind1 = v2 + (((i+1)*(i+1))/2);
		int ind2 = v1 + (((i+1)*(i+1))/2);
		m[ind2] = min(m[ind2], m[ind1] + m[pos2] );
		temp1[i^1] = m[ind2];
	}

	for(int j = 0; j < v1; j++){
		//int ind3 = opt_matpos2((2*k)^1,j);
		int ind3 = j + vvi;
		//int ind4 = opt_matpos2( 2*k,j);
		int ind4 = j + vi;
		//result[n*((2*k)^1) + j] = min(result[n*((2*k)^1) + j], result[n*((2*k)^1) + 2*k] + result[n*(2*k) + j]);
		m[ind3] = min(m[ind3], m[pos2] + m[ind4]);
	}
	for(int j = 0; j < v1; j++){
		//int ind3 = opt_matpos2((2*k)^1,j);
		int ind3 = j + vvi;
		//int ind4 = opt_matpos2(2*k,j);
		int ind4 = j + vi;
		//result[n*2*k + j] = min(result[n*2*k + j], result[n*2*k + ((2*k)^1)] + result[n*((2*k)^1) + j]);
		m[ind4] = min(m[ind4], m[pos1] + m[ind3]);
	}


	//for(unsigned k = 2*v; k < 2*v + 2; k++){
		/******
			 incremental Floyd-Warshall : v in pivot position.
			 This is vectorized. 
		******/
		for(unsigned i = 0; i < v1; i++){
			int i2 = (i|1);
			int br = i2;
			int j;
			int ind_ik, ind_ikk;
			//ind_ik = (i^1) + ((((k^1) + 1)*((k^1) + 1))/2);
			ind_ik = (i^1) + vvi;
			//ind_ikk = (i^1) + ((((kk^1) + 1)*((kk^1) + 1))/2);
			ind_ikk = (i^1) + vi;
			//double ik = m[n*i + k];
			v_double_type ik = v_set1_double(m[ind_ik]);
			v_double_type ikk = v_set1_double(m[ind_ikk]);
			double ikd = m[ind_ik];
			double ikkd = m[ind_ikk];
			for(j = 0; j < br/v_length; j++){
				//int ind_kj = j + (((k + 1)*(k + 1))/2);
				int ind_kj = j*v_length + vi;
				//double kj = m[n*k + j];
				//double kj = m[ind_kj];
				v_double_type kj = v_load_double(m + ind_kj);
				//int ind_kkj = j + (((kk + 1)*(kk + 1))/2);
				int ind_kkj = j*v_length + vvi;
				//double kkj = m[ind_kkj];
				v_double_type kkj = v_load_double(m + ind_kkj);
				int ind_ij = j*v_length + (((i + 1)*(i + 1))/2);
				//m[n*i + j] = min(m[n*i + j], ik + kj);
				v_double_type ij = v_load_double(m + ind_ij);
				v_double_type op1 = v_add_double(ik,kj);
				v_double_type op2 = v_add_double(ikk,kkj);
				v_double_type op3 = v_min_double(op1,op2);
				v_double_type res = v_min_double(ij,op3);
				//m[ind_ij] = min(m[ind_ij], ik + kj);
				//m[ind_ij] = min(m[ind_ij], ikk + kkj);
				v_store_double(m + ind_ij,res);
			}
			for(j = (br/v_length)*v_length; j <= br; j++){
				//int ind_kj = j + (((k + 1)*(k + 1))/2);
				int ind_kj = j + vi;
				//double kj = m[n*k + j];
				double kj = m[ind_kj];
				//int ind_kkj = j + (((kk + 1)*(kk + 1))/2);
				int ind_kkj = j + vvi;
				double kkj = m[ind_kkj];
				int ind_ij = j + (((i + 1)*(i + 1))/2);
				//m[n*i + j] = min(m[n*i + j], ik + kj);
				m[ind_ij] = min(m[ind_ij], ikd + kj);
				m[ind_ij] = min(m[ind_ij], ikkd + kkj);
			}
		}


		for(unsigned i = 2*v + 2; i < n; i++){
			int i2 = (i|1);
			int br = v1;
			int j;
			int ind_ik, ind_ikk;
			ind_ik = v1 + (((i + 1)*(i + 1))/2);				
			ind_ikk = v2 + (((i + 1)*(i + 1))/2);				
			//double ik = m[n*i + k];
			v_double_type ik = v_set1_double(m[ind_ik]);
			v_double_type ikk = v_set1_double(m[ind_ikk]);
			double ikd = m[ind_ik];
			double ikkd = m[ind_ikk];
			int b = min(l,i2);
			for(j = 0; j < br/v_length; j++){
				//int ind_kj = j + (((k + 1)*(k + 1))/2);
				int ind_kj = j*v_length + vi;
				//double kj = m[n*k + j];
				//double kj = m[ind_kj];
				v_double_type kj = v_load_double(m + ind_kj);
				//int ind_kkj = j + (((kk + 1)*(kk + 1))/2);
				int ind_kkj = j*v_length + vvi;
				//double kkj = m[ind_kkj];
				v_double_type kkj = v_load_double(m + ind_kkj);
				int ind_ij = j*v_length + (((i + 1)*(i + 1))/2);
				//m[n*i + j] = min(m[n*i + j], ik + kj);
				v_double_type ij = v_load_double(m + ind_ij);
				v_double_type op1 = v_add_double(ik,kj);
				v_double_type op2 = v_add_double(ikk,kkj);
				v_double_type op3 = v_min_double(op1,op2);
				v_double_type res = v_min_double(ij,op3);
				//m[ind_ij] = min(m[ind_ij], ik + kj);
				//m[ind_ij] = min(m[ind_ij], ikk + kkj);
				v_store_double(m + ind_ij,res);
			}
			for(j = (br/v_length)*v_length; j <= br; j++){
				//int ind_kj = j + (((k + 1)*(k + 1))/2);
				int ind_kj = j + vi;
				//double kj = m[n*k + j];
				double kj = m[ind_kj];
				//int ind_kkj = j + (((kk + 1)*(kk + 1))/2);
				int ind_kkj = j + vvi;
				double kkj = m[ind_kkj];
				int ind_ij = j + (((i + 1)*(i + 1))/2);
				//m[n*i + j] = min(m[n*i + j], ik + kj);
				m[ind_ij] = min(m[ind_ij], ikd + kj);
				m[ind_ij] = min(m[ind_ij], ikkd + kkj);
			}

			for(j = v1 + 2; j <= b; j++){
				//double kj = m[n*k + j];
				double kj = temp2[j];
				//int ind_kkj = (kk^1) + ((((j^1) + 1)*((j^1) + 1))/2);
				//int ind_kkj = temp1[j];
				double kkj = temp1[j];
				//m[n*i + j] = min(m[n*i + j], ik + kj);
				int ind_ij = j + (((i + 1)*(i + 1))/2);
				m[ind_ij] = min(m[ind_ij], ikd + kj);
				m[ind_ij] = min(m[ind_ij], ikkd + kkj);
			}
			if(b < i2){
				for(j = b/v_length; j < i2/v_length; j++){
					v_double_type kj = v_load_double(temp2 + j*v_length);
					v_double_type kkj = v_load_double(temp1 + j*v_length);
					int ind_ij = j*v_length + (((i + 1)*(i + 1))/2);
					v_double_type ij = v_load_double(m + ind_ij);
					v_double_type op1 = v_add_double(ik,kj);
					v_double_type op2 = v_add_double(ikk,kkj);
					v_double_type op3 = v_min_double(op1,op2);
					v_double_type res = v_min_double(ij,op3);
					v_store_double(m + ind_ij,res);
				}
				for(j = (i2/v_length)*v_length; j <=i2; j++){
					//int ind_kj = (k^1) + ((((j^1) + 1)*((j^1) + 1))/2);
					//int ind_kj = v2 + ((((j^1) + 1)*((j^1) + 1))/2);
					//double kj = m[n*k + j];
					double kj = temp2[j];
					//int ind_kkj = (kk^1) + ((((j^1) + 1)*((j^1) + 1))/2);
					//int ind_kkj = v1 + ((((j^1) + 1)*((j^1) + 1))/2);
					double kkj = temp1[j];
					//m[n*i + j] = min(m[n*i + j], ik + kj);
					int ind_ij = j + (((i + 1)*(i + 1))/2);
					m[ind_ij] = min(m[ind_ij], ikd + kj);
					m[ind_ij] = min(m[ind_ij], ikkd + kkj);
				}
			}
		}
	
	oo->nni = 2*dim*(dim+1);
	bool res;	
	if(is_int){
		res = strengthning_int_dense(oo, temp1, n);
	}
	else{
        	res = strengthning_dense(oo, temp1, n);
	}
    free(temp1);
    free(temp2);
	return res;
}




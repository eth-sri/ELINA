/*
	Copyright 2015 Software Reliability Lab, ETH Zurich

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

		http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
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
			__m256d vik = _mm256_set1_pd(m[ind_ik]);
			__m256d vikk = _mm256_set1_pd(m[ind_ikk]);
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
			for(j = 0; j <br1/4; j++){
				//double kj = m[n*k + j];
				int ind_kj = j*4 + (((k + 1)*(k + 1))/2);
				__m256d vkj = _mm256_loadu_pd(m+ind_kj);
				int ind_kkj = j*4 + (((kk + 1)*(kk + 1))/2);
				__m256d vkkj = _mm256_loadu_pd(m+ind_kkj);
				//double jk = m[n*j + k];
				//int ind_jk = k + (((j + 1)*(j + 1))/2);
				//double jk = m[ind_jk];
				//m[n*i + j] = min(m[n*i + j], ik + kj);
				int ind_ij = j*4 + (((i + 1)*(i + 1))/2);
				__m256d vij = _mm256_loadu_pd(m+ind_ij);
				__m256d op1 = _mm256_add_pd(vik,vkj);
				__m256d op2 = _mm256_add_pd(vikk,vkkj);
				__m256d op3 = _mm256_min_pd(op1, op2);
				__m256d res = _mm256_min_pd(vij,op3);
				_mm256_storeu_pd(m+ind_ij,res);
				//m[n*j + i] = min(m[n*j + i], jk + ki);
			}
			for(j = (br1/4)*4; j < br1; j++){
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
	int mod = l%4;
	if(mod){
		l = l + (4 - mod);
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
			__m256d ik = _mm256_set1_pd(m[ind_ik]);
			__m256d ikk = _mm256_set1_pd(m[ind_ikk]);
			double ikd = m[ind_ik];
			double ikkd = m[ind_ikk];
			for(j = 0; j < br/4; j++){
				//int ind_kj = j + (((k + 1)*(k + 1))/2);
				int ind_kj = j*4 + vi;
				//double kj = m[n*k + j];
				//double kj = m[ind_kj];
				__m256d kj = _mm256_loadu_pd(m + ind_kj);
				//int ind_kkj = j + (((kk + 1)*(kk + 1))/2);
				int ind_kkj = j*4 + vvi;
				//double kkj = m[ind_kkj];
				__m256d kkj = _mm256_loadu_pd(m + ind_kkj);
				int ind_ij = j*4 + (((i + 1)*(i + 1))/2);
				//m[n*i + j] = min(m[n*i + j], ik + kj);
				__m256d ij = _mm256_loadu_pd(m + ind_ij);
				__m256d op1 = _mm256_add_pd(ik,kj);
				__m256d op2 = _mm256_add_pd(ikk,kkj);
				__m256d op3 = _mm256_min_pd(op1,op2);
				__m256d res = _mm256_min_pd(ij,op3);
				//m[ind_ij] = min(m[ind_ij], ik + kj);
				//m[ind_ij] = min(m[ind_ij], ikk + kkj);
				_mm256_storeu_pd(m + ind_ij,res);
			}
			for(j = (br/4)*4; j <= br; j++){
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
			__m256d ik = _mm256_set1_pd(m[ind_ik]);
			__m256d ikk = _mm256_set1_pd(m[ind_ikk]);
			double ikd = m[ind_ik];
			double ikkd = m[ind_ikk];
			int b = min(l,i2);
			for(j = 0; j < br/4; j++){
				//int ind_kj = j + (((k + 1)*(k + 1))/2);
				int ind_kj = j*4 + vi;
				//double kj = m[n*k + j];
				//double kj = m[ind_kj];
				__m256d kj = _mm256_loadu_pd(m + ind_kj);
				//int ind_kkj = j + (((kk + 1)*(kk + 1))/2);
				int ind_kkj = j*4 + vvi;
				//double kkj = m[ind_kkj];
				__m256d kkj = _mm256_loadu_pd(m + ind_kkj);
				int ind_ij = j*4 + (((i + 1)*(i + 1))/2);
				//m[n*i + j] = min(m[n*i + j], ik + kj);
				__m256d ij = _mm256_loadu_pd(m + ind_ij);
				__m256d op1 = _mm256_add_pd(ik,kj);
				__m256d op2 = _mm256_add_pd(ikk,kkj);
				__m256d op3 = _mm256_min_pd(op1,op2);
				__m256d res = _mm256_min_pd(ij,op3);
				//m[ind_ij] = min(m[ind_ij], ik + kj);
				//m[ind_ij] = min(m[ind_ij], ikk + kkj);
				_mm256_storeu_pd(m + ind_ij,res);
			}
			for(j = (br/4)*4; j <= br; j++){
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
				for(j = b/4; j < i2/4; j++){
					__m256d kj = _mm256_loadu_pd(temp2 + j*4);
					__m256d kkj = _mm256_loadu_pd(temp1 + j*4);
					int ind_ij = j*4 + (((i + 1)*(i + 1))/2);
					__m256d ij = _mm256_loadu_pd(m + ind_ij);
					__m256d op1 = _mm256_add_pd(ik,kj);
					__m256d op2 = _mm256_add_pd(ikk,kkj);
					__m256d op3 = _mm256_min_pd(op1,op2);
					__m256d res = _mm256_min_pd(ij,op3);
					_mm256_storeu_pd(m + ind_ij,res);
				}
				for(j = (i2/4)*4; j <=i2; j++){
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
	
	return res;
}




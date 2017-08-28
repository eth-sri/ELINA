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

#include "opt_oct_closure_dense_scalar.h"


bool strengthning_int_dense_scalar(opt_oct_mat_t * oo, double *temp, int n){
	
	double *m = oo->mat;
	/*****
		Store strengthening operands in an array
	******/
	for(int i = 0; i < n; i++){
		int ind1 = opt_matpos2(i^1, i);
		temp[i] = ceil(m[ind1]/2);
	}
	
	/******
		Apply scalar strengthening
	*******/
	for(int i = 0; i < n; i++){
		for(int j = 0; j <= (i|1); j++){
			int ind = j + (((i+1)*(i+1))/2);
			m[ind] = min(m[ind], temp[i^1] + temp[j]);
		}
			
	}

	/******
		Check for negative cycle
	*******/
	for(int i = 0; i < n; i++){
		int ind = i + (((i+1)*(i+1))/2);
		if(m[ind] < 0){
			return true;
		}
		else{
			m[ind] = 0;
		}
	}
	return false;
}

bool strengthning_dense_scalar(opt_oct_mat_t * oo, double *temp, int n){
	double *m = oo->mat;
	/*****
		Store strengthening operands in an array
	******/
	for(int i = 0; i < n; i++){
		int ind1 = opt_matpos2(i^1,i);
		temp[i] = m[ind1];
	}
	
	/******
		Apply scalar strengthening
	*******/
	for(int i = 0; i < n; i++){
		for(int j = 0; j <= (i|1); j++){
			int ind = j + (((i+1)*(i+1))/2);
			m[ind] = min(m[ind], (temp[i^1] + temp[j])/2);
		}
			
	}
	
	/******
		Check for negative cycle
	*******/
	for(int i = 0; i < n; i++){
		int ind = i + (((i+1)*(i+1))/2);
		if(m[ind] < 0){
			return true;
		}
		else{
			m[ind] = 0;
		}
	}
	return false;
}


bool floyd_warshall_dense_scalar(opt_oct_mat_t *oo, double *temp1, double *temp2, int dim, bool is_int){
    double *m = oo->mat;
    int size = 4 * dim * dim;
    int n = 2*dim; 
    double count = 0;
    /******
		Floyd Warshall step
    *******/
    for(int k = 0; k < dim; k++){
	//int k2 = k==0 ? k + 1 : k;
	//int pos1 = matpos2(2*k, (2*k)^1);
	int pos1 = ((2*k)^1) + ((((2*k) + 1)*((2*k) + 1))/2);
	//int pos2 = matpos2((2*k)^1, 2*k);
	int pos2 = (2*k) + (((((2*k)^1) + 1)*(((2*k)^1) + 1))/2);
	
	/******
		For k-th iteration update 2k and (2k+1)-th
		row and column first. 
	*******/
	for(int i = 2*k + 2; i < n;i++){
		//int ind1 = matpos2(i,((2*k)^1));
		int ind1 = ((2*k)^1) + (((i+1)*(i+1))/2);
		//int ind2 = matpos2(i,2*k);
		int ind2 = (2*k) + (((i+1)*(i+1))/2);
		//int ind2 = n*i + ((2*k)^1);
		//int ind1 = n*i + (2*k);
		//m[ind2] = std::min(m[ind2], m[ind1] + m[n*(2*k) + ((2*k)^1)] );
		m[ind1] = min(m[ind1], m[ind2] + m[pos1] );
		temp2[i^1] = m[ind1];
	}


	
	for(int i = 2*k + 2; i < n; i++){
		//int ind1 = matpos2(i,((2*k)^1));
		int ind1 = ((2*k)^1) + (((i+1)*(i+1))/2);
		//int ind2 = matpos2(i,2*k);
		int ind2 = (2*k) + (((i+1)*(i+1))/2);
		//int ind2 = n*i + ((2*k)^1);
		//int ind1 = n*i + (2*k);
		//m[ind1] = std::min(m[ind1], m[ind2] + m[n*((2*k)^1) + (2*k)] );
		m[ind2] = min(m[ind2], m[ind1] + m[pos2] );
		temp1[i^1] = m[ind2];
	}

	for(int j = 0; j < (2*k); j++){
		//int ind3 = matpos2((2*k)^1,j);
		int ind3 = j + (((((2*k)^1)+1)*(((2*k)^1)+1))/2);
		//int ind4 = matpos2( 2*k,j);
		int ind4 = j + ((((2*k)+1)*((2*k)+1))/2);
		//m[n*((2*k)^1) + j] = std::min(m[n*((2*k)^1) + j], m[n*((2*k)^1) + 2*k] + m[n*(2*k) + j]);
		m[ind3] = min(m[ind3], m[pos2] + m[ind4]);
	}
	for(int j = 0; j < (2*k); j++){
		//int ind3 = matpos2((2*k)^1,j);
		int ind3 = j + (((((2*k)^1)+1)*(((2*k)^1)+1))/2);
		//int ind4 = matpos2(2*k,j);
		int ind4 = j + ((((2*k)+1)*((2*k)+1))/2);
		//m[n*2*k + j] = std::min(m[n*2*k + j], m[n*2*k + ((2*k)^1)] + m[n*((2*k)^1) + j]);
		m[ind4] = min(m[ind4], m[pos1] + m[ind3]);
	}

	/*******
		This is the main loop. Apply k-th iteration on
		rest of elements, which is equivalent to 2k and (2k+1)-th
		iteration in APRON strong closure algorithm.
	********/
	for(int i = 0; i < 2*k; i++){
		int i2 = (i%2==0) ? (i + 1): i;
		int br = i2 < 2*k ? i2 : 2*k - 1;
		//int ind1 = matpos2(i,2*k);
		int ind1 = (i^1) + (((((2*k)^1) + 1)*(((2*k)^1) + 1))/2);
		//int ind2 = matpos2(i, ((2*k)^1));
		int ind2 = (i^1) + ((((2*k) + 1)*((2*k) + 1))/2);
		//double t1 = m[n*(2*k) + (i^1)];
		//double t2 = m[n*((2*k)^1) + (i^1)];
		double t1 = m[ind2];
		double t2 = m[ind1];
		//int j2 = (j/2)*2;
			for(int j = 0; j <=br; j++){
				//int ind3 = matpos2((2*k)^1,j);
				int ind3 = j + (((((2*k)^1)+1)*(((2*k)^1)+1))/2);
				//int ind4 = matpos2(2*k,j);
				int ind4 = j + ((((2*k)+1)*((2*k)+1))/2);
				//int ind5 = matpos2(i,j);
				int ind5 = j + (((i+1)*(i+1))/2);
				//int ind2 = matpos2(k,j);
				//double op1 = t1 + m[n*((2*k)^1) + j];
				//double op2 = t2 + m[n*(2*k) + j];
				double op1 = t1 + m[ind3];
				double op2 = t2 + m[ind4];
				double op3 = min(op1, op2);
				m[ind5] = min(m[ind5],op3);
				count = count + 4;
			}
			for(int j = (2*k) + 2; j <=i2; j++){
				//int ind3 = matpos2(n,(2*k)^1,j);
				//int ind4 = matpos2(n, 2*k,j);
				//int ind2 = matpos2(k,j);
				//double op1 = t1 + m[n*(j^1) + 2*k];
				//double op2 = t2 + m[n*(j^1) + ((2*k)^1)];
				//int ind5 = matpos2(i,j);
				int ind5 = j + (((i+1)*(i+1))/2);
				double op1 = t1 + temp1[j];
				double op2 = t2 + temp2[j];
				double op3 = min(op1, op2);
				m[ind5] = min(m[ind5],op3 );
				count = count + 4;
			}
		//}
	}

	for(int i = 2*k + 2; i < n; i++){
		int i2 = (i%2==0) ? (i + 1): i;
		int br = i2 < 2*k ? i2 : 2*k - 1;
		//int ind1 = matpos2(i,(2*k)^1);
		int ind1 = ((2*k)^1) + (((i+1)*(i+1))/2);
		//int ind2 = matpos2(i,2*k);
		int ind2 = (2*k) + (((i+1)*(i+1))/2);
		//double t1 = m[n*i + ((2*k)^1)];
		//double t2 = m[n*i + 2*k];
		double t1 = m[ind1];
		double t2 = m[ind2];
		//int j2 = (j/2)*2;
			for(int j = 0; j <= br; j++){
				//int ind3 = matpos2((2*k)^1,j);
				int ind3 = j + (((((2*k)^1)+1)*(((2*k)^1)+1))/2);
				//int ind4 = matpos2( 2*k,j);
				int ind4 = j + ((((2*k)+1)*((2*k)+1))/2);
				//int ind5 = matpos2(i,j);
				int ind5 = j + (((i+1)*(i+1))/2);
				//int ind2 = matpos2(k,j);
				//double op1 = t1 + m[n*((2*k)^1) + j];
				//double op2 = t2 + m[n*(2*k) + j];
				double op1 = t1 + m[ind3];
				double op2 = t2 + m[ind4];
				double op3 = min(op1, op2);
				m[ind5] = min(m[ind5],op3 );
				count = count + 4;
			}
			for(int j = (2*k) + 2; j <= i2; j++){
				//int ind3 = matpos2(n,(2*k)^1,j);
				//int ind4 = matpos2(n, 2*k,j);
				//int ind2 = matpos2(k,j);
				//double op1 = t1 + m[n*(j^1) + 2*k];
				//double op2 = t2 + m[n*(j^1) + ((2*k)^1)];
				//int ind5 = matpos2(i,j);
				int ind5 = j + (((i+1)*(i+1))/2);
				double op1 = t1 + temp1[j];
				double op2 = t2 + temp2[j];
				double op3 = min(op1, op2);
				m[ind5] = min(m[ind5],op3 );
				count = count + 4;
			}
		//}
	}
	
    }
}

bool strong_closure_dense_scalar(opt_oct_mat_t *oo, double *temp1, double *temp2, int dim, bool is_int){
    floyd_warshall_dense_scalar(oo,temp1,temp2,dim,is_int);
    int n = 2*dim;
    oo->nni = 2*dim*(dim+1);
    if(is_int){
	return strengthning_int_dense_scalar(oo,temp1,n);
    }
    else{
    	return strengthning_dense_scalar(oo,temp1,n);
    }
}


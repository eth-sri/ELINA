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


#include "opt_oct_closure_dense.h"
#include <assert.h>

#define U_i 32
#define U_j 32



double strong_closure_calc_perf_dense(double cycles, int dim){
  double n = 2*dim;
  return  (n*n*n)/cycles;
}

bool strengthning_int_dense(opt_oct_mat_t * oo, double *temp, int n){
        double *m = oo->mat;
	/*****
		Store strengthening operands in an array
	******/
	for(int i = 0; i < n; i++){
		//int ind1 = matpos2(i^1, i);
		int ind1 = i + ((((i^1) + 1)*((i^1) + 1))/2);
		temp[i] = ceil(m[ind1]/2);
	}
	
	/******
		Apply vectorized strengthening
	*******/
	for(int i = 0; i < n; i++){
		int i2 = (i|1);
		v_double_type t1 = v_set1_double(temp[i^1]);
		double *p = m + (((i+1)*(i+1))/2);
		for(int j = 0; j < (i2/v_length); j++){
			//int ind = j*8 + (((i+1)*(i+1))/2);
			v_double_type t2 = v_load_double(temp + j*v_length);
			v_double_type op1 = v_add_double(t1,t2);
			v_double_type op2 = v_load_double(p + j*v_length);
			v_double_type res = v_min_double(op1, op2);
			//m[ind] = min(m[ind], temp[i^1] + temp[j]);
			v_store_double(m + (((i+1)*(i+1))/2) + j*v_length, res);
		}
		for(int j = (i2/v_length)*v_length; j<= i2; j++){
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

bool strengthning_dense(opt_oct_mat_t * oo, double *temp, int n){
	/*****
		Store strengthening operands in an array
	******/
	double *m = oo->mat;
	for(int i = 0; i < n; i++){
		//int ind1 = matpos2(i^1,i);
		int ind1 = i + ((((i^1) + 1)*((i^1) + 1))/2);
		temp[i] = m[ind1]/2;
	}
	
	/******
		Apply vectorized strengthening
	*******/
	for(int i = 0; i < n; i++){
		int i2 = (i|1);
		v_double_type t1 = v_set1_double(temp[i^1]);
		double *p = m + (((i+1)*(i+1))/2);
		for(int j = 0; j < (i2/v_length); j++){
			//int ind = j*8 + (((i+1)*(i+1))/2);
			v_double_type t2 = v_load_double(temp + j*v_length);
			v_double_type op1 = v_add_double(t1,t2);
			v_double_type op2 = v_load_double(p + j*v_length);
			v_double_type res = v_min_double(op1, op2);
			//m[ind] = min(m[ind], temp[i^1] + temp[j]);
			v_store_double(m + (((i+1)*(i+1))/2) + j*v_length, res);
		}
		for(int j = (i2/v_length)*v_length; j<= i2; j++){
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

bool floyd_warshall_dense(opt_oct_mat_t *oo, double *temp1, double *temp2, int dim, bool is_int){
    double *m = oo->mat;
    int size = 4 * dim * dim;
    int n = 2*dim; 
    double count = 0;
    /******
		Floyd Warshall step
    *******/
    for(int k = 0; k < dim; k++){
	//int pos1 = matpos2(2*k, (2*k)^1);
	int ki = ((((2*k) + 1)*((2*k) + 1))/2);
	int kki = (((((2*k)^1) + 1)*(((2*k)^1) + 1))/2);
	int pos1 = ((2*k)^1) + ki;
	//int pos2 = matpos2((2*k)^1, 2*k);
	int pos2 = (2*k) + kki;
	/******
		For k-th iteration update 2k and (2k+1)-th
		row and column first. This part is non-vectorized.
	*******/
	for(int i = 2*k + 2; i < n;i++){
		//int ind1 = matpos2(i,((2*k)^1));
		int ind1 = ((2*k)^1) + (((i+1)*(i+1))/2);
		//int ind2 = matpos2(i,2*k);
		int ind2 = (2*k) + (((i+1)*(i+1))/2);
		//int ind2 = n*i + ((2*k)^1);
		//int ind1 = n*i + (2*k);
		//m[ind2] = min(m[ind2], m[ind1] + m[n*(2*k) + ((2*k)^1)] );
		m[ind1] = min(m[ind1], m[ind2] + m[pos1] );
                count = count + 2;
		temp2[i^1] = m[ind1];
	}


	for(int i = 2*k + 2; i < n; i++){
		//int ind1 = matpos2(i,((2*k)^1));
		int ind1 = ((2*k)^1) + (((i+1)*(i+1))/2);
		//int ind2 = matpos2(i,2*k);
		int ind2 = (2*k) + (((i+1)*(i+1))/2);
		//int ind2 = n*i + ((2*k)^1);
		//int ind1 = n*i + (2*k);
		//m[ind1] = min(m[ind1], m[ind2] + m[n*((2*k)^1) + (2*k)] );
		m[ind2] = min(m[ind2], m[ind1] + m[pos2] );
		temp1[i^1] = m[ind2];
		count = count + 2;
	}

	for(int j = 0; j < (2*k); j++){
		//int ind3 = matpos2((2*k)^1,j);
		int ind3 = j + kki;
		//int ind4 = matpos2( 2*k,j);
		int ind4 = j + ki;
		//m[n*((2*k)^1) + j] = min(m[n*((2*k)^1) + j], m[n*((2*k)^1) + 2*k] + m[n*(2*k) + j]);
		m[ind3] = min(m[ind3], m[pos2] + m[ind4]);
		count = count + 2;
	}
	for(int j = 0; j < (2*k); j++){
		//int ind3 = matpos2((2*k)^1,j);
		int ind3 = j + kki;
		//int ind4 = matpos2(2*k,j);
		int ind4 = j + ki;
		//m[n*2*k + j] = min(m[n*2*k + j], m[n*2*k + ((2*k)^1)] + m[n*((2*k)^1) + j]);
		m[ind4] = min(m[ind4], m[pos1] + m[ind3]);
		count = count + 2;
	}
	double *p1 = m + kki;
	double *p2 = m + ki;
	int l = (2*k + 2);
	int mod = l%v_length;
	if(mod){
		l = l + (v_length - mod);
	}
	/*******
		This is the vectorized main loop. Apply k-th iteration on
		rest of elements, which is equivalent to 2k and (2k+1)-th
		iteration in APRON strong closure algorithm.
	********/
	for(int i = 0; i < 2*k; i++){
		int i2 = (i%2==0) ? (i + 1): i;
		int br = i2 < 2*k ? i2 : 2*k - 1;
		//int ind1 = matpos2(i,2*k);
		int ind1 = (i^1) + kki;
		//int ind2 = matpos2(i, ((2*k)^1));
		int ind2 = (i^1) + ki;
		//double t1 = m[n*(2*k) + (i^1)];
		//double t2 = m[n*((2*k)^1) + (i^1)];
		double ft1 = m[ind2];
		double ft2 = m[ind1];
		v_double_type t1 = v_set1_double(m[ind2]);
		v_double_type t2 = v_set1_double(m[ind1]);
		
		int b = min(l,i2);
		double *p = m + (((i+1)*(i+1))/2);
		for(int j = 0; j <br/v_length; j++){
				
			v_double_type t3 = v_load_double(p1 + j*v_length);
			v_double_type op1 = v_add_double(t1,t3);
			//double op2 = t2 + m[ind4];
			v_double_type t4 = v_load_double(p2 + j*v_length);
			v_double_type op2 = v_add_double(t2,t4);
			//double op3 = min(op1, op2);
			v_double_type op3 = v_min_double(op1,op2);
			v_double_type op4 = v_load_double(p + j*v_length);
			v_double_type res = v_min_double(op3, op4);
			v_store_double(p + j*v_length, res);
			//m[ind5] = min(m[ind5],op3 );
			count = count + 4;
		}
			
		for(int j = (br/v_length)*v_length; j<=br;j++){
			int ind3 = j + kki;
			int ind4 = j + ki;
			int ind5 = j + (((i+1)*(i+1))/2);
			double op1 = ft1 + m[ind3];
			double op2 = ft2 + m[ind4];
			double op3 = min(op1, op2);
			m[ind5] = min(m[ind5],op3 );
		}
			
		for(int j = 2*k + 2; j <= b; j++){
			int ind5 = j + (((i+1)*(i+1))/2);
			double op1 = ft1 + temp1[j];
			double op2 = ft2 + temp2[j];
			double op3 = min(op1, op2);
			m[ind5] = min(m[ind5],op3 );
		}

		if(b < i2){
			for(int j = b/v_length; j <i2/v_length; j++){
				v_double_type t3 = v_load_double(temp1 + j*v_length);
				v_double_type op1 = v_add_double(t1,t3);
				//double op2 = t2 + temp2[j^1];
				v_double_type t4 = v_load_double(temp2 + j*v_length);
				v_double_type op2 = v_add_double(t2,t4);
				//double op3 = min(op1, op2);
				v_double_type op3 = v_min_double(op1, op2);
				v_double_type op4 = v_load_double(p + j*v_length);
				//m[ind5] = min(m[ind5],op3 );
				v_double_type res = v_min_double(op3, op4);
				v_store_double(p + j*v_length, res);
				count = count + 4;
			}
			for(int j = (i2/v_length)*v_length; j<=i2; j++){
				int ind5 = j + (((i+1)*(i+1))/2);
				double op1 = ft1 + temp1[j];
				double op2 = ft2 + temp2[j];
				double op3 = min(op1, op2);
				m[ind5] = min(m[ind5],op3 );
			}
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
		double ft1 = m[ind1];
		double ft2 = m[ind2];
		
		v_double_type t1 = v_set1_double(m[ind1]);
		v_double_type t2 = v_set1_double(m[ind2]);
		int b = min(l,i2);
		
		double *p = m + (((i+1)*(i+1))/2);
	
		for(int j = 0; j < br/v_length; j++){
			//int ind3 = j + kki;
			//int ind4 = matpos2( 2*k,j);
			//int ind4 = j + ki;
			//int ind5 = matpos2(i,j);
			//int ind5 = j + (((i+1)*(i+1))/2);
			//int ind2 = matpos2(k,j);
			//double op1 = t1 + m[n*((2*k)^1) + j];
			//double op2 = t2 + m[n*(2*k) + j];
			//double op1 = t1 + m[ind3];
			v_double_type t3 = v_load_double(p1 + j*v_length);
			v_double_type op1 = v_add_double(t1,t3);
			//double op2 = t2 + m[ind4];
			v_double_type t4 = v_load_double(p2 + j*v_length);
			v_double_type op2 = v_add_double(t2,t4);
			//double op3 = min(op1, op2);
			v_double_type op3 = v_min_double(op1,op2);
			v_double_type op4 = v_load_double(p + j*v_length);
			//m[ind5] = min(m[ind5],op3 );
			v_double_type res = v_min_double(op3,op4);
			v_store_double(p + j*v_length,res);
			count = count + 4;
		}
		for(int j = (br/v_length)*v_length; j<=br; j++){
			int ind3 = j + kki;
			int ind4 = j + ki;
			int ind5 = j + (((i+1)*(i+1))/2);
			double op1 = ft1 + m[ind3];
			double op2 = ft2 + m[ind4];
			double op3 = min(op1, op2);
			m[ind5] = min(m[ind5],op3 );
			count = count + 4;
		}
		for(int j = 2*k + 2; j <= b; j++){
			int ind5 = j + (((i+1)*(i+1))/2);
			double op1 = ft1 + temp1[j];
			double op2 = ft2 + temp2[j];
			double op3 = min(op1, op2);
			m[ind5] = min(m[ind5],op3 );
		}
		if(b < i2){
			
			for(int j = b/v_length; j < i2/v_length; j++){
				//int ind3 = matpos2(n,(2*k)^1,j);
				//int ind4 = matpos2(n, 2*k,j);
				//int ind2 = matpos2(k,j);
				//double op1 = t1 + m[n*(j^1) + 2*k];
				//double op2 = t2 + m[n*(j^1) + ((2*k)^1)];
				//int ind5 = matpos2(i,j);
				//int ind5 = j + (((i+1)*(i+1))/2);
				//double op1 = t1 + temp1[j]		
				v_double_type t3 = v_load_double(temp1 + j*v_length);
				v_double_type op1 = v_add_double(t1,t3);
				//double op2 = t2 + temp2[j];
				v_double_type t4 = v_load_double(temp2 + j*v_length);
				v_double_type op2 = v_add_double(t2,t4);
				//double op3 = min(op1, op2);
				v_double_type op3 = v_min_double(op1,op2);
				v_double_type op4 = v_load_double(p + j*v_length);
				//m[ind5] = min(m[ind5],op3 );
				v_double_type res = v_min_double(op3,op4);
				v_store_double(p + j*v_length, res);
				count = count + 4;
			}
			for(int j = (i2/v_length)*v_length; j <=i2; j++){
				int ind5 = j + (((i+1)*(i+1))/2);
				double op1 = ft1 + temp1[j];
				double op2 = ft2 + temp2[j];
				double op3 = min(op1, op2);
				m[ind5] = min(m[ind5],op3 );
			}
		}
	}
	
    }
}

bool strong_closure_dense(opt_oct_mat_t *oo, double *temp1, double *temp2, int dim, bool is_int){
    floyd_warshall_dense(oo,temp1,temp2,dim,is_int);
    int n = 2*dim;
    oo->nni = 2*dim*(dim+1);
    if(is_int){
	return strengthning_int_dense(oo,temp1,n);
    }
    else{
    	return strengthning_dense(oo,temp1,n);
    }
}


void print_dense(double *m, int dim){
   int n = 2*dim;
    for (int i = 0; i < 2*dim; ++i){
	int ni = n*i;
        for (int j = 0; j < 2*dim; ++j){
            printf("%.15f \t", m[ni + j]);
        }
        printf("\n");
    }
    printf("\n\n");
}



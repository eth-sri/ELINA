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
		__m256d t1 = _mm256_set1_pd(temp[i^1]);
		double *p = m + (((i+1)*(i+1))/2);
		for(int j = 0; j < (i2/4); j++){
			//int ind = j*8 + (((i+1)*(i+1))/2);
			__m256d t2 = _mm256_loadu_pd(temp + j*4);
			__m256d op1 = _mm256_add_pd(t1,t2);
			__m256d op2 = _mm256_loadu_pd(p + j*4);
			__m256d res = _mm256_min_pd(op1, op2);
			//m[ind] = min(m[ind], temp[i^1] + temp[j]);
			_mm256_storeu_pd(m + (((i+1)*(i+1))/2) + j*4, res);
		}
		for(int j = (i2/4)*4; j<= i2; j++){
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
		__m256d t1 = _mm256_set1_pd(temp[i^1]);
		double *p = m + (((i+1)*(i+1))/2);
		for(int j = 0; j < (i2/4); j++){
			//int ind = j*8 + (((i+1)*(i+1))/2);
			__m256d t2 = _mm256_loadu_pd(temp + j*4);
			__m256d op1 = _mm256_add_pd(t1,t2);
			__m256d op2 = _mm256_loadu_pd(p + j*4);
			__m256d res = _mm256_min_pd(op1, op2);
			//m[ind] = min(m[ind], temp[i^1] + temp[j]);
			_mm256_storeu_pd(m + (((i+1)*(i+1))/2) + j*4, res);
		}
		for(int j = (i2/4)*4; j<= i2; j++){
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
	int mod = l%4;
	if(mod){
		l = l + (4 - mod);
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
		__m256d t1 = _mm256_set1_pd(m[ind2]);
		__m256d t2 = _mm256_set1_pd(m[ind1]);
		
		int b = min(l,i2);
		double *p = m + (((i+1)*(i+1))/2);
		for(int j = 0; j <br/4; j++){
				
			__m256d t3 = _mm256_loadu_pd(p1 + j*4);
			__m256d op1 = _mm256_add_pd(t1,t3);
			//double op2 = t2 + m[ind4];
			__m256d t4 = _mm256_loadu_pd(p2 + j*4);
			__m256d op2 = _mm256_add_pd(t2,t4);
			//double op3 = min(op1, op2);
			__m256d op3 = _mm256_min_pd(op1,op2);
			__m256d op4 = _mm256_loadu_pd(p + j*4);
			__m256d res = _mm256_min_pd(op3, op4);
			_mm256_storeu_pd(p + j*4, res);
			//m[ind5] = min(m[ind5],op3 );
			count = count + 4;
		}
			
		for(int j = (br/4)*4; j<=br;j++){
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
			for(int j = b/4; j <i2/4; j++){
				__m256d t3 = _mm256_loadu_pd(temp1 + j*4);
				__m256d op1 = _mm256_add_pd(t1,t3);
				//double op2 = t2 + temp2[j^1];
				__m256d t4 = _mm256_loadu_pd(temp2 + j*4);
				__m256d op2 = _mm256_add_pd(t2,t4);
				//double op3 = min(op1, op2);
				__m256d op3 = _mm256_min_pd(op1, op2);
				__m256d op4 = _mm256_loadu_pd(p + j*4);
				//m[ind5] = min(m[ind5],op3 );
				__m256d res = _mm256_min_pd(op3, op4);
				_mm256_storeu_pd(p + j*4, res);
				count = count + 4;
			}
			for(int j = (i2/4)*4; j<=i2; j++){
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
		
		__m256d t1 = _mm256_set1_pd(m[ind1]);
		__m256d t2 = _mm256_set1_pd(m[ind2]);
		int b = min(l,i2);
		
		double *p = m + (((i+1)*(i+1))/2);
	
		for(int j = 0; j < br/4; j++){
			//int ind3 = j + kki;
			//int ind4 = matpos2( 2*k,j);
			//int ind4 = j + ki;
			//int ind5 = matpos2(i,j);
			//int ind5 = j + (((i+1)*(i+1))/2);
			//int ind2 = matpos2(k,j);
			//double op1 = t1 + m[n*((2*k)^1) + j];
			//double op2 = t2 + m[n*(2*k) + j];
			//double op1 = t1 + m[ind3];
			__m256d t3 = _mm256_loadu_pd(p1 + j*4);
			__m256d op1 = _mm256_add_pd(t1,t3);
			//double op2 = t2 + m[ind4];
			__m256d t4 = _mm256_loadu_pd(p2 + j*4);
			__m256d op2 = _mm256_add_pd(t2,t4);
			//double op3 = min(op1, op2);
			__m256d op3 = _mm256_min_pd(op1,op2);
			__m256d op4 = _mm256_loadu_pd(p + j*4);
			//m[ind5] = min(m[ind5],op3 );
			__m256d res = _mm256_min_pd(op3,op4);
			_mm256_storeu_pd(p + j*4,res);
			count = count + 4;
		}
		for(int j = (br/4)*4; j<=br; j++){
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
			
			for(int j = b/4; j < i2/4; j++){
				//int ind3 = matpos2(n,(2*k)^1,j);
				//int ind4 = matpos2(n, 2*k,j);
				//int ind2 = matpos2(k,j);
				//double op1 = t1 + m[n*(j^1) + 2*k];
				//double op2 = t2 + m[n*(j^1) + ((2*k)^1)];
				//int ind5 = matpos2(i,j);
				//int ind5 = j + (((i+1)*(i+1))/2);
				//double op1 = t1 + temp1[j]		
				__m256d t3 = _mm256_loadu_pd(temp1 + j*4);
				__m256d op1 = _mm256_add_pd(t1,t3);
				//double op2 = t2 + temp2[j];
				__m256d t4 = _mm256_loadu_pd(temp2 + j*4);
				__m256d op2 = _mm256_add_pd(t2,t4);
				//double op3 = min(op1, op2);
				__m256d op3 = _mm256_min_pd(op1,op2);
				__m256d op4 = _mm256_loadu_pd(p + j*4);
				//m[ind5] = min(m[ind5],op3 );
				__m256d res = _mm256_min_pd(op3,op4);
				_mm256_storeu_pd(p + j*4, res);
				count = count + 4;
			}
			for(int j = (i2/4)*4; j <=i2; j++){
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



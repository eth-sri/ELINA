#include "maxpool_convex_hull.h"


void populate_matrix_maxpool(dd_MatrixPtr M, double *lb, double *ub,size_t pool_size, size_t val){
	size_t *arr = (size_t *)malloc((pool_size-1)*sizeof(size_t));
	size_t i,j=0;
	for(i=0; i < pool_size; i++){
		if(i!=val){
			arr[j] = i;
			j++;
		}
	}

	size_t rowsize = 3*pool_size + 1;
	size_t colsize = pool_size + 2;
	for(i=0; i < rowsize; i++){
		for(j=0; j < colsize; j++){
			dd_set_d(M->matrix[i][j],0.0);
		}
	}

	for(i=0; i < pool_size; i++){
		dd_set_d(M->matrix[2*i][0],-lb[i]);
		dd_set_d(M->matrix[2*i][i+1],1);
		dd_set_d(M->matrix[2*i+1][0], ub[i]);
		dd_set_d(M->matrix[2*i+1][i+1], -1);
	}
	size_t start = 2*pool_size;
	for(i=0; i < pool_size-1; i++){
		dd_set_d(M->matrix[start][val+1],1);
		dd_set_d(M->matrix[start][arr[i]+1],-1);
		start++;
	}

        dd_set_d(M->matrix[start][val+1],1); 
	dd_set_d(M->matrix[start][pool_size+1],-1);
	start++;
	dd_set_d(M->matrix[start][val+1],-1); 
	dd_set_d(M->matrix[start][pool_size+1],1);
	
	M->representation=dd_Inequality;
	free(arr);
}

dd_MatrixPtr maxpool_deeppoly_approx(double *lb, double *ub,size_t pool_size){
	dd_set_global_constants();
	dd_MatrixPtr *H = (dd_MatrixPtr *)malloc(pool_size*sizeof(dd_MatrixPtr)); 
	dd_MatrixPtr A = NULL;
	dd_MatrixPtr *V = (dd_MatrixPtr *)malloc(pool_size*sizeof(dd_MatrixPtr));
	dd_PolyhedraPtr *Poly = (dd_PolyhedraPtr*)malloc(pool_size*sizeof(dd_PolyhedraPtr));
	dd_PolyhedraPtr P = NULL;
	dd_ErrorType *errs = (dd_ErrorType*)malloc(pool_size*sizeof(dd_ErrorType));
	dd_ErrorType err;
	dd_rowrange m = 3*pool_size + 1;
	dd_colrange d = pool_size+2;
	size_t i ;
	for(i=0; i < pool_size; i++){
		H[i] = dd_CreateMatrix(m,d);
		populate_matrix_maxpool(H[i],lb,ub,pool_size,i);
		Poly[i] = dd_DDMatrix2Poly(H[i],errs+i);
		if(errs[i]!=dd_NoError){
			goto _L99;
		}
		V[i] = dd_CopyGenerators(Poly[i]);
	}


	//dd_WriteMatrix(stdout, M1);

	P = dd_DDMatrix2Poly(V[0],&err);
	if (err!=dd_NoError) goto _L99;
	for(i=1; i < pool_size; i++){
		dd_DDInputAppend(&P,V[i],&err);
		if(err!=dd_NoError){
			goto _L99;
		}
	}
	
	A = dd_CopyInequalities(P);

_L99:
	for(i=0; i < pool_size; i++){
		if(H[i]!=NULL){
			dd_FreeMatrix(H[i]);
		}
		if(V[i]!=NULL){
			dd_FreeMatrix(V[i]);
		}
		if(Poly[i]!=NULL){
			dd_FreePolyhedra(Poly[i]);
		}
	}
	if(P!=NULL){
		dd_FreePolyhedra(P);
	}
	free(H);
	free(V);
	free(Poly);
	free(errs);
	dd_free_global_constants();
	return A;
}

//int main(){
//	size_t pool_size = 10;
//	double lb[pool_size], ub[pool_size];
//	size_t i;
//	for(i=0; i < pool_size; i++){
//		lb[i] = -1;
//		ub[i] = 1;
//	}
//	dd_MatrixPtr res = maxpool_deeppoly_approx(lb,ub,pool_size);
//	dd_WriteMatrix(stdout,res);
//	dd_FreeMatrix(res);
//}

#include "expr.h"

void elina_double_interval_add_expr_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup, double inf, double sup, double inf_expr, double sup_expr){
	*res_inf = inf + inf_expr;
	*res_sup = sup + sup_expr;
	double maxA = fmax(fabs(inf_expr),fabs(sup_expr));
	double tmp1, tmp2;
	elina_double_interval_mul(&tmp1,&tmp2, inf, sup, maxA*pr->ulp, maxA*pr->ulp);
	*res_inf += tmp1;
	*res_sup += tmp2;
}


void elina_double_interval_add_cst_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup, double inf, double sup, double inf_expr, double sup_expr){
	elina_double_interval_add_expr_coeff(pr, res_inf, res_sup, inf, sup, inf_expr, sup_expr);
	*res_inf += pr->min_denormal;
	*res_sup += pr->min_denormal;	
}


void elina_double_interval_mul_expr_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup, double inf, double sup, double inf_expr, double sup_expr){
	elina_double_interval_mul(res_inf,res_sup,inf,sup,inf_expr,sup_expr);
	double maxA = fmax(fabs(inf_expr),fabs(sup_expr));
	double tmp1, tmp2;
	elina_double_interval_mul(&tmp1,&tmp2, inf, sup, maxA*pr->ulp, maxA*pr->ulp);
	*res_inf += tmp1;
	*res_sup += tmp2;
}

void elina_double_interval_mul_cst_coeff(fppoly_internal_t *pr, double * res_inf, double *res_sup, double inf, double sup, double inf_expr, double sup_expr){
	elina_double_interval_mul_expr_coeff(pr, res_inf, res_sup, inf, sup, inf_expr, sup_expr);
	*res_inf += pr->min_denormal;
	*res_sup += pr->min_denormal;	
}


void expr_fprint(FILE * stream, expr_t *expr){
	if((expr->inf_coeff==NULL) || (expr->sup_coeff==NULL)){
		fprintf(stdout,"+ [%g, %g]\n",-expr->inf_cst,expr->sup_cst);
		return;
	}
	size_t size = expr->size;
	size_t i;
	for(i=0; i < size; i++){
		if(i==0){
			if(expr->type==DENSE){
				fprintf(stream, "[%g, %g]x0 ", -expr->inf_coeff[0],expr->sup_coeff[0]);
			}
			else{
				fprintf(stream, "[%g, %g]x%zu ", -expr->inf_coeff[0],expr->sup_coeff[0],expr->dim[0]);
			}
		}
		
		else{
			if(expr->type==DENSE){
				fprintf(stream,"+ [%g, %g]x%zu ",-expr->inf_coeff[i],expr->sup_coeff[i],i);
			}
			else{
				fprintf(stream,"+ [%g, %g]x%zu ",-expr->inf_coeff[i],expr->sup_coeff[i],expr->dim[i]);
			}
		}
	}
	
	fprintf(stdout,"+ [%g, %g]\n",-expr->inf_cst,expr->sup_cst);
	
}


void expr_print(expr_t * expr){
	expr_fprint(stdout, expr);	
}

expr_t * alloc_expr(void){
	expr_t *expr = (expr_t *)malloc(sizeof(expr_t));
	expr->inf_coeff = NULL;
	expr->sup_coeff = NULL;
	expr->dim = NULL;
	return expr;
}

expr_t * create_dense_expr(double *coeff, double cst, size_t size){
	expr_t *expr = (expr_t *)malloc(sizeof(expr_t));
	expr->inf_coeff = (double *)malloc(size*sizeof(double));
	expr->sup_coeff = (double *)malloc(size*sizeof(double));
	expr->dim= NULL;
	size_t i;
	expr->size = size;
	expr->inf_cst = -cst;
	expr->sup_cst = cst;
	expr->type = DENSE;
	for(i=0; i < size; i++){
		expr->inf_coeff[i] = -coeff[i];
		expr->sup_coeff[i] = coeff[i];
	}
	return expr;
}


expr_t * create_cst_expr(double l, double u){
	expr_t *expr = (expr_t*)malloc(sizeof(expr_t));
	expr->inf_coeff = NULL;
	expr->sup_coeff = NULL;
	expr->dim = NULL;
	expr->type = SPARSE;
	expr->size = 0;
	expr->inf_cst = l;
	expr->sup_cst = u;
	return expr;
}

expr_t * create_sparse_expr(double *coeff, double cst, size_t *dim, size_t size){
	expr_t *expr = (expr_t *)malloc(sizeof(expr_t));
	if(size>0){
		expr->inf_coeff = (double *)malloc(size*sizeof(double));
		expr->sup_coeff = (double *)malloc(size*sizeof(double));
		expr->dim = (size_t *)malloc(size*sizeof(size_t));
	}
	else{
		expr->inf_coeff = NULL;
		expr->sup_coeff = NULL;
		expr->dim = NULL;
	}
	size_t i;
	expr->size = size;
	expr->inf_cst = -cst;
	expr->sup_cst = cst;
	expr->type = SPARSE;
	for(i=0; i < size; i++){
		expr->inf_coeff[i] = -coeff[i];
		expr->sup_coeff[i] = coeff[i];
		expr->dim[i] = dim[i];
	}
	return expr;
}


void free_expr(expr_t *expr){
	if(expr->inf_coeff){
		free(expr->inf_coeff);
		expr->inf_coeff = NULL;
	}
	if(expr->sup_coeff){
		free(expr->sup_coeff);
		expr->sup_coeff = NULL;
	}
	if(expr->type==SPARSE && expr->dim){
		free(expr->dim);
	}
	expr->dim = NULL;
	free(expr);
	expr = NULL;  
}

expr_t * copy_cst_expr(expr_t *src){
	expr_t *dst = (expr_t *)malloc(sizeof(expr_t));
	dst->inf_coeff = NULL;
	dst->sup_coeff = NULL;
	dst->inf_cst = src->inf_cst;
	dst->sup_cst = src->sup_cst; 
	dst->type = src->type;
	dst->dim = NULL;
	dst->size = src->size; 
	return dst;
}



expr_t * copy_expr(expr_t *src){
	expr_t *dst = (expr_t *)malloc(sizeof(expr_t));
	dst->inf_coeff = (double *)malloc(src->size*sizeof(double));
	dst->sup_coeff = (double *)malloc(src->size*sizeof(double));
	
	size_t i;
	dst->inf_cst = src->inf_cst;
	dst->sup_cst = src->sup_cst; 
	dst->type = src->type;
	for(i=0; i < src->size; i++){
		dst->inf_coeff[i] = src->inf_coeff[i];
		dst->sup_coeff[i] = src->sup_coeff[i];
	}
	if(src->type==SPARSE){
		dst->dim = (size_t *)malloc(src->size*sizeof(size_t));
		for(i=0; i < src->size; i++){
			dst->dim[i] = src->dim[i];
		}
	}
	dst->size = src->size; 
	return dst;
}

expr_t* concretize_dense_sub_expr(fppoly_internal_t *pr, expr_t * expr, double *inf, double *sup, size_t start, size_t size){
	expr_t * res = (expr_t *)malloc(sizeof(expr_t));
	res->inf_coeff = (double *)malloc(start*sizeof(double));
	res->sup_coeff = (double *)malloc(start*sizeof(double));
	size_t i;
	res->inf_cst = expr->inf_cst;
	res->sup_cst = expr->sup_cst;
	res->type = expr->type;
	for(i=0; i < start; i++){
		res->inf_coeff[i] = expr->inf_coeff[i];
		res->sup_coeff[i] = expr->sup_coeff[i];
	}
	for(i=start; i< size;i++){
		double tmp1,tmp2;
		elina_double_interval_mul_expr_coeff(pr,&tmp1,&tmp2,inf[i-start],sup[i-start],expr->inf_coeff[i],expr->sup_coeff[i]);
		res->inf_cst += tmp1;
		res->sup_cst += tmp2;
	}
	res->size = start;
	return res;
}


void merge_sparse_expr(expr_t *expr, size_t l, size_t m, size_t r) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    /* create temp arrays */
    size_t *L = (size_t *)malloc(n1*sizeof(size_t));
    size_t *R = (size_t *)malloc(n2*sizeof(size_t));
    double *L2 = (double *)malloc(n1*sizeof(double));
    double *R2 = (double *)malloc(n2*sizeof(double));
    double *L3 = (double *)malloc(n1*sizeof(double));
    double *R3 = (double *)malloc(n2*sizeof(double));
    
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++) {
        L[i] = expr->dim[l + i];
        L2[i] = expr->inf_coeff[l + i];
	L3[i] = expr->sup_coeff[l + i];
    }
    for (j = 0; j < n2; j++) {
        R[j] = expr->dim[m + 1 + j];
        R2[j] = expr->inf_coeff[m + 1 + j];
	R3[j] = expr->sup_coeff[m + 1 + j];
    }

    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            expr->dim[k] = L[i];
            expr->inf_coeff[k] = L2[i];
	    expr->sup_coeff[k] = L3[i];
            i++;
        } else {
            expr->dim[k] = R[j];
            expr->inf_coeff[k] = R2[j];
	    expr->sup_coeff[k] = R3[j];
            j++;
        }
        k++;
    }

    /* Copy the remaining elements of L[], if there
       are any */
    while (i < n1) {
        expr->dim[k] = L[i];
        expr->inf_coeff[k] = L2[i];
	expr->sup_coeff[k] = L3[i];
        i++;
        k++;
    }

    /* Copy the remaining elements of R[], if there
       are any */
    while (j < n2) {
        expr->dim[k] = R[j];
        expr->inf_coeff[k] = R2[j];
	expr->sup_coeff[k] = R3[j];
        j++;
        k++;
    }
    free(L);
    free(R);
    free(L2);
    free(R2);
    free(L3);
    free(R3);
}


/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void merge_sort_sparse_expr(expr_t *expr, size_t l, size_t r) {
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        size_t m = l + (r - l) / 2;

        // Sort first and second halves
        merge_sort_sparse_expr(expr, l, m);
        merge_sort_sparse_expr(expr, m + 1, r);

        merge_sparse_expr(expr, l, m, r);
    }
}

void sort_sparse_expr(expr_t *expr){
	merge_sort_sparse_expr(expr,0,expr->size-1);
}


expr_t * multiply_expr(fppoly_internal_t *pr, expr_t *expr, double mul_inf, double mul_sup){
	expr_t * res = alloc_expr();
	if(expr->size > 0){
		res->inf_coeff = malloc(expr->size*sizeof(double));
		res->sup_coeff = malloc(expr->size*sizeof(double));
	}
	else{
		res->inf_coeff = NULL;		
		res->sup_coeff = NULL;
	}
	res->type = expr->type;
	size_t i;
	for(i=0; i < expr->size; i++){
		//res->coeff[i] = mul_coeff*expr->coeff[i];
		elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],mul_inf,mul_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
		
	}
	if(expr->type==SPARSE){
		if(expr->size>0){
			res->dim = (size_t*)malloc(expr->size*sizeof(size_t));
			for(i=0; i < expr->size; i++){
				res->dim[i] = expr->dim[i];
			}
		}
		else{
			res->dim = NULL;
		}
	}
	res->size = expr->size;
	
	elina_double_interval_mul_cst_coeff(pr,&res->inf_cst,&res->sup_cst,mul_inf,mul_sup,expr->inf_cst,expr->sup_cst);
	
	//res->cst = mul_coeff*expr->cst;
	return res;
}


expr_t * multiply_cst_expr(fppoly_internal_t *pr, expr_t *expr, double mul_inf, double mul_sup){
	expr_t * res = alloc_expr();
	res->inf_coeff = NULL;		
	res->sup_coeff = NULL;
	res->dim = NULL;
	res->type = expr->type;
	res->size = expr->size;
	elina_double_interval_mul_cst_coeff(pr,&res->inf_cst,&res->sup_cst,mul_inf,mul_sup,expr->inf_cst,expr->sup_cst);
	//res->cst = mul_coeff*expr->cst;
	return res;
}






void add_cst_expr(fppoly_internal_t *pr, expr_t * exprA, expr_t *exprB){
	double maxA = fmax(fabs(exprA->inf_cst),fabs(exprA->sup_cst));
	double maxB = fmax(fabs(exprB->inf_cst),fabs(exprB->sup_cst));
	exprA->inf_cst = exprA->inf_cst + exprB->inf_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal; 
	exprA->sup_cst = exprA->sup_cst + exprB->sup_cst + (maxA + maxB)*pr->ulp + pr->min_denormal; 
	return;
}

//A = A + B
void add_expr(fppoly_internal_t *pr,expr_t * exprA, expr_t * exprB){
	//
	size_t sizeB = exprB->size;
	if(sizeB==0){
		double maxA = fmax(fabs(exprA->inf_cst),fabs(exprA->sup_cst));
		double maxB = fmax(fabs(exprB->inf_cst),fabs(exprB->sup_cst));
		exprA->inf_cst = exprA->inf_cst + exprB->inf_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal; 
		exprA->sup_cst = exprA->sup_cst + exprB->sup_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal; 
		return;
	}
	size_t i;
	if(exprA->size==0){
		
		exprA->size = exprB->size;
		double maxA = fmax(fabs(exprA->inf_cst),fabs(exprA->sup_cst));
		double maxB = fmax(fabs(exprB->inf_cst),fabs(exprB->sup_cst));
		exprA->inf_cst += exprB->inf_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal;
		exprA->sup_cst += exprB->sup_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal;
		exprA->inf_coeff = (double *)malloc(sizeB*sizeof(double));
		exprA->sup_coeff = (double *)malloc(sizeB*sizeof(double));
		for(i=0; i < sizeB; i++){
			exprA->inf_coeff[i] = exprB->inf_coeff[i];
			exprA->sup_coeff[i] = exprB->sup_coeff[i];
		} 
		exprA->type = exprB->type;
		if(exprA->type==SPARSE){
			exprA->dim = (size_t *)malloc(sizeB*sizeof(size_t));
			for(i=0; i < sizeB; i++){
				exprA->dim[i] = exprB->dim[i];
			}
		} 
		
		return;
	}
	else{
		size_t sizeA = exprA->size;
		assert(sizeA==sizeB);
		double maxA = fmax(fabs(exprA->inf_cst),fabs(exprA->sup_cst));
		double maxB = fmax(fabs(exprB->inf_cst),fabs(exprB->sup_cst));
		exprA->inf_cst += exprB->inf_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal;
		exprA->sup_cst += exprB->sup_cst  + (maxA + maxB)*pr->ulp + pr->min_denormal;
		if(exprA->type==DENSE){
			if(exprB->type==DENSE){
				//printf("AA\n");
				//fflush(stdout);
				for(i=0; i < sizeB; i++){
					maxA = fmax(fabs(exprA->inf_coeff[i]),fabs(exprA->sup_coeff[i]));
					maxB = fmax(fabs(exprB->inf_coeff[i]),fabs(exprB->sup_coeff[i]));
					exprA->inf_coeff[i] = exprA->inf_coeff[i] + exprB->inf_coeff[i] + (maxA + maxB)*pr->ulp;
					exprA->sup_coeff[i] = exprA->sup_coeff[i] + exprB->sup_coeff[i] + (maxA + maxB)*pr->ulp;
				}
			}
			else{
				//printf("AB\n");	
				//fflush(stdout);
				size_t k = 0;
				for(i=0; i < sizeA; i++){
					if(k < sizeB && exprB->dim[k]==i){
						maxA = fmax(fabs(exprA->inf_coeff[i]),fabs(exprA->sup_coeff[i]));
						maxB = fmax(fabs(exprB->inf_coeff[k]),fabs(exprB->sup_coeff[k]));
						exprA->inf_coeff[i] = exprA->inf_coeff[i] + exprB->inf_coeff[k] + (maxA + maxB)*pr->ulp ;
						exprA->sup_coeff[i] = exprA->sup_coeff[i] + exprB->sup_coeff[k] + (maxA + maxB)*pr->ulp;
						k++;
					}
				}
			}
		}
		else{
			size_t sizeB = exprB->size;
			size_t k;
			double * new_inf_coeff;
			double * new_sup_coeff;
			if(exprB->type==DENSE){
				//printf("BB\n");
				//fflush(stdout);
				i=0;
				new_inf_coeff = (double *)malloc(sizeB*sizeof(double));
				new_sup_coeff = (double *)malloc(sizeB*sizeof(double));				
				for(k=0; k < sizeB; k++){
					if(i < sizeA && exprA->dim[i] == k){
						maxA = fmax(fabs(exprA->inf_coeff[i]),fabs(exprA->sup_coeff[i]));
						maxB = fmax(fabs(exprB->inf_coeff[k]),fabs(exprB->sup_coeff[k]));
						new_inf_coeff[k] = exprA->inf_coeff[i] + exprB->inf_coeff[k] + (maxA + maxB)*pr->ulp;
						new_sup_coeff[k] = exprA->sup_coeff[i] + exprB->sup_coeff[k] + (maxA + maxB)*pr->ulp;
						i++;
					}
					else{
						new_inf_coeff[k] = exprB->inf_coeff[k];
						new_sup_coeff[k] = exprB->sup_coeff[k];
					}
				}
				exprA->type = DENSE;
				exprA->size = sizeB;
				free(exprA->dim);
				exprA->dim = NULL;
			}
			else{
				//printf("BC %zu %zu\n",exprA->size, exprB->size);
				//fflush(stdout);
				i=0;
				k=0;
				size_t l = 0;
				//printf("before sort\n");
				//expr_print(exprA);
				//fflush(stdout);
				//if(exprA->size>0){
				//	sort_sparse_expr(exprA);
				//}
				//printf("after sort\n");
				//expr_print(exprA);
				//fflush(stdout);
				
				new_inf_coeff = (double *)malloc((sizeA+sizeB)*sizeof(double));
				new_sup_coeff = (double *)malloc((sizeA+sizeB)*sizeof(double));
				size_t * new_dim = (size_t *)malloc((sizeA+sizeB)*sizeof(size_t));
				while(i < sizeA && k < sizeB){
					if(exprA->dim[i] < exprB->dim[k]){
						new_inf_coeff[l] = exprA->inf_coeff[i];
						new_sup_coeff[l] = exprA->sup_coeff[i];
						new_dim[l] = exprA->dim[i];
						i++;
						
					}
					else if(exprB->dim[k] < exprA->dim[i]){
						new_inf_coeff[l] = exprB->inf_coeff[k];
						new_sup_coeff[l] = exprB->sup_coeff[k];
						new_dim[l] = exprB->dim[k];
						k++;
					}
					else{
						maxA = fmax(fabs(exprA->inf_coeff[i]),fabs(exprA->sup_coeff[i]));
						maxB = fmax(fabs(exprB->inf_coeff[k]),fabs(exprB->sup_coeff[k]));
						new_inf_coeff[l] = exprA->inf_coeff[i] + exprB->inf_coeff[k] + (maxA + maxB)*pr->ulp;
						new_sup_coeff[l] = exprA->sup_coeff[i] + exprB->sup_coeff[k] + (maxA + maxB)*pr->ulp;
						new_dim[l] = exprA->dim[i];
						i++;
						k++;
					}
					l++;
				}
				while(i < sizeA){
					new_inf_coeff[l] = exprA->inf_coeff[i];
					new_sup_coeff[l] = exprA->sup_coeff[i];
					new_dim[l] = exprA->dim[i];
					i++;
					l++;
				}
				while(k < sizeB){
					new_inf_coeff[l] = exprB->inf_coeff[k];
					new_sup_coeff[l] = exprB->sup_coeff[k];
					new_dim[l] = exprB->dim[k];
					k++;
					l++;
				}
				
				new_inf_coeff = (double*)realloc(new_inf_coeff,l*sizeof(double));
				new_sup_coeff = (double*)realloc(new_sup_coeff,l*sizeof(double));
				free(exprA->dim);
				exprA->dim = NULL;
				new_dim = (size_t *)realloc(new_dim,l*sizeof(size_t));
				exprA->dim = new_dim;
				exprA->size = l;
			}
			if(exprA->inf_coeff){
				free(exprA->inf_coeff);
				exprA->inf_coeff = NULL;
			}
			if(exprA->sup_coeff){
				free(exprA->sup_coeff);
				exprA->sup_coeff = NULL;
			}
			exprA->inf_coeff = new_inf_coeff;
			exprA->sup_coeff = new_sup_coeff;
			
		}
	}
}

expr_t * expr_replace_bounds_affine(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons, bool is_lower){
	if(expr->size==0){
		return copy_cst_expr(expr);
	}	
	if(expr->inf_coeff==NULL || expr->sup_coeff==NULL ){
		return alloc_expr();
	}
	
	size_t num_neurons = expr->size;
	size_t i,k;
	expr_t * res;
	if(expr->type==DENSE){
		k = 0;
	}
	else{
		k = expr->dim[0];		
	}
		
	expr_t * mul_expr = NULL;
	neuron_t * neuron_k = neurons[k];
	if(is_lower){
		if(expr->sup_coeff[0] < 0){
			mul_expr = neuron_k->uexpr;
		}
		else if(expr->inf_coeff[0]<0){
			mul_expr = neuron_k->lexpr;
		}
	}
	else{
		if(expr->sup_coeff[0] < 0){
			mul_expr = neuron_k->lexpr;
		}
		else if(expr->inf_coeff[0]<0){
			mul_expr = neuron_k->uexpr;
		}
	}	
	
	if(mul_expr==NULL){
		
		double tmp1=0.0, tmp2=0.0;
		if(expr->inf_coeff[0]!=0 || expr->sup_coeff[0]!=0){
			elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,neuron_k->lb,neuron_k->ub,expr->inf_coeff[0],expr->sup_coeff[0]);
		}
		double coeff[1];
		size_t dim[1];
		coeff[0] = 0;
		dim[0] = 0;
		if(is_lower){
			res = create_sparse_expr(coeff,-tmp1,dim,1);
		}
		else{
			res = create_sparse_expr(coeff,tmp2,dim,1);
		}
	}
	else if(mul_expr->size==0){
			
		res = multiply_cst_expr(pr,mul_expr,expr->inf_coeff[0],expr->sup_coeff[0]);
			
	}
	else{
			
		res = multiply_expr(pr,mul_expr,expr->inf_coeff[0],expr->sup_coeff[0]);
	}
    	
	for(i=1; i < num_neurons; i++){
		if(expr->type==DENSE){
			k = i;
		}
		else{
			k = expr->dim[i];
		}
		neuron_k = neurons[k];
		mul_expr = NULL;
		if(is_lower){
			if(expr->sup_coeff[i] < 0){
				mul_expr = neuron_k->uexpr;
			}
			else if(expr->inf_coeff[i]<0){
				mul_expr = neuron_k->lexpr;
			}
		}
		else{
			if(expr->sup_coeff[i] < 0){
				mul_expr = neuron_k->lexpr;
			}
			else if(expr->inf_coeff[i]<0){
				mul_expr = neuron_k->uexpr;
			}
		}
		if(expr->sup_coeff[i]==0 && expr->inf_coeff[i]==0){
			continue;
		}
		expr_t * tmp_mul_expr = NULL;
		if(expr->sup_coeff[i] < 0 || expr->inf_coeff[i]<0){
			if(mul_expr->size==0){
				tmp_mul_expr = multiply_cst_expr(pr,mul_expr, expr->inf_coeff[i],expr->sup_coeff[i]);
				add_cst_expr(pr,res,tmp_mul_expr);
				free_expr(tmp_mul_expr);
			}
			else{
				tmp_mul_expr = multiply_expr(pr, mul_expr, expr->inf_coeff[i],expr->sup_coeff[i]);
				
				add_expr(pr,res,tmp_mul_expr);
				free_expr(tmp_mul_expr);
				
			}
		}
		else{
			//printf("WTF2 %g %g\n",expr->inf_coeff[i],expr->sup_coeff[i]);
			//fflush(stdout);
			double tmp1, tmp2;
			elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,neuron_k->lb,neuron_k->ub,expr->inf_coeff[i],expr->sup_coeff[i]);
			if(is_lower){
				res->inf_cst = res->inf_cst + tmp1;
				res->sup_cst = res->sup_cst - tmp1;
			}
			else{
				res->inf_cst = res->inf_cst - tmp2;
				res->sup_cst = res->sup_cst + tmp2;
			}
		}
		
	}
	res->inf_cst = res->inf_cst + expr->inf_cst; 
	res->sup_cst = res->sup_cst + expr->sup_cst; 
	return res;
}

expr_t * lexpr_replace_bounds_affine(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
	return expr_replace_bounds_affine(pr, expr, neurons, true);
}

expr_t * uexpr_replace_bounds_affine(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
	return expr_replace_bounds_affine(pr, expr, neurons, false);
}

expr_t * expr_replace_bounds_activation(fppoly_internal_t * pr, expr_t * expr, neuron_t ** neurons, bool is_lower){
	size_t num_neurons = expr->size;
	size_t i,k;
	expr_t * res = alloc_expr();  
	res->inf_coeff = (double *)malloc(num_neurons*sizeof(double));
	res->sup_coeff = (double *)malloc(num_neurons*sizeof(double));
	res->inf_cst = expr->inf_cst;
	res->sup_cst = expr->sup_cst;
	res->type = expr->type;
	res->size = num_neurons;  

	for(i = 0; i < num_neurons; i++){
		if(expr->type==DENSE){
			k = i;
		}
		else{
			k = expr->dim[i];
		}
		neuron_t *neuron_k = neurons[k];
		
		if((expr->sup_coeff[i]==0) && (expr->inf_coeff[i]==0)){
			res->inf_coeff[i] = 0.0;
			res->sup_coeff[i] = 0.0;
			continue;
		}
		expr_t * mul_expr = NULL;
		if(is_lower){
			if(expr->sup_coeff[i] < 0){
				mul_expr = neuron_k->uexpr;
			}
			else if(expr->inf_coeff[i]<0){
				mul_expr = neuron_k->lexpr;
			}
		}
		else{
			if(expr->sup_coeff[i] < 0){
				mul_expr = neuron_k->lexpr;
			}
			else if(expr->inf_coeff[i]<0){
				mul_expr = neuron_k->uexpr;
			}
		}
		if(expr->sup_coeff[i]<0 || expr->inf_coeff[i] < 0){
			double lambda_inf = mul_expr->inf_coeff[0];
			double lambda_sup = mul_expr->sup_coeff[0];
			double mu_inf = mul_expr->inf_cst;
			double mu_sup = mul_expr->sup_cst;
			//res->coeff[i] = lambda*expr->coeff[i];
			//res->cst = res->cst + expr->coeff[i]*mu;
			elina_double_interval_mul_expr_coeff(pr,&res->inf_coeff[i],&res->sup_coeff[i],lambda_inf,lambda_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
			double tmp1, tmp2;
			elina_double_interval_mul_cst_coeff(pr,&tmp1,&tmp2,mu_inf,mu_sup,expr->inf_coeff[i],expr->sup_coeff[i]);
			res->inf_cst = res->inf_cst + tmp1 + pr->min_denormal;
			res->sup_cst = res->sup_cst + tmp2 + pr->min_denormal;
			
		}
		else{
			res->inf_coeff[i] = 0.0;
			res->sup_coeff[i] = 0.0;
			double tmp1, tmp2;
			elina_double_interval_mul_expr_coeff(pr, &tmp1,&tmp2, neuron_k->lb, neuron_k->ub, expr->inf_coeff[i],expr->sup_coeff[i]);
			if(is_lower){
				res->inf_cst = res->inf_cst + tmp1;
				res->sup_cst = res->sup_cst - tmp1;
			}
			else{
				res->inf_cst = res->inf_cst - tmp2;
				res->sup_cst = res->sup_cst + tmp2;
			}
			
		}
	}
	if(expr->type==SPARSE){
		res->dim = (size_t*)malloc(num_neurons*sizeof(size_t));
		for(i=0; i < num_neurons; i++){
			res->dim[i] = expr->dim[i];
		}
	}
	
	return res;
}

expr_t * lexpr_replace_bounds_activation(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
	return expr_replace_bounds_activation(pr, expr, neurons, true);
}

expr_t * uexpr_replace_bounds_activation(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
	return expr_replace_bounds_activation(pr, expr, neurons, false);
}


expr_t * lexpr_replace_bounds(fppoly_internal_t * pr, expr_t * expr, neuron_t ** neurons, bool is_activation){
  //	if(is_activation){
  //		return lexpr_replace_bounds_activation(pr, expr, neurons);
  //	}
  //	else{
  return lexpr_replace_bounds_affine(pr, expr, neurons);
  //	}
}

expr_t * uexpr_replace_bounds(fppoly_internal_t * pr, expr_t * expr, neuron_t ** neurons, bool is_activation){
  //	if(is_activation){
  //		return uexpr_replace_bounds_activation(pr, expr, neurons);
  //	}
  //	else{
  return uexpr_replace_bounds_affine(pr, expr, neurons);
  //	}
}


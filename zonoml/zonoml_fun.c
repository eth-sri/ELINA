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

#include "zonoml_fun.h"


/* cretae zonotope from perturbed network inputs of size intdim+realdim */
elina_abstract0_t* zonotope_from_network_input(elina_manager_t* man, size_t intdim, size_t realdim, double* inf_array, double* sup_array)
{
    start_timing();
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_OF_BOX);
    zonotope_t* res = zonotope_alloc(man,intdim,realdim);
    size_t i = 0;
    for (i=0; i<intdim+realdim; i++) {
	//elina_interval_fprint(stdout,tinterval[i]);
        //printf("%g %g\n", inf_array[i],sup_array[i]);
	res->box_inf[i] = -inf_array[i];
	res->box_sup[i] = sup_array[i];
	res->paf[i] = zonotope_aff_alloc_init(pr);
	res->paf[i]->itv_inf = INFINITY;
	res->paf[i]->itv_sup = INFINITY;
	if (inf_array[i] > sup_array[i]){
		
		 res->paf[i] = pr->bot;
	}
	else if ((inf_array[i]==-INFINITY) && (sup_array[i]==INFINITY)){ 
		res->paf[i] = pr->top;
	}
	else if (isfinite(inf_array[i]) && (inf_array[i] == sup_array[i])) {
		res->paf[i]->c_inf = -inf_array[i];
		res->paf[i]->c_sup = sup_array[i];
	}
	else if (isfinite(inf_array[i]) && isfinite(sup_array[i])){
		 zonotope_aff_add_itv(pr, res->paf[i], -inf_array[i], sup_array[i], IN);
	}
	else{
	 	res->paf[i]->c_inf = -inf_array[i];
		res->paf[i]->c_sup = sup_array[i];
	}
        //elina_interval_set(res->paf[i]->itv,tinterval[i]);
	res->paf[i]->pby++;
    }

    man->result.flag_best = true;
    man->result.flag_exact = true;
    record_timing(zonoml_network_input_time);
    return abstract0_of_zonotope(man,res);
}



zonotope_aff_t* zonotope_aff_mul_weight(zonotope_internal_t* pr, zonotope_aff_t* src, double lambda)
{
   
    if ((lambda==0)|| zonotope_aff_is_known_to_be_zero(pr, src)) {
        return zonotope_aff_alloc_init(pr);
    } else if (zonotope_aff_is_bottom(pr, src)) {
	
        return zonotope_aff_bottom(pr);
    } else if (zonotope_aff_is_top(pr, src) || !isfinite(lambda)) {
        return zonotope_aff_top(pr);
    } else {
         
        zonotope_aff_t* dst = NULL;
        zonotope_aaterm_t *p,*q;
        q = NULL;
        dst = zonotope_aff_alloc_init(pr);
        
        elina_double_interval_mul(&dst->c_inf,&dst->c_sup, -lambda, lambda, src->c_inf,src->c_sup);
        double maxA =  fmax(fabs(src->c_inf),fabs(src->c_sup));
        double tmp1, tmp2;
        elina_double_interval_mul(&tmp1,&tmp2, -lambda, lambda, maxA*pr->ulp, maxA*pr->ulp);
        double fp_err_inf = tmp1+ pr->min_denormal;
        double fp_err_sup = tmp2 + pr->min_denormal;
        dst->c_inf+= tmp1 + pr->min_denormal;
        dst->c_sup+= tmp2 + pr->min_denormal;
        if (src->q) {
            
            dst->q = q = zonotope_aaterm_alloc_init();
            for (p=src->q; p; p=p->n) {
		
                double tmp_inf = 0.0;
		double tmp_sup = 0.0;
                elina_double_interval_mul(&tmp_inf, &tmp_sup, -lambda, lambda, p->inf, p->sup);
                maxA =  fmax(fabs(p->inf),fabs(p->sup));
                elina_double_interval_mul(&tmp1,&tmp2, -lambda, lambda, maxA*pr->ulp, maxA*pr->ulp);
                fp_err_inf += tmp1;
                fp_err_sup += tmp2;
		q->inf = tmp_inf+tmp1;
		q->sup = tmp_sup+tmp2;
                q->pnsym = p->pnsym;
                if (p->n) {
                    /* continue */
                    q->n = zonotope_aaterm_alloc_init();
                    q = q->n;
                } else {
                    /* the last iteration */
                    dst->end = q; 
                }
            }
            
        }
        
        dst->l = src->l;
	
        elina_double_interval_mul(&dst->itv_inf, &dst->itv_sup,  -lambda, lambda,src->itv_inf, src->itv_sup);
        dst->itv_inf += fp_err_inf;
        dst->itv_sup += fp_err_sup;
        return dst;
    }
}

zonotope_aff_t * zonotope_aff_from_dense_weights(zonotope_internal_t * pr, double * weights, size_t offset, size_t size, zonotope_t *z){
    zonotope_aff_t *res = zonotope_aff_alloc_init(pr);
    size_t i;	
    
    for(i=0; i < size; i++){
	zonotope_aff_t *tmp;
	zonotope_aff_t *aff = z->paf[offset+i];
	//printf("%.30f\n", weights[i]);
	//zonotope_aff_fprint(pr,stdout,aff);
        //printf("%.30f\n",aff->itv_sup);
	//printf("\n");
	tmp = zonotope_aff_mul_weight(pr,aff,weights[i]);
	//printf("after mul : %.15f %.15f\n", -tmp->itv_inf ,tmp->itv_sup);
	//zonotope_aff_fprint(pr,stdout,tmp);
	//printf("\n");
	//fflush(stdout);
	zonotope_aff_t *tmp1 = res;
	res = zonotope_aff_add(pr,tmp1,tmp,z);
	//printf("after add: %.15f %.15f\n", -res->itv_inf,res->itv_sup);
	//zonotope_aff_fprint(pr,stdout,tmp);
	//printf("\n");
	//fflush(stdout);	
        zonotope_aff_free(pr,tmp);
        zonotope_aff_free(pr,tmp1);
		
    }
    //printf("\n");
    return res;
}

zonotope_aff_t * zonotope_aff_from_dense_weights_bias(zonotope_internal_t* pr, double * weights, double bias, size_t offset, size_t size, zonotope_t *z){
    zonotope_aff_t * res = zonotope_aff_from_dense_weights(pr, weights, offset, size, z);
    res->c_inf += -bias;
    res->c_sup += bias;
    res->itv_inf += -bias;
    res->itv_sup += bias;	
    return res;
}


zonotope_aff_t * zonotope_aff_from_sparse_weights_bias(zonotope_internal_t* pr, double * weights, double bias, elina_dim_t *dim, size_t size, zonotope_t *z){
    zonotope_aff_t *res = zonotope_aff_alloc_init(pr);
    res->c_inf = -bias;
    res->c_sup = bias;
    res->itv_inf = -bias;
    res->itv_sup = bias;
    size_t i;
	
    for(i=0; i < size; i++){
	
	zonotope_aff_t *tmp;
        elina_dim_t var = dim[i];
	
	zonotope_aff_t *aff = z->paf[var];
	
	tmp = zonotope_aff_mul_weight(pr,aff,weights[i]);
	zonotope_aff_t *tmp1 = res;
	res = zonotope_aff_add(pr,tmp1,tmp,z);	
        zonotope_aff_free(pr,tmp);
        zonotope_aff_free(pr,tmp1);
		
    }
   
    return res;
}


void * handle_ffn_matmult_zono_parallel(void *args){
	zonoml_ffn_matmult_thread_t * data = (zonoml_ffn_matmult_thread_t *)args;
	zonotope_internal_t * pr = data->pr;
	zonotope_t * z = data->z;
	size_t start_offset = data->start_offset;
	size_t idx_start = data->start;
	size_t idx_end = data->end;
	double **weights = data->weights;
	double *bias = data->bias;
	bool has_bias = data->has_bias;
	size_t expr_offset = data->expr_offset;
	size_t expr_size = data->expr_size;

	size_t offset = start_offset + idx_start;
	size_t i;
	for (i=idx_start; i< idx_end; i++) {
		zonotope_aff_check_free(pr, z->paf[offset]);
	
		z->paf[offset] = has_bias ? zonotope_aff_from_dense_weights_bias(pr, weights[i], bias[i], expr_offset, expr_size, z) : 
					    zonotope_aff_from_dense_weights(pr, weights[i], expr_offset, expr_size, z);
		if (zonotope_aff_is_top(pr, z->paf[offset])) {
	    	     zonotope_aff_check_free(pr, z->paf[offset]);
	    	     z->paf[offset] = pr->top;
		} 
		else if (zonotope_aff_is_bottom(pr, z->paf[offset])) {
	    	     zonotope_aff_check_free(pr, z->paf[offset]);
	    	     z->paf[offset] = pr->bot;
		}
		z->box_inf[offset] = z->paf[offset]->itv_inf;
		z->box_sup[offset] = z->paf[offset]->itv_sup;
		z->paf[offset]->pby++;
		offset++;
    	}
	return NULL; 
}

// assumes that the variables for the FFN matmult have already been added, the first assignment is from start_offset
elina_abstract0_t* ffn_matmult_zono(elina_manager_t * man, bool destructive, elina_abstract0_t* element, elina_dim_t start_offset,
			       double **weights, double * bias,  size_t num_var, size_t expr_offset, size_t expr_size){
   start_timing();
   zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
   zonotope_t *z = zonotope_of_abstract0(element);
   zonotope_t* res = zonotope_copy(man, z);
   //printf("input\n");
   //zonotope_fprint(stdout,man,res,NULL);
   size_t i = 0;
   for (i=0; i<res->dims; i++) {
	res->paf[i]->itv_inf = res->box_inf[i];
	res->paf[i]->itv_sup = res->box_sup[i];
    }
    ffn_matmult_zono_parallel(pr, res, start_offset, weights, bias,  num_var, expr_offset, expr_size, handle_ffn_matmult_zono_parallel, true);
    man->result.flag_best = false;
    man->result.flag_exact = false;
    record_timing(zonoml_ffn_matmult_time);
    //printf("Output\n");
    //zonotope_fprint(stdout,man,res,NULL);
    return abstract0_of_zonotope(man,res);
}


elina_abstract0_t* ffn_matmult_without_bias_zono(elina_manager_t * man, bool destructive, elina_abstract0_t* element, elina_dim_t start_offset,
			       double **weights, size_t num_var, size_t expr_offset, size_t expr_size){
   start_timing();
   zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
   zonotope_t *z = zonotope_of_abstract0(element);
   zonotope_t* res = zonotope_copy(man, z);
   size_t i = 0;
   for (i=0; i<res->dims; i++) {
	res->paf[i]->itv_inf = res->box_inf[i];
	res->paf[i]->itv_sup = res->box_sup[i];
    }
    ffn_matmult_zono_parallel(pr, res, start_offset, weights, NULL,  num_var, expr_offset, expr_size, handle_ffn_matmult_zono_parallel, false);
    man->result.flag_best = false;
    man->result.flag_exact = false;
    record_timing(zonoml_ffn_matmult_time);
    return abstract0_of_zonotope(man,res);
}

typedef enum bias_op{
	ADD,
        SUB1,
	SUB2,
        MUL,
}bias_op;

void create_linexpr_array_with_bias(elina_linexpr0_t ** expr_array, elina_dim_t * tdim_array, size_t offset, double *bias, size_t num_var, bias_op OP){
	size_t i;
	for (i=0; i<num_var; i++) {
		double cst, coeff;
		if(OP==MUL){
			cst = 0;
			coeff = bias[i];
		}
		else{
			if(OP==SUB1){
				coeff = 1;
				cst = -bias[i];
                        } else if (OP == SUB2) {
                          coeff = -1;
                          cst = bias[i];
                        } else {
                          coeff = 1;
                          cst = bias[i];
                        }
                }
		expr_array[i] = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
		elina_coeff_set_scalar_double(&expr_array[i]->cst,cst);
		elina_linterm_t *lterm = &expr_array[i]->p.linterm[0];
		lterm->dim = offset;
		elina_coeff_set_interval_double(&lterm->coeff,coeff,coeff); 
		tdim_array[i] = offset;
		offset++;
   	}
}

elina_abstract0_t * ffn_unary_bias_zono(elina_manager_t * man, bool destructive, elina_abstract0_t* element, elina_dim_t start_offset,
			        double * bias, size_t num_var, bias_op OP){
   start_timing();
   zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
   zonotope_t *z = zonotope_of_abstract0(element);
   zonotope_t* res = zonotope_copy(man, z);
   size_t i = 0;
   elina_linexpr0_t ** expr_array = (elina_linexpr0_t **)malloc(num_var*sizeof(elina_linexpr0_t *));
   elina_dim_t * tdim_array = (elina_dim_t *)malloc(num_var*sizeof(elina_dim_t));
   create_linexpr_array_with_bias(expr_array, tdim_array, start_offset, bias, num_var, OP);
   res = zonotope_assign_linexpr_array(man,true,res,tdim_array,expr_array,num_var,NULL);
   for(i=0; i < num_var; i++){
	elina_linexpr0_free(expr_array[i]);
   } 
    free(expr_array);
    free(tdim_array);
    man->result.flag_best = false;
    man->result.flag_exact = false;
    record_timing(zonoml_ffn_matmult_time);
    return abstract0_of_zonotope(man,res);
}

elina_abstract0_t* ffn_add_bias_zono(elina_manager_t * man, bool destructive, elina_abstract0_t* element, elina_dim_t start_offset,
			        double * bias, size_t num_var){
   return ffn_unary_bias_zono(man, destructive, element, start_offset, bias, num_var, ADD);
}

elina_abstract0_t* ffn_sub_bias_zono(elina_manager_t * man, bool destructive, elina_abstract0_t* element, elina_dim_t start_offset,
			        double * bias, bool is_minuend, size_t num_var){
   if(is_minuend==true){
	return ffn_unary_bias_zono(man, destructive, element, start_offset, bias, num_var, SUB1);
   }
   else{
	return ffn_unary_bias_zono(man, destructive, element, start_offset, bias, num_var, SUB2);
   }
}



elina_abstract0_t* ffn_mul_bias_zono(elina_manager_t * man, bool destructive, elina_abstract0_t* element, elina_dim_t start_offset,
			        double * bias, size_t num_var){
   return ffn_unary_bias_zono(man, destructive, element, start_offset, bias, num_var, MUL);
}





void * handle_conv_matmult_zono_parallel(void *args){
	zonoml_conv_matmult_thread_t * data = (zonoml_conv_matmult_thread_t *)args;
	zonotope_internal_t * pr = data->pr;
	zonotope_t * z = data->z;
	size_t start_offset = data->start_offset;
	size_t idx_start = data->start;
	size_t idx_end = data->end;
	double *filter_weights = data->filter_weights;
	double *filter_bias = data->filter_bias;

	size_t expr_offset = data->expr_offset;
	size_t *input_size = data->input_size;
	size_t *filter_size = data->filter_size;
	size_t num_filters = data->num_filters;
	size_t *strides = data->strides;
	size_t *output_size = data->output_size;		
	long int pad_top = data->pad_top;
	long int pad_left = data->pad_left;  
	bool has_bias = data->has_bias;
	
	size_t out_x, out_y, out_z, inp_z;
	long int x_shift, y_shift;
     
        size_t num_pixels = input_size[0]*input_size[1]*input_size[2];
	
	size_t mat_x;
	size_t o12 = output_size[1]*output_size[2];
	for (mat_x=idx_start; mat_x< idx_end; mat_x++) {
	     size_t out_x = mat_x/o12;
	     size_t out_y = (mat_x - out_x*o12) / output_size[2];
	     size_t out_z =  mat_x -out_x*o12 - out_y*output_size[2];
	     size_t num_coeff = input_size[2]*filter_size[0]*filter_size[1];
	     size_t actual_coeff = 0;
	     double *coeff = (double *)malloc(num_coeff*sizeof(double));
	     elina_dim_t *dim = (elina_dim_t *)malloc(num_coeff*sizeof(elina_dim_t));
	     size_t i=0;
	     for(inp_z=0; inp_z <input_size[2]; inp_z++) {
		  for(x_shift = 0; x_shift < (long int)filter_size[0]; x_shift++) {
		      for(y_shift =0; y_shift < (long int)filter_size[1]; y_shift++) {
			  long int x_val = out_x*strides[0]+x_shift-pad_top;	
			  long int y_val = out_y*strides[1]+y_shift-pad_left;
			  if(y_val<0 || y_val >= (long int)input_size[1]){
			     continue;
			  }
				     
			  if(x_val<0 || x_val >= (long int)input_size[0]){
			     continue;
			  }
				     
			  size_t mat_offset = x_val*input_size[1]*input_size[2] + y_val*input_size[2] + inp_z;
				     
			  size_t mat_y = expr_offset + mat_offset;
		          size_t filter_index = x_shift*filter_size[1]*input_size[2]*output_size[2] + y_shift*input_size[2]*output_size[2] + inp_z*output_size[2] + out_z;	
			  if(mat_offset>=num_pixels){		 
			     continue;
		          }    
			  coeff[i] = filter_weights[filter_index];
			  dim[i] = mat_y;
			  i++;
			  actual_coeff++;
		      }
		   }
	     }
			
	     double cst = has_bias ? filter_bias[out_z] : 0.0;

             z->paf[start_offset+mat_x] = zonotope_aff_from_sparse_weights_bias(pr, coeff, cst, dim, actual_coeff, z);
			
             if (zonotope_aff_is_top(pr, z->paf[start_offset+mat_x])) {
	    	 zonotope_aff_check_free(pr, z->paf[start_offset+mat_x]);
	    	 z->paf[start_offset+mat_x] = pr->top;
	     } 
             else if (zonotope_aff_is_bottom(pr, z->paf[start_offset+mat_x])) {
	    	 zonotope_aff_check_free(pr, z->paf[start_offset+mat_x]);
	    	 z->paf[start_offset+mat_x] = pr->bot;
	     }
	     z->box_inf[start_offset+mat_x] = z->paf[start_offset+mat_x]->itv_inf;
	     z->box_sup[start_offset+mat_x] = z->paf[start_offset+mat_x]->itv_sup;
	     z->paf[start_offset+mat_x]->pby++;
	     free(coeff);
	     free(dim);
		
    	}
	return NULL; 
}


// assumes that the variables for the convolutional matmult have already been added, two dimensional filter, the first assignment is from start_offset
elina_abstract0_t* conv_matmult_zono(elina_manager_t* man, bool destructive, elina_abstract0_t* element, elina_dim_t start_offset, double *filter_weights, double * filter_bias,  
				      size_t * input_size, size_t expr_offset, size_t *filter_size, size_t num_filters, size_t *strides, size_t *output_size, size_t pad_top, size_t pad_left, bool has_bias){
	start_timing();
	zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
	
	zonotope_t *z = zonotope_of_abstract0(element);
   	zonotope_t* res = zonotope_copy(man, z);
	elina_dimension_t dims = zonotope_dimension(man,res);
	size_t i, j;
	size_t num_pixels = input_size[0]*input_size[1]*input_size[2];
	
	size_t num_out_neurons = output_size[0]*output_size[1]*output_size[2];
	

	conv_matmult_zono_parallel(pr, res, start_offset, filter_weights, filter_bias, num_out_neurons,
				   expr_offset, input_size, filter_size, num_filters, strides, output_size, 
				   pad_top, pad_left, handle_conv_matmult_zono_parallel, has_bias);
	
    	man->result.flag_best = false;
    	man->result.flag_exact = false;
    	record_timing(zonoml_conv_matmult_time);
    	return abstract0_of_zonotope(man,res);
}

elina_abstract0_t *handle_gather_layer(elina_manager_t* man, bool destructive, elina_abstract0_t * abs, size_t *indexes){
        elina_dimension_t dimension = elina_abstract0_dimension(man,abs);
        size_t size = dimension.intdim + dimension.realdim;
	elina_dimperm_t * dimperm = elina_dimperm_alloc(size);
	size_t i;
	for(i=0; i < size; i++){
               
		dimperm->dim[indexes[i]] = i;
	}
	elina_abstract0_t *res = elina_abstract0_permute_dimensions(man,destructive, abs,dimperm);
	elina_dimperm_free(dimperm);
	return res;
}

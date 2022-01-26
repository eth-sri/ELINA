/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
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


#ifndef _ZONOML_INTERNAL_H_
#define _ZONOML_INTERNAL_H_

#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>
#include  <fenv.h>
#include <errno.h>
#include "zonotope.h"
#include "zonoml.h"
#include "rdtsc.h"

#include <pthread.h>
//#include <sys/sysinfo.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(TIMING)
#define start_timing()				\
tsc_counter start, end;				\
double cycles;					\
CPUID();					\
RDTSC(start)
    
#define record_timing(counter)		\
RDTSC(end);				\
CPUID();				\
cycles = (double)(COUNTER_DIFF(end, start));	\
counter += cycles
    
    extern double zonoml_relu_time;
    extern double zonoml_sigmoid_time;
    extern double zonoml_tanh_time;
    extern double zonoml_pool_time;
    extern double zonoml_network_input_time;
    extern double zonoml_conv_matmult_time;
    extern double zonoml_ffn_matmult_time;

#else
    #define start_timing()
    #define record_timing(X)
#endif


 


typedef struct zonoml_ffn_matmult_thread_t{
	size_t start;
	size_t end;
	zonotope_internal_t* pr;
	zonotope_t *z;
	elina_dim_t start_offset;
	double ** weights;
	double * bias;
	size_t expr_offset;
	size_t expr_size;
	bool has_bias;	
}zonoml_ffn_matmult_thread_t;


typedef struct zonoml_relu_thread_t{
	size_t start;
	size_t end;
	size_t start_offset;
    	zonotope_t *z;
	elina_manager_t *man;
	bool create_new_noise_symbol;
	zonotope_noise_symbol_t ** epsilon_map;
	zonotope_noise_symbol_t ** epsilon_map_extra;
}zonoml_relu_thread_t;

typedef struct zonoml_s_curve_thread_t{
	size_t start;
	size_t end;
	size_t start_offset;
    	zonotope_t *z;
	elina_manager_t *man;
	zonotope_noise_symbol_t ** epsilon_map;
	bool is_sigmoid;
	char * is_used;
}zonoml_s_curve_thread_t;

typedef struct zonoml_conv_matmult_thread_t{
	size_t start;
	size_t end;
	zonotope_internal_t* pr;
	zonotope_t *z;
	elina_dim_t start_offset;
	double * filter_weights;
	double * filter_bias;
	size_t * input_size;
	size_t expr_offset;
	size_t *filter_size;
	size_t num_filters;
	size_t *strides;
	long int pad_top;
	long int pad_left;
	size_t *output_size;
	bool has_bias;	
}zonoml_conv_matmult_thread_t;


typedef struct zonoml_pool_thread_t{
	size_t start;
	size_t end;
	zonotope_internal_t* pr;
	zonotope_t *z;
	size_t src_offset;
	size_t * pool_size;
	size_t * input_size;
	size_t dst_offset;
	size_t *strides;
	long int pad_top;
	long int pad_left;
	size_t *output_size;
	zonotope_noise_symbol_t ** epsilon_map;
	char * is_used;
	bool is_maxpool;
}zonoml_pool_thread_t;


static inline zonotope_t* zonotope_of_abstract0(elina_abstract0_t* a)
{
  return (zonotope_t*)a->value;
}

static inline elina_abstract0_t* abstract0_of_zonotope(elina_manager_t* man, zonotope_t* z)
{
  elina_abstract0_t* r = malloc(sizeof(elina_abstract0_t));
  assert(r);
  r->value = z;
  r->man = elina_manager_copy(man);
  return r;
}

static inline zonotope_noise_symbol_t* zonotope_noise_symbol_add_with_index(zonotope_internal_t *pr, noise_symbol_t type, size_t index)
    /* increment the global index of used noise symbols and add the noise symbol in pr->eps */
{
    uint_t dim = pr->dim;
    zonotope_noise_symbol_t* res;
	
    /* resize epsilon array */
    if (dim % 1024 == 0){ pr->epsilon = (zonotope_noise_symbol_t**)realloc(pr->epsilon, (dim+1024)*sizeof(zonotope_noise_symbol_t*));
}
    res = pr->epsilon[dim] = (zonotope_noise_symbol_t*)malloc(sizeof(zonotope_noise_symbol_t));
    if (type == IN) {
	if (pr->epssize % 1024 == 0) pr->inputns = (uint_t*)realloc(pr->inputns, (pr->epssize+1024)*sizeof(uint_t));
	pr->inputns[pr->epssize] = dim; pr->epssize++;
    }
    res->index = dim;
    res->type = type;
    pr->dim++;
    return res;
}

static inline void ffn_matmult_zono_parallel(zonotope_internal_t* pr, zonotope_t *z, elina_dim_t start_offset,
			       			    double **weights, double * bias,  size_t num_out_neurons,
						    size_t expr_offset, size_t expr_size, void *(*function)(void *), bool has_bias){
	//int num_threads = get_nprocs();
	int num_threads = sysconf(_SC_NPROCESSORS_ONLN);
	if(num_threads <1){
		num_threads = 1;
	}
	
	zonoml_ffn_matmult_thread_t args[num_threads];
	pthread_t threads[num_threads];
	int i;
	//printf("start %zu\n",num_out_neurons);
	//fflush(stdout);
	if((int)num_out_neurons < num_threads){
		for (i = 0; i < (int)num_out_neurons; i++){
	    		args[i].start = i; 
	    		args[i].end = i+1;
			args[i].pr = pr;
			args[i].z = z;
			args[i].start_offset = start_offset;
			args[i].weights = weights;
			args[i].bias = bias;
			args[i].expr_offset = expr_offset;
			args[i].expr_size = expr_size;   
			args[i].has_bias = has_bias;			
	    		pthread_create(&threads[i], NULL,function, (void*)&args[i]);
			
	  	}
		for (i = 0; i < (int)num_out_neurons; i = i + 1){
			pthread_join(threads[i], NULL);
		}
	}
	else{
		size_t idx_start = 0;
		size_t idx_n = num_out_neurons / num_threads;
		size_t idx_end = idx_start + idx_n;
		
		
	  	for (i = 0; i < num_threads; i++){
	    		args[i].start = idx_start; 
	    		args[i].end = idx_end;
			args[i].pr = pr;
			args[i].z = z;
			args[i].start_offset = start_offset;
			args[i].weights = weights;
			args[i].bias = bias;
			args[i].expr_offset = expr_offset;
			args[i].expr_size = expr_size;   
			args[i].has_bias = has_bias;
	    		pthread_create(&threads[i], NULL,function, (void*)&args[i]);
			idx_start = idx_end;
			idx_end = idx_start + idx_n;
	    		if(idx_end>num_out_neurons){
				idx_end = num_out_neurons;
			}
			if((i==num_threads-2)){
				idx_end = num_out_neurons;
				
			}
	  	}

		for (i = 0; i < num_threads; i = i + 1){
			pthread_join(threads[i], NULL);
		}
	}
	
}

static inline void relu_zono_parallel(elina_manager_t* man, zonotope_t *z, elina_dim_t start_offset, elina_dim_t num_out_neurons, bool create_new_noise_symbol, void *(*function)(void *)){
	//int num_threads = get_nprocs();
	int num_threads = sysconf(_SC_NPROCESSORS_ONLN);
	if(num_threads <1){
		num_threads = 1;
	}
	int i;
	
	zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
	size_t offset = start_offset;
	zonotope_noise_symbol_t ** epsilon_map = NULL;
	zonotope_noise_symbol_t ** epsilon_map_extra = NULL;

	if(create_new_noise_symbol){
		epsilon_map = (zonotope_noise_symbol_t **)malloc(num_out_neurons*sizeof(zonotope_noise_symbol_t*));
		epsilon_map_extra = (zonotope_noise_symbol_t **)malloc(num_out_neurons*sizeof(zonotope_noise_symbol_t*));
	}

	for(i=0; i < (int)num_out_neurons; i++){
		elina_interval_t *itv = zonotope_bound_dimension(man,z,offset);
		double inf = itv->inf->val.dbl;
		double sup = itv->sup->val.dbl;
		double inf_l = -inf;
		double inf_u = inf;
		double sup_l = -sup;
		double sup_u = sup;
		
			
		double width_l = sup_l + inf_u;
		double width_u = sup_u + inf_l;
		
		double bound_l, bound_u;
		elina_double_interval_mul(&bound_l, &bound_u, inf_l, inf_u, sup_l, sup_u);
		double tmp = bound_l;
		bound_l = bound_u;
		bound_u = tmp;
		elina_double_interval_div(&bound_l, &bound_u, bound_l, bound_u, width_l, width_u);
		if(create_new_noise_symbol){
			if((inf<0) && (sup>0)){
				epsilon_map[i] = zonotope_noise_symbol_add(pr, IN);
				if((-inf>sup) && (-bound_l>1)){
					epsilon_map_extra[i] = zonotope_noise_symbol_add(pr, IN);
				}	
			}
			else{
				epsilon_map[i] = NULL;
			}
		}
		elina_interval_free(itv);
		offset++;
	}
	
	zonoml_relu_thread_t args[num_threads];
	pthread_t threads[num_threads];
	
	//printf("start %zu\n",num_out_neurons);
	//fflush(stdout);
	if((int)num_out_neurons < num_threads){
		for (i = 0; i < (int)num_out_neurons; i++){
	    		args[i].start = i; 
	    		args[i].end = i+1;
			args[i].man = man;
			args[i].z = z;
            		args[i].start_offset = start_offset;
			args[i].epsilon_map = epsilon_map;	
			args[i].epsilon_map_extra = epsilon_map_extra;
			args[i].create_new_noise_symbol = create_new_noise_symbol;
	    		pthread_create(&threads[i], NULL,function, (void*)&args[i]);
			
	  	}
		for (i = 0; i < (int)num_out_neurons; i = i + 1){
			pthread_join(threads[i], NULL);
		}
	}
	else{
		size_t idx_start = 0;
		size_t idx_n = num_out_neurons / num_threads;
		size_t idx_end = idx_start + idx_n;
		
		
	  	for (i = 0; i < num_threads; i++){
	    		args[i].start = idx_start; 
	    		args[i].end = idx_end;
			args[i].man = man;
			args[i].z = z;
			args[i].start_offset = start_offset;
			args[i].epsilon_map = epsilon_map;
			args[i].epsilon_map_extra = epsilon_map_extra;
			args[i].create_new_noise_symbol = create_new_noise_symbol;
	    		pthread_create(&threads[i], NULL,function, (void*)&args[i]);
			idx_start = idx_end;
			idx_end = idx_start + idx_n;
	    		if(idx_end>num_out_neurons){
				idx_end = num_out_neurons;
			}
			if((i==num_threads-2)){
				idx_end = num_out_neurons;
				
			}
	  	}

		for (i = 0; i < num_threads; i = i + 1){
			pthread_join(threads[i], NULL);
		}
	}
	if(create_new_noise_symbol){
		free(epsilon_map);
		free(epsilon_map_extra);
	}
	
}
	

static inline void s_curve_zono_parallel(elina_manager_t* man, zonotope_t *z, elina_dim_t start_offset, elina_dim_t num_out_neurons, void *(*function)(void *), bool is_sigmoid){
	//int num_threads = get_nprocs();
	int num_threads = sysconf(_SC_NPROCESSORS_ONLN);
	if(num_threads <1){
		num_threads = 1;
	}
	int i;
	
	zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
	
	zonotope_noise_symbol_t ** epsilon_map = (zonotope_noise_symbol_t **)malloc(num_out_neurons*sizeof(zonotope_noise_symbol_t*));
	char *is_used = (char*)calloc(num_out_neurons,sizeof(char));
	for(i=0; i < (int)num_out_neurons; i++){
	    epsilon_map[i] = zonotope_noise_symbol_add(pr, IN);
	}
	
	zonoml_s_curve_thread_t args[num_threads];
	pthread_t threads[num_threads];
	
	
	if((int)num_out_neurons < num_threads){
		for (i = 0; i < (int)num_out_neurons; i++){
	    		args[i].start = i; 
	    		args[i].end = i+1;
			args[i].man = man;
			args[i].z = z;
            		args[i].start_offset = start_offset;
			args[i].epsilon_map = epsilon_map;	
			args[i].is_sigmoid = is_sigmoid;
			args[i].is_used = is_used;
	    		pthread_create(&threads[i], NULL,function, (void*)&args[i]);
			
	  	}
		for (i = 0; i < (int)num_out_neurons; i = i + 1){
			pthread_join(threads[i], NULL);
		}
	}
	else{
		size_t idx_start = 0;
		size_t idx_n = num_out_neurons / num_threads;
		size_t idx_end = idx_start + idx_n;
		
		
	  	for (i = 0; i < num_threads; i++){
	    		args[i].start = idx_start; 
	    		args[i].end = idx_end;
			args[i].man = man;
			args[i].z = z;
			args[i].start_offset = start_offset;
			args[i].epsilon_map = epsilon_map;
			args[i].is_sigmoid = is_sigmoid;
			args[i].is_used = is_used;
	    		pthread_create(&threads[i], NULL,function, (void*)&args[i]);
			idx_start = idx_end;
			idx_end = idx_start + idx_n;
	    		if(idx_end>num_out_neurons){
				idx_end = num_out_neurons;
			}
			if((i==num_threads-2)){
				idx_end = num_out_neurons;
				
			}
	  	}

		for (i = 0; i < num_threads; i = i + 1){
			pthread_join(threads[i], NULL);
		}
	}
	for(i=0; i < (int)num_out_neurons; i++){
		if(!is_used[i]){
			zonotope_delete_constrained_noise_symbol(pr, epsilon_map[i]->index, z);
		}
	}
	free(epsilon_map);
	free(is_used);
	
}



static inline void conv_matmult_zono_parallel(zonotope_internal_t* pr, zonotope_t *z, elina_dim_t start_offset,
			       			    double * filter_weights, double * filter_bias,  size_t num_out_neurons,
						    size_t expr_offset, size_t *input_size, size_t *filter_size, size_t num_filters,
						    size_t *strides, size_t *output_size, long int pad_top,
						    long int pad_left, void *(*function)(void *), bool has_bias){
	//int num_threads = get_nprocs();
	int num_threads = sysconf(_SC_NPROCESSORS_ONLN);
	if(num_threads <1){
		num_threads = 1;
	}
	
	zonoml_conv_matmult_thread_t args[num_threads];
	pthread_t threads[num_threads];
	int i;
	//printf("start %zu\n",num_out_neurons);
	//fflush(stdout);
	if((int)num_out_neurons < num_threads){
		for (i = 0; i < (int)num_out_neurons; i++){
	    		args[i].start = i; 
	    		args[i].end = i+1;
			args[i].pr = pr;
			args[i].z = z;
			args[i].start_offset = start_offset;
			args[i].filter_weights = filter_weights;
			args[i].filter_bias = filter_bias;
			args[i].expr_offset = expr_offset;
			args[i].input_size = input_size;
			args[i].filter_size = filter_size;
			args[i].num_filters = num_filters;
			args[i].strides = strides;
			args[i].output_size = output_size;
			args[i].pad_top = pad_top;
			args[i].pad_left = pad_left;   
			args[i].has_bias = has_bias;			
	    		pthread_create(&threads[i], NULL,function, (void*)&args[i]);
			
	  	}
		for (i = 0; i < (int)num_out_neurons; i = i + 1){
			pthread_join(threads[i], NULL);
		}
	}
	else{
		size_t idx_start = 0;
		size_t idx_n = num_out_neurons / num_threads;
		size_t idx_end = idx_start + idx_n;
		
		
	  	for (i = 0; i < num_threads; i++){
	    		args[i].start = idx_start; 
	    		args[i].end = idx_end;
			args[i].pr = pr;
			args[i].z = z;
			args[i].start_offset = start_offset;
			args[i].filter_weights = filter_weights;
			args[i].filter_bias = filter_bias;
			args[i].expr_offset = expr_offset;
			args[i].input_size = input_size;
			args[i].filter_size = filter_size;
			args[i].num_filters = num_filters;
			args[i].strides = strides;
			args[i].output_size = output_size;
			args[i].pad_top = pad_top;
			args[i].pad_left = pad_left;   
			args[i].has_bias = has_bias;
	    		pthread_create(&threads[i], NULL,function, (void*)&args[i]);
			idx_start = idx_end;
			idx_end = idx_start + idx_n;
	    		if(idx_end>num_out_neurons){
				idx_end = num_out_neurons;
			}
			if((i==num_threads-2)){
				idx_end = num_out_neurons;
				
			}
	  	}

		for (i = 0; i < num_threads; i = i + 1){
			pthread_join(threads[i], NULL);
		}
	}
	
}


static inline void pool_zono_parallel(zonotope_internal_t* pr, zonotope_t *z, size_t src_offset, size_t *pool_size,
			       	           size_t num_out_neurons, size_t dst_offset, size_t *input_size, size_t *strides,
					  size_t *output_size, long int pad_top, long int pad_left, bool is_maxpool, void *(*function)(void *)){
	//int num_threads = get_nprocs();
	int num_threads = sysconf(_SC_NPROCESSORS_ONLN);
	if(num_threads <1){
		num_threads = 1;
	}
	
	zonoml_pool_thread_t args[num_threads];
	pthread_t threads[num_threads];
	int i;
	zonotope_noise_symbol_t ** epsilon_map = (zonotope_noise_symbol_t **)malloc(num_out_neurons*sizeof(zonotope_noise_symbol_t*));
	elina_dimchange_t dimadd;
    	elina_dimchange_init(&dimadd, 0, num_out_neurons);
	//elina_dimension_t dims = zonotope_dimension(pr->man,z);
        //elina_dim_t num_var = dims.intdim + dims.realdim;
	char * is_used = (char*)calloc(num_out_neurons,sizeof(char));
	for(i=0; i < (int)num_out_neurons; i++){
		epsilon_map[i] = zonotope_noise_symbol_add(pr, IN);
		dimadd.dim[i] = dst_offset;		
		
	}
	z = zonotope_add_dimensions(pr->man,true,z,&dimadd,false);
	
	//printf("start %zu\n",num_out_neurons);
	//fflush(stdout);
	if((int)num_out_neurons < num_threads){
		for (i = 0; i < (int)num_out_neurons; i++){
	    		args[i].start = i; 
	    		args[i].end = i+1;
			args[i].pr = pr;
			args[i].z = z;
			args[i].src_offset = src_offset;
			args[i].pool_size = pool_size;
			args[i].dst_offset = dst_offset;
			args[i].input_size = input_size;
			args[i].strides = strides;
			args[i].output_size = output_size;
			args[i].pad_top = pad_top;
			args[i].pad_left = pad_left;
			args[i].epsilon_map = epsilon_map;
			args[i].is_used = is_used;   
			args[i].is_maxpool = is_maxpool;
	    		pthread_create(&threads[i], NULL,function, (void*)&args[i]);
			
	  	}
		for (i = 0; i < (int)num_out_neurons; i = i + 1){
			pthread_join(threads[i], NULL);
		}
	}
	else{
		size_t idx_start = 0;
		size_t idx_n = num_out_neurons / num_threads;
		size_t idx_end = idx_start + idx_n;
		
		
	  	for (i = 0; i < num_threads; i++){
	    		args[i].start = idx_start; 
	    		args[i].end = idx_end;
			args[i].pr = pr;
			args[i].z = z;
			args[i].src_offset = src_offset;
			args[i].pool_size = pool_size;
			args[i].dst_offset = dst_offset;
			args[i].input_size = input_size;
			args[i].strides = strides;
			args[i].output_size = output_size;
			args[i].pad_top = pad_top;
			args[i].pad_left = pad_left;  
			args[i].epsilon_map = epsilon_map;
			args[i].is_used = is_used;
		        args[i].is_maxpool = is_maxpool;	
	    		pthread_create(&threads[i], NULL,function, (void*)&args[i]);
			idx_start = idx_end;
			idx_end = idx_start + idx_n;
	    		if(idx_end>num_out_neurons){
				idx_end = num_out_neurons;
			}
			if((i==num_threads-2)){
				idx_end = num_out_neurons;
				
			}
	  	}

		for (i = 0; i < num_threads; i = i + 1){
			pthread_join(threads[i], NULL);
		}
	}
	for(i=0; i < (int)num_out_neurons; i++){
		
		if(!is_used[i]){
			
			zonotope_delete_constrained_noise_symbol(pr, epsilon_map[i]->index, z);
			
		}
	}
	free(epsilon_map);
	free(is_used);
	
}

#ifdef __cplusplus
}
#endif

#endif

#include "maxpool_approx.h"

size_t handle_maxpool_layer(elina_manager_t *man, elina_abstract0_t *element, 
			   size_t *pool_size, size_t *input_size, size_t *predecessors){
	//assert(dimensionality==3);
	//printf("maxpool start\n");
	//fflush(stdout);
	//printf("maxpool start\n");
	//fflush(stdout);
	assert(pool_size[0]==2 && pool_size[1]==2 && pool_size[2]==1);
	//assert(stride[0]==2 && stride[1]==2 && stride[2]==1);
	
	size_t i,j,k;
	size_t * output_size = (size_t *)malloc(3*sizeof(size_t));
	for(i=0; i < 3; i++){
		output_size[i] = input_size[i]/pool_size[i];
	}


	size_t num_input_neurons = input_size[0]*input_size[1]*input_size[2];
	size_t num_out_neurons = output_size[0]*output_size[1]*output_size[2];

    	size_t o12 = output_size[1]*output_size[2];
   	size_t i12 = input_size[1]*input_size[2];
    	size_t p01 = pool_size[0]*pool_size[1];
	
	fppoly_t * fp = fppoly_of_abstract0(element);
	size_t numlayers = fp->numlayers;
	fppoly_add_new_layer(fp,num_out_neurons, MAXPOOL, NONE);
	size_t out_pos;
	double * inf = (double *) calloc(p01,sizeof(double));
	double * sup = (double *) calloc(p01,sizeof(double));
	size_t * pool_map = (size_t *)calloc(p01,sizeof(size_t));
	neuron_t ** out_neurons = fp->layers[numlayers]->neurons;
	fp->layers[numlayers]->predecessors = predecessors;
	size_t count = 0;
	for(out_pos=0; out_pos<num_out_neurons; out_pos++){
		size_t out_x = out_pos / o12;
		size_t out_y = (out_pos-out_x*o12) / output_size[2];
		size_t out_z = out_pos-out_x*o12 - out_y*output_size[2];
		size_t inp_x = out_x*pool_size[0];
		size_t inp_y = out_y*pool_size[1];
		size_t inp_z = out_z;
		size_t inp_pos = inp_x*i12 + inp_y*input_size[2] + inp_z;
		size_t pool_start_dim = out_pos*pool_size[0]*pool_size[1];
		//printf("inpXYZ: %zu, %zu, %zu %zu %zu\n", inp_x, inp_y, inp_z, out_pos, num_out_neurons);
        	//printf("outXYZ: %zu, %zu, %zu\n", out_x, out_y, out_z);
		//fflush(stdout);
		size_t x_shift, y_shift, l = 0;
		double sum_u = 0.0;
		double sum_l = 0.0;
		double max_u = -INFINITY;
		double max_l = -INFINITY;
		
		size_t max_l_var = 0.0; 
		size_t max_u_var = 0.0;
		size_t min_width_var = 0.0;
		double min_width = INFINITY;
		for(x_shift = 0; x_shift < pool_size[0]; x_shift++){
			for(y_shift = 0; y_shift < pool_size[1]; y_shift++){
				size_t pool_cur_dim = inp_pos + x_shift*i12 + y_shift*input_size[2];
				//printf("pool_cur_dim %zu %zu %zu\n",pool_cur_dim,fp->layers[numlayers-1]->dims,numlayers);
				//fflush(stdout);
				pool_map[l] = pool_cur_dim;
				// use the ReLU bounds from the previous layer
				double lb = -fp->layers[numlayers-1]->neurons[pool_cur_dim]->lb;
				double ub = fp->layers[numlayers-1]->neurons[pool_cur_dim]->ub;
				if(ub<=0){
					inf[l] = 0.0;
					sup[l] = 0.0;
				}
				else if(lb>0){
					inf[l] = lb;
					sup[l] = ub;
				}
				else{
					inf[l] = 0;
					sup[l] = ub;
				}
				//printf("inf: %g %g\n",inf[l],sup[l]);
				//fflush(stdout);
				sum_u = sum_u + sup[l];
				sum_l = sum_l + inf[l];
				if(sup[l]>max_u){
					max_u = sup[l];
					max_u_var = pool_map[l];
				}
				if(inf[l] > max_l){
					max_l = inf[l];
					max_l_var = pool_map[l];
				}
				if((ub-lb) < min_width){
					min_width = ub - lb;
					min_width_var = pool_map[l];
				}
				l++;
			}
		}
		
		bool flag = false;
		size_t var = 0;
		for(j=0; j < p01; j++){
			bool is_greater = true;
			for(k = 0;  k < p01; k++){
				if(k==j)continue;
				if((inf[k]==sup[k]) && (inf[j]>=sup[k])){
					continue;
				}
				else if((inf[j]==inf[k]) && (sup[j]==sup[k]) && (inf[j]==sup[j])){
					continue;
				}
				else if(inf[j]<=sup[k]){
					is_greater = false;
					break;
				}
			}
			if(is_greater){
				flag = true;
				var =pool_map[j];
				break;
			}
		}
		//printf("max_l: %gmax_u: %g\n",max_l,max_u);
		//fflush(stdout);
		if(flag){
		//if(0){
			//x_new = x_var
			count++;
			//printf("out_pos: %zu\n",out_pos);
			//fflush(stdout);
			double coeff[1];
			size_t dim[1];
			coeff[0] = 1;
			dim[0] = var;
			out_neurons[out_pos]->lexpr = create_sparse_expr(coeff,0,dim,1);
			out_neurons[out_pos]->uexpr = create_sparse_expr(coeff,0,dim,1);
			//out_neurons[out_pos]->expr = create_sparse_expr(coeff,0,dim,1);
		}
		else{
			//max_l	<= x_new <= max_u
			double lcoeff[1];
			size_t ldim[1];
			lcoeff[0] = 1;
			ldim[0] = max_l_var;
			//lcoeff[0] = 0;
			//ldim[0] = 0;
			//printf("max_l: %gmax_u: %g\n",max_l,max_u);
			//fflush(stdout);
			//expr_t * expr = (expr_t *)malloc(sizeof(expr_t));
			//expr->inf_coeff= expr->sup_coeff = expr->dim = NULL;
			//expr->size = 0;
			//expr->inf_cst = -max_l;
			
			//out_neurons[out_pos]->lexpr = NULL;//create_sparse_expr(NULL,max_l,NULL,0);
			out_neurons[out_pos]->lexpr = create_sparse_expr(lcoeff,0,ldim,1);
			//double *ucoeff = (double *)malloc(p01*sizeof(double));
			//size_t *udim = (size_t *)malloc(p01*sizeof(size_t));
			//for(j=0; j < p01; j++){
			//	ucoeff[j] = 1.0;
			//	udim[j] = pool_map[j];
			//}
			double ucoeff[1];
			size_t udim[1];
			ucoeff[0] = 0;
			udim[0] = 0;
			out_neurons[out_pos]->uexpr = create_sparse_expr(ucoeff,max_u,udim,1);
			//out_neurons[out_pos]->uexpr = create_sparse_expr(ucoeff,max_l-sum_l,udim,p01);
			//sort_sparse_expr(out_neurons[out_pos]->uexpr);
			//free(ucoeff);
			//free(udim);
			//out_neurons[out_pos]->lexpr = create_cst_expr(-max_l,max_l);
			//out_neurons[out_pos]->uexpr = create_cst_expr(-max_u,max_u);
		}
		out_neurons[out_pos]->lb = -max_l;
		out_neurons[out_pos]->ub = max_u;		
	}
	//update_state_using_previous_layers_parallel(man,fp,numlayers);
        free(inf);
	free(sup);
	free(pool_map);
	free(output_size);
        //printf("count: %zu\n",count);
	//fflush(stdout);
	//printf("return here2\n");
	//fppoly_fprint(stdout,man,fp,NULL);
	//fflush(stdout);
	return num_out_neurons;
}



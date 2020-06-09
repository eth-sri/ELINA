#include "lstm_approx.h"


expr_t * lexpr_unroll_lstm_layer(fppoly_internal_t *pr, expr_t * expr, neuron_t ** neurons){
	return NULL;
}


void create_lstm_layer(elina_manager_t *man, elina_abstract0_t *abs, size_t h, size_t *predecessors, size_t num_predecessors){
	fppoly_t *fp = fppoly_of_abstract0(abs);
	size_t numlayers = fp->numlayers;
	fppoly_add_new_layer(fp,h, predecessors, num_predecessors, false);
	fp->lstm_index = numlayers;
}

void handle_lstm_layer(elina_manager_t *man, elina_abstract0_t *abs, double **weights,  double *bias, size_t d, size_t h, size_t * predecessors, size_t num_predecessors){
	fppoly_t *fp = fppoly_of_abstract0(abs);
	fppoly_internal_t *pr = fppoly_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
	size_t lstm_index = fp->lstm_index;
	layer_t *layer = fp->layers[lstm_index];
	neuron_t **out_neurons = fp->layers[lstm_index]->neurons;
	fp->layers[lstm_index]->predecessors = predecessors;
	size_t i;
	neuron_t * neuron = neuron_alloc();
	bool first_time_step = (layer->h_t_inf==NULL && layer->h_t_sup==NULL);
	size_t k = h + d;
	if(first_time_step){
		layer->h_t_inf = (double*)malloc(h*sizeof(double));
		layer->h_t_sup = (double*)malloc(h*sizeof(double));
		layer->c_t_inf = (double*)malloc(h*sizeof(double));
		layer->c_t_sup = (double*)malloc(h*sizeof(double));
	}

	// TODO: Fix, debug: 	for(i=0; i< h; i++){
	for(i=0; i< 1; i++){
	  //printf("i = %d\n",(int)i);
		expr_t *f_t_lexpr, *i_t_lexpr, *o_t_lexpr, *c_t_lexpr;
		if(first_time_step){
			i_t_lexpr =  create_dense_expr(weights[i],bias[i],d);
			c_t_lexpr =  create_dense_expr(weights[h+i],bias[h+i],d);	
			f_t_lexpr =  create_dense_expr(weights[2*h+i],bias[2*h+i],d);
			o_t_lexpr =  create_dense_expr(weights[3*h+i],bias[3*h+i],d);
		}		
		else{
			expr_t * tmp1 = create_dense_expr(weights[i],bias[i],d+h);
			expr_t * tmp2 = create_dense_expr(weights[h+i],bias[h+i],d+h);
			expr_t * tmp3 = create_dense_expr(weights[2*h+i],bias[2*h+i],d+h);
			expr_t * tmp4 = create_dense_expr(weights[3*h+i],bias[3*h+i],d+h);
			i_t_lexpr = concretize_dense_sub_expr(pr, tmp1, layer->h_t_inf, layer->h_t_sup, d, d+h);
			c_t_lexpr = concretize_dense_sub_expr(pr, tmp2, layer->h_t_inf, layer->h_t_sup, d, d+h);
			f_t_lexpr = concretize_dense_sub_expr(pr, tmp3, layer->h_t_inf, layer->h_t_sup, d, d+h);
			o_t_lexpr = concretize_dense_sub_expr(pr, tmp4, layer->h_t_inf, layer->h_t_sup, d, d+h);
			free_expr(tmp1);	
			free_expr(tmp2);
			free_expr(tmp3);
			free_expr(tmp4);
		}

		//expr_print(f_t_lexpr);

		//printf("computing forget...\n");
		expr_t *f_t_uexpr = copy_expr(f_t_lexpr);
		expr_t *tmp_f_t_lexpr = copy_expr(f_t_lexpr);
		expr_t *tmp_f_t_uexpr = copy_expr(f_t_uexpr);
		double lb_f_t = get_lb_using_previous_layers(man, fp, tmp_f_t_lexpr, lstm_index);
		double ub_f_t = get_ub_using_previous_layers(man, fp, tmp_f_t_uexpr, lstm_index);
		/* free_expr(tmp_f_t_lexpr); */
		/* free_expr(tmp_f_t_uexpr); */

		neuron->lb = lb_f_t;
		neuron->ub = ub_f_t;
		//printf("forget gate before sigmoid: lb = %lf, ub = %lf\n",neuron->lb, neuron->ub);
		//expr_print(f_t_lexpr);
		//expr_print(f_t_uexpr);
		lb_f_t = apply_sigmoid_lexpr(pr, &f_t_lexpr, neuron);
		ub_f_t = apply_sigmoid_uexpr(pr, &f_t_uexpr, neuron);
		//printf("forget gate after sigmoid: lb_f_t = %lf, ub_f_t = %lf\n",lb_f_t,ub_f_t);
		//expr_print(f_t_lexpr);
		//expr_print(f_t_uexpr);
		//printf("forget gate done\n\n");

		//printf("computing input...\n");
		expr_t *i_t_uexpr = copy_expr(i_t_lexpr);
		expr_t *tmp_i_t_lexpr = copy_expr(i_t_lexpr);
		expr_t *tmp_i_t_uexpr = copy_expr(i_t_uexpr);
		double lb_i_t = get_lb_using_previous_layers(man, fp, tmp_i_t_lexpr,lstm_index);
		double ub_i_t = get_ub_using_previous_layers(man, fp, tmp_i_t_uexpr, lstm_index);	
		/* free_expr(tmp_i_t_lexpr); */
		/* free_expr(tmp_i_t_uexpr); */
		neuron->lb = lb_i_t;
		neuron->ub = ub_i_t;
		//printf("input gate before sigmoid: lb = %lf, ub = %lf\n",neuron->lb, neuron->ub);
		//expr_print(i_t_uexpr);
		lb_i_t = apply_sigmoid_lexpr(pr, &i_t_lexpr, neuron);
		ub_i_t = apply_sigmoid_uexpr(pr, &i_t_uexpr, neuron);
		//expr_print(i_t_uexpr);
		//printf("input gate after sigmoid: lb_i_t = %lf, ub_i_t = %lf\n",lb_i_t,ub_i_t);
		//printf("input gate done\n\n");

		//printf("computing output...\n");
		expr_t *o_t_uexpr = copy_expr(o_t_lexpr);
		expr_t *tmp_o_t_lexpr = copy_expr(o_t_lexpr);
		expr_t *tmp_o_t_uexpr = copy_expr(o_t_uexpr);
		double lb_o_t = get_lb_using_previous_layers(man, fp, tmp_o_t_lexpr, lstm_index);
		double ub_o_t = get_ub_using_previous_layers(man, fp, tmp_o_t_uexpr, lstm_index);
		/* free_expr(tmp_o_t_lexpr); */
		/* free_expr(tmp_o_t_uexpr); */

		neuron->lb = lb_o_t;
		neuron->ub = ub_o_t;		
		//printf("output gate before sigmoid: lb = %lf, ub = %lf\n",neuron->lb, neuron->ub);
		lb_o_t = apply_sigmoid_lexpr(pr, &o_t_lexpr, neuron);
		ub_o_t = apply_sigmoid_uexpr(pr, &o_t_uexpr, neuron);
		//printf("output gate after sigmoid: lb = %lf, ub = %lf\n",lb_o_t,ub_o_t);
		out_neurons[i]->lb = lb_o_t;
		out_neurons[i]->ub = ub_o_t;
		out_neurons[i]->lexpr = o_t_lexpr;
		out_neurons[i]->uexpr = o_t_uexpr;
		//printf("output gate done\n\n");

		//printf("computing control state...\n");
		//printf("control expression:\n");
		//expr_print(c_t_lexpr);
		//printf("...\n");
		expr_t *c_t_uexpr = copy_expr(c_t_lexpr);
		expr_t *tmp_c_t_lexpr = copy_expr(c_t_lexpr);
		expr_t *tmp_c_t_uexpr = copy_expr(c_t_uexpr);
		double lb_c_t = get_lb_using_previous_layers(man, fp, tmp_c_t_lexpr, lstm_index);
		double ub_c_t = get_ub_using_previous_layers(man, fp, tmp_c_t_uexpr, lstm_index);
		neuron->lb = lb_c_t;
		neuron->ub = ub_c_t;
		//expr_print(c_t_lexpr);
		//expr_print(c_t_uexpr);
		//printf("control before tanh: lb = %lf, ub = %lf\n",neuron->lb,neuron->ub);
		lb_c_t = apply_tanh_lexpr(pr,&c_t_lexpr, neuron);
		ub_c_t = apply_tanh_uexpr(pr,&c_t_uexpr, neuron);			
		//printf("control after tanh: lb = %lf, ub = %lf\n",lb_c_t,ub_c_t);
		//printf("control expression:\n");
		//expr_print(c_t_lexpr);
		//expr_print(c_t_uexpr);

		//printf("=======================\n");

		//printf("multiplying control by input:\n");
		expr_t *tmp_l, *tmp_u;
		double width1 = ub_i_t + lb_i_t;
		double width2 = ub_c_t + lb_c_t;
		tmp_l = c_t_lexpr;
		tmp_u = c_t_uexpr;
		//printf("control: [%lf %lf], input: [%lf %lf]\n",lb_c_t,ub_c_t,lb_i_t,ub_i_t);
		//printf("control before multiplying by input:\n");
		//expr_print(c_t_lexpr);
		//expr_print(c_t_uexpr);
		if(width1 < width2){
		  //printf("concretize input\n");
			c_t_lexpr = multiply_expr(pr,c_t_lexpr,lb_i_t,ub_i_t);
			c_t_uexpr = multiply_expr(pr,c_t_uexpr,lb_i_t,ub_i_t);
		}
		else{
		  //printf("concretize control\n");
            if(lb_c_t<0){
                c_t_lexpr = multiply_expr(pr,i_t_lexpr,lb_c_t,ub_c_t);
                c_t_uexpr = multiply_expr(pr,i_t_uexpr,lb_c_t,ub_c_t);
            }
            else if(ub_c_t<0){
                c_t_lexpr = multiply_expr(pr,i_t_uexpr,lb_c_t,ub_c_t);
                c_t_uexpr = multiply_expr(pr,i_t_lexpr,lb_c_t,ub_c_t);
            }
            else{
                c_t_lexpr = multiply_expr(pr,i_t_lexpr,0,0);
                c_t_uexpr = multiply_expr(pr,i_t_uexpr,0,0);
                double tmp1, tmp2;
                elina_double_interval_mul_expr_coeff(pr,&tmp1,&tmp2,lb_i_t,ub_i_t,lb_c_t,ub_c_t);
                c_t_lexpr->inf_cst += tmp1;
                c_t_lexpr->sup_cst += tmp2;
                c_t_uexpr->inf_cst += tmp1;
                c_t_uexpr->sup_cst += tmp2;
            }
		}

		//printf("control after multiplying by input:\n");
		//expr_print(c_t_lexpr);
		//expr_print(c_t_uexpr);

		free_expr(tmp_l);
		free_expr(tmp_u);

		//printf("here\n\n\n");
		//printf("====================================\n");
		
		if(!first_time_step){
            if(layer->c_t_inf[i]<0){
                tmp_l = multiply_expr(pr,f_t_lexpr,layer->c_t_inf[i],layer->c_t_sup[i]);
                tmp_u = multiply_expr(pr,f_t_uexpr,layer->c_t_inf[i],layer->c_t_sup[i]);
            }
            else if(layer->c_t_sup[i]<0){
                tmp_l = multiply_expr(pr,f_t_uexpr,layer->c_t_inf[i],layer->c_t_sup[i]);
                tmp_u = multiply_expr(pr,f_t_lexpr,layer->c_t_inf[i],layer->c_t_sup[i]);
            }
            else{
                tmp_l = multiply_expr(pr,f_t_lexpr,0,0);
                tmp_u = multiply_expr(pr,f_t_uexpr,0,0);
                double tmp1, tmp2;
                elina_double_interval_mul_expr_coeff(pr,&tmp1,&tmp2,lb_f_t,ub_f_t,layer->c_t_inf[i],layer->c_t_sup[i]);
                tmp_l->inf_cst += tmp1;
                tmp_l->sup_cst += tmp2;
                tmp_u->inf_cst += tmp1;
                tmp_u->sup_cst += tmp2;
            }
			add_expr(pr,c_t_lexpr,tmp_l);
			add_expr(pr,c_t_uexpr,tmp_u);
			free_expr(tmp_l);
			free_expr(tmp_u);
		}
		layer->c_t_inf[i] = get_lb_using_previous_layers(man, fp, c_t_lexpr, lstm_index);
		layer->c_t_sup[i] = get_ub_using_previous_layers(man, fp, c_t_uexpr, lstm_index);

		neuron->lb = layer->c_t_inf[i];
		neuron->ub = layer->c_t_sup[i];

		//printf("c_t ---> lb = %lf, ub = %lf\n", neuron->lb, neuron->ub);

		lb_c_t = apply_tanh_lexpr(pr,&c_t_lexpr, neuron);
		ub_c_t = apply_tanh_uexpr(pr,&c_t_uexpr, neuron);
		
		width1 = ub_o_t + lb_o_t;
		width2 = ub_c_t + lb_c_t; 

		expr_t * h_t_lexpr, *h_t_uexpr;
		if(width1<width2){
			h_t_lexpr = multiply_expr(pr,c_t_lexpr,lb_o_t,ub_o_t);
			h_t_uexpr = multiply_expr(pr,c_t_uexpr,lb_o_t,ub_o_t);
		}
		else{
			h_t_lexpr =  multiply_expr(pr,o_t_lexpr,lb_c_t,ub_c_t);
			h_t_uexpr =  multiply_expr(pr,o_t_uexpr,lb_c_t,ub_c_t);
		}

		layer->h_t_inf[i] = get_lb_using_previous_layers(man, fp, h_t_lexpr, lstm_index);
		layer->h_t_sup[i] = get_ub_using_previous_layers(man, fp, h_t_uexpr, lstm_index);

		free_expr(f_t_lexpr);
		free_expr(f_t_uexpr);
		free_expr(i_t_lexpr);
		free_expr(i_t_uexpr);
		free_expr(c_t_lexpr);
		free_expr(c_t_uexpr);
		free_expr(h_t_lexpr);
		free_expr(h_t_uexpr);
	}
	free_neuron(neuron);
	//update_state_using_previous_layers_parallel(man,fp,numlayers);
	return;
}


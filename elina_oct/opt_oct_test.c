/*
	Copyright 2016 Software Reliability Lab, ETH Zurich

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
#include "oct_internal.h"
#include "opt_oct_hmat.h"

#define num_comp 5

int dim, e_dim;
opt_oct_internal_t *pr1;
oct_internal_t *pr2;
ap_manager_t *moo, *mo;
array_comp_list_t * acl1,*acl2;
ap_dim_t  *tdim;

double abs_d(double num){
	if(num < 0){
		return -num;
	}
	return num;
}

void print_oct(bound_t *b, int dim){
	for(int i = 0; i < 2*dim; i++){
		for(int j = 0; j < 2*dim; j++){
			bound_print(b[matpos2(i,j)]);
			fprintf(stdout,"\t");
		}
		fprintf(stdout,"\n");
	}
	fprintf(stdout,"\n");
}


void random_oct(oct_t **o, opt_oct_t **oo, array_comp_list_t * acl, int dim)
{
  opt_oct_mat_t * oom;
  double *d;
  bound_t* b;
  int i,j;
  int x,y;
  num_t n;
  *oo = opt_oct_alloc_internal(pr1, dim,0);
  *o = oct_alloc_internal(pr2,dim,0);
  (*oo)->m = opt_hmat_alloc_top(dim);
  free_array_comp_list((*oo)->m->acl);
  (*o)->m = hmat_alloc_top(pr2,dim);
  oom = (*oo)->m;
  oom->acl = acl;
  //oom->is_dense = false;
  d = oom->mat;
  b = (*o)->m;
  num_init(n);
  comp_list_t * cl = acl->head;
  while(cl != NULL){
	  unsigned short int * ca = to_sorted_array(cl,dim);
	  unsigned short int comp_size = cl->size;
	  for (i=0;i<2*comp_size;i++){
	    int i1 = (i%2==0)? 2*ca[i/2] : 2*ca[i/2] + 1;
	    for (j=0;j< 2*comp_size;j++) {
	      int j1 = (j%2==0)? 2*ca[j/2]: 2*ca[j/2]+1;
	      if(j1 > (i1|1)){
		    break;
	      }
	    
	     // if (i1==j1) continue;
	        int ind = j1 + (((i1 + 1)*(i1 + 1))/2);
		y = rand()%4+1;
		x = rand()%20+2;      
		num_set_double(n,x/y);
		bound_set_num(b[ind],n);
		d[ind] = x/y;
	    }
	}
	free(ca);
	cl = cl->next;
  }
  num_clear(n);
}


void random_array_comp_list(array_comp_list_t * acl){
	for(int i = 0; i < num_comp; i++){
		comp_list_t *cl = create_comp_list();
		int comp_size = dim/num_comp;
		int i1 = comp_size*i;
		for(int j = 0; j < comp_size;j++){
			int j1 = rand()%comp_size;
			int num = i1 + j1;
			//printf("%d\n",num);
			if(!contains_comp(cl,num)){
				insert_comp(cl,num);
			}
			//insert_comp(cl,i*3+2);
			
		}
		if(cl->size > 0){
			insert_comp_list(acl,cl);
		}
		else{
			free_comp_list(cl);
		}
	}
	
}

void generate_input_unary(oct_t ** o, opt_oct_t **oo){
	array_comp_list_t * acl = create_array_comp_list();
	random_array_comp_list(acl);
	random_oct(o,oo,acl,dim);
}

void generate_input_binary(oct_t ** o1, oct_t **o2,opt_oct_t **oo1, opt_oct_t **oo2){
	generate_input_unary(o1,oo1);
	generate_input_unary(o2,oo2);

}

bool verify(oct_t *o, opt_oct_t *oo){
	opt_oct_mat_t *oom = oo->closed==NULL ? oo->m : oo->closed; 
	if(!oom->is_dense){
		convert_to_dense_mat(oom,o->dim,false);
	}
       //oom->is_dense = true;
       double *d = oom->mat;
       bound_t *b = o->closed==NULL ? o->m : o->closed;
       if(!d && !b){
		//printf("Output Verified\n");
		return 1;
	}
	for(unsigned short int i = 0; i < 2*o->dim; i++){
		for(unsigned short int j = 0; j <= (i|1);j++,b++,d++){
			if(i==j)continue;
			double d1 = *d;
			double d2;
			double_set_num(&d2, *b);
			if(abs_d(d1-d2) > 10e-6){
				return 0;
			}
			
		}
	}
	
	return 1;
}


/********
	Test for Join
********/
void test_join(void){
	oct_t *o1,*o2,*o;
	opt_oct_t *oo1,*oo2,*oo;
	generate_input_binary(&o1,&o2,&oo1,&oo2);
	
	o = oct_join(mo,false,o1,o2);
	oo = opt_oct_join(moo,false,oo1,oo2);
	if(verify(o,oo)){
		printf("Join Verified\n");
	}
	else{
		printf("Join Not Verified\n");
	}
	opt_oct_free(moo,oo);
	opt_oct_free(moo,oo1);	
	opt_oct_free(moo,oo2);	
	oct_free(mo,o); 
	oct_free(mo,o1); 
	oct_free(mo,o2); 
	
}

/********
	Test for Meet
********/
void test_meet(void){
	oct_t *o1,*o2,*o;
	opt_oct_t *oo1,*oo2,*oo;
	generate_input_binary(&o1,&o2,&oo1,&oo2);
	
	o = oct_meet(mo,false,o1,o2);
	oo = opt_oct_meet(moo,false,oo1,oo2);
	if(verify(o,oo)){
		printf("Meet Verified\n");
	}
	else{
		printf("Meet Not Verified\n");
	}
	opt_oct_free(moo,oo);
	opt_oct_free(moo,oo1);	
	opt_oct_free(moo,oo2);	
	oct_free(mo,o); 
	oct_free(mo,o1); 
	oct_free(mo,o2); 
	
}


/*********
	Test for Join Array
*********/
void test_join_array(){
	oct_t *oa[10], *o;
	opt_oct_t *ooa[10], *oo;
	for(int i = 0; i < 10; i++){
		generate_input_unary(&oa[i],&ooa[i]);
	}
	o = oct_join_array(mo ,oa, 10);
	oo = opt_oct_join_array(moo, ooa, 10);
	if(verify(o,oo)){
		printf("Join Array verified\n");
	}	
	else{
		printf("Join Array Not Verified\n");
	}
	for(int i = 0; i < 10; i++){
		oct_free(mo,oa[i]);
		opt_oct_free(mo,ooa[i]);
	}
	oct_free(mo,o);
	opt_oct_free(moo,oo);
}


/*********
	Test for Meet Array
*********/
void test_meet_array(){
	oct_t *oa[10], *o;
	opt_oct_t *ooa[10], *oo;
	for(int i = 0; i < 10; i++){
		generate_input_unary(&oa[i],&ooa[i]);
	}
	o = oct_meet_array(mo ,oa, 10);
	oo = opt_oct_meet_array(moo, ooa, 10);
	if(verify(o,oo)){
		printf("Meet Array verified\n");
	}	
	else{
		printf("Meet Array Not Verified\n");
	}
	for(int i = 0; i < 10; i++){
		oct_free(mo,oa[i]);
		opt_oct_free(mo,ooa[i]);
	}
	oct_free(mo,o);
	opt_oct_free(moo,oo);
}


/********
	Test for Widening
********/
void test_widening(void){
	oct_t *o1,*o2,*o;
	opt_oct_t *oo1,*oo2,*oo;
	generate_input_binary(&o1,&o2,&oo1,&oo2);
	
	o = oct_widening(mo,o1,o2);
	oo = opt_oct_widening(moo,oo1,oo2);
	if(verify(o,oo)){
		printf("Widening Verified \n");
	}
	else{
		printf("Widening Not verified\n");
	}
	opt_oct_free(moo,oo);
	opt_oct_free(moo,oo1);	
	opt_oct_free(moo,oo2);	
	oct_free(mo,o); 
	oct_free(mo,o1); 
	oct_free(mo,o2); 
	
}

/*********
	Test for Widening Thresholds
*********/
void test_widening_thresholds(void){
	oct_t *o1, *o2, *o;
	opt_oct_t *oo1, *oo2, *oo;
	generate_input_binary(&o1,&o2,&oo1,&oo2);
	ap_scalar_t* t[10];
    	for (int i=0;i<10;i++) t[i]=ap_scalar_alloc_set_double((rand()%30-15)*0.25);
	o = oct_widening_thresholds(mo,o1,o2,t,10);
	oo = opt_oct_widening_thresholds(moo,oo1,oo2,t,10);
	if(verify(o,oo)){
		printf("Widening Threshold Verified\n");
	}
	else{
		printf("Widening Threshold Not Verified\n");
	}
	for(int i = 0; i < 10; i++)
	ap_scalar_free(t[i]);
	oct_free(mo,o);
	oct_free(mo,o1);
	oct_free(mo,o2);
	opt_oct_free(mo,oo);
	opt_oct_free(mo,oo1);
	opt_oct_free(mo,oo2);
}

/********
	Test for Narrowing
********/
void test_narrowing(void){
	opt_oct_t *oo1,*oo2,*oo;
	oct_t *o1,*o2,*o;
	
	generate_input_binary(&o1,&o2,&oo1,&oo2);
	
	o = oct_narrowing(mo,o1,o2);
	oo = opt_oct_narrowing(moo,oo1,oo2);
	if(verify(o,oo)){
		printf("Narrowing Verified \n");
	}
	else{
		printf("Narrowing Not verified\n");
	}
	opt_oct_free(moo,oo);
	opt_oct_free(moo,oo1);	
	opt_oct_free(moo,oo2);	
	oct_free(mo,o); 
	oct_free(mo,o1); 
	oct_free(mo,o2); 
	
}

/********
	Test for Bound Dim
********/
void test_bound_dim(void){
	oct_t *o;
	opt_oct_t *oo;
	generate_input_unary(&o, &oo);
	bool flag = 0;
	for(int i = 0; i < 10;i++){
		int pos = rand()%dim;
		//printf("Pos is %d\n",pos);
		ap_interval_t * inv1 = opt_oct_bound_dimension(moo,oo,pos);
		ap_interval_t * inv2 = oct_bound_dimension(mo,o,pos);
		//verify(oo1, o1);
		if(!ap_interval_equal(inv1,inv2)){
			fprintf(stdout,"Bound Dim not verified\n");
			flag = 1;
			ap_interval_free(inv1);
			ap_interval_free(inv2);
			break;
		}
		ap_interval_free(inv1);
		ap_interval_free(inv2);
	}
	if(!flag){
		fprintf(stdout, "Bound Dim Verified\n");
	}
	opt_oct_free(moo,oo);
	oct_free(mo,o); 
	
}


/********
	Test for Add Epsilon
********/
void test_add_epsilon(void){
	oct_t *o, *o1;
	opt_oct_t *oo, *oo1;
	generate_input_unary(&o, &oo);
	bool flag = 0;
	for(int i = 0; i < 10;i++){
		int pos = rand()%dim;
		ap_scalar_t *r = ap_scalar_alloc();
		ap_scalar_set_double(r,pos);
		oo1 = opt_oct_add_epsilon(moo,oo,r);
		o1 = oct_add_epsilon(mo,o,r);
		//oo1 = opt_oct_add_epsilon_bin(moo,oo,oo2,r);
		//o1 = oct_add_epsilon_bin(mo,o,o2,r);
		if(!verify(o1,oo1)){
			printf("Add Epsilong not verified\n");
			opt_oct_free(moo,oo1);
			oct_free(mo,o1);
			flag = 1;
			break;
		}
		ap_scalar_free(r);
		opt_oct_free(moo,oo1);
		oct_free(mo,o1);
	}
	if(!flag){
		printf("Add Epsilon Verified\n");
	}
	opt_oct_free(moo,oo);
	oct_free(mo,o); 
}


/********
	Test for Add Epsilon Binary
********/
void test_add_epsilon_bin(void){
	oct_t *o1, *o2, *o;
	opt_oct_t *oo1, *oo2, *oo;
	generate_input_binary(&o1, &o2, &oo1, &oo2);
	bool flag = 0;
	for(int i = 0; i < 10;i++){
		int pos = rand()%dim;
		ap_scalar_t *r = ap_scalar_alloc();
		ap_scalar_set_double(r,pos);
		//oo1 = opt_oct_add_epsilon(moo,oo,r);
		//o1 = oct_add_epsilon(mo,o,r);
		oo = opt_oct_add_epsilon_bin(moo,oo1,oo2,r);
		o = oct_add_epsilon_bin(mo,o1,o2,r);
		if(!verify(o, oo)){
			printf("Add Epsilon Bin Not Verified\n");
			opt_oct_free(moo,oo);
			oct_free(mo,o);
			flag = 1;
			break;
		}
		ap_scalar_free(r);
		opt_oct_free(moo,oo);
		oct_free(mo,o);
	}
	if(!flag){
		printf("Add Epsilon Bin Verified\n");
	}
	opt_oct_free(moo,oo1);
	oct_free(mo,o1); 
	opt_oct_free(moo,oo2);	
	oct_free(mo,o2); 
	
}


int main(int argc, char **argv){
	dim = atoi(argv[1]);
	int size = 2*dim*(dim+1);
	moo = opt_oct_manager_alloc();
	if(!moo)return 1;
	mo = oct_manager_alloc();
	if(!mo)return 1;
	for (int i=0;i<AP_FUNID_SIZE;i++) {
		moo->option.funopt[i].flag_exact_wanted = true;
    		moo->option.funopt[i].flag_best_wanted = true;
		moo->option.funopt[i].algorithm = 0;
    		mo->option.funopt[i].flag_exact_wanted = true;
    		mo->option.funopt[i].flag_best_wanted = true;
		mo->option.funopt[i].algorithm = 0;
  	}
  	for (int i=0;i<AP_EXC_SIZE;i++){
		moo->option.abort_if_exception[i] = true;
    		mo->option.abort_if_exception[i] = true;
  	}
	
	pr1 = opt_oct_init_from_manager(moo, 0, dim);
	pr2 = oct_init_from_manager(mo,0,dim);
	tdim = (ap_dim_t *)malloc(e_dim*sizeof(ap_dim_t));
	test_join();
	test_meet();
	test_join_array();
	test_meet_array();
	test_widening();
	test_widening_thresholds();
	test_narrowing();
	test_bound_dim();
	test_add_epsilon();
	test_add_epsilon_bin();
	ap_manager_free(moo);
	ap_manager_free(mo);
	//free(tdim);
	return 0;
}

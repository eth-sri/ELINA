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

#include "comp_list.h"

void unite_comp_lists(comp_list_t *cl1,comp_list_t *cl2,char *map,int i, int j,unsigned short int n){
	comp_t *c2 = cl2->head;
	while(c2!=NULL){
		unsigned short int num = c2->num;
		if(!map[n*i + num]){
			insert_comp(cl1,num);
			map[n*i + num] = 1;
		}
		c2 = c2->next;
	}
}

array_comp_list_t * extract(double *m, unsigned short int n){
	
	comp_list_t **comp = (comp_list_t **)calloc(n,sizeof(comp_list_t *));
	char *map = (char *)calloc(n*n,sizeof(char));
	unsigned short int * owner = (unsigned short int *)calloc(n,sizeof(unsigned short int));
	for(int i = 0; i < n; i++){
		comp[i]= create_comp_list();
		owner[i] = i;
	}
	for(int i =0; i < n; i++){
		for(int j = 0; j <=i; j++){
			int flag = 0;
			for(int i_1 = 0; i_1 < 2; i_1++){
				for(int j_1 = 0;j_1 < 2; j_1++){
					int i_2 = i*2 + i_1;
					int j_2 = j*2 + j_1;
					int ind = j_2 + (((i_2+1)*(i_2+1))/2);
					if((i_2 != j_2) && (m[ind] != INFINITY)){
						flag = 1;
					}
				}
			}
			if(flag){
				unsigned short int oi,oj;
			
				//printf("%d\t%d\n",i,j);
				oi = owner[i];
				//int s1 = comp[oi].size;
				if(!map[n*oi + i]){
					insert_comp(comp[oi],i);
					map[n*oi + i] = 1;
				}				
				
				oj  = owner[j];
				if(!map[n*oj + j]){
					insert_comp(comp[oj],j);
					map[n*oj + j] = 1;
				}
				unite_comp_lists(comp[oi],comp[oj],map,oi,oj,n);
				comp_t *c = comp[oi]->head;
				for(int k = 0; k < comp[oi]->size; k++){
					unsigned short int num = c->num;
					if(num!=oi ){
						owner[num] = owner[oi];
						//free((blk + k)->ind);
					}
					c = c->next;
				}
			}
		}
	}
	free(map);
	int size = 0;
	//double * map = (double *)calloc(n,sizeof(double));
	array_comp_list_t * res = create_array_comp_list();
	for(int i = 0; i < n; i++){
		if((owner[i]==i) && (comp[i]->size > 0)){
			comp_list_t * dst = copy_comp_list(comp[i]);
			insert_comp_list(res,dst);
		}
	}
	for(int i = 0; i < n; i++){
		free_comp_list(comp[i]);
	}
	free(comp);
	free(owner);
	return res;
}


array_comp_list_t * extract_comps(char *m, unsigned short int n){
	comp_list_t **comp = (comp_list_t **)calloc(n,sizeof(comp_list_t *));
	char *map = (char *)calloc(n*n,sizeof(char));
	unsigned short int * owner = (unsigned short int *)calloc(n,sizeof(unsigned short int));
	for(int i = 0; i < n; i++){
		//comp_list_t * cl =  create_comp_list();
		//comp + i = create_comp_list();
		//comp[i] = *cl;
		comp[i]= create_comp_list();
		//comp[i]->size = 0;
		owner[i] = i;
		//comp[i].owner = i;
	}
	for(int i =0; i < n; i++){
		for(int j = 0; j < n; j++){
			if((i!=j) && (m[n*i+j])){
				unsigned short int oi,oj;
			
				//printf("%d\t%d\n",i,j);
				oi = owner[i];
				//int s1 = comp[oi].size;
				if(!map[n*oi + i]){
					insert_comp(comp[oi],i);
					map[n*oi + i] = 1;
				}				
				
				oj  = owner[j];
				if(!map[n*oj + j]){
					insert_comp(comp[oj],j);
					map[n*oj + j] = 1;
				}
				unite_comp_lists(comp[oi],comp[oj],map,oi,oj,n);
				comp_t *c = comp[oi]->head;
				for(int k = 0; k < comp[oi]->size; k++){
					unsigned short int num = c->num;
					if(num!=oi ){
						owner[num] = owner[oi];
						//free((blk + k)->ind);
					}
					c = c->next;
				}
			}
		}
	}
	free(map);
	int size = 0;
	//double * map = (double *)calloc(n,sizeof(double));
	array_comp_list_t * res = create_array_comp_list();
	for(int i = 0; i < n; i++){
		if((owner[i]==i) && (comp[i]->size > 0)){
			comp_list_t * dst = copy_comp_list(comp[i]);
			insert_comp_list(res,dst);
		}
		//else{
			//free(comp + i);
		//}
	}
	//for(int i = 0; i < n; i++){
		//free_comp_list(comp + i);
	//}
	for(int i = 0; i < n; i++){
		free_comp_list(comp[i]);
	}
	free(comp);
	free(owner);
	return res;
}


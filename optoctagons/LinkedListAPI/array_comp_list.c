/*
	Copyright 2015 Department of Computer Science, ETH Zurich

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


#include "comp_list.h"

void print_array_comp_list(array_comp_list_t *acl, unsigned short int n){
	fprintf(stdout,"%d\n",acl->size);
	if(!acl || !acl->size){
		return;
	}
	comp_list_t * cl = acl->head;
	
	while(cl!=NULL){
		fprintf(stdout,"c\n");
		print_comp_list(cl,n);
		cl = cl->next;
	}
	fprintf(stdout,"\n");
	fflush(stdout);
}

void free_array_comp_list(array_comp_list_t * acl){
	if(acl==NULL){
		return;
	}
	
	comp_list_t *cl = acl->head;
	
	while(cl!=NULL){
		comp_list_t *temp = cl;
		cl = cl->next;
		free_comp_list(temp);
		temp = NULL;
	}
	free(acl);
	acl = NULL;
}

array_comp_list_t * create_array_comp_list(){
	array_comp_list_t * acl = (array_comp_list_t *)malloc(sizeof(array_comp_list_t));
	acl -> head = NULL; 
	//acl->tail = NULL;
	acl->size = 0;
	return acl;
}

array_comp_list_t * copy_array_comp_list(array_comp_list_t *src){
	if(!src){
		return NULL;
	}
	array_comp_list_t * dst = create_array_comp_list();
	comp_list_t * s = src->head;
	while(s!=NULL){
		comp_list_t * d = copy_comp_list(s);
		insert_comp_list(dst,d);
		s = s->next;
	}
	return dst;
}

void insert_comp_list(array_comp_list_t *acl, comp_list_t * cl){
	comp_list_t *h = acl->head;
	if(!h){
		acl->head = cl;
		acl->head->next=NULL;
		//acl->tail = cl;
		acl->size++;
		return;
	}
	//acl->tail->next = cl;
	//acl->tail = acl->tail->next;
	cl->next = acl->head;
	acl->head = cl;
	acl->size++;
	return;
}

void insert_comp_list_with_union(array_comp_list_t * acl, comp_list_t *cl,unsigned short int n){
	comp_list_t * cl1 = acl->head;
	char *map = create_map(cl,n);
	while(cl1!=NULL){
		if(is_disjoint(cl,cl1,n)){
			cl1 = cl1->next;
			continue;
		}
		else if(is_included(cl1,cl,n)){
			remove_comp_list(acl,cl1);
		}
		else if(is_included(cl,cl1,n)){
			free(map);
			free_comp_list(cl);
			return;
		}
		
		else{
			union_comp_list(cl,cl1,map);
			remove_comp_list(acl,cl1);
			
		}
		cl1 = cl1->next;
	}
	insert_comp_list(acl,cl);
	free(map);
}

comp_list_t * find(array_comp_list_t *acl,unsigned short int num){
	unsigned short int num_comp = acl->size;
	comp_list_t *cl = acl->head;
	
	for(int l = 0; l < num_comp; l++){
		if(contains_comp(cl,num)){
			return cl;
		}
		cl = cl->next;
	}
	return NULL;
}

void remove_comp_list(array_comp_list_t *acl, comp_list_t *clr){
	if(!acl || !clr){
		return;
	}
	
	unsigned short int num_comp = acl->size;
	if(!num_comp){
		return;
	}
	comp_list_t * cl = acl->head;
	if(!cl){
		return;
	}
	if(cl==clr){
		acl->head = acl->head->next;
		acl->size--;
		free_comp_list(cl);
		return;
	}
	if(num_comp==1){
		return;
	}
	//cl = cl->next;
	for(int l = 1; l < num_comp; l++){
		if(cl->next==clr){
			comp_list_t *temp = cl->next;
			cl->next = cl->next->next;
			free_comp_list(temp);
			acl->size--;
			return;
		}
		cl = cl->next;
	}
}

int is_equal_array_comp_list(array_comp_list_t *acl1, array_comp_list_t *acl2, unsigned short int n){
	int num_comp1 = acl1->size;
	int num_comp2 = acl2->size;
	if(num_comp1 != num_comp2){
		return 0;
	}
	
	char *map2 = (char *)calloc(n,sizeof(char));
	comp_list_t * cl1 = acl1->head;
	while(cl1 != NULL){
		char *map1 = create_map(cl1,n);
		comp_list_t * cl2 = acl2->head;
		int flag = 0;
		while(cl2!=NULL){
			char *map2 = create_map(cl2,n);
			if(is_equal_map(map1,map2,n)){
				flag = 1;
				free(map2);
				break;	
			}
			cl2 = cl2->next;
			free(map2);
		}
		if(!flag){
			free(map1);
			return 0;
		}
		free(map1);
		cl1 = cl1->next;
	}
	
	return 1;
}


/******
	Check if acl1 is greater than acl2
******/
int is_lequal_array_comp_list(array_comp_list_t * acl1, array_comp_list_t * acl2, unsigned short int n){
	char * map1 = (char *)calloc(n,sizeof(char));
	char * map2 = (char *)calloc(n,sizeof(char));
	comp_list_t * cl1 = acl1->head;
	while(cl1 != NULL){
		comp_t * c1 = cl1->head;
		while(c1!=NULL){
			unsigned short int num = c1->num;
			map1[num] = 1;
			c1 = c1->next;
		}
		cl1 = cl1->next;
	}
	comp_list_t * cl2 = acl2->head;
	while(cl2 != NULL){
		comp_t * c2 = cl2->head;
		while(c2 != NULL){
			unsigned short int num = c2->num;
			map2[num] = 1;
			c2 = c2->next;
		}
		cl2 = cl2->next;
	}
	for(unsigned short int i = 0; i < n; i++){
		if(!map2[i]){
			if(map1[i]){
				return 0;
			}
		}
		
	}
	return 1;
}

void clear_array_comp_list(array_comp_list_t *acl){
	comp_list_t * cl = acl->head;
	while(cl!=NULL){
		remove_comp_list(acl,cl);
		cl = cl->next;
	}
}

int is_connected(array_comp_list_t *acl, unsigned short int i, unsigned short int j){
	comp_list_t * li = find(acl,i);
	comp_list_t * lj = find(acl,j);
	if((li==NULL) && (lj==NULL)){
		return 0;
	}
	if(li!=lj){
		return 0;
	}
	return 1;
}

char * create_array_map(array_comp_list_t * acl, unsigned short int n){
	char * map = (char *)calloc(n,sizeof(char));
	comp_list_t * cl = acl->head;
	unsigned short int l = 1;
	while(cl!=NULL){
		comp_t * c = cl->head;
		while(c!=NULL){
			unsigned short int num = c->num;
			map[num] = l;
			c = c->next;
		}
		l++;
		cl = cl->next;
	}
	return map;
}


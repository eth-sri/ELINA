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

/***

***/

int is_included(comp_list_t *cl1, comp_list_t *cl2, unsigned short int n){
	char *map = (char *)calloc(n,sizeof(char));
	//int s = c2->size;
	comp_t *c = cl2->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		map[num] = 1;
		c = c->next;
	}
	//s = c1->size;
	c = cl1->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		if(!map[num]){
			free(map);
			return 0;
		}
		c = c->next;
	}
	free(map);
	return 1;
}

int is_disjoint(comp_list_t * cl1, comp_list_t *cl2, unsigned short int n){
	char * map = (char *)calloc(n,sizeof(char));
	//int s = c1->size;
	comp_t *c = cl1->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		map[num] = 1;
		c = c->next;
	}
	//s = c2->size;
	c = cl2->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		if(map[num]){
			free(map);
			return 0;
		}
		c=c->next;
	}
	free(map);
	return 1;
}


char * create_map(comp_list_t *cl, unsigned short int n){
	char *map = (char *)calloc(n,sizeof(char));
	comp_t *c = cl->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		map[num] = 1;
		c = c->next;
	}
	return map;
}

void union_comp_list(comp_list_t * cl1,comp_list_t * cl2, char *map){
	comp_t *c2 = cl2->head;
	while(c2!=NULL){
		unsigned short int num = c2->num;
		if(!map[num]){
			insert_comp(cl1,num);
			map[num] = 1;
		}
		c2 = c2->next;
	}
}

void union_comp_list_direct(comp_list_t *cl1, comp_list_t *cl2, unsigned short int n){
	char *map = create_map(cl1,n);
	union_comp_list(cl1,cl2,map);
	free(map);
}


comp_list_t * comp_array_union_direct(array_comp_list_t * acl1, array_comp_list_t *acl2, unsigned short int *om1, unsigned short int *om2,  unsigned short int n){
	char* map = (char *)calloc(n,sizeof(char));
	int nc1 = acl1->size;
	int nc2 = acl2->size;
	comp_list_t * res = create_comp_list();
	comp_list_t *cl1 = acl1->head;
	for(int i = 0; i < nc1; i++){
		if(!om1[i]){
			cl1 = cl1->next;
			continue;
		}
		comp_t *c1 = cl1->head;
		while(c1!=NULL){
			unsigned short int num = c1->num;
			if(!map[num]){
				map[num] = 1;
				insert_comp(res,num);
			}	
			c1 = c1->next;
		}
		cl1 = cl1->next;
	}
	comp_list_t * cl2 = acl2->head;
	for(int i = 0; i < nc2; i++){
		if(!om2[i]){
			cl2 = cl2->next;
			continue;
		}
		comp_t *c2 = cl2->head;
		while(c2!=NULL){
			unsigned short int num = c2->num;
			if(!map[num]){
				map[num] = 1;
				insert_comp(res,num);
			}
			c2 = c2->next;
		}
		cl2 = cl2->next;
	}
	free(map);
	return res;
}




array_comp_list_t * union_array_comp_list(array_comp_list_t *acl1, array_comp_list_t *acl2, unsigned short int n){
	int s1 = acl1->size;
	if(!s1){
		return copy_array_comp_list(acl2);
	}
	int s2 = acl2->size;
	if(!s2){
		return copy_array_comp_list(acl1);
	}
	unsigned short int * dis_map1 = (unsigned short int *)calloc(s1,sizeof(unsigned short int));
	unsigned short int * dis_map2 = (unsigned short int *)calloc(s2,sizeof(unsigned short int));
	int tnc = s1+s2;
	char *overlap_map = (char *)calloc(tnc*tnc,sizeof(char));
	comp_list_t * cl1 = acl1->head;
	for(int i = 0; i < s1; i++){
		comp_list_t *cl2 = acl2->head;
		for(int j = 0; j < s2; j++){
			if(is_included(cl1,cl2,n)){
				dis_map2[j]++;
				//continue;
			}
			else if(is_included(cl2,cl1,n)){
				dis_map1[i]++;
				//continue;
			}
			else if(is_disjoint(cl1,cl2,n)){
				dis_map1[i]++;
				dis_map2[j]++;
				//continue;
			}
			else{
				overlap_map[tnc*i + (s1+j)] = 1;
				overlap_map[tnc*(s1+j) + i] = 1;
			}
			cl2 = cl2->next;
		}
		cl1 = cl1->next;
	}

	array_comp_list_t *res = create_array_comp_list();
	cl1 = acl1->head;
	for(int i = 0; i < s1; i++){
		if(dis_map1[i]==s2){
			//comp_list_t *dst = create_comp_list();
			comp_list_t *dst = copy_comp_list(cl1);
			insert_comp_list(res,dst);
		}
		cl1 = cl1->next;
	}
	
	comp_list_t *cl2 = acl2->head;
	for(int i = 0; i < s2; i++){
		if(dis_map2[i]==s1){
			//comp_list_t *dst = create_comp_list();
			comp_list_t *dst = copy_comp_list(cl2);
			insert_comp_list(res,dst);
		}
		cl2 = cl2->next;
	}
	
	array_comp_list_t * acl = extract_comps(overlap_map,tnc);
	
	comp_list_t *cl = acl->head;
	for(int i = 0; i < acl->size; i++){
		comp_t *c = cl->head; 
		int s = cl->size;
		unsigned short int * om1 = (unsigned short int *)calloc(s1,sizeof(unsigned short int));
		unsigned short int * om2 = (unsigned short int *)calloc(s2,sizeof(unsigned short int));
		for(int j = 0; j < s; j++){
			unsigned short int num = c->num;
			if(num>=s1){
				om2[num-s1] = 1;
			}
			else{
				om1[num] = 1;
			}
			c = c->next;
		}
		comp_list_t * res1 = comp_array_union_direct(acl1,acl2,om1,om2,n);
		
		//comp_list_t *dst = create_comp_list();
		//comp_list_t * dst = copy_comp_list(res1);
		free(om1);
		free(om2);
		insert_comp_list(res,res1);
		cl = cl->next;	
	}
	free(overlap_map);
	free(dis_map1);
	free(dis_map2);
	
	free_array_comp_list(acl);
	return res;
}

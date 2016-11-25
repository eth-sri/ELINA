/*
	Copyright 2015 Software Reliability Lab, ETH Zurich

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

comp_list_t * intersection_comp_list(comp_list_t *c1, comp_list_t *c2, unsigned short int n){
	unsigned short int *map = (unsigned short int *)calloc(n,sizeof(unsigned short int));
	comp_list_t *res = create_comp_list();
	int l = 0;
	comp_t *c = c1->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		map[num] = 1;
		c = c->next;
	}
	c = c2->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		if(map[num]){
			insert_comp(res,num);
		}
		c = c->next;
	}
	free(map);
	return res;
}



comp_list_t * compute_diff(char *map, comp_list_t * cl, unsigned short int n){
	comp_t * c = cl->head;
	comp_list_t * res = create_comp_list();
	int flag = 0;
	while(c!=NULL){
		unsigned short int num = c ->num;
		if(!map[num]){
			insert_comp(res,num);
		}	
		else{
			flag = 1;
		}
		c = c->next;
	}
	if(flag){
		return res;
	}
	else{
		free_comp_list(res);
		return create_comp_list();
	}
}

/***
 c1 - c2 as well as c2 - c1
**/
array_comp_list_t * intersection_comp_list_compute_diff_both(comp_list_t *c1, comp_list_t *c2, unsigned short int n){
	unsigned short int *map1 = (unsigned short int *)calloc(n,sizeof(unsigned short int));
	unsigned short int *map2 = (unsigned short int *)calloc(n,sizeof(unsigned short int));
        array_comp_list_t * res = create_array_comp_list();
	comp_list_t *cl1 = create_comp_list();
        comp_list_t *cl2 = create_comp_list();
	comp_list_t *cl3 = create_comp_list();
	int l = 0;
	comp_t *c = c2->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		map2[num] = 1;
		c = c->next;
	}
	c = c1->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		map1[num] = 1;
		c = c->next;
	}
	c = c1->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		if(map2[num]){
			insert_comp(cl1,num);
		}
		else{
			insert_comp(cl2,num);
		}
		c = c->next;	
	}
	c = c2->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		if(!map1[num]){
			insert_comp(cl3,num);
		}
		c = c->next;
	}
	free(map1);
	free(map2);
        insert_comp_list(res,cl1);
        insert_comp_list(res,cl2);
	insert_comp_list(res,cl3);
	return res;
}



/***
 c1 - c2
**/
array_comp_list_t * intersection_comp_list_compute_diff(comp_list_t *c1, comp_list_t *c2, unsigned short int n){
	unsigned short int *map = (unsigned short int *)calloc(n,sizeof(unsigned short int));
        array_comp_list_t * res = create_array_comp_list();
	comp_list_t *cl1 = create_comp_list();
        comp_list_t *cl2 = create_comp_list();
	int l = 0;
	comp_t *c = c2->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		map[num] = 1;
		c = c->next;
	}
	c = c1->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		if(map[num]){
			insert_comp(cl1,num);
		}
		else{
			insert_comp(cl2,num);
		}
		c = c->next;
	}
	free(map);
        insert_comp_list(res,cl1);
        insert_comp_list(res,cl2);
	return res;
}


array_comp_list_t * intersection_array_comp_list(array_comp_list_t *acl1, array_comp_list_t *acl2, unsigned short int n){
	int s1 = acl1->size;
	int s2 = acl2->size;
	array_comp_list_t * res = create_array_comp_list();
	//int l = 0;
	if(!s1 || !s2){
		return res;
	}
	comp_list_t * cl1 = acl1->head;
	while(cl1!=NULL){
		comp_list_t * cl2 = acl2->head;
		while(cl2!=NULL){
			comp_list_t * src = intersection_comp_list(cl1,cl2,n);
			if(src->size>0){
				insert_comp_list(res,src);
			}
			else{
				free_comp_list(src);
			}
			cl2 = cl2->next;
		}
		cl1 = cl1->next;
	}
	return res;
}




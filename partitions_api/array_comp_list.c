/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
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

void print_array_comp_list(array_comp_list_t *acl, unsigned short int n){
	if(!acl || !acl->size){
		return;
	}
	fprintf(stdout,"%d\n",acl->size);
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

void insert_comp_list_tail(array_comp_list_t *acl, comp_list_t *cl){
	comp_list_t *h = acl->head;
	if(!h){
		acl->head = cl;
		acl->head->next=NULL;
		//acl->tail = cl;
		acl->size++;
		return;
	}
	while(h->next!=NULL){
		h = h->next;
	}
	h->next = cl;
	cl->next = NULL;
	acl->size++;
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
			comp_list_t * tmp = cl1;
			cl1 = cl1->next;
			remove_comp_list(acl,tmp);
		}
		else if(is_included(cl,cl1,n)){
			free(map);
			free_comp_list(cl);
			return;
		}
		
		else{
			union_comp_list(cl,cl1,map);
			comp_list_t *tmp = cl1;
			cl1 = cl1->next;
			remove_comp_list(acl,tmp);	
		}
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

short int find_index(array_comp_list_t *acl,comp_list_t * src){
	unsigned short int num_comp = acl->size;
	comp_list_t *cl = acl->head;
	
	for(int l = 0; l < num_comp; l++){
		if(cl == src){
			return l;
		}
		cl = cl->next;
	}
	return -1;
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
	
	//char *map2 = (char *)calloc(n,sizeof(char));
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

unsigned short int * create_array_map(array_comp_list_t * acl, unsigned short int n){
	unsigned short int * map = (unsigned short int *)calloc(n,sizeof(unsigned short int));
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

short int is_comp_list_included(array_comp_list_t *acl, comp_list_t *clr, unsigned short int n){
	comp_list_t * cl = acl->head;
	short int l = 0;
	while(cl!=NULL){
		if(is_included(clr,cl,n)){
			return l;
		}
		l++;
		cl = cl->next;
	}
	return -1;
}

short int is_covered(array_comp_list_t *acl, comp_list_t *clr, unsigned short int n){
	char *map = (char *)calloc(n,sizeof(unsigned short int));
	comp_list_t *cl = acl->head;
	while(cl!=NULL){
		comp_t * c = cl->head;
		while(c != NULL){
			unsigned short int num = c->num;
			map[num] = 1;
			c = c->next;
		}
		cl = cl->next;
	}
	
	comp_t * c = clr->head;
	while(c != NULL){
		unsigned short int num = c->num;
		if(!map[num]){
			free(map);
			return 0;
		}
		c = c->next;	
	}
	return 1;
}


char * create_intersection_map(array_comp_list_t * acl1, array_comp_list_t * acl2, unsigned short int num){
	char * map1 = (char *)calloc(num,sizeof(char));
	char * map2 = (char *)calloc(num,sizeof(char));
	char * map = (char *)calloc(num,sizeof(char));
	comp_list_t * cl1 = acl1->head;
	while(cl1!=NULL){
		comp_t * c1 = cl1->head;
		while(c1 != NULL){
			unsigned short int num1 = c1->num;
			map1[num1] = 1;
			c1 = c1->next;
		}
		cl1 = cl1->next;
	} 
	comp_list_t * cl2 = acl2->head;
	while(cl2!=NULL){
		comp_t * c2 = cl2->head;
		while(c2 != NULL){
			unsigned short int num2 = c2->num;
			map2[num2] = 1;
			c2 = c2->next;
		}
		cl2 = cl2->next;
	}
	unsigned short int i;
	for(i=0; i < num; i++){
		map[i] = map1[i] & map2[i];
	}
	free(map1);
	free(map2);
	return map;
}

size_t array_comp_list_serialize_common(void* dst, array_comp_list_t* acl, int dry_run){
  size_t idx = 0;

  if(!dry_run){
    *(unsigned short int*)(dst + idx) = acl->size;
  }
  idx += sizeof(unsigned short int);

  idx += comp_list_serialize_common(dst + idx, acl->head, acl->size, dry_run);
  return idx;
}

array_comp_list_t* array_comp_list_deserialize(void* p, size_t* size){
  array_comp_list_t* acl = (array_comp_list_t*)malloc(sizeof(array_comp_list_t));
  size_t idx = 0;

  acl->size = *(unsigned short int*)(p + idx);
  idx += sizeof(unsigned short int);

  size_t head_size;
  acl->head = comp_list_deserialize(p + idx, acl->size, &head_size);
  idx += head_size;

  *size = idx;
  return acl;
}

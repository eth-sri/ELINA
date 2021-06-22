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


#include "comp_list.h"

void fprint_comp_list(FILE* stream, comp_list_t *cl, unsigned short int n){
	if(!cl || !cl->size){
		return;
	}
	unsigned short int * ca = to_sorted_array(cl,n);
	unsigned short int comp_size = cl->size;
	for(unsigned short int i = 0; i < comp_size; i++){
		fprintf(stream,"%d ",ca[i]);
	}
	free(ca);
	fprintf(stream,"\n");
	fflush(stream);
}

void print_comp_list(comp_list_t * cl, unsigned short int n){
	fprint_comp_list(stdout,cl,n);
}


comp_list_t * create_comp_list(){
	comp_list_t * cl = (comp_list_t *)malloc(sizeof(comp_list_t));
	cl->head = NULL;
	//cl->tail = NULL;
	cl->next = NULL;
	cl->size = 0;
	return cl;
}

comp_list_t * copy_comp_list(comp_list_t *src){
	if(!src){
		return NULL;
	}
	comp_list_t * dst = create_comp_list();
	comp_t *s = src->head;
	while(s!=NULL){
		unsigned short int num = s->num;
		insert_comp(dst,num);
		s = s->next;
	}
	return dst;
}

void free_comp_list(comp_list_t *cl){
	if(cl==NULL){
		return;
	}
	comp_t * c = cl->head;
	while(c!=NULL){
		comp_t *temp = c;
		c = c->next;
		free(temp);
		temp = NULL;
	}
	//free(cl->head);
	//free(cl->tail);
	free(cl);
	//cl->head = NULL;
	cl= NULL;
	
}

void insert_comp(comp_list_t *cl, unsigned short int num){
	comp_t *h = cl->head;
	comp_t * c = (comp_t *)malloc(sizeof(comp_t));
	c->num = num;
	c->next = NULL;
	if(!h){
		cl->head = c;
		//cl->tail = c;
		cl->size++;
		return;
	}
	//cl->tail->next = c;
	//cl->tail = cl->tail->next;
	c->next = cl->head;
	cl->head = c;
	cl->size++;
	return;
}


int contains_comp(comp_list_t *cl, unsigned short int num){
	comp_t * c = cl->head;
	while(c!=NULL){
		if(c->num==num){
			return 1;
		}
		c = c->next;
	}
	return 0;
}

unsigned short int comp_list_size(comp_list_t *cl){
	//comp_t *h = cl->head;
	//comp_t *t = cl->tail;
	return cl->size;
}

void remove_comp(comp_list_t *cl, unsigned short int num){
	if(!cl){
		return;
	}
	comp_t * c = cl->head;
	if(!c){
		return;
	}
	if(cl->head->num==num){
		cl->head = cl->head->next;
		cl->size--;
		free(c);
		return;
	}
	if(cl->size==1){
		return;
	}
	while(c->next!=NULL){
		if(c->next->num==num){
			comp_t *temp = c->next;
			/*if(c->next==cl->tail){
				cl->tail = c;
				free(temp);
				c->next = NULL;
			}
			else{*/
				c->next = c->next->next;
				free(temp);
			//}
			cl->size--;
			return;
		}
		c = c->next;
	}
	return;
}

unsigned short int * to_sorted_array(comp_list_t * cl,unsigned short int n){
	unsigned short int comp_size = cl->size;
	char * map = (char *)calloc(n,sizeof(char));
	unsigned short int * res = (unsigned short int *)malloc(comp_size*sizeof(unsigned short int));
	comp_t *c = cl->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		map[num] = 1;
		c = c->next;
	}
	int l = 0;
	for(int i = 0; i < n; i++){
		if(map[i]){
			res[l] = i;
			l++;
		}
	}
	free(map);
	return res;
}

int is_equal_map(char *map1, char *map2, unsigned short int n){
	for(int i = 0; i < n; i++){
		if(map1[i]!=map2[i]){
			return 0;
		}
	}
	return 1;
}

void create_comp_list_map(comp_list_t *cl, unsigned short int n,char *map){
	if(cl==NULL){
		return;
	}
	comp_t * c = cl->head;
	while(c!=NULL){
		unsigned short int num = c->num;
		map[num] = 1;
		c = c->next;
	}
	
} 

unsigned short int * map_index(unsigned short int * dst, unsigned short int * src, unsigned short int comp_size){
	unsigned short int * res = (unsigned short int *)calloc(comp_size, sizeof(unsigned short int));
	unsigned short int i, j = 0;
	for(i = 0; i < comp_size; i++){
		unsigned short int num = dst[i];
		while(src[j] != num){
			j++;
		}
		res[i] = j;
		//printf("i: %d j: %d num: %d\n",i,j,num);
		//fflush(stdout);
		j++;
	}
	return res;
}

size_t comp_serialize_common(void* dst, comp_t* head, unsigned short int length, int dry_run){
  size_t idx = 0;
  comp_t* c = head;
  for(unsigned short int i = 0; i < length; i++) {
    if(!dry_run){
      *(unsigned short int*)(dst + idx) = c->num;
    }
    idx += sizeof(unsigned short int);
    c = c->next;
  }
  return idx;
}

size_t comp_list_serialize_common(void* dst, comp_list_t* head, unsigned short int length, int dry_run){
  size_t idx = 0;
  comp_list_t* cl = head;
  for(unsigned short int i = 0; i < length; i++){
    if(!dry_run){
      *(unsigned short int*)(dst + idx) = cl->size;
    }
    idx += sizeof(unsigned short int);
    idx += comp_serialize_common(dst + idx, cl->head, cl->size, dry_run);
    cl = cl->next;
  }
  return idx;
}

comp_t* comp_deserialize(void* p, unsigned short int length, size_t* size){
  if(length == 0){
    *size = 0;
    return NULL;
  }else{
    comp_t** buf = (comp_t**)malloc(sizeof(comp_t*) * length);
    size_t idx = 0;

    for(unsigned short int i = 0; i < length; i++){
      buf[i] = (comp_t*)malloc(sizeof(comp_t) * length);

      buf[i]->num = *(unsigned short int*)(p + idx);
      idx += sizeof(unsigned short int);
    }
    for(unsigned short int i = 0; i < length - 1; i++){
      buf[i]->next = buf[i + 1];
    }
    buf[length - 1]->next = NULL;

    *size = idx;
    comp_t* head = buf[0];
    free(buf);
    return head;
  }
}

comp_list_t* comp_list_deserialize(void* p, unsigned short int length, size_t* size){
  if(length == 0){
    *size = 0;
    return NULL;
  }else{
    comp_list_t** buf = (comp_list_t**)malloc(sizeof(comp_list_t*) * length);
    size_t idx = 0;

    for(size_t i = 0; i < length; i++){
      buf[i] = (comp_list_t*)malloc(sizeof(comp_list_t));

      buf[i]->size = *(unsigned short int*)(p + idx);
      idx += sizeof(unsigned short int);

      size_t comp_size;
      buf[i]->head = comp_deserialize(p + idx, buf[i]->size, &comp_size);
      idx += comp_size;
    }
    for(size_t i = 0; i < length - 1; i++){
      buf[i]->next = buf[i + 1];
    }
    buf[length - 1]->next = NULL;

    *size = idx;
    comp_list_t* head = buf[0];
    free(buf);
    return head;
  }
}

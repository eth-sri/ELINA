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




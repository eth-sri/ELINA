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

#ifndef __COMP_LIST_H
#define __COMP_LIST_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct comp_t{
	unsigned short int num;
	struct comp_t * next;
}comp_t;

typedef struct comp_list_t{
	comp_t *head;
	//comp_t *tail;
	struct comp_list_t * next;
	unsigned short int size;
}comp_list_t;

typedef struct array_comp_list_t{
	comp_list_t *head;
	//comp_list_t *tail;
	unsigned short int size;	
}array_comp_list_t; 

/****
Basic Linked List Insert, Delete and Find functions for Component List
*****/
comp_list_t * create_comp_list(void);
comp_list_t * copy_comp_list(comp_list_t *src);
void free_comp_list(comp_list_t *cl);
unsigned short int comp_list_size(comp_list_t *cl);
void insert_comp(comp_list_t *cl, unsigned short int num);
int contains_comp(comp_list_t *cl, unsigned short int num);
void remove_comp(comp_list_t *cl, unsigned short int num);
void fprint_comp_list(FILE* stream, comp_list_t *cl, unsigned short int n);
void print_comp_list(comp_list_t *cl, unsigned short int n);
unsigned short int * to_sorted_array(comp_list_t * cl,unsigned short int n);
unsigned short int * map_index(unsigned short int * dst, unsigned short int * src, unsigned short int comp_size);

/*****
Set Operations on two Component Lists
*****/

comp_list_t * intersection_comp_list(comp_list_t *c1, comp_list_t *c2, unsigned short int n);
int is_disjoint(comp_list_t * cl1, comp_list_t *cl2, unsigned short int n);
int is_disjoint_with_map(comp_list_t *cl, char *map);
short int is_map_disjoint(unsigned short int * cmap, char *map2, unsigned short int n);
int is_included(comp_list_t *cl1, comp_list_t *cl2, unsigned short int n);
void union_comp_list(comp_list_t * cl1,comp_list_t * cl2, char *map);
void union_comp_list_direct(comp_list_t *cl1, comp_list_t *cl2, unsigned short int n);
void unite_comp_lists(comp_list_t *cl1,comp_list_t *cl2,char *map,int i, int j,unsigned short int n);
char * create_map(comp_list_t *cl, unsigned short int n);
int is_equal_map(char *map1, char *map2, unsigned short int n);

/****
Basic Linked List Insert, Delete and Find functions for List of Component List
*****/

array_comp_list_t * create_array_comp_list(void);
array_comp_list_t * copy_array_comp_list(array_comp_list_t *src);
void free_array_comp_list(array_comp_list_t * acl);
void insert_comp_list(array_comp_list_t *acl, comp_list_t * cl);
void insert_comp_list_tail(array_comp_list_t *acl, comp_list_t * cl);
void insert_comp_list_with_union(array_comp_list_t * acl, comp_list_t *cl,unsigned short int n);
comp_list_t * find(array_comp_list_t *acl,unsigned short int num);
short int find_index(array_comp_list_t *acl,comp_list_t *cl);
void remove_comp_list(array_comp_list_t *acl, comp_list_t *cl);
void print_array_comp_list(array_comp_list_t *acl,unsigned short int n);
void clear_array_comp_list(array_comp_list_t *acl);
unsigned short int * create_array_map(array_comp_list_t * acl, unsigned short int n);
comp_list_t * compute_diff(char *map, comp_list_t * cl, unsigned short int n);

/*****
Intersection and Union of two Component Lists
*****/

array_comp_list_t * intersection_array_comp_list(array_comp_list_t *acl1, array_comp_list_t *acl2,unsigned short int n);
array_comp_list_t * intersection_comp_list_compute_diff(comp_list_t *c1, comp_list_t *c2, unsigned short int n);
array_comp_list_t * intersection_comp_list_compute_diff_both(comp_list_t *c1, comp_list_t *c2, unsigned short int n);
comp_list_t * comp_array_union_direct(array_comp_list_t * acl1, array_comp_list_t *acl2, unsigned short int *om1, unsigned short int *om2,  unsigned short int n);
array_comp_list_t * union_array_comp_list(array_comp_list_t *acl1, array_comp_list_t *acl2, unsigned short int n);
int is_equal_array_comp_list(array_comp_list_t *acl1, array_comp_list_t *acl2, unsigned short int n);
int is_lequal_array_comp_list(array_comp_list_t * acl1, array_comp_list_t * acl2, unsigned short int n);
int is_connected(array_comp_list_t *acl, unsigned short int i, unsigned short int j);
short int is_comp_list_included(array_comp_list_t *acl, comp_list_t *cl, unsigned short int n);
short int is_covered(array_comp_list_t *acl, comp_list_t *clr, unsigned short int n);
char * create_intersection_map(array_comp_list_t * acl1, array_comp_list_t * acl2, unsigned short int num);
/***
Extracting Components
**/
array_comp_list_t * extract(double *m, unsigned short int n);
array_comp_list_t * extract_comps(char *m, unsigned short int n);

/*****
Serialization
*****/
size_t comp_list_serialize_common(void* dst, comp_list_t* cl, unsigned short int length, int dry_run);
size_t array_comp_list_serialize_common(void* dst, array_comp_list_t* acl, int dry_run);

comp_list_t* comp_list_deserialize(void* p, unsigned short int length, size_t* size);
array_comp_list_t* array_comp_list_deserialize(void* p, size_t* size);

/*static int matpos(int i, int j){
	return j + ((i+1)*(i+1))/2;
}

static int matpos2(int i , int j){
	 if (j>i) return matpos(j^1,i^1);
  	 else return matpos(i,j);
}*/

#endif

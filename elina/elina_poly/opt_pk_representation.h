/* ********************************************************************** */
/* pk_representation.h: General management of polyhedra  */
/* ********************************************************************** */

/* This file is part of the APRON Library, released under LGPL license.  Please
   read the COPYING file packaged in the distribution */

#ifndef _OPT_PK_REPRESENTATION_H_
#define _OPT_PK_REPRESENTATION_H_

#include "opt_pk.h"
#include "opt_pk_matrix.h"
#include "opt_pk_cherni.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ********************************************************************** */
/* I. Memory */
/* ********************************************************************** */

/* Allocates a polyedron and fills its structure with null values, which
   corresponds to a bottom polyhedron. */
opt_pk_array_t * opt_pk_array_alloc(opt_pk_t ** poly, array_comp_list_t *acl, unsigned short int maxcols);

opt_pk_t* opt_poly_alloc(size_t intdim, size_t realdim);

/* Free all the members of the polyhedron structure (GMP semantics) */
void opt_poly_clear(opt_pk_t* op);

void opt_poly_array_clear(opt_pk_internal_t *opk, opt_pk_array_t * op);


void opt_poly_set(opt_pk_t* oa, opt_pk_t* ob);

void opt_pk_fprint(FILE* stream, elina_manager_t *man, opt_pk_t* op,
	       char** name_of_dim);
/* ********************************************************************** */
/* II. Control of internal representation */
/* ********************************************************************** */

/* Minimization function, in the sense of minimized dual representation This
   function minimizes if not already done the given polyhedron.
   Raise an exception, but still transmit it (pk->exn not reseted).
*/
void opt_poly_chernikova(elina_manager_t* man, opt_pk_t* poly, char* msg);

void opt_pk_convert(elina_manager_t *man, opt_pk_array_t * op, char *msg);
/* Same as poly_chernikova, but in addition ensure normalized epsilon
   constraints. */
void opt_poly_minimize(elina_manager_t* man, opt_pk_t* op);

/* Same as poly_chernikova2, but in addition normalize matrices by Gauss
   elimination and sorting */
//void poly_chernikova3(elina_manager_t* man, pk_t* poly, char* msg);

/* Put the polyhedron with minimized constraints and frames.  If in addition
   the integer man->option->canonicalize.algorithm is strictly positive,
   normalize equalities and lines, and also strict constraints */
//void poly_canonicalize(elina_manager_t* man, pk_t* poly);

/* Has the polyhedron normalized strict constraints ? */
//static inline bool poly_is_conseps(pk_internal_t* pk, pk_t* po);

/* Is the polyhedron under (strong) normal form ? */
//bool pk_is_canonical(elina_manager_t* man, pk_t* po);

/* Make available the matrix of constraints (resp. frames). The matrix will
   remain unavailable iff the polyhedron appears to be empty */
//static inline void poly_obtain_C(elina_manager_t* man, pk_t* po, char* msg);
static inline void opt_poly_obtain_F(elina_manager_t* man, opt_pk_t* po, char* msg);

/* Assuming the the matrix of constraints (resp. frames) is available, sort it,
   and take care of the saturation matrices. */
void opt_poly_obtain_sorted_C(opt_pk_internal_t* opk, opt_pk_t* op);
//void poly_obtain_sorted_F(pk_internal_t* pk, pk_t* poly);

/* Assuming one of the saturation matrix is available, make satC (resp. satF)
   available. */
static inline void opt_poly_obtain_satC(opt_pk_t* poly);
static inline void opt_poly_obtain_satF(opt_pk_t* poly);


static inline void opt_poly_obtain_satC(opt_pk_t* poly)
{
  if (!poly->satC){
    assert(poly->F && poly->satF);
    poly->satC = opt_satmat_transpose(poly->satF,poly->F->nbrows);
  }
}

static inline void opt_poly_obtain_satF(opt_pk_t* poly)
{
  if (!poly->satF){
    assert(poly->C && poly->satC);
    poly->satF = opt_satmat_transpose(poly->satC,poly->C->nbrows);
  }
}

//swap two components
void poly_swap(opt_pk_t *poly1, opt_pk_t *poly2);

/* Exchange C and F, sat C and satF, nbeq and nbline */
static inline void opt_poly_dual(opt_pk_t* poly)
{
  void* ptr;
  size_t nb;
  ptr = poly->C; poly->C = poly->F; poly->F = ptr;
  ptr = poly->satC; poly->satC = poly->satF; poly->satF = ptr;
  nb = poly->nbeq; poly->nbeq = poly->nbline; poly->nbline = nb;
}

static inline void opt_poly_array_dual(opt_pk_array_t* op)
{
  array_comp_list_t * acl = op->acl;
  if(op->is_bottom || !acl){
	return;
  }
  unsigned short int num_comp = acl->size;
  unsigned short int k;
  opt_pk_t ** poly = op->poly;
  for(k=0; k < num_comp; k++){
	opt_poly_dual(poly[k]);
  }
}

char meet_cons_one_comp(opt_pk_internal_t *opk, opt_pk_t **poly_a, array_comp_list_t *acla, 
			unsigned short int **ca_arr, opt_pk_t *poly, unsigned short int *ca,  char * map);

void meet_cons_with_map(opt_pk_internal_t *opk, opt_pk_array_t *oa, opt_pk_t **poly, unsigned short int *rmapa, 
		unsigned short int **ca_arr, size_t *counterC,  char * map, char * exclusion_map);

void meet_cons(opt_pk_internal_t *opk, opt_pk_array_t *oa, opt_pk_t **poly, unsigned short int *rmapa, 
		unsigned short int **ca_arr, size_t *counterC, char * map);



size_t cartesian_product_vertices_one_comp(opt_pk_t **poly_a, array_comp_list_t *acl,
					 unsigned short int ** ca_arr, size_t * num_vertex_a,
					 opt_pk_t * poly, size_t counter, unsigned short int  *ca, char * map);

void cartesian_product_vertices_with_map(opt_pk_array_t *oa, opt_pk_t ** poly, 
				 unsigned short int *rmapa, unsigned short int ** ca_arr, 
				  size_t * num_vertex_a, size_t *num_vertex, size_t *counterF,
				  char * exclusion_map);


void cartesian_product_vertices(opt_pk_array_t *oa, opt_pk_t ** poly, 
				 unsigned short int *rmapa, unsigned short int ** ca_arr, 
				  size_t * num_vertex_a, size_t *num_vertex, size_t *counterF);

void meet_rays_one_comp(opt_pk_t **poly_a, array_comp_list_t *acla, unsigned short int **ca_arr,
	       		size_t * nblinemap, opt_pk_t *poly, size_t * num_vertex_a,
	        	unsigned short int  * ca, size_t start, char * map);

void meet_rays_with_map(opt_pk_array_t *oa, opt_pk_t **poly, unsigned short int *rmapa, 
	       unsigned short int **ca_arr, size_t * num_vertex_a, size_t * counterF,
		char *exclsuion_map);

void meet_rays(opt_pk_array_t *oa, opt_pk_t **poly, unsigned short int *rmapa, 
	       unsigned short int **ca_arr, size_t * num_vertex_a, size_t * counterF);

void combine_satmat(opt_pk_internal_t *opk, opt_pk_t *poly, unsigned short int comp_size, size_t end, bool con_to_ray);
/* Like poly_chernikova, etc, but if usual is false, take into account the fact that
   poly_dual has been applied */ 
//void poly_chernikova_dual(elina_manager_t* man,
//			  pk_t* poly,
//			  char* msg,
//			  bool usual);
//static inline void poly_obtain_C_dual(elina_manager_t* man,pk_t* po, char* msg, bool usual);
//static inline void poly_obtain_F_dual(elina_manager_t* man,pk_t* po, char* msg, bool usual);

static inline void opt_poly_obtain_C(elina_manager_t* man, opt_pk_t* op, char* msg)
{
  if (!op->C) opt_poly_chernikova(man,op,msg);
}

static inline void opt_poly_obtain_F(elina_manager_t* man, opt_pk_t* op, char* msg)
{
  if (!op->F) opt_poly_chernikova(man,op,msg);
}

/* ********************************************************************** */
/* V. Checking */
/* ********************************************************************** */

/* Checks gauss elimination */
//bool matrix_check_gauss(matrix_t* mat, size_t nbeq);

/* Perform rather detailed and costly checks on a polyhedron, to detect
   inconsistencies */
//bool poly_check(pk_internal_t* pk, pk_t* po);
/* Same as previous function, but usual==false mean that poly_dual has been
   applied */
//bool poly_check_dual(pk_internal_t* pk, pk_t* po, bool usual);

/* ********************************************************************** */
/* ********************************************************************** */
/* INLINE DEFINITIONS */
/* ********************************************************************** */

//static inline bool poly_is_conseps(pk_internal_t* pk, pk_t* po)
//{
//  return (!pk->strict || 
//	  po->realdim==0 || 
//	  (po->status & pk_status_conseps));;
//}
//  static inline bool poly_is_minimaleps(pk_internal_t* pk, pk_t* po)
//{
//  return (!pk->strict ||
//	  po->realdim==0 ||
//	  (po->status & pk_status_minimaleps));
//}


double abs_diff(double a, opt_numint_t b);

void quasi_removal(opt_pk_internal_t *opk, opt_pk_t * oc);


void replace_ineq_with_eq(opt_pk_internal_t *opk, opt_pk_t * op);

void opt_poly_copy(opt_pk_t *dst, opt_pk_t *src);

#ifdef __cplusplus
}
#endif

#endif

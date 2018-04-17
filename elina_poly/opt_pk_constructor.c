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

/* ********************************************************************** */
/* opt_pk_constructor.c: constructors and accessors */
/* ********************************************************************** */

#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"
#include "opt_pk_representation.h"
#include "opt_pk_constructor.h"
#include "opt_pk_user.h"



/* ********************************************************************** */
/* I. Constructors */
/* ********************************************************************** */

/* ====================================================================== */
/* Empty polyhedron */
/* ====================================================================== */


void opt_poly_set_bottom(opt_pk_internal_t* opk, opt_pk_array_t* op)
{
  opt_poly_array_clear(opk,op);
  op->is_bottom = true;
}

/*
The empty polyhedron is just defined by the absence of both
constraints matrix and frames matrix.
*/

opt_pk_array_t* opt_pk_bottom(elina_manager_t* man, size_t intdim, size_t realdim)
{
  #if defined(TIMING)
 	 	start_timing();
  #endif
  
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_BOTTOM);
  opt_pk_array_t* op = opt_pk_array_alloc(NULL,NULL,intdim+realdim+opk->dec);
  op->is_bottom = true;
  opt_pk_internal_realloc_lazy(opk,intdim+realdim);
  //op->status = opt_pk_status_conseps | opt_pk_status_minimaleps;
  man->result.flag_exact = man->result.flag_best = true;
  #if defined(TIMING)
 	 	record_timing(bottom_time);
  #endif
  return op;
}

/* ====================================================================== */
/* Universe polyhedron */
/* ====================================================================== */

void opt_matrix_fill_constraint_top(opt_pk_internal_t* opk, opt_matrix_t* oc, size_t start)
{
    opt_vector_clear(oc->p[start+0],oc->nbcolumns);
    oc->p[start+0][0] = 1;
    oc->p[start+0][opt_polka_cst] = 1;
}

void opt_poly_set_top(opt_pk_internal_t* opk, opt_pk_array_t* op)
{
  size_t i;
  size_t dim;
  opt_poly_array_clear(opk,op);
  array_comp_list_t * acl = create_array_comp_list();
  op->acl = acl;
  op->is_bottom = false;
  opt_pk_t ** poly = (opt_pk_t **)malloc(sizeof(opt_pk_t *));
  op->poly = poly;
}

opt_pk_array_t* opt_pk_top(elina_manager_t* man, size_t intdim, size_t realdim)
{
  
  #if defined(TIMING)
 	 	start_timing();
  #endif
  opt_pk_array_t* op;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_TOP);
  opt_pk_internal_realloc_lazy(opk,intdim+realdim);

  op = opt_pk_array_alloc(NULL,NULL,intdim + realdim+opk->dec);
  opt_poly_set_top(opk,op);
  man->result.flag_exact = man->result.flag_best = true;
  #if defined(TIMING)
 	 	record_timing(top_time);
  #endif
  return op;
}

/* ====================================================================== */
/* Hypercube polyhedron */
/* ====================================================================== */


/* The matrix is supposed to be big enough */
static 
int opt_matrix_fill_constraint_interval(opt_pk_internal_t* opk,
				opt_matrix_t* oc, size_t start,
				elina_interval_t* interval,
				size_t intdim, size_t realdim,
				bool integer)
{
  size_t k;
  bool ok;
  k = start;
  int res1 = elina_scalar_infty(interval->inf);
  int res2 = elina_scalar_infty(interval->sup);
  bool res3 = elina_scalar_equal(interval->inf,interval->sup);
  if (!res1 && !res2 && res3){
      ok = opt_vector_set_dim_bound(opk,oc->p[k],
				 0, interval->sup, 0,
				 intdim,realdim,
				 integer);
      if (!ok){
	return -1;
      }
      k++;
   }
   else {
      /* inferior bound */
      if (!res1){
	opt_vector_set_dim_bound(opk,oc->p[k],
			     0, interval->inf, -1,
			     intdim,realdim,
			     integer);
	k++;
      }
      /* superior bound */
      if (!res2){
	opt_vector_set_dim_bound(opk,oc->p[k],
			     0, interval->sup, 1,
			     intdim,realdim,
			     integer);
	k++;
      }
    }
  return (int)k;
}

/* Abstract an hypercube defined by the array of intervals of size
   intdim+realdim.  */

opt_pk_array_t* opt_pk_of_box(elina_manager_t* man,
		size_t intdim, size_t realdim,
		elina_interval_t** array)
{
  elina_dim_t k;
  int res;
  size_t dim;
  opt_pk_array_t* op;

  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_OF_BOX);
  opt_pk_internal_realloc_lazy(opk,intdim+realdim);

  dim = intdim + realdim;
  opt_pk_t ** poly = (opt_pk_t **)malloc(dim*sizeof(opt_pk_t *));
  array_comp_list_t *acl = create_array_comp_list();
  for(k=0; k < dim; k++){
	unsigned short int k1 = dim - 1 - k;
  	poly[k1] = opt_poly_alloc(1,0);
	poly[k1]->C = opt_matrix_alloc(3,3,false);
	poly[k1]->status = opt_pk_status_conseps;
	poly[k1]->C->p[0][0] = 1;
	poly[k1]->C->p[0][1] = 1;
	comp_list_t *cl = create_comp_list();
	insert_comp(cl,k+opk->dec);
	insert_comp_list(acl,cl);
	res = opt_matrix_fill_constraint_interval(opk,poly[k1]->C,opk->dec-1,array[k],1,0,true);
	if (res==-1){
		unsigned short int j;
    		for(j=0; j < dim; j++){
			if(poly[k1]->C){
				opt_matrix_free(poly[k1]->C);
			}
			free(poly[k1]);
		}
		free_array_comp_list(acl);
    		return opt_pk_bottom(man,intdim,realdim);
  	}
  	poly[k1]->C->nbrows = (size_t)res;
  }
  op = opt_pk_array_alloc(poly,acl,intdim+realdim+opk->dec);
  man->result.flag_exact = man->result.flag_best = true;
  return op;
}

/* ********************************************************************** */
/* II. Accessors */
/* ********************************************************************** */

/* Return the dimensions of the polyhedra */
elina_dimension_t opt_pk_dimension(elina_manager_t* man, opt_pk_array_t* op){
  elina_dimension_t res;
  //if(op->maxcols==1){
  //printf("maxcols: %d\n",op->maxcols);
  
  res.intdim = op->maxcols - 2;
  res.realdim = 0;
  return res;
}


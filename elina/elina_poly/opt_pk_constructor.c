/*
	Copyright 2016 Software Reliability Lab, ETH Zurich

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
int opt_matrix_fill_constraint_box(opt_pk_internal_t* opk,
				opt_matrix_t* oc, size_t start,
				elina_interval_t** box,
				size_t intdim, size_t realdim,
				bool integer)
{
  size_t k;
  elina_dim_t i;
  bool ok;
  itv_t itv;
  k = start;

  itv_init(itv);
  for (i=0; i<intdim+realdim; i++){
    itv_set_elina_interval(opk->itv,itv,box[i]);
    if (itv_is_point(opk->itv,itv)){
      ok = opt_vector_set_dim_bound(opk,oc->p[k],
				 (elina_dim_t)i, bound_numref(itv->sup), 0,
				 intdim,realdim,
				 integer);
      if (!ok){
	itv_clear(itv);
	return -1;
      }
      k++;
    }
    else {
      /* inferior bound */
      if (!bound_infty(itv->inf)){
	opt_vector_set_dim_bound(opk,oc->p[k],
			     (elina_dim_t)i, bound_numref(itv->inf), -1,
			     intdim,realdim,
			     integer);
	k++;
      }
      /* superior bound */
      if (!bound_infty(itv->sup)){
	opt_vector_set_dim_bound(opk,oc->p[k],
			     (elina_dim_t)i, bound_numref(itv->sup), 1,
			     intdim,realdim,
			     integer);
	k++;
      }
    }
  }
  itv_clear(itv);
  return (int)k;
}

/* Abstract an hypercube defined by the array of intervals of size
   intdim+realdim.  */

opt_pk_t* opt_pk_of_box(elina_manager_t* man,
		size_t intdim, size_t realdim,
		elina_interval_t** array)
{
  int k;
  size_t dim;
  opt_pk_t* op;

  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_OF_BOX);
  opt_pk_internal_realloc_lazy(opk,intdim+realdim);

  dim = intdim + realdim;
  op = opt_poly_alloc(intdim,realdim);
  op->status = opt_pk_status_conseps;

  dim = intdim + realdim;
  op->C = opt_matrix_alloc(opk->dec-1 + 2*dim, opk->dec + dim, false);

  /* constraints */
  opt_matrix_fill_constraint_top(opk,op->C,0);
  k = opt_matrix_fill_constraint_box(opk,op->C,opk->dec-1,array,intdim,realdim,true);
  if (k==-1){
    opt_matrix_free(op->C);
    op->C = NULL;
    return op;
  }
  op->C->nbrows = (size_t)k;
  
  //assert(poly_check(pk,po));
  man->result.flag_exact = man->result.flag_best = true;
  return op;
}

/* ********************************************************************** */
/* II. Accessors */
/* ********************************************************************** */

/* Return the dimensions of the polyhedra */
elina_dimension_t opt_pk_dimension(elina_manager_t* man, opt_pk_array_t* op){
  elina_dimension_t res;
  res.intdim = op->maxcols - 2;
  res.realdim = 0;
  return res;
}


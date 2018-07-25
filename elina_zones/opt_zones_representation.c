/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
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


#include "opt_zones_internal.h"
#include "opt_zones.h"
#include "opt_mat.h"

// Internal allocation function
opt_zones_t * opt_zones_alloc_internal(opt_zones_internal_t *pr, unsigned short int dim, unsigned short int intdim){
	opt_zones_t *o = (opt_zones_t *)malloc(sizeof(opt_zones_t));
	o->intdim = intdim;
	o->dim = dim;
	o->closed = NULL;
	o->m = NULL;
	return o;
}


// Internal copy function
opt_zones_t * opt_zones_copy_internal(opt_zones_internal_t *pr, opt_zones_t *o){
	opt_zones_t *r = opt_zones_alloc_internal(pr,o->dim, o->intdim);
	r->m = opt_zones_mat_copy(o->m,o->dim);
	r->closed = opt_zones_mat_copy(o->closed,o->dim);
	return r;	
}

// Internal function to free up memory
void opt_zones_free_internal(opt_zones_internal_t *pr, opt_zones_t *o){
	if(o->m){
		opt_zones_mat_free(o->m);
	}
	if(o->closed){
		opt_zones_mat_free(o->closed);
	}
	o->m = NULL;
	o->closed = NULL;
	free(o);
}

/********* Allocate a new zones abstract element by setting the appropriate matrix **************/
opt_zones_t* opt_zones_set_mat(opt_zones_internal_t* pr, opt_zones_t* o, 
			       opt_zones_mat_t* m, opt_zones_mat_t* closed, 
			       bool destructive)
{
  opt_zones_t* r;
  if (destructive) {
    /* free non-aliased matrices */
    if (o->m && o->m!=m && o->m!=closed){
      opt_zones_mat_free(o->m);
      o->m = NULL;
    }
    if (o->closed && o->closed!=m && o->closed!=closed){
      opt_zones_mat_free(o->closed);
      o->closed = NULL;
    }
    r = o;
  }
  else {
    /* copy aliased matrices */
    r = opt_zones_alloc_internal(pr,o->dim,o->intdim);
    if (m && (o->m==m || o->closed==m))
    {    
	 m = opt_zones_mat_copy(m,o->dim);
    }
    if (closed && (o->m==closed || o->closed==closed)){
      
      closed = opt_zones_mat_copy(closed,o->dim);
    }
  }
  r->m = m;
  r->closed = closed;
  return r;
}


/******** Copy Operator *****************/
opt_zones_t* opt_zones_copy(elina_manager_t* man, opt_zones_t* o)
{
  
  
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_COPY,0);
  opt_zones_t *res = opt_zones_copy_internal(pr,o);
  
  return res;
}


/********* Free Operator ***************/
void opt_zones_free(elina_manager_t* man, opt_zones_t* o)
{
   
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_FREE,0);
  opt_zones_free_internal(pr,o);
  
}

/********* Allocate a bottom abstract element **************/
opt_zones_t* opt_zones_bottom(elina_manager_t* man, unsigned short int intdim, unsigned short int realdim)
{
  
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_BOTTOM,0);
  opt_zones_t* o = opt_zones_alloc_internal(pr,intdim+realdim,intdim);
  
  return o;
}

/********** Allocate a top abstract element ****************/
opt_zones_t* opt_zones_top(elina_manager_t* man, unsigned short int intdim, unsigned short int realdim)
{
  
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_TOP,0);
  opt_zones_t* r = opt_zones_alloc_internal(pr,intdim+realdim,intdim);
  //printf("top %p %d %d %d\n",r,r->dim,intdim, realdim);
  //fflush(stdout);
  r->closed = opt_zones_mat_alloc_top(r->dim);
  
  return r;
}

/********** Dimension Operator ************************/
elina_dimension_t opt_zones_dimension(elina_manager_t* man, opt_zones_t* o)
{
  
  elina_dimension_t r;
  r.intdim = o->intdim;
  r.realdim = o->dim-o->intdim;
   
  return r;
}


/* Measuring the size of a zone abstract element */
int opt_zones_size(elina_manager_t* man, opt_zones_t* o)
{
  if (!o->m) return 1;
  int size = (o->dim+1)*(o->dim+1);
  return size;
}


/******
 No need to compute closure if o->closed is not NULL
****/
void opt_zones_sparse_weak_closure(opt_zones_internal_t *pr,  opt_zones_t *o){
	if(o->closed || !o->m){
		return;
	}
	unsigned short int dim = o->dim;
	o->closed = opt_zones_mat_copy(o->m,dim);
	
	#if defined(TIMING)
		start_timing();
	#endif
	
	closure_comp_sparse(o->closed,dim);
	if(strengthening_intra_comp_zones(o->closed,dim)){
		opt_zones_mat_free(o->closed);
		opt_zones_mat_free(o->m);
		o->closed = NULL;
		o->m = NULL;
	}
	
        #if defined(TIMING)
		record_timing(zones_closure_time);
	#endif
}



void opt_zones_cache_closure(opt_zones_internal_t *pr, opt_zones_t *o){
	if(o->closed || !o->m){
		return;
	}
	o->closed = opt_zones_mat_copy(o->m,o->dim);
	if(opt_zones_mat_closure(o->closed,o->dim)){
		opt_zones_mat_free(o->closed);
		opt_zones_mat_free(o->m);
		o->closed = NULL;
		o->m = NULL;
	}
}


/*We throw an exception for minimize*/
void opt_zones_minimize(elina_manager_t* man, opt_zones_t* o)
{
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_MINIMIZE,0);
  elina_manager_raise_exception(man,ELINA_EXC_NOT_IMPLEMENTED,pr->funid,
			     "not implemented");
}


/* We throw an exception for canoncalize */
void opt_zones_canonicalize(elina_manager_t* man, opt_zones_t* o)
{
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_CANONICALIZE,0);
  elina_manager_raise_exception(man,ELINA_EXC_NOT_IMPLEMENTED,pr->funid,
			     "not implemented");
}

/* We compute hash  */
int opt_zones_hash(elina_manager_t* man, opt_zones_t* o)
{
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_HASH,0);
  if (pr->funopt->algorithm>=0) opt_zones_cache_closure(pr,o);
  if (o->closed || o->m) {
    int r = 0;
    opt_zones_mat_t *oo = o->closed ? o->closed : o->m;
    double *m = oo->mat;
    unsigned short int i,j;
    for (i=0;i<o->dim+1;i++)
      for (j=0;j<o->dim+1;j++,m++)
	r = r*37 + *m;
    return r;
  }
  return 0;
}

/* We throw an exception just like APRON */
void opt_zones_approximate(elina_manager_t* man, opt_zones_t* o, int algorithm)
{
  opt_zones_internal_t* pr = opt_zones_init_from_manager(man,ELINA_FUNID_APPROXIMATE,0);
  elina_manager_raise_exception(man,ELINA_EXC_NOT_IMPLEMENTED,pr->funid,
			     "not implemented");
}
/****

Topological closure
****/
opt_zones_t* opt_zones_closure(elina_manager_t *man, bool destructive, opt_zones_t *o){
	opt_zones_internal_t *pr = (opt_zones_internal_t *)opt_zones_init_from_manager(man, ELINA_FUNID_CLOSURE,0);
	if(destructive)return o;
	return opt_zones_copy_internal(pr,o);
}


/*****
Print Timing Information
*****/

void opt_zones_fprint(FILE* stream, elina_manager_t* man, opt_zones_t * o,char** name_of_dim){
	#if defined(TIMING)
		fprintf(stdout,"Times are in CPU Cycles\n");
		fprintf(stdout,"Top: %g\n",zones_top_time);
		fprintf(stdout,"Free: %g\n",zones_free_time);
		fprintf(stdout,"Copy: %g\n",zones_copy_time);
		fprintf(stdout,"Closure: %g\n",zones_closure_time);
		fprintf(stdout,"Is_Equal: %g\n",zones_is_equal_time);
		fprintf(stdout,"Is_Lequal: %g\n",zones_is_lequal_time);
		fprintf(stdout,"Meet_Abstract: %g\n",zones_meet_time);
		fprintf(stdout,"Join: %g\n",zones_join_time);
		fprintf(stdout,"Widening: %g\n",zones_widening_time);
		fprintf(stdout,"Add_dimension: %g\n",zones_add_dimension_time);
		fprintf(stdout,"Permute_dimension: %g\n",zones_permute_dimension_time);
		fprintf(stdout,"Meet_Lincons_Array: %g\n",zones_meet_lincons_time);
		fprintf(stdout,"Forget_Array %g\n",zones_forget_array_time);
		fprintf(stdout,"Oct_to_Box: %g\n",zones_to_box_time);
		fprintf(stdout,"Alloc: %g\n",zones_alloc_time);
		fprintf(stdout,"Is_Top: %g\n",zones_is_top_time);
		fprintf(stdout,"Expand: %g\n",zones_expand_time);
		fprintf(stdout,"Fold: %g\n",zones_fold_time);
		fprintf(stdout,"Sat_Lincons: %g\n",zones_sat_lincons_time);
		fprintf(stdout,"Assign Linexpr: %g\n",zones_assign_linexpr_time);
        	fprintf(stdout,"Narrowing Time: %g\n",zones_narrowing_time);
		fprintf(stdout,"Is Unconstrained time: %g\n",zones_is_unconstrained_time);
		double total_time = zones_top_time + zones_free_time + zones_copy_time + zones_closure_time + zones_is_equal_time + 
				    zones_is_lequal_time + zones_meet_time + zones_join_time + zones_widening_time + zones_add_dimension_time
				    + zones_permute_dimension_time + zones_meet_lincons_time + zones_forget_array_time + zones_to_box_time + zones_alloc_time + 
				    zones_is_top_time + zones_expand_time + zones_fold_time + zones_sat_lincons_time + zones_assign_linexpr_time + zones_narrowing_time + zones_is_unconstrained_time;
		fprintf(stdout,"Total Zones Analysis: %g\n",total_time);
		fflush(stdout);
	#endif
}


/*****
Manager specific
***/

void opt_zones_internal_free(opt_zones_internal_t *pr){
	free(pr->tmp);
	free(pr->tmp2);
	pr->tmp = NULL;
	pr->tmp2 = NULL;
	free(pr);
}


elina_manager_t* opt_zones_manager_alloc(void)
{
  size_t i;
  elina_manager_t* man;
  opt_zones_internal_t* pr;

  if (!elina_fpu_init()) {
   fprintf(stderr,"opt_zones_manager_alloc cannot change the FPU rounding mode\n");
  }

  pr = (opt_zones_internal_t*)malloc(sizeof(opt_zones_internal_t));
  assert(pr);
  pr->tmp_size = 10;
  pr->tmp = (double *)calloc(pr->tmp_size,sizeof(double));
  assert(pr->tmp);
  init_array_zones(pr->tmp,pr->tmp_size);
  pr->tmp2 = calloc(pr->tmp_size,sizeof(long));
  assert(pr->tmp2);
  
  man = elina_manager_alloc("opt_zones","1.1 with double", pr,
			 (void (*)(void*))opt_zones_internal_free);

  pr->man = man;

  man->funptr[ELINA_FUNID_COPY] = &opt_zones_copy;
  man->funptr[ELINA_FUNID_FREE] = &opt_zones_free;
  man->funptr[ELINA_FUNID_ASIZE] = &opt_zones_size;
  man->funptr[ELINA_FUNID_MINIMIZE] = &opt_zones_minimize;
  man->funptr[ELINA_FUNID_CANONICALIZE] = &opt_zones_canonicalize;
  man->funptr[ELINA_FUNID_HASH] = &opt_zones_hash;
  man->funptr[ELINA_FUNID_APPROXIMATE] = &opt_zones_approximate;
  man->funptr[ELINA_FUNID_FPRINT] = &opt_zones_fprint;
  //man->funptr[ELINA_FUNID_FPRINTDIFF] = &opt_zones_fprintdiff;
  //man->funptr[ELINA_FUNID_FDUMP] = &opt_zones_fdump;
  //man->funptr[ELINA_FUNID_SERIALIZE_RAW] = &opt_zones_serialize_raw;
  //man->funptr[ELINA_FUNID_DESERIALIZE_RAW] = &opt_zones_deserialize_raw;
  man->funptr[ELINA_FUNID_BOTTOM] = &opt_zones_bottom;
  man->funptr[ELINA_FUNID_TOP] = &opt_zones_top;
  //man->funptr[ELINA_FUNID_OF_BOX] = &opt_zones_of_box;
  man->funptr[ELINA_FUNID_DIMENSION] = &opt_zones_dimension;
  man->funptr[ELINA_FUNID_IS_BOTTOM] = &opt_zones_is_bottom;
  man->funptr[ELINA_FUNID_IS_TOP] = &opt_zones_is_top;
  man->funptr[ELINA_FUNID_IS_LEQ] = &opt_zones_is_leq;
  man->funptr[ELINA_FUNID_IS_EQ] = &opt_zones_is_eq;
  man->funptr[ELINA_FUNID_IS_DIMENSION_UNCONSTRAINED] = &opt_zones_is_dimension_unconstrained;
  //man->funptr[ELINA_FUNID_SAT_INTERVAL] = &opt_zones_sat_interval;
  //man->funptr[ELINA_FUNID_SAT_LINCONS] = &opt_zones_sat_lincons_timing;
  //man->funptr[ELINA_FUNID_SAT_TCONS] = &opt_zones_sat_tcons;
  man->funptr[ELINA_FUNID_BOUND_DIMENSION] = &opt_zones_bound_dimension;
  //man->funptr[ELINA_FUNID_BOUND_LINEXPR] = &opt_zones_bound_linexpr;
  //man->funptr[ELINA_FUNID_BOUND_TEXPR] = &opt_zones_bound_texpr;
  man->funptr[ELINA_FUNID_TO_BOX] = &opt_zones_to_box;
  man->funptr[ELINA_FUNID_TO_LINCONS_ARRAY] = &opt_zones_to_lincons_array;
  man->funptr[ELINA_FUNID_TO_TCONS_ARRAY] = &opt_zones_to_tcons_array;
  //man->funptr[ELINA_FUNID_TO_GENERATOR_ARRAY] = &opt_zones_to_generator_array;
  man->funptr[ELINA_FUNID_MEET] = &opt_zones_meet;
  //man->funptr[ELINA_FUNID_MEET_ARRAY] = &opt_zones_meet_array;
  man->funptr[ELINA_FUNID_MEET_LINCONS_ARRAY] = &opt_zones_meet_lincons_array;
  man->funptr[ELINA_FUNID_MEET_TCONS_ARRAY] = &opt_zones_meet_tcons_array;
  man->funptr[ELINA_FUNID_JOIN] = &opt_zones_join;
  //man->funptr[ELINA_FUNID_JOIN_ARRAY] = &opt_zones_join_array;
  //man->funptr[ELINA_FUNID_ADD_RAY_ARRAY] = &opt_zones_add_ray_array;
  //man->funptr[ELINA_FUNID_ASSIGN_LINEXPR_ARRAY] = &opt_zones_assign_linexpr_array;
  //man->funptr[ELINA_FUNID_SUBSTITUTE_LINEXPR_ARRAY] = &opt_zones_substitute_linexpr_array;
  man->funptr[ELINA_FUNID_ASSIGN_TEXPR_ARRAY] = &opt_zones_assign_texpr_array;
  //man->funptr[ELINA_FUNID_SUBSTITUTE_TEXPR_ARRAY] = &opt_zones_substitute_texpr_array;
  man->funptr[ELINA_FUNID_ADD_DIMENSIONS] = &opt_zones_add_dimensions;
  man->funptr[ELINA_FUNID_REMOVE_DIMENSIONS] = &opt_zones_remove_dimensions;
  man->funptr[ELINA_FUNID_PERMUTE_DIMENSIONS] = &opt_zones_permute_dimensions;
  man->funptr[ELINA_FUNID_FORGET_ARRAY] = &opt_zones_forget_array;
  man->funptr[ELINA_FUNID_EXPAND] = &opt_zones_expand;
  //man->funptr[ELINA_FUNID_FOLD] = &opt_zones_fold;
  man->funptr[ELINA_FUNID_WIDENING] = &opt_zones_widening;
  man->funptr[ELINA_FUNID_CLOSURE] = &opt_zones_closure;
    
	
  for (i=0;i<ELINA_EXC_SIZE;i++)
    elina_manager_set_abort_if_exception(man,i,false);
  
  return man;
}


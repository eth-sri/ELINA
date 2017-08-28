/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
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

#include "opt_oct_hmat.h"

opt_oct_t* opt_oct_assign_linexpr(elina_manager_t* man,
			  bool destructive, opt_oct_t* o,
			  elina_dim_t d, elina_linexpr0_t* lexpr,
			  opt_oct_t* dest)
{
  opt_oct_internal_t* pr =
    opt_oct_init_from_manager(man,ELINA_FUNID_ASSIGN_LINEXPR_ARRAY,2*(o->dim+1+5));
  opt_uexpr u = opt_oct_uexpr_of_linexpr(pr,pr->tmp,lexpr,o->intdim,o->dim);
  opt_oct_mat_t* src;
  bool respect_closure;
  if (d >= o->dim) {
    return NULL;
  }

  if (dest && !dest->closed && !dest->m)
    /* definitively empty due to dest*/
    return opt_oct_set_mat(pr,o,NULL,NULL,destructive);

  if (u.type==OPT_EMPTY)
    /* definitively empty due to empty expression */
    return opt_oct_set_mat(pr,o,NULL,NULL,destructive);

  /* useful to close only for non-invertible assignments */
  if ((u.type!=OPT_UNARY || u.i!=d) && pr->funopt->algorithm>=0)
    opt_oct_cache_closure(pr,o);
  src = o->closed ? o->closed : o->m;
  if (!src) return opt_oct_set_mat(pr,o,NULL,NULL,destructive); /* empty */

  /* can / should we try to respect the closure */
  respect_closure = (src==o->closed) && (pr->funopt->algorithm>=0) && (!dest);

  if (!destructive) src = opt_hmat_copy(src,o->dim);

  /* go */
  #if defined(TIMING)
  	start_timing();
  #endif
  opt_hmat_assign(pr,u,src,o->dim,d,&respect_closure);
  #if defined(TIMING)
  	record_timing(assign_linexpr_time);
  #endif

  /* exact on Q if zeroary or unary, closed arg and no conv error */
  if (u.type==OPT_BINARY || u.type==OPT_OTHER) flag_incomplete;
  else if (num_incomplete || o->intdim) flag_incomplete;
  else if (!o->closed) flag_algo;
  else if (pr->conv) flag_conv;

  /* intersect with dest */
  if (dest) {
    opt_oct_mat_t* src2 = dest->closed ? dest->closed : dest->m;
    meet_half(src,src,src2,o->dim,true);
  }
  
  if (respect_closure) return opt_oct_set_mat(pr,o,NULL,src,destructive);
  else return opt_oct_set_mat(pr,o,src,NULL,destructive);
}



opt_oct_t* opt_oct_assign_linexpr_array(elina_manager_t* man,
				bool destructive, opt_oct_t* o,
				elina_dim_t* tdim,
				elina_linexpr0_t** lexpr,
				size_t size,
				opt_oct_t* dest)
{
  if (size==1)
    return opt_oct_assign_linexpr(man,destructive,o,tdim[0],lexpr[0],dest);

  opt_oct_internal_t* pr =
    opt_oct_init_from_manager(man,ELINA_FUNID_ASSIGN_LINEXPR_ARRAY,2*(o->dim+size+5));
  elina_dim_t* d = (elina_dim_t*) pr->tmp2;
  opt_oct_mat_t *src, *dst;
  size_t i;
  elina_dim_t p = o->dim;
  int inexact = 0;
  bool respect_closure = false; /* TODO */
  int src_size = 2*(o->dim)*(o->dim+1);
  
  /* checks */
  if(size<=0){
    return NULL;
  }
  for (i = 0; i < o->dim; i++) {
    d[i] = 0;
  }
  for (i=0;i<size;i++) {
    if (tdim[i] >= o->dim) {
      return NULL;
    }
    if(d[tdim[i]]){
	return NULL;
    }			 /* tdim has duplicate */
    d[tdim[i]] = 1;
  }

  if (dest && !dest->closed && !dest->m)
    /* definitively empty due to dest*/
    return opt_oct_set_mat(pr,o,NULL,NULL,destructive);
  if (pr->funopt->algorithm>=0) opt_oct_cache_closure(pr,o);
  src = o->closed ? o->closed : o->m;
  if (!src) return opt_oct_set_mat(pr,o,NULL,NULL,destructive); /* empty */

  /* add temporary dimensions to hold destination variables */
  dst = opt_hmat_alloc_top(o->dim+size);
  #if defined(TIMING)
  	start_timing();
  #endif
  opt_hmat_set_array(dst->mat,src->mat,src_size);

  /* coefs in expr for temporary dimensions are set to 0 */
  for (i=0;i<2*size;i++){
    pr->tmp[2*o->dim+i+2] = 0;
  }
  /* perform assignments */
  for (i=0;i<size;i++) {

    opt_uexpr u = opt_oct_uexpr_of_linexpr(pr,pr->tmp,lexpr[i],o->intdim,o->dim);

    if (u.type==OPT_EMPTY) {
      opt_hmat_free(dst);
      return opt_oct_set_mat(pr,o,NULL,NULL,destructive);
    }

    if (u.type==OPT_BINARY || u.type==OPT_OTHER) inexact = 1;
    
    opt_hmat_assign(pr,u,dst,o->dim+size,o->dim+i,&respect_closure);
    
  }
  #if defined(TIMING)
  	record_timing(assign_linexpr_time);
  #endif
  /* now close & remove temporary variables */
  if (pr->funopt->algorithm>=0) {
    if (opt_hmat_strong_closure(dst,o->dim+size)) {
      /* empty */
      opt_hmat_free(dst);
      return opt_oct_set_mat(pr,o,NULL,NULL,destructive);
    }
  }
  else flag_algo;
  if (!destructive) src = opt_hmat_alloc(src_size);
  for (i = 0; i < o->dim; i++) {
    d[i] = i;
  }
  for (i=0;i<size;i++) {
    d[o->dim+i] = tdim[i];
    d[tdim[i]] = o->dim;
  }
  opt_hmat_permute(src,dst,o->dim,o->dim+size,d);
  opt_hmat_free(dst);

  /* intersect with dest */
  if (dest) {
    opt_oct_mat_t * src2 = dest->closed ? dest->closed : dest->m;
    meet_half(src,src,src2,o->dim,true);
  }

  if (inexact || o->intdim) flag_incomplete;
  else if (!o->closed) flag_algo;
  else if (pr->conv) flag_conv;

  return opt_oct_set_mat(pr,o,src,NULL,destructive);
}


opt_oct_t* opt_oct_meet_lincons_array(elina_manager_t* man,
			      bool destructive, opt_oct_t* o,
			      elina_lincons0_array_t* array)
{
  printf(".");
  fflush(stdout);
  //printf("meet lincons INPUT\n");
  //elina_lincons0_array_t arr1 = opt_oct_to_lincons_array(man,o);
  //elina_lincons0_array_fprint(stdout,&arr1,NULL);
  //elina_lincons0_array_clear(&arr1);
  //elina_lincons0_array_fprint(stdout,array,NULL);
  //fflush(stdout);
  opt_oct_internal_t* pr =
    opt_oct_init_from_manager(man,ELINA_FUNID_MEET_LINCONS_ARRAY,2*(o->dim+8));
  if (!o->closed && !o->m)
    /* definitively empty */
    return opt_oct_set_mat(pr,o,NULL,NULL,destructive);
  else {
    bool exact, respect_closure;
    int i;
    opt_oct_mat_t * oo = o->closed ? o->closed : o->m;
    /* can / should we try to respect closure */
    respect_closure = (oo==o->closed) && (pr->funopt->algorithm>=0);
    int size = 2 * (o->dim) * (o->dim + 1);
    if (!destructive) oo = opt_hmat_copy(oo,o->dim);

    /* go */
   
    #if defined(TIMING)
  	start_timing();
    #endif
    bool res = opt_hmat_add_lincons(pr,oo,o->intdim,o->dim,array,&exact,&respect_closure);
	
    #if defined(TIMING)
	record_timing(meet_lincons_time);
    #endif
    if (res) {
      /* empty */
      if (!destructive) {
	opt_hmat_free(oo);
	oo = NULL;
      }
      return opt_oct_set_mat(pr,o,NULL,NULL,destructive);
    }
    else {
      /* exact if octagonal constraints & no conversion error */
      if (num_incomplete || !exact) flag_incomplete;
      else if (pr->conv) flag_conv;
	opt_oct_t * r;
      if (respect_closure) r = opt_oct_set_mat(pr,o,NULL,oo,destructive);
      else r = opt_oct_set_mat(pr,o,oo,NULL,destructive);
	/*printf("meet lincons OUT\n");
  elina_lincons0_array_t arr2 = opt_oct_to_lincons_array(man,r);
  elina_lincons0_array_fprint(stdout,&arr2,NULL);
  elina_lincons0_array_clear(&arr2);*/
	return r;
    }
  }
}


opt_oct_t* opt_oct_meet_tcons_array(elina_manager_t* man,
			    bool destructive, opt_oct_t* o,
			    elina_tcons0_array_t* array)
{
  
  return elina_generic_meet_intlinearize_tcons_array(man,destructive,o,array,
						  ELINA_SCALAR_DOUBLE,
						  ELINA_LINEXPR_INTLINEAR,
						  &opt_oct_meet_lincons_array);
}

opt_oct_t* opt_oct_assign_texpr_array(elina_manager_t* man,
			      bool destructive, opt_oct_t* o,
			      elina_dim_t* tdim,
			      elina_texpr0_t** texpr,
			      int size,
			      opt_oct_t* dest)
{
  opt_oct_t *b = elina_generic_assign_texpr_array(man,destructive,o,tdim,texpr,size,dest);
  return b;
  //return elina_generic_assign_texpr_array(man,destructive,o,tdim,texpr,size,dest);
}

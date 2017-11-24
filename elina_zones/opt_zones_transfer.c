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

#include "opt_mat.h"


/****************************
	Assign with linear/non-linear expression
****************************/

opt_zones_t* opt_zones_assign_texpr_array(elina_manager_t* man,
			      bool destructive, opt_zones_t* o,
			      elina_dim_t* tdim,
			      elina_texpr0_t** texpr,
			      size_t size,
			      opt_zones_t* dest)
{    
  
   opt_zones_mat_t * oz = o->closed ? o->closed : o->m;
	comp_list_t *clb=NULL;
	bool need_refine = false;
	if(oz&&!oz->is_dense){
		
		array_comp_list_t * acla = oz->acl;
		
		
		elina_texpr0_t * expr = texpr[0];
		clb = create_comp_list();
		texpr0_to_comp_list_zones(clb,expr);
	
		comp_list_t *cld = find(acla,tdim[0]);
	
		if(!contains_comp(clb,tdim[0])){
			if(cld!=NULL && is_disjoint(cld,clb,o->dim) && (cld->size>1)){
				need_refine = true;
			}
			insert_comp(clb,tdim[0]);
		}
	
		
	}
	
   opt_zones_t *r = elina_generic_assign_texpr_array(man,destructive,o,tdim,texpr,size,dest);
	
   opt_zones_mat_t * or = r->closed ? r->closed : r->m;
	if(or && !or->is_dense){
          if (need_refine) {

            array_comp_list_t *acl = or->acl;
            comp_list_t *clb_copy = create_comp_list();
            comp_list_t *clv = find(acl, tdim[0]);
            comp_t *cb = clb->head;
            while (cb != NULL) {
              if (contains_comp(clv, cb->num)) {
                remove_comp(clv, cb->num);
                insert_comp(clb_copy, cb->num);
              }
              cb = cb->next;
            }
            insert_comp_list_tail(acl, clb_copy);
          }
                free_comp_list(clb);
	}
   return r;
}

/*******************************
	Meet with a linear constraint
********************************/
opt_zones_t* opt_zones_meet_lincons_array(elina_manager_t* man,
			      bool destructive, opt_zones_t* o,
			      elina_lincons0_array_t* array)
{
  printf(".");

  fflush(stdout);
  opt_zones_internal_t* pr =
    opt_zones_init_from_manager(man,ELINA_FUNID_MEET_LINCONS_ARRAY,2*(o->dim+8));
  if (!o->closed && !o->m)
    /* definitively empty */
    return opt_zones_set_mat(pr,o,NULL,NULL,destructive);
  else {
    bool exact, respect_closure;
    size_t i;
    opt_zones_mat_t* oz = o->closed ? o->closed : o->m;
	
    /* can / should we try to respect closure */
    respect_closure = (oz==o->closed) && (pr->funopt->algorithm>=0);

    if (!destructive) oz = opt_zones_mat_copy(oz,o->dim);

    /* go */
    
    #if defined (TIMING)
	start_timing();
    #endif
    bool res = opt_zones_mat_add_lincons(pr,oz,o->intdim,o->dim,array,&exact,&respect_closure);
    // printf("res: %d\n",res);
    // print_mat(oz,o->dim);
#if defined(TIMING)
		record_timing(zones_meet_lincons_time);
     #endif
    if (res) {
      /* empty */
      if (!destructive) opt_zones_mat_free(oz);
      return opt_zones_set_mat(pr,o,NULL,NULL,destructive);
    }
    else {
      /* exact if zonal constraints & no conversion error */
      if (zone_num_incomplete || !exact) zones_flag_incomplete;
      else if (pr->conv) zone_flag_conv;
      opt_zones_t * r;
      if (respect_closure) r= opt_zones_set_mat(pr,o,NULL,oz,destructive);
      else r = opt_zones_set_mat(pr,o,oz,NULL,destructive);
     
      return r;
    }
  }
}


/****************************
	Meet with linear/non-linear constraint
****************************/

opt_zones_t* opt_zones_meet_tcons_array(elina_manager_t* man,
			    bool destructive, opt_zones_t* o,
			    elina_tcons0_array_t* array)
{
  
  opt_zones_t *r = elina_generic_meet_intlinearize_tcons_array(man,destructive,o,array,
						  ELINA_SCALAR_DOUBLE,
						  ELINA_LINEXPR_INTLINEAR,
						  &opt_zones_meet_lincons_array);
  
  return r;
}


/*
	Copyright 2015 Software Reliability Lab, ETH Zurich

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


#include "opt_oct_hmat.h"


bool opt_oct_is_bottom(ap_manager_t* man, opt_oct_t* o)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_IS_BOTTOM,0);
  if (pr->funopt->algorithm>=0){ opt_oct_cache_closure(pr,o);}
  //m = o->closed ? o->closed : o->m;
  if (o->closed) {
    /* definitively non empty on Q */
    if (num_incomplete || o->intdim) { flag_incomplete; }
    return false;
  }
  else if (!o->m){
    /* definitively empty */
    return true;
  }
  else {
    /* no closure => we don't know */
    flag_algo;
    return false;
  }
}

bool opt_oct_is_top(ap_manager_t* man, opt_oct_t* o)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_IS_TOP,0);
  int i,j;
  opt_oct_mat_t* m = o->m ? o->m : o->closed;
  if (!m) return false;
  return is_top_half(m,o->dim);
}

bool opt_oct_is_leq(ap_manager_t* man, opt_oct_t* o1, opt_oct_t* o2)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_IS_LEQ,0);
   if((o1->dim != o2->dim) || (o1->intdim != o2->intdim))return false;
  if (pr->funopt->algorithm>=0){ opt_oct_cache_closure(pr,o1);}
  if (!o1->closed && !o1->m) {
    /* a1 definitively empty */
    return true;
  }
  else if (!o2->closed && !o2->m) {
    /* a2 definitively empty */
    if (o1->closed) {
      /* a1 not empty on Q */
      if (num_incomplete || o1->intdim) { flag_incomplete; }
      return false;
    }
    else { flag_algo; return false; }
  }
  else {
    opt_oct_mat_t *oo1 = o1->closed ? o1->closed : o1->m;
    opt_oct_mat_t *oo2 = o2->closed ? o2->closed : o2->m;
    bool res= is_lequal_half(oo1, oo2, o1->dim);
    //if(res){
	//opt_oct_fprint(stdout,man,oo1,NULL);
    //}
    return res;
  }
}


bool opt_oct_is_eq(ap_manager_t* man, opt_oct_t* o1, opt_oct_t* o2)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_IS_EQ,0);
   if((o1->dim != o2->dim) || (o1->intdim != o2->intdim))return false;
  if (pr->funopt->algorithm>=0) {
    opt_oct_cache_closure(pr,o1);
    opt_oct_cache_closure(pr,o2);
  }
  if (!o1->closed && !o1->m) {
    if (!o2->closed && !o2->m) {
      /* both are empty */
      return true;
    }
    else if (o2->closed) {
      /* a1 empty, e2 not empty on Q */
      if (num_incomplete || o1->intdim) { flag_incomplete; }
      return false;
    }
    else { flag_algo; return false; }
  }
  else if (!o2->closed && !o2->m) {
    if (o1->closed) {
      /* a2 empty, e1 not empty on Q */
      if (num_incomplete || o1->intdim) { flag_incomplete; }
      return false;
    }
    else { flag_algo; return false; }
  }
  else {
    opt_oct_mat_t *oo1 = o1->closed ? o1->closed : o1->m;
    opt_oct_mat_t *oo2 = o2->closed ? o2->closed : o2->m;
    bool res = is_equal_half(oo1,oo2,o1->dim);
    
    
    //if(eq_count==10262){
	//print_hmat(oo1,o1->dim);
	//print_hmat(oo2,o2->dim);
    //}
    return res;
  }
}

ap_tcons0_array_t opt_oct_to_tcons_array(ap_manager_t* man, opt_oct_t* o)
{
  return ap_generic_to_tcons_array(man,o);
}


ap_interval_t** opt_oct_to_box(ap_manager_t* man, opt_oct_t* o)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_TO_BOX,0);
  ap_interval_t** in = ap_interval_array_alloc(o->dim);
  size_t i;
  if (pr->funopt->algorithm>=0) {
	opt_oct_cache_closure(pr,o);
  }
  #if defined(TIMING)
	start_timing();
  #endif
  if (!o->closed && !o->m) {
    /* definitively empty */
    for (i=0;i<o->dim;i++)
      ap_interval_set_bottom(in[i]);
  }
  else {
    /* put variable bounds */
    
    
    opt_oct_mat_t* oo = o->closed ? o->closed : o->m;
    double *m = oo->mat;
    if(!oo->is_dense){
	    array_comp_list_t * acl = oo->acl;
	     for (i=0;i<o->dim;i++){
	     	ap_interval_set_top(in[i]);
	     }
	    comp_list_t * cl = acl->head;
	    while(cl!=NULL){
		    int comp_size = cl->size;
		    comp_t * c = cl->head;
		    for (i=0;i<comp_size;i++){
		      int i1 = c->num;
		      opt_interval_of_bounds(pr,in[i1],
					 m[opt_matpos(2*i1,2*i1+1)],m[opt_matpos(2*i1+1,2*i1)],true);
			 c = c->next;
		    }
		   cl = cl->next;
	    }
    }
    else{
	 for (i=0;i<o->dim;i++){
     		 opt_interval_of_bounds(pr,in[i],
			 m[opt_matpos(2*i,2*i+1)],m[opt_matpos(2*i+1,2*i)],true);
    	}
    }
    man->result.flag_exact = false;
    if (!o->closed) flag_algo;
    else if (num_incomplete || o->intdim) flag_incomplete;
    else if (pr->conv) flag_conv;
   
  }
  #if defined(TIMING)
	record_timing(oct_to_box_time);
  #endif
  return in;
}

ap_interval_t* opt_oct_bound_texpr(ap_manager_t* man,
			       opt_oct_t* o, ap_texpr0_t* expr)
{
  return ap_generic_bound_texpr(man,o,expr,NUM_AP_SCALAR,false);
}


ap_interval_t* opt_oct_bound_dimension(ap_manager_t* man,
				   opt_oct_t* o, ap_dim_t dim)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_BOUND_DIMENSION,0);
  ap_interval_t* r = ap_interval_alloc();
  if(dim>=o->dim){
	ap_interval_free(r);
	return NULL;
  }
  if (pr->funopt->algorithm>=0) opt_oct_cache_closure(pr,o);
  if (!o->closed && !o->m) {
    /* really empty */
    ap_interval_set_bottom(r);
  }
  else if (o->closed) {
    /* optimal in Q */
    opt_oct_mat_t * oo = o->closed;
    double *mm = oo->mat;
    if(!oo->is_dense){
	if(find(oo->acl,dim)==NULL){
		ap_interval_set_top(r);
	}
	else{
		opt_interval_of_bounds(pr,r,
		       mm[opt_matpos(2*dim,2*dim+1)],
		       mm[opt_matpos(2*dim+1,2*dim)],true);
	}
    }
    else{
    	opt_interval_of_bounds(pr,r,
		       mm[opt_matpos(2*dim,2*dim+1)],
		       mm[opt_matpos(2*dim+1,2*dim)],true);
    }
    if (num_incomplete || o->intdim) flag_incomplete;
    else if (pr->conv) flag_conv;
  }
  else {
    /* not optimal */
    opt_oct_mat_t * oo = o->m;
    double *mm = oo->mat;
    if(!oo->is_dense){
	if(find(oo->acl,dim)==NULL){
		ap_interval_set_top(r);
	}
	else{
		opt_interval_of_bounds(pr,r,
		       mm[opt_matpos(2*dim,2*dim+1)],
		       mm[opt_matpos(2*dim+1,2*dim)],true);
	}
    }
    else{
    	opt_interval_of_bounds(pr,r,
		       mm[opt_matpos(2*dim,2*dim+1)],mm[opt_matpos(2*dim+1,2*dim)],
		       true);
    }
    flag_algo;
  }
  return r;
}

ap_lincons0_array_t opt_oct_to_lincons_array(ap_manager_t* man, opt_oct_t* o)
{
  ap_lincons0_array_t ar;
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_TO_LINCONS_ARRAY,0);
  
  if (!o->closed && !o->m) {
    /* definitively empty */
    ar = ap_lincons0_array_make(1);
    ar.p[0] = ap_lincons0_make_unsat();
  }
  else {
    /* put non-oo constraint bounds only */
    opt_oct_mat_t* oo = o->closed ? o->closed : o->m;
    
    double *m = oo->mat;
    int i,j,n=0;
    int size = 2*(o->dim)*(o->dim + 1);
    ar = ap_lincons0_array_make(size);
    if(!oo->is_dense){
	    array_comp_list_t * acl = oo->acl;
	    char * map = (char *)calloc(o->dim,sizeof(char));
	    unsigned short int * cm = (unsigned short int *)calloc(o->dim,sizeof(unsigned short int));
	    comp_list_t * cl = acl->head;
	    int l = 0;
	    while(cl!=NULL){
		comp_t *c = cl->head;
		while(c != NULL){
			unsigned short int num = c->num;
			map[num] = 1;
			cm[num] = l;
			c = c->next;
		}
		cl = cl->next;
		l++;
	    }
	    for (i=0;i<2*o->dim;i++){
	      if(!map[i/2]){
		continue;
	      }
	      for (j=0;j<=(i|1);j++) {
		if(!map[j/2])continue;
		if(cm[j/2] != cm[i/2])continue;
		if ((i==j) ||(m[opt_matpos2(i,j)]==INFINITY)) continue;
		ar.p[n] = opt_lincons_of_bound(pr,i,j,m[opt_matpos2(i,j)]);
		n++;
	      }
	   }
    	free(map);
    	free(cm);
    }
    else{
	for (i=0;i<2*o->dim;i++){
      		for (j=0;j<=(i|1);j++,m++) {
			if ((i==j) || (*m==INFINITY)) continue;

			ar.p[n] = opt_lincons_of_bound(pr,i,j,*m);
			n++;
      		}
	}
   }
    ar.size = n;
   // m = o->closed ? o->closed : o->m;
     
    if (pr->conv) flag_conv;
  }
  return ar;
}

/***********************************
	Check if the bound of a variable is within the given interval
***********************************/
bool opt_oct_sat_interval(ap_manager_t* man, opt_oct_t* o,
		      ap_dim_t dim, ap_interval_t* i)
{
  opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_SAT_INTERVAL,0);
  if(dim >= o->dim){
	return false;
  }
  if (pr->funopt->algorithm>=0) opt_oct_cache_closure(pr,o);
  if (!o->closed && !o->m) {
    /* really empty */
    return true;
  }
  else {
    opt_oct_mat_t* oo = o->closed ? o->closed : o->m;
    double *m = oo->mat;
    ap_interval_t* b = ap_interval_alloc();
    bool r;
    if(!oo->is_dense){
	array_comp_list_t *acl = oo->acl;
	if(find(acl,dim)==NULL){
		ap_interval_set_top(b);
	}
	else{
		/* get (possibly approximated) bounds */
		opt_interval_of_bounds(pr,b,
		       m[opt_matpos(2*dim,2*dim+1)],m[opt_matpos(2*dim+1,2*dim)],true);
	}
    }
    else{
	/* get (possibly approximated) bounds */
	opt_interval_of_bounds(pr,b,
		 m[opt_matpos(2*dim,2*dim+1)],m[opt_matpos(2*dim+1,2*dim)],true);
    }
    
    
    /* compare with i */
    r = (ap_scalar_cmp(b->inf,i->inf)>=0) && (ap_scalar_cmp(b->sup,i->sup)<=0);
    ap_interval_free(b);
    if (r) return true; /* definitively saturates */
    else
      /* definitely does not saturate on Q if closed & no conv error */
      if (num_incomplete || o->intdim) { flag_incomplete; return false; }
      else if (!o->closed) { flag_algo; return false; }
      else if (pr->conv) { flag_conv; return false; }
      else return false;
  }
}

/******************************
 Is dimension unconstrained
*****************************/
bool opt_oct_is_dimension_unconstrained(ap_manager_t* man, opt_oct_t* o,
				    ap_dim_t dim)
{
  opt_oct_internal_t* pr =
    opt_oct_init_from_manager(man,AP_FUNID_IS_DIMENSION_UNCONSTRAINED,0);
  if(dim<o->dim){
	return false;
  }
  if (!o->closed && !o->m)
    /* definitively empty */
    return false;
  else {
    opt_oct_mat_t * oo = o->closed ? o->closed : o->m;
    double * m = oo->mat;
    size_t i, d2=2*dim;
    if(!oo->is_dense){
	array_comp_list_t *acl = oo->acl;
	comp_list_t * cl = acl->head;
	while(cl!=NULL){
		unsigned short int *ca = to_sorted_array(cl,dim);
		unsigned short int comp_size = cl->size;
		unsigned short int j;
		for(j=0; j< comp_size; j++){
			unsigned short int j1 = ca[j];
			if(j1==dim){
				if(m[opt_matpos2(d2,d2+1)]!=INFINITY || m[opt_matpos2(d2+1,d2)]!=INFINITY){
					return false;
				}
			}
			else{
				if(m[opt_matpos2(2*j1,d2)]!=INFINITY){
					return false;
				}
				if(m[opt_matpos2(2*j1+1,d2)]!=INFINITY){
					return false;
				}
				if(m[opt_matpos2(2*j1,d2+1)]!=INFINITY){
					return false;
				}
				if(m[opt_matpos2(2*j1+1,d2+1)]!=INFINITY){
					return false;
				}
			}
		} 
		free(ca);
		cl = cl->next;
	}
    }
    else{
	for (i=0;i<2*o->dim;i++) {
      		if ((m[opt_matpos2(i,d2)]!=INFINITY) && (i!=d2)) return false;
      		if ((m[opt_matpos2(i,d2+1)]!=INFINITY) && (i!=d2+1)) return false;
    	}
    }
    
    return true;
  }
}



/****

SAT Constraints
****/



bool opt_oct_sat_lincons(ap_manager_t *man,opt_oct_internal_t* pr,opt_oct_t* o,
		     ap_lincons0_t* lincons)
{
  
 
    opt_oct_mat_t * oo = o->closed ? o->closed : o->m;
    double *b = oo->mat;   
    size_t i, ui, uj;
    ap_constyp_t c = lincons->constyp;
    opt_uexpr u;
    bool r;

    switch (c) {

      /* skipped */
    case AP_CONS_EQMOD:
    case AP_CONS_DISEQ:
      return false;

      /* handled */
    case AP_CONS_EQ:
    case AP_CONS_SUPEQ:
    case AP_CONS_SUP:
      break;

      /* error */
    default:
      assert(0);
    }

    u = opt_oct_uexpr_of_linexpr(pr,pr->tmp,lincons->linexpr0,o->intdim,o->dim);

    switch (u.type) {

    case OPT_EMPTY:
      /* the empty set has all properties */
     {
      return true;
     }

    case OPT_ZERO:
      if ((c==AP_CONS_SUPEQ && (pr->tmp[0]<=0)) ||
	  /* [-a,b] >= 0 <=> a <= 0 */
	  (c==AP_CONS_SUP && (pr->tmp[0]<0)) ||
	  /* [-a,b] > 0 <=> a < 0 */
	  (c==AP_CONS_EQ && (pr->tmp[0]==0) && (pr->tmp[1]==0))
	  /* [-a,b] = 0 <=> a=b=0 */
	  )
       {
	return true;
	} /* always saturates */
      else {
	/* does not always saturate on Q, if closed and no conv error */
	if (num_incomplete || o->intdim) { flag_incomplete; return false; }
	else if (!o->closed) { flag_algo; return false; }
	else if (pr->conv) { flag_conv; return false; }
	return false;
      }

   case OPT_UNARY:
      if (u.coef_i==1) ui = 2*u.i; else ui = 2*u.i+1;
      //bound_mul_2(pr->tmp[0],pr->tmp[0]);
      pr->tmp[0] = 2*pr->tmp[0];
      //bound_mul_2(pr->tmp[1],pr->tmp[1]);
      pr->tmp[1] = 2*pr->tmp[1];
      //bound_badd(pr->tmp[0],b[matpos(ui,ui^1)]);
      if(!(oo->is_dense) && (find(oo->acl,u.i)==NULL)){
	      pr->tmp[0] += INFINITY;
	      pr->tmp[1] += INFINITY;
      }
      else{
	      pr->tmp[0] += b[opt_matpos(ui,ui^1)];
	      //bound_badd(pr->tmp[1],b[matpos(ui^1,ui)]);
	      pr->tmp[1] += b[opt_matpos(ui^1,ui)];
      }
      if ((pr->tmp[0] <=0) &&
	  /* c_i X_i + [-a,b] >= 0 <=> -c_i X_i + a <= 0 */
	  (c!=AP_CONS_SUP || (pr->tmp[0]<0)) &&
	  /* c_i X_i + [-a,b] >  0 <=> -c_i X_i + a <  0 */
	  (c!=AP_CONS_EQ || (pr->tmp[1]<=0))
	  /* c_i X_i + [-a,b] <= 0 <=>  c_i X_i + b <= 0 */
	  ){
	return true;
	} /* always saturates */
      else {
	/* does not always saturate on Q, if closed and no conv error */
	if (num_incomplete || o->intdim) { flag_incomplete; return false; }
	else if (!o->closed) { flag_algo; return false; }
	else if (pr->conv) { flag_conv; return false; }
	return false;
      }

    case OPT_BINARY:
      if ( u.coef_i==1) ui = 2*u.i; else ui = 2*u.i+1;
      if ( u.coef_j==1) uj = 2*u.j; else uj = 2*u.j+1;
      if(!(oo->is_dense) && !(is_connected(oo->acl,u.i,u.j))){
		pr->tmp[0] += INFINITY;
		pr->tmp[1] += INFINITY;
      }
      else{
	      //bound_badd(pr->tmp[0],b[matpos(uj,ui^1)]);
	      pr->tmp[0] += b[opt_matpos2(uj,ui^1)];
	      //bound_badd(pr->tmp[1],b[matpos(uj^1,ui)]);
	      pr->tmp[1] += b[opt_matpos2(uj^1,ui)];
      }
      if ((pr->tmp[0]<=0) &&
	  /* c_i X_i + c_j X_j + [-a,b] >= 0 <=> -c_i X_i - c_j X_j + a <= 0 */
	  (c!=AP_CONS_SUP || (pr->tmp[0]<0)) &&
	  /* c_i X_i + c_j X_j + [-a,b] >  0 <=> -c_i X_i - c_j X_j + a <  0 */
	  (c!=AP_CONS_EQ || (pr->tmp[1]<=0))
	  /* c_i X_i + c_j X_j + [-a,b] <= 0 <=>  c_i X_i + c_j X_j + b <= 0 */
	  ){
		return true;
	   }
      else {
	/* does not saturate on Q, when closed and no conv error */
	if (num_incomplete || o->intdim) { flag_incomplete; return false; }
	else if (!o->closed) { flag_algo; return false; }
	else if (pr->conv) { flag_conv; return false; }
	return false;
      }

    case OPT_OTHER:
      /* no clue */
      flag_incomplete;
      return false;

    default:
      assert(0);
      return false; /* unreachable */
    }
}

bool opt_oct_sat_lincons_timing(ap_manager_t* man, opt_oct_t* o,
		     ap_lincons0_t* lincons){

   opt_oct_internal_t* pr = opt_oct_init_from_manager(man,AP_FUNID_SAT_LINCONS,
					     2*(o->dim+1));
   if (pr->funopt->algorithm>=0) opt_oct_cache_closure(pr,o);
   if (!o->closed && !o->m) {
    	/* really empty */
   	 return true;
   }
   else{
	#if defined(TIMING)
		start_timing();
  	#endif
	bool res = opt_oct_sat_lincons(man,pr,o,lincons);
	#if defined(TIMING)
		record_timing(sat_lincons_time);
  	#endif
	return res;
   }
}



bool opt_oct_sat_tcons(ap_manager_t* man, opt_oct_t* o,
		   ap_tcons0_t* cons)
{
  return ap_generic_sat_tcons(man,o,cons,NUM_AP_SCALAR,false);
}


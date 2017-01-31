/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
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

/* ********************************************************************** */
/* opt_pk_assign.c: Assignements and Substitutions */
/* ********************************************************************** */

#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"
#include "opt_pk_representation.h"
#include "opt_pk_user.h"
#include "opt_pk_constructor.h"
#include "opt_pk_meetjoin.h"
#include "opt_pk_assign.h"
#include "opt_pk_project.h"


/* ====================================================================== */
/* Matrix transformations: several variables and expressions */
/* ====================================================================== */

/* The list of pair (variable,expr) is given by an array of type
   equation_t.

   IMPRTANT: the array tdim should be sorted in ascending order.
*/

/* insertion sort for sorting the array tdim */
static
void opt_pk_asssub_isort(elina_dim_t* tdim, opt_numint_t** tvec, size_t size)
{
  size_t i,j;

  for (i=1; i<size; i++){
    elina_dim_t dim = tdim[i];
    opt_numint_t* vec = tvec[i];
    for (j=i; j>0; j--){
      if (tdim[j-1]>dim){
	tdim[j] = tdim[j-1];
	tvec[j] = tvec[j-1];
      }
      else
	break;
    }
    tdim[j]=dim;
    tvec[j]=vec;
  }
}


/* ====================================================================== */
/* Inversion of a (deterministic) linear expression */
/* ====================================================================== */
static
void opt_vector_invert_expr(opt_pk_internal_t* opk,
			opt_numint_t* ntab,
			elina_dim_t dim,
			opt_numint_t* tab,
			size_t size)
{
  size_t i;
  size_t var = opk->dec+dim;
  int sgn = opt_numint_sgn(tab[var]);

  assert(sgn!=0);
  if (sgn>0){
    ntab[0] = tab[var];
    ntab[var] = tab[0];
    for (i=1; i<size; i++){
      if (i!=var)
	ntab[i] = -tab[i];
    }
  } else {
    ntab[0] = -tab[var];
    ntab[var] = -tab[0];
    for (i=1; i<size; i++){
      if (i!=var)
	ntab[i] = tab[i];
    }
  }
  opt_vector_normalize_expr(opk,ntab,size);
  return;
}



/* ********************************************************************** */
/* IV. Assignement/Substitution of a single dimension */
/* ********************************************************************** */

/* ====================================================================== */
/* Assignement/Substitution by a *deterministic* linear expression */
/* ====================================================================== */

opt_pk_array_t* opt_poly_asssub_linexpr_det(bool assign, elina_manager_t* man,
			      bool destructive,
			      opt_pk_array_t* oa,
			      elina_dim_t dim, elina_linexpr0_t* linexpr0)
{
  
  int sgn;
  opt_pk_array_t* op;
  opt_pk_internal_t* opk = (opt_pk_internal_t*)man->internal;
  unsigned short int maxcols = oa->maxcols;
  op = destructive ? oa : opt_pk_array_alloc(NULL,NULL,maxcols);
 
  /**************************
		Handle Independent components
  ***************************/
  array_comp_list_t * acla = oa->acl;
  unsigned short int var = dim + opk->dec;
  comp_list_t *clb = linexpr0_to_comp_list(opk,linexpr0);
   
  if(!contains_comp(clb,var)){
	insert_comp(clb,var);
  }
  array_comp_list_t * aclb = create_array_comp_list();
  insert_comp_list(aclb,clb);
  array_comp_list_t * acl = union_array_comp_list(acla,aclb,maxcols);
  
  unsigned short int num_comp = acl->size;
  
  /*************************
		Transform the Linexpr
  *************************/
  short int res = is_comp_list_included(acl,clb,maxcols);
  unsigned short int j,k = 0,comp_size=0;
  comp_list_t * cl = acl->head;
  elina_linexpr0_t * dst = NULL;
  unsigned short int * ca = NULL;
  while(cl!=NULL){
	if(k==res){
		comp_size = cl->size;
		ca = to_sorted_array(cl,maxcols);
		dst = copy_linexpr0_with_comp_list(opk,linexpr0,ca, comp_size);
		break;
	}
	k++;
	cl=cl->next;
  }

  /* Convert linear expression */
  
  opt_vector_set_elina_linexpr0(opk,
			 opk->poly_numintp,
			 dst,
			 comp_size,1);
  
  unsigned short int nvar = 0;
  for(k=0; k < comp_size;k++){
	if(ca[k]==var){
		nvar = k;
		break;
	}
  }
  
   
  
  sgn = opt_numint_sgn(opk->poly_numintp[nvar+2]);
  unsigned short int num_compa = acla->size;
  
  comp_list_t * cla = acla->head;
  //elina_lincons0_array_t arr1 = opt_pk_to_lincons_array(man,op);
  //elina_lincons0_array_fprint(stdout,&arr1,NULL);
  unsigned short int * rmapa = (unsigned short int *)calloc(num_compa, sizeof(unsigned short int));
  size_t * nbmapa = (size_t *)calloc(num_comp,sizeof(size_t));
  size_t * nbeqmapa = (size_t *)calloc(num_comp,sizeof(size_t));
  char * disjoint_map = (char *)calloc(num_compa,sizeof(char));
  opt_pk_t ** poly_a = oa->poly;
  
  for(k=0; k < num_compa; k++){
	short int res_a = is_comp_list_included(acl,cla,maxcols);
	rmapa[k] = res_a;
	nbmapa[res_a] = nbmapa[res_a] + poly_a[k]->C->nbrows;
        nbeqmapa[res_a] = nbeqmapa[res_a] + poly_a[k]->nbeq;
	cla = cla->next;
  }
  
  opt_pk_t ** poly = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
  cl = acl->head;
  
  for(k=0; k < num_comp; k++){
	unsigned short int comp_size = cl->size;
	poly[k] = opt_poly_alloc(comp_size,0);
	cl = cl->next;
  }
	opt_matrix_t *mat = opt_matrix_alloc(nbmapa[res]+1, poly[res]->intdim+2,false);
	poly[res]->C = mat;
	mat->p[0][0] = 1;
	mat->p[0][1] = 1;
	size_t countC = 1;
	cla = acla->head;
	for(k=0; k < num_compa; k++){
		if(rmapa[k]==res){
			unsigned short int * ca_a = to_sorted_array(cla,maxcols);
			unsigned short int comp_size_a = cla->size;
			opt_pk_t * oal = poly_a[k];
			opt_matrix_t * src = oal->C;
			size_t nbconsa = src->nbrows;
			opt_numint_t ** sa = src->p;
			opt_matrix_t *dst_mat = mat;
			size_t count = countC;
			opt_numint_t ** da = dst_mat->p;
			size_t i;
			for(i=0; i < nbconsa; i++){
				opt_numint_t * sai = sa[i];
				opt_numint_t * dai = da[count];
				dai[0] = sai[0];
				dai[1] = sai[1];
				unsigned short int j1 = 0;
				for(j=0; j < comp_size_a; j++){
					while(ca[j1]!=ca_a[j]){
						j1++;
					}
					dai[j1+2] = sai[j+2];
					j1++; 
				}
				count++;
			} 
			countC = count;
			free(ca_a);
		}
		cla = cla->next;
	}
	//if(k==res){
		if (!sgn){ /* Expression is not invertible */
			//opt_pk_t * tmp = opt_poly_alloc(comp_size,0);
			poly[res]->nbeq = nbeqmapa[res];
			elina_dim_t tdim = nvar;
			/*#if defined (CONVERT)
				
				tmp->F = matF;
				if(opk->exn){
					opk->exn = ELINA_EXC_NONE;
					opt_poly_set_top(opk,op);
					man->result.flag_best = man->result.flag_exact = false;
					goto  opt_poly_asssub_linexpr_det_exit;
				}
				poly[k]->F = assign ?
      					     opt_matrix_assign_variable(opk, destructive, tmp->F, tdim, opk->poly_numintp) :
      					     opt_matrix_substitute_variable(opk, destructive, tmp->F, tdim, opk->poly_numintp);
				
				poly[k]->nbline = nblinemapa[k];
				opt_poly_chernikova(man,poly[k],"gen to cons");
				if(poly[k]->satC)
				opt_satmat_fprint(stdout,poly[k]->satC);
			#else*/
				
				//opt_matrix_t * cons = opt_matrix_alloc(1,comp_size+2,true);
				opt_poly_projectforget_array(false,
						  man,poly[res],poly[res],&tdim,1,true);
			    	size_t nbcons = poly[res]->C->nbrows;
			    	opt_matrix_resize_rows(poly[res]->C,nbcons+1);
				opt_numint_t * ov = opk->poly_numintp;
				opt_numint_t * dpi = poly[res]->C->p[nbcons];
				opt_vector_copy(dpi,ov,mat->nbcolumns);
				dpi[0] = 0;
				dpi[nvar+opk->dec] = -1; 
				poly[res]->nbeq++;			
				poly[res]->is_minimized = false;
				opt_poly_chernikova(man,poly[res],"non invertible assign");
			//#endif
	  	}
	  
	  	else { /* Expression is invertible and we have constraints */
	    		/* Invert the expression in opk->poly_numintp2 */
	    		opt_vector_invert_expr(opk,
					       opk->poly_numintp2,
					       nvar, opk->poly_numintp,
					       mat->nbcolumns);
	    		poly[res]->C = assign ? 
				     opt_matrix_substitute_variable(opk,true,mat, nvar, opk->poly_numintp2) :
				     opt_matrix_assign_variable(opk,true,mat, nvar, opk->poly_numintp2);
			poly[res]->nbeq = nbeqmapa[res];
			poly[res]->is_minimized = true;
			opt_poly_chernikova(man,poly[res],"gen to cons");
	  	}
		
	//}
	//else{
	  for(k=0; k < num_compa; k++){
		unsigned short int ind = rmapa[k];
		if(ind==res){	
			continue;
		}
		disjoint_map[k] = 1;
		if(destructive){
			poly[ind]->C = poly_a[k]->C;
			poly[ind]->nbeq = poly_a[k]->nbeq;
			poly[ind]->F = poly_a[k]->F;
			poly[ind]->satF = poly_a[k]->satF;
			poly[ind]->satC = poly_a[k]->satC;
			poly[ind]->nbline = poly_a[k]->nbline;
		}
		else{
			poly[ind]->C = opt_matrix_copy(poly_a[k]->C);
			poly[ind]->nbeq = poly_a[k]->nbeq;
			poly[ind]->F = poly_a[k]->F ? opt_matrix_copy(poly_a[k]->F) : NULL; 
			poly[ind]->satF = poly_a[k]->satF ? opt_satmat_copy(poly_a[k]->satF) : NULL; 
			poly[ind]->satC = poly_a[k]->satC ? opt_satmat_copy(poly_a[k]->satC) : NULL; 
			poly[ind]->nbline = poly_a[k]->nbline;
		}	
	 }
	//cl = cl->next;
    //}	
   	
    opt_poly_asssub_linexpr_det_exit:
	    if(destructive){
		for(k=0; k < num_compa;k++){
			if(!disjoint_map[k]){
				opt_poly_clear(poly_a[k]);
			}
			free(poly_a[k]);
		}
		free_array_comp_list(acla);
		free(poly_a);	 
	    }
	    free(rmapa);
	    free(nbmapa);
	    free(nbeqmapa);
	    free(disjoint_map);
		free(ca);
    elina_linexpr0_free(dst);
    free_array_comp_list(aclb);
    op->poly = poly;
    op->acl = acl;
   
    return op;
}



/* ====================================================================== */
/* Assignement/Substitution by a linear expression */
/* ====================================================================== */
static
opt_pk_array_t* opt_poly_asssub_linexpr(bool assign, bool lazy,
			  elina_manager_t* man,
			  bool destructive,
			  opt_pk_array_t* oa,
			  elina_dim_t dim, elina_linexpr0_t* linexpr,
			  opt_pk_array_t* ob)
{
  
  opt_pk_array_t* op;
  opt_pk_internal_t* opk = (opt_pk_internal_t*)man->internal;
  opt_pk_internal_realloc_lazy(opk,oa->maxcols-1);
  array_comp_list_t * acla = oa->acl;
  unsigned short int maxcols = oa->maxcols;
  unsigned short int intdim = maxcols - 2;
  /* Return empty if empty */
  if (oa->is_bottom || !acla){
    man->result.flag_best = man->result.flag_exact = true;
    return destructive ? oa : opt_pk_bottom(man,intdim, 0);
  }

  unsigned short int num_compa = acla->size;
  unsigned short int k;
  opt_pk_t ** poly_a = oa->poly;
  /* Minimize the argument if option say so */
  if (!lazy){
    for(k=0; k < num_compa;k++){
	opt_pk_t * oak = poly_a[k];
	if(!lazy){
	   opt_poly_chernikova(man,oak,"of the argument");
	}
	else{
	     opt_poly_obtain_C(man,oak,"of the argument");
	}
	if (opk->exn){
	    opk->exn = ELINA_EXC_NONE;
	    man->result.flag_best = man->result.flag_exact = false;
	    if (destructive){
		opt_poly_set_top(opk,oa);
		return oa;
	    } 
	    else {
		return opt_pk_top(man,intdim,0);
	    }
	}
	if(!oak->C && !oak->F){
	   man->result.flag_best = man->result.flag_exact = true;
    	   return destructive ? oa : opt_pk_bottom(man,intdim,0);
	}
     }
  }
   
  /* Choose the right technique */
  if (elina_linexpr0_is_linear(linexpr)){
    op = opt_poly_asssub_linexpr_det(assign, man,destructive,oa,dim,linexpr);
    
    if (ob){
      opt_poly_meet(true,lazy,man,op,op,ob);
	 
    }
  }
  else {
    op = elina_generic_asssub_linexpr_array(true,man,destructive,oa,&dim,&linexpr,1,ob);
  }


  /* Minimize the result if option say so */
  opt_pk_t ** poly = op->poly;
  array_comp_list_t * acl = op->acl;
  unsigned short int num_comp = acl->size;
  if ( !lazy){
    for(k=0; k < num_comp; k++){
	opt_pk_t * opp = poly[k];
	opt_poly_chernikova(man,opp,"assign result");    	
	if (opk->exn){
	    opk->exn = ELINA_EXC_NONE;
	    man->result.flag_best = man->result.flag_exact = false;
	    if (ob){
	        op = opt_pk_copy(man,ob);
	    } 
	    else{
		opt_poly_set_top(opk,op);
	    }
	    return op;
	}
    }
  }
  /* Is the result exact or best ? */
  if (opk->funopt->flag_best_wanted || opk->funopt->flag_exact_wanted){
    man->result.flag_best = man->result.flag_exact =
      (dim < (intdim) || !elina_linexpr0_is_real(linexpr, intdim)) ?
      false :
      true;
  }
  else {
    man->result.flag_best = man->result.flag_exact = (intdim==0);
  }
  return op;
}


/* ********************************************************************** */
/* III. Assignement/Substitution of several dimensions */
/* ********************************************************************** */

/* ====================================================================== */
/* Assignement/Substitution by several *deterministic* linear expressions */
/* ====================================================================== */

opt_pk_array_t* opt_poly_asssub_linexpr_array_det(elina_manager_t* man,
				    bool destructive,
				    opt_pk_array_t* oa,
				    elina_dim_t* tdim, elina_linexpr0_t** texpr,
				    size_t size)
{
  size_t i;
  elina_dim_t* tdim2;
  opt_numint_t** tvec;
  size_t nbcols;
  opt_matrix_t* mat;
  opt_pk_array_t* op;
  opt_pk_internal_t* opk = (opt_pk_internal_t*)man->internal;

  op = destructive ? oa : opt_pk_array_alloc(NULL,NULL,oa->maxcols);
  array_comp_list_t * acla = oa->acl;
  unsigned short int maxcols = oa->maxcols;
  unsigned short int intdim = maxcols - 2;
  /* Return empty if empty */
  if (oa->is_bottom || !acla){
    man->result.flag_best = man->result.flag_exact = true;
    opt_poly_set_bottom(opk,op);
    return op;
  }

  
  /* Convert linear expressions */
  nbcols = oa->maxcols;
  tvec = (opt_numint_t**)malloc(size*sizeof(opt_numint_t*));
  for (i=0; i<size; i++){
    tvec[i] = opt_vector_alloc(nbcols);
    opt_vector_set_elina_linexpr0(opk,
			   tvec[i],
			   texpr[i],
			   intdim,1);
  }
  /* Copy tdim because of sorting */
  tdim2 = (elina_dim_t*)malloc(size*sizeof(elina_dim_t));
  memcpy(tdim2,tdim,size*sizeof(elina_dim_t));
  opt_pk_asssub_isort(tdim2,tvec,size);
  /* Perform the operation */
 // mat =
   // opt_matrix_assign_variables(opk, oa->C, tdim2, tvec, size) :
  /* Free allocated stuff */
  for (i=0; i<size; i++){
    opt_vector_free(tvec[i],nbcols);
  }
  free(tvec);
  free(tdim2);

  /* Update polyhedra */
  if (destructive){
    opt_poly_array_clear(opk,op);
  }
  //op->status = 0;
  return op;
}


/* ====================================================================== */
/* Assignement/Substitution by an array of linear expressions */
/* ====================================================================== */
static
opt_pk_array_t* opt_poly_asssub_linexpr_array(bool lazy,
				elina_manager_t* man,
				bool destructive,
				opt_pk_array_t* oa,
				elina_dim_t* tdim, elina_linexpr0_t** texpr, size_t size,
				opt_pk_array_t* ob)
{
  opt_pk_array_t* op;
  opt_pk_internal_t* opk = (opt_pk_internal_t*)man->internal;
  array_comp_list_t * acla = oa->acl;
  unsigned short int maxcols = oa->maxcols;
  unsigned short int intdim = maxcols - 2;

  /* Return empty if empty */
  if (oa->is_bottom || !acla){
    man->result.flag_best = man->result.flag_exact = true;
    return destructive ? oa : opt_pk_bottom(man,intdim,0);
  }

  unsigned short int num_compa = acla->size;
  unsigned short int k;
  opt_pk_t ** poly_a = oa->poly;
  /* Minimize the argument if option say so */
  if ( !lazy){
    for(k=0; k < num_compa;k++){
	opt_pk_t * oak = poly_a[k];
    	opt_poly_minimize(man,oak);
    	if (opk->exn){
      		opk->exn = ELINA_EXC_NONE;
      		man->result.flag_best = man->result.flag_exact = false;
      		if (destructive){
			opt_poly_set_top(opk,oa);
			return oa;
      		} else {
			return opt_pk_top(man,intdim,0);
      		}
    	}
     }
  }
  
  /* Choose the right technique */
  if (elina_linexpr0_array_is_linear(texpr,size)){
    op = opt_poly_asssub_linexpr_array_det(man,destructive,oa,tdim,texpr,size);
    if (ob){
      opt_poly_meet(true,lazy,man,op,op,ob);
    }
  }
  else {
    op = elina_generic_asssub_linexpr_array(true,man,destructive,oa,tdim,texpr,size,ob);
  }

  /* Minimize the result if option say so */
  opt_pk_t ** poly = op->poly;
  array_comp_list_t * acl = op->acl;
  unsigned short int num_comp = acl->size;
  if ( !lazy){
    for(k=0; k < num_comp; k++){
	opt_pk_t * opp = poly[k];
    	opt_poly_minimize(man,opp);
    	if (opk->exn){
      		opk->exn = ELINA_EXC_NONE;
      		man->result.flag_best = man->result.flag_exact = false;
      		if (ob){
		  op = opt_pk_copy(man,ob);
		} else{
		  opt_poly_set_top(opk,op);
		}
      		return op;
    	}
    }
   
  }

  /* Is the result exact or best ? */
  if (opk->funopt->flag_best_wanted || opk->funopt->flag_exact_wanted){
    size_t i;
    man->result.flag_best = true;
    for (i=0;i<size;i++){
      if (tdim[i] < intdim || !elina_linexpr0_is_real(texpr[i], intdim)){
	man->result.flag_best = false;
	break;
      }
    }
    man->result.flag_exact = man->result.flag_best;
  }
  else {
    man->result.flag_best = man->result.flag_exact = (intdim==0);
  }

  return op;
}


/* ====================================================================== */
/* Assignement/Substitution by an array of tree expressions */
/* ====================================================================== */
opt_pk_array_t* opt_poly_asssub_texpr_array(bool assign,
			      bool lazy,
			      elina_manager_t* man,
			      bool destructive,
			      opt_pk_array_t* oa,
			      elina_dim_t* tdim, elina_texpr0_t** texpr, size_t size,
			      opt_pk_array_t* ob)
{
  
  opt_pk_array_t* op;
  opt_pk_internal_t* opk = (opt_pk_internal_t*)man->internal;
  array_comp_list_t * acla = oa->acl;
  op = destructive ? oa : opt_pk_array_alloc(NULL,NULL,oa->maxcols);
  /* Return empty if empty */
  if(oa->is_bottom || !acla){
	 man->result.flag_best = man->result.flag_exact = true;
	 opt_poly_set_bottom(opk,op);
    	 return op;
  }

  unsigned short int num_compa = acla->size;
  unsigned short int k;
  opt_pk_t ** poly_a = oa->poly; 
  unsigned short int intdim = oa->maxcols - 2;
  /* Minimize the argument if option say so */
  //if ( !lazy){
    for(k=0; k < num_compa; k++){
	opt_pk_t * oak = poly_a[k];
	opt_poly_chernikova(man,oak,"asssub texpr array");
	if(opk->exn){
	   opk->exn = ELINA_EXC_NONE;
	   man->result.flag_best = man->result.flag_exact = false;
	   if (destructive){
	       opt_poly_set_top(opk,oa);
	       return oa;
	   } 
	   else {
	        return opt_pk_top(man,intdim,0);
	   }
        }
	if(!oak->C && !oak->F){
		man->result.flag_best = man->result.flag_exact = true;
    		return destructive ? oa : opt_pk_bottom(man,intdim,0);
	}
	
    }
    
  //}
  
  /* Choose the right technique */
  if (elina_texpr0_array_is_scalar(texpr,size) &&
      elina_texpr0_array_is_interval_linear(texpr,size)){
    elina_abstract0_t abs;
    abs.value = oa;
    abs.man = man;
    
    if (size==1){
      elina_linexpr0_t* linexpr0 =
	elina_intlinearize_texpr0(man,&abs,texpr[0],NULL,
			       ELINA_SCALAR_MPQ,false);
      op = opt_poly_asssub_linexpr_det(assign,man,destructive,
				   oa, tdim[0], linexpr0);
      elina_linexpr0_free(linexpr0);
	
    }
    else {
      elina_linexpr0_t** tlinexpr =
	elina_intlinearize_texpr0_array(man,&abs,texpr,size,NULL,
				     ELINA_SCALAR_MPQ,false);
      op = opt_poly_asssub_linexpr_array_det(man,destructive,
					 oa,tdim,tlinexpr,size);
      
      elina_linexpr0_array_free(tlinexpr,size);
    }
    if (ob){
      opt_poly_meet(true,lazy,man,op,op,ob);
    }
  }
  else {
    op = elina_generic_asssub_texpr_array(assign,man,destructive,oa,tdim,texpr,size,ob);
  }
  
  /* Minimize the result if option say so */
  array_comp_list_t * acl = op->acl;
  unsigned short int num_comp = acl->size;
  opt_pk_t ** poly = op->poly;
  if ( !lazy){
    for(k=0; k < num_comp; k++){
	opt_pk_t * opp = poly[k];
	opt_poly_chernikova(man,opp,"asssub texpr array");
	if(opk->exn){
	   opk->exn = ELINA_EXC_NONE;
	   man->result.flag_best = man->result.flag_exact = false;
      	   if (ob){ 
		if(destructive){
			opt_poly_array_clear(opk,op);
			free(op);
		}
		else{
			free(op);
		}
	       op = opt_pk_copy(man,ob);
	   } else {opt_poly_set_top(opk,op);}
      			return op;
		}
    }
   
  }
  /* Is the result exact or best ? */
  if (opk->funopt->flag_best_wanted || opk->funopt->flag_exact_wanted){
    man->result.flag_best = elina_texpr0_array_is_interval_linear(texpr,size);
    man->result.flag_exact = man->result.flag_best;
  }
  else {
    man->result.flag_best = man->result.flag_exact = false;
  }
  return op;
}


/* ********************************************************************** */
/* V. Assignement/Substitution: interface */
/* ********************************************************************** */

opt_pk_array_t* opt_pk_assign_linexpr_array(elina_manager_t* man,
			      bool destructive, opt_pk_array_t* oa,
			      elina_dim_t* tdim, elina_linexpr0_t** texpr,
			      size_t size,
			      opt_pk_array_t* ob)
{
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
  
  opt_pk_array_t* op;
  #if defined(TIMING)
 	 start_timing();
  #endif

  op = size==1 ? opt_poly_asssub_linexpr(true, opk->funopt->algorithm<=0, man,destructive,oa,tdim[0],texpr[0],ob) :
       		 opt_poly_asssub_linexpr_array( opk->funopt->algorithm<=0,
			      man,destructive,oa,tdim,texpr,size,ob);
  //assert(poly_check(pk,po));
  #if defined(TIMING)
 	 	record_timing(assign_linexpr_time);
  #endif
 
  return op;
}

opt_pk_array_t* opt_pk_assign_texpr_array(elina_manager_t* man,
			    bool destructive, opt_pk_array_t* oa,
			    elina_dim_t* tdim,
			    elina_texpr0_t** texpr,
			    size_t size,
			    opt_pk_array_t* dest)
{
  
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_ASSIGN_TEXPR_ARRAY);
  opt_pk_array_t* op;
  #if defined(TIMING)
 	 start_timing();
  #endif
  op = opt_poly_asssub_texpr_array(true,
			       opk->funopt->algorithm<=0,
			       man,destructive,oa,tdim,texpr,size,dest);
  #if defined(TIMING)
 	 	record_timing(assign_linexpr_time);
  #endif
  return op;
}


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

/* ********************************************************************** */
/* opt_pk_assign.c: Assignements and Substitutions */
/* ********************************************************************** */

#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"
#include "opt_pk_representation.h"
#include "opt_pk_user.h"
#include "opt_pk_test.h"
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
  //printf("ASSIGN INPUT\n");
  //elina_lincons0_array_t arr = opt_pk_to_lincons_array(man,oa);
  //elina_lincons0_array_fprint(stdout,&arr,NULL);
  //printf("x%d:= ",dim);
  //elina_linexpr0_fprint(stdout,linexpr0,NULL);
  //printf("\n");
  //elina_lincons0_array_clear(&arr);
  //fflush(stdout);
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
  comp_list_t * cli = find(acla,var);
  bool need_refine = false; 
  comp_list_t * clb_copy = NULL;
  if(!contains_comp(clb,var)){
	
	if(cli!=NULL && is_disjoint(cli,clb,maxcols) && (cli->size>1)){
        	need_refine = true;
	}
	insert_comp(clb,var);
	if(need_refine){
		clb_copy = copy_comp_list(clb);
	}
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
  size_t * num_vertex_a = (size_t *)calloc(num_compa,sizeof(size_t));
  size_t  nbvertex = 0;
  size_t  nbline = 0;
  size_t  nbcons = 0;
  size_t  nbeq = 0;
  char * disjoint_map = (char *)calloc(num_compa,sizeof(char));
  char * line_map = (char *)calloc(comp_size, sizeof(char));
  opt_pk_t ** poly_a = oa->poly;
  // whether to use generators or constraints
  bool flag1 = !sgn && assign;
  bool flag2 = sgn && ! assign;
  for(k=0; k < num_compa; k++){
	short int res_a = is_comp_list_included(acl,cla,maxcols);
	rmapa[k] = res_a;
	opt_pk_t *oak = poly_a[k];
	if(res_a==res){
		disjoint_map[k] = 1;
		unsigned short int *ca_a = to_sorted_array(cla,maxcols);
		if(flag1){
			num_vertex_a[k] = opt_generator_rearrange(oak->F,oak->satF);
			if(!nbvertex){
	 			nbvertex = num_vertex_a[k];
      			}
      			else{
	 			nbvertex = nbvertex*num_vertex_a[k];
      			}
      			nbline += oak->F->nbrows - num_vertex_a[k];
			nbcons +=  poly_a[k]->C->nbrows;
			nbeq +=  poly_a[k]->nbeq;
		}
		else if(flag2){
			num_vertex_a[k] = opt_generator_rearrange(oak->F,oak->satF);
			if(!nbvertex){
	 			nbvertex = num_vertex_a[k];
      			}
      			else{
	 			nbvertex = nbvertex*num_vertex_a[k];
      			}
      			nbline += oak->F->nbrows - num_vertex_a[k];
		}
		else{
			nbcons +=  poly_a[k]->C->nbrows;
			nbeq +=  poly_a[k]->nbeq;
		}
		
		unsigned short int k2=0, k3 = 0;
		for(k2=0; k2< cla->size; k2++){
			while(ca[k3]!=ca_a[k2]){
				k3++;
			}
			line_map[k3] = 1;
		}
		free(ca_a);
	}
	cla = cla->next;
  }
  opt_pk_t ** poly;
  if(need_refine && assign){
	poly = (opt_pk_t **)malloc((num_comp+1)*sizeof(opt_pk_t *));
  }
  else{
	poly = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
  }
 
  cl = acl->head;
  char * exc_map = (char *)calloc(num_comp,sizeof(char));
  for(k=0; k < num_comp; k++){
	unsigned short int comp_size = cl->size;
	poly[k] = opt_poly_alloc(comp_size,0);
	cl = cl->next;
  }
	size_t num_uncons_var = 0;
	opt_matrix_t *matC = NULL;
	opt_matrix_t *matF = NULL; 
	if(flag1){
		for(k=0; k < comp_size; k++){
			if(!line_map[k]){
				nbline++;
			}
		}
		matF =  opt_matrix_alloc(nbvertex+nbline+1, poly[res]->intdim+2,false);
		matC = opt_matrix_alloc(nbcons+1, poly[res]->intdim+2,false);
	}
	else if(flag2){
		matF =  opt_matrix_alloc(nbvertex+nbline+1, poly[res]->intdim+2,false);
	}
	else{
		matC = opt_matrix_alloc(nbcons+1, poly[res]->intdim+2,false);
	}
	if(flag1 || flag2){
		fuse_generators_intersecting_blocks(matF,poly_a,acla,ca,num_vertex_a,disjoint_map,maxcols);
		size_t nbrows = matF->nbrows;
		for(k=0; k < comp_size; k++){
			if(!line_map[k]){
				matF->p[nbrows][k+opk->dec] = 1;
				nbrows++;
			}
		}
		if(!nbvertex){
			matF->p[nbrows][0] = 1;
			matF->p[nbrows][1] = 1;
			nbrows++;
		}
		matF->nbrows=nbrows;
		poly[res]->F = matF;
	}
	if(!flag2){
		//matC->p[0][0] = 1;
		//matC->p[0][1] = 1;
		size_t countC = 0;
		cla = acla->head;
		bool pos_flag = true;
		for(k=0; k < num_compa; k++){
			if(rmapa[k]==res){
				unsigned short int * ca_a = to_sorted_array(cla,maxcols);
				unsigned short int comp_size_a = cla->size;
				opt_pk_t * oal = poly_a[k];
				opt_matrix_t * src = oal->C;
				size_t nbconsa = src->nbrows;
				opt_numint_t ** sa = src->p;
				opt_matrix_t *dst_mat = matC;
				size_t count = countC;
				opt_numint_t ** da = dst_mat->p;
				size_t i;
				bool flag = false;
				for(i=0; i < nbconsa; i++){
					opt_numint_t * sai = sa[i];
					opt_numint_t * dai = da[count];
					if(opt_vector_is_positivity_constraint(opk,sai,src->nbcolumns)){
						flag = true;
					}
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
				if(!flag){
					pos_flag = false;
				}
				countC = count;
				free(ca_a);
			}
			cla = cla->next;
		}
		if(pos_flag){
			matC->p[countC][0] = 1;
			matC->p[countC][1] = 1;
			countC++;			
		}
		matC->nbrows = countC;
		poly[res]->C = matC;
	}
	//     
	
		if (!sgn){ /* Expression is not invertible */
			//opt_pk_t * tmp = opt_poly_alloc(comp_size,0);
			//
			//elina_dim_t tdim = nvar;
			if(assign){
				//poly[res]->nbline = nbline;
				
				poly[res]->F = opt_matrix_alloc(matF->nbrows,matF->nbcolumns,false);
				poly[res]->nbline = opt_matrix_assign_variable(opk, poly[res]->F, matF, nvar, opk->poly_numintp);
				
				poly[res]->nbeq = nbeq;
				
				elina_dim_t *tdim=(elina_dim_t *)malloc(sizeof(elina_dim_t));
				tdim[0] = nvar;
				opt_poly_projectforget_array(false,
						  man,poly[res],poly[res],tdim,1,true);
			    	size_t nbcons = poly[res]->C->nbrows;
			    	opt_matrix_resize_rows(poly[res]->C,nbcons+1);
				opt_numint_t * ov = opk->poly_numintp;
				opt_numint_t * dpi = poly[res]->C->p[nbcons];
				opt_vector_copy(dpi,ov,matC->nbcolumns);
				dpi[0] = 0;
				dpi[nvar+opk->dec] = -ov[0]; 
				poly[res]->nbeq++;			
				poly[res]->is_minimized = false;
				if(!opk->exn){
					/*if( need_refine){
					
						comp_list_t * clv = find(acl,var);
					
					
						comp_t * cb = clb_copy->head;
						while(cb!=NULL){
							remove_comp(clv,cb->num);
							cb = cb->next;
						}
						unsigned short int *ca1 = to_sorted_array(clv,maxcols);
						unsigned short int *ca2 = to_sorted_array(clb_copy,maxcols);
						unsigned short int * ind_map_a = map_index(ca1,ca,clv->size);
						unsigned short int * ind_map_b = map_index(ca2,ca,clb_copy->size);
						opt_pk_t * tmp = poly[res];
						opt_matrix_t * F = tmp->F;
						opt_matrix_t * C = tmp->C;
						bool is_pos = false;
						poly[res] = opt_poly_alloc(clv->size,0);
						poly[res]->C = opt_matrix_alloc(C->nbrows+1,clv->size+opk->dec,false);
						poly[res]->F = opt_matrix_alloc(F->nbrows,clv->size+opk->dec,false);
						poly[res]->nbeq = split_matrix(opk,poly[res]->C,C,ind_map_a,clv->size, &is_pos);
						//if(!is_pos){
						//	size_t nbrows = poly[res]->C->nbrows;
						//	poly[res]->C->p[nbrows][0] = 1;
						//	poly[res]->C->p[nbrows][1] = 1;
						//	poly[res]->C->nbrows++;
						//}
						poly[res]->nbline = split_matrix(opk,poly[res]->F,F,ind_map_a,clv->size,&is_pos); 

						is_pos = false;
						poly[num_comp] = opt_poly_alloc(clb_copy->size,0);
						poly[num_comp]->C = opt_matrix_alloc(C->nbrows+1,clb_copy->size+opk->dec,false);
						poly[num_comp]->F = opt_matrix_alloc(F->nbrows,clb_copy->size+opk->dec,false); 
						poly[num_comp]->nbeq = split_matrix(opk,poly[num_comp]->C,C,ind_map_b,clb_copy->size, &is_pos);
						//if(!is_pos){
						//	size_t nbrows = poly[num_comp]->C->nbrows;
						//	poly[num_comp]->C->p[nbrows][0] = 1;
						//	poly[num_comp]->C->p[nbrows][1] = 1;
						//	poly[num_comp]->C->nbrows++;
						//}
						poly[num_comp]->nbline = split_matrix(opk,poly[num_comp]->F,F,ind_map_b,clb_copy->size, &is_pos); 
					
					
						poly[res]->satC = opt_satmat_alloc(poly[res]->F->nbrows,opt_bitindex_size(poly[res]->C->nbrows));
						combine_satmat(opk,poly[res],clv->size,poly[res]->C->nbrows,true);
						poly[num_comp]->satC = opt_satmat_alloc(poly[num_comp]->F->nbrows,opt_bitindex_size(poly[num_comp]->C->nbrows));
						combine_satmat(opk,poly[num_comp],clb_copy->size,poly[num_comp]->C->nbrows,true);

						 
						insert_comp_list_tail(acl,clb_copy);
						 
					
						free(ca1);
						free(ca2);
						free(ind_map_a);
						free(ind_map_b);
					
						opt_matrix_free(C);
						opt_matrix_free(F);
						free(tmp);
					
					}
					else {*/
						poly[res]->satC = opt_satmat_alloc(poly[res]->F->nbrows,opt_bitindex_size(poly[res]->C->nbrows));
						combine_satmat(opk,poly[res],matC->nbcolumns - opk->dec,poly[res]->C->nbrows,true);
					//}
				}
				else{
					opk->exn = ELINA_EXC_NONE;
					exc_map[res] = 1;
				}
				opt_matrix_free(matF);
				free(tdim);
			}
			else{
				poly[res]->C = opt_matrix_substitute_variable(opk,true, matC, nvar, opk->poly_numintp);
				if(!opk->exn){
					opt_poly_chernikova(man,poly[res],"non invertible assign");
					
				}
				else{
					opk->exn = ELINA_EXC_NONE;
					exc_map[res] = 1;
				}
			}
				
				//opt_matrix_t * cons = opt_matrix_alloc(1,comp_size+2,true);
				
				
	  	}
	  
	  	else { /* Expression is invertible and we have constraints */
	    		/* Invert the expression in opk->poly_numintp2 */
	    		opt_vector_invert_expr(opk,
					       opk->poly_numintp2,
					       nvar, opk->poly_numintp,
					       matC? matC->nbcolumns: matF->nbcolumns);
			
			if(assign){
				
				poly[res]->C = opt_matrix_substitute_variable(opk,true,matC, nvar, opk->poly_numintp2);
				if(opk->exn){
					opk->exn = ELINA_EXC_NONE;
					exc_map[res] = 1;
				}
			}
			else{
				poly[res]->F = opt_matrix_alloc(matF->nbrows,matF->nbcolumns,false);
				poly[res]->nbline =  opt_matrix_assign_variable(opk,poly[res]->F,matF, nvar, opk->poly_numintp2);
				if(opk->exn){
					opk->exn = ELINA_EXC_NONE;
					exc_map[res] = 1;
				}
				opt_matrix_free(matF);
			}

			//poly[res]->nbeq = nbeqmapa[res];
			//poly[res]->is_minimized = true;
			
			opt_poly_chernikova(man,poly[res],"gen to cons");
			if(opk->exn){
				opk->exn = ELINA_EXC_NONE;
				exc_map[res] = 1;
			}
	  	}
		
	//}
	//else{
	
	  for(k=0; k < num_compa; k++){
		disjoint_map[k] = 0;
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
    free(num_vertex_a);
    free(rmapa);
    free(disjoint_map);
    free(ca);
    free(line_map);
    elina_linexpr0_free(dst);
    free_array_comp_list(aclb);
    k=0;
    cl = acl->head;
    while(k < num_comp){
	opt_pk_t *oak = poly[k];
	if(exc_map[k]){
		comp_list_t * tmp = cl;
		cl = cl->next;
		remove_comp_list(acl,tmp);
		unsigned short int k1;
		for(k1=k; k1 < num_comp - 1; k1++){
			poly[k1] = poly[k1+1];
		}
		opt_poly_clear(oak);
		num_comp--; 
	}
	else{
		k++;
		cl=cl->next;
	}
    }
    op->poly = poly;
    op->acl = acl;
    free(exc_map);
      //printf("ASSIGN OUTPUT\n");
	
	//elina_lincons0_array_t arr1 = opt_pk_to_lincons_array(man,op);
	//elina_lincons0_array_fprint(stdout,&arr1,NULL);
	//elina_lincons0_array_clear(&arr1);
	//fflush(stdout);
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
    op = elina_generic_asssub_linexpr_array(assign,man,destructive,oa,&dim,&linexpr,1,ob);
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

opt_pk_array_t* opt_poly_asssub_linexpr_array_det(bool assign, elina_manager_t* man,
				    bool destructive,
				    opt_pk_array_t* oa,
				    elina_dim_t* tdim, elina_linexpr0_t** texpr,
				    size_t size)
{
  size_t i;
  elina_dim_t* tdim2;
  opt_numint_t** tvec;
  //size_t nbcols;
  opt_matrix_t* mat;
  opt_pk_array_t* op;
  opt_pk_internal_t* opk = (opt_pk_internal_t*)man->internal;
  unsigned short int k;
  //op = destructive ? oa : opt_pk_array_alloc(NULL,NULL,oa->maxcols);
  opt_pk_t ** poly_a = oa->poly;
  array_comp_list_t * acla = oa->acl;
  unsigned short int num_compa = acla->size;
  unsigned short int maxcols = oa->maxcols;
  //unsigned short int intdim = maxcols - 2;
  /* Return empty if empty */
  if (oa->is_bottom || !acla){
    man->result.flag_best = man->result.flag_exact = true;
    
    return opt_pk_bottom(man,oa->maxcols-opk->dec,0);
  }

  
  /* Convert linear expressions */
  //nbcols = oa->maxcols;
  
  array_comp_list_t * aclb = create_array_comp_list();
  comp_list_t ** clb_arr = (comp_list_t **)malloc(size*sizeof(comp_list_t *));
  for (i=0; i<size; i++){
    
    comp_list_t *clb = linexpr0_to_comp_list(opk,texpr[i]);
     if(!contains_comp(clb,tdim[i]+opk->dec)){
		insert_comp(clb,tdim[i]+opk->dec);
     }
     clb_arr[i] = copy_comp_list(clb);
     insert_comp_list_with_union(aclb,clb,maxcols);	
  }

  /**********************************
	Compute LUB of partitions
  **********************************/
  array_comp_list_t *acl = union_array_comp_list(acla,aclb,maxcols);
  unsigned short int num_comp = acl->size;
  unsigned short int ** ca_arr = (unsigned short int **)malloc(num_comp*sizeof(unsigned short int *));
  unsigned short int * comp_size_map = (unsigned short int *)calloc(num_comp,sizeof(unsigned short int));
  comp_list_t * cl = acl->head;
  for(k=0; k < num_comp; k++){
	unsigned short int comp_size = cl->size; 
	ca_arr[k] = to_sorted_array(cl,maxcols);
	comp_size_map[k] = comp_size;
	cl = cl->next;
  }
  /**********************************
	Factor assignment statement according to LUB
  ***********************************/

  unsigned short int * rmapb = (unsigned short int *)calloc(size, sizeof(unsigned short int));
  size_t * nbmapb = (size_t *)calloc(num_comp,sizeof(size_t));
  elina_linexpr0_t ** expr_array = (elina_linexpr0_t **)malloc(size*sizeof(elina_linexpr0_t *));
  for(i=0;i < size; i++){
	unsigned short int ind = is_comp_list_included(acl,clb_arr[i],maxcols);
	rmapb[i] = ind;
	expr_array[i] = copy_linexpr0_with_comp_list(opk,texpr[i],ca_arr[ind],comp_size_map[ind]);
	nbmapb[ind]++;
  }
 

  /*********************************
	Factor A according to LUB
  **********************************/
  unsigned short int * rmapa = (unsigned short int *)calloc(num_compa, sizeof(unsigned short int));
  char * disjoint_map = (char *)calloc(num_compa,sizeof(char));
  size_t * nbgenmapa = (size_t *)calloc(num_comp,sizeof(size_t));
  size_t * nblinemapa = (size_t *)calloc(num_comp,sizeof(size_t));
  size_t * num_vertex_a = (size_t *)calloc(num_compa,sizeof(size_t));
  size_t * num_vertex = (size_t *)calloc(num_comp,sizeof(size_t));
  unsigned short int * array_map_a = create_array_map(acla,maxcols);
  comp_list_t * cla = acla->head;
  for(k=0; k < num_compa; k++){
	opt_pk_t * oak = poly_a[k];
	short int res = is_comp_list_included(acl,cla,maxcols);
	rmapa[k] = res;
	if(assign){
		opt_poly_obtain_satF(oak);
		num_vertex_a[k] = opt_generator_rearrange(oak->F,oak->satF);
		if(num_vertex_a[k]){
			if(!num_vertex[res]){
				num_vertex[res] = num_vertex_a[k];
			}
			else{
				num_vertex[res] = num_vertex[res] * num_vertex_a[k];
			}
		}
		nbgenmapa[res] = nbgenmapa[res] + oak->F->nbrows;
		nblinemapa[res] = nblinemapa[res] + oak->nbline;
        }
        else{
		nbgenmapa[res] = nbgenmapa[res] + oak->C->nbrows;
		nblinemapa[res] = nblinemapa[res] + oak->nbeq;
	}
	cla = cla->next;
  }
  
  opt_pk_t ** poly = (opt_pk_t **)malloc(num_comp*sizeof(opt_pk_t *));
  size_t * counterF = (size_t *)calloc(num_comp,sizeof(size_t));
  char * pos_con_map = (char *)calloc(num_comp,sizeof(char));
  cl = acl->head;
  for(k=0; k < num_comp; k++){
	unsigned short int comp_size = cl->size;
	poly[k] = opt_poly_alloc(comp_size,0);
	pos_con_map[k] = 1;
	cl = cl->next;
  }

  cl = acl->head;  	
  for(k=0; k < num_comp; k++){
	unsigned short int l,j;
	if(!nbmapb[k]){
		for(l=0; l < num_compa;l++){
			if(rmapa[l]==k){
				break;
			}
		}
		disjoint_map[l] = 1;
		if(destructive){
			poly[k]->C = poly_a[l]->C;
			poly[k]->nbeq = poly_a[l]->nbeq;
			poly[k]->F = poly_a[l]->F;
			poly[k]->satF = poly_a[l]->satF;
			poly[k]->satC = poly_a[l]->satC;
			poly[k]->nbline = poly_a[l]->nbline;
		}
		else{
			
			poly[k]->C =poly_a[l]->C ? opt_matrix_copy(poly_a[l]->C) : NULL;
			poly[k]->nbeq = poly_a[l]->nbeq;
			poly[k]->F = poly_a[l]->F ? opt_matrix_copy(poly_a[l]->F) : NULL; 
			poly[k]->satF = poly_a[l]->satF ? opt_satmat_copy(poly_a[l]->satF) : NULL; 
			poly[k]->satC = poly_a[l]->satC ? opt_satmat_copy(poly_a[l]->satC) : NULL; 
			poly[k]->nbline = poly_a[l]->nbline;
		}
	}
	else{
		unsigned short int comp_size = comp_size_map[k];
		unsigned short int k1;
		unsigned short int * ca = ca_arr[k];
		if(assign){
			
			//poly[k]->C = opt_matrix_alloc(nbmapa[k]+1, comp_size+2,false);
			//poly[k]->nbeq = nbeqmapa[k];		
			unsigned short int nblines = 0;
			for(k1=0; k1 < comp_size; k1++){
				unsigned short int var = ca[k1];
				if(!array_map_a[var]){
					nblines++;
				}
			}
			poly[k]->F = opt_matrix_alloc(nbgenmapa[k]+2*num_vertex[k]+nblines+1, comp_size+2,false);
			num_vertex[k] = 0;
			nblines = 0;
			for(k1=0; k1 < comp_size; k1++){
				unsigned short int var = ca[k1];
				if(!array_map_a[var]){
					poly[k]->F->p[nblines][k1+2] = 1;
					nblines++;
				}
			}
			poly[k]->nbline = nblinemapa[k] + nblines;
			if(nblines==comp_size){
				poly[k]->F->p[nblines][0] = 1;
				poly[k]->F->p[nblines][1] = 1;
				nblines++;
			}
			counterF[k] = nblines;
			poly[k]->F->nbrows = nblines;
		}
		else{
			poly[k]->C = opt_matrix_alloc(nbgenmapa[k]+1,comp_size+2,false);
			poly[k]->nbeq = nblinemapa[k]; 
			counterF[k] = 0;
		}
	}
	cl = cl->next;
  }
  
  free(array_map_a);
	
  if(assign){
	  // cartesian product of vertices
	  cartesian_product_vertices_with_map(oa, poly, rmapa, ca_arr, num_vertex_a, num_vertex, counterF, disjoint_map);

	  // meet of rays
	  meet_rays_with_map(oa, poly, rmapa, ca_arr, num_vertex_a, counterF, disjoint_map);
  }

  else{
		meet_cons_with_map(opk,oa,poly,rmapa,ca_arr,counterF,pos_con_map,disjoint_map);	
		for(k=0; k < num_comp; k++){
		    size_t count = counterF[k];
	            if(pos_con_map[k] && nbmapb[k]){
		       poly[k]->C->p[count][0] = 1;
		       poly[k]->C->p[count][1] = 1;
		       counterF[k]++;
	 	       poly[k]->C->nbrows = counterF[k];
	            }
  		}
  }
  char * exc_map = (char *)calloc(num_comp,sizeof(char));
  /* Copy tdim because of sorting */
  for(k=0; k < num_comp; k++){
        if(nbmapb[k]){
		tvec = (opt_numint_t **)malloc(nbmapb[k]*sizeof(opt_numint_t *));
		tdim2 = (elina_dim_t*)malloc(nbmapb[k]*sizeof(elina_dim_t));
		unsigned short int l = 0;
		unsigned short int comp_size = comp_size_map[k];
		unsigned short int * ca = ca_arr[k];
		for(i=0; i <size; i++){
			if(rmapb[i]==k){
				tvec[l] = opt_vector_alloc(comp_size+2);
				opt_vector_set_elina_linexpr0(opk,tvec[l],expr_array[i],comp_size,1);
				unsigned short int k1 = 0;
				unsigned short int var = tdim[i] + opk->dec;
				while(ca[k1]!=var){
					k1++;
				}
				tdim2[l] = k1;
				l++;
			}
			
		}
	  	opt_pk_asssub_isort(tdim2,tvec,nbmapb[k]);
		opt_pk_t * oak = poly[k];
	   	/* Perform the assignment operation */
		opt_matrix_t * tmp = assign ? oak->F : oak->C;
		if(assign){
	 		poly[k]->F = opt_matrix_assign_variables(opk, oak->F, tdim2, tvec, nbmapb[k]);
			if(opk->exn){
				opk->exn = ELINA_EXC_NONE;
				exc_map[k] = 1;
			}
		}
		else{
			poly[k]->C = opt_matrix_substitute_variables(opk,oak->C, tdim2, tvec, nbmapb[k]);
			if(opk->exn){
				opk->exn = ELINA_EXC_NONE;
				exc_map[k] = 1;
			}
		}
		opt_matrix_free(tmp);
	  	/* Free allocated stuff */
		for(i=0; i < nbmapb[k]; i++){
			opt_vector_free(tvec[i],comp_size+2);
		}
		free(tvec);
		free(tdim2);
	}
	free(ca_arr[k]);
  }
  
  for (i=0; i<size; i++){
    	free_comp_list(clb_arr[i]);
        elina_linexpr0_free(expr_array[i]);
  }

  /* Free up the auxilliary structures*/
  if (destructive){
    free_array_comp_list(oa->acl);
    for(k=0; k < num_compa; k++){
	if(!disjoint_map[k]){
		opt_poly_clear(poly_a[k]);
	}
	free(poly_a[k]);
    }
    free(poly_a);
    //opt_poly_array_clear(opk,op);
  }
  free(clb_arr);
  free(rmapb);
  free(nbmapb);
  free(rmapa);
  free(disjoint_map);
  free(nbgenmapa);
  free(nblinemapa);
  free(num_vertex_a);
  free(num_vertex);
  free(expr_array);
  free_array_comp_list(aclb);
  free(comp_size_map);
  free(ca_arr);
  free(pos_con_map);
  free(counterF);
  op = destructive ? oa : opt_pk_array_alloc(NULL,NULL,maxcols);
  
  unsigned short int k1=0;
  unsigned short int bound = num_comp;
  cl = acl->head;
  for(k=0; k < num_comp; k++){
	opt_pk_t * oak = poly[k1];
	if(exc_map[k]){
		comp_list_t * tmp = cl;
		cl = cl->next;
		remove_comp_list(acl,tmp);
		unsigned short int k2;
		for(k2=k1; k2 < bound - 1; k2++){
			poly[k2] = poly[k2+1];
		}
		opt_poly_clear(oak);
		bound--; 
	}
	else{
		k1++;
		cl=cl->next;
	}
  }
  op->poly = poly;
  op->acl = acl;
  free(exc_map);
  //op->status = 0;
  return op;
}


/* ====================================================================== */
/* Assignement/Substitution by an array of linear expressions */
/* ====================================================================== */
static
opt_pk_array_t* opt_poly_asssub_linexpr_array(bool assign, bool lazy,
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
  /* Minimize the argument  */
  //if ( !lazy){
    for(k=0; k < num_compa;k++){
    	opt_poly_chernikova(man,poly_a[k],"assign linexpr array input");
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
  //}
  
  /* Choose the right technique */
  if (elina_linexpr0_array_is_linear(texpr,size)){
    op = opt_poly_asssub_linexpr_array_det(assign,man,destructive,oa,tdim,texpr,size);
    if (ob){
      opt_poly_meet(true,lazy,man,op,op,ob);
    }
  }
  else {
    op = elina_generic_asssub_linexpr_array(assign,man,destructive,oa,tdim,texpr,size,ob);
  }

  /* Minimize the result if option say so */
  opt_pk_t ** poly = op->poly;
  array_comp_list_t * acl = op->acl;
  unsigned short int num_comp = acl->size;
  if ( !lazy){
    for(k=0; k < num_comp; k++){
	//opt_pk_t * opp = poly[k];
    	opt_poly_chernikova(man,poly[k],"assign linexpr array output");
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
      op = opt_poly_asssub_linexpr_array_det(assign,man,destructive,
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
       		 opt_poly_asssub_linexpr_array(true, opk->funopt->algorithm<=0,
			      man,destructive,oa,tdim,texpr,size,ob);
  //assert(poly_check(pk,po));
  #if defined(TIMING)
 	 	record_timing(assign_linexpr_time);
  #endif
 
  return op;
}

opt_pk_array_t* opt_pk_substitute_linexpr_array(elina_manager_t* man,
				  bool destructive, opt_pk_array_t* oa,
				  elina_dim_t* tdim, elina_linexpr0_t** texpr,
				  size_t size,
				  opt_pk_array_t* ob)
{
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_SUBSTITUTE_LINEXPR_ARRAY);
  opt_pk_array_t* op;
  #if defined(TIMING)
 	 start_timing();
  #endif
  op =
    size==1 ?
    opt_poly_asssub_linexpr(false,
			opk->funopt->algorithm<=0,
			man,destructive,oa,tdim[0],texpr[0],ob)
    :
    opt_poly_asssub_linexpr_array(false,
			      opk->funopt->algorithm<=0,
			      man,destructive,oa,tdim,texpr,size,ob);
  #if defined(TIMING)
 	 	record_timing(substitute_linexpr_time);
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


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

/* ********************************************************************** */
/* opt_pk_expandfold.c: expanding and folding dimensions */
/* ********************************************************************** */


#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"
#include "opt_pk_representation.h"
#include "opt_pk_user.h"
#include "opt_pk_constructor.h"
#include "opt_pk_meetjoin.h"


/* ********************************************************************** */
/* I. Expand */
/* ********************************************************************** */

/* ---------------------------------------------------------------------- */
/* Matrix */
/* ---------------------------------------------------------------------- */

/* Expand the dimension dim of the matrix into (dimsup+1)
   dimensions, with dimsup new dimensions inserted just before
   offset. */
static
opt_matrix_t* opt_matrix_expand(opt_pk_internal_t* opk,
			bool destructive,
			opt_matrix_t* C,
			elina_dim_t dim,
			size_t offset,
			size_t dimsup)
{
  elina_dimchange_t* dimchange;
  size_t i,j,row,nb;
  unsigned short int col;
  size_t nbrows;
  unsigned short int nbcols;
  opt_numint_t** p;
  opt_matrix_t* nC;

  if (dimsup==0){
    return destructive ? C : opt_matrix_copy(C);
  }
  nbrows = C->nbrows;
  nbcols = C->nbcolumns;
  col = opk->dec + dim;
  /* Count the number of constraints to duplicate */
  nb=0;
  p = C->p;
  for (i=0; i<nbrows; i++){
    if (p[i][col]!=0)
      nb++;
  }
  /* Redimension matrix */
  dimchange = elina_dimchange_alloc(dimsup,0);
  for (i=0;i<dimsup;i++){
    dimchange->dim[i]=offset;
  }
  nC = opt_matrix_add_dimensions(opk,destructive,C,dimchange, false);
  elina_dimchange_free(dimchange);
  opt_matrix_resize_rows(nC,nbrows+nb*dimsup);
  if (nb==0)
    return nC;

  /* Duplicate constraints */
  p = nC->p;
  row = nbrows;
  for (i=0; i<nbrows; i++){
    if (p[i][col]!=0){
      for (j=offset;j < offset+dimsup; j++){
	
	opt_vector_copy(p[row],
		    p[i],
		    nbcols+dimsup);
	p[row][opk->dec+j] = p[row][col];
	p[row][col] = 0;
	row++;
      }
    }
  }
  nC->_sorted = false;
  return nC;
}

/* ---------------------------------------------------------------------- */
/* Polyhedra Expand */
/* ---------------------------------------------------------------------- */

opt_pk_array_t* opt_pk_expand(elina_manager_t* man,
		bool destructive, opt_pk_array_t* oa,
		elina_dim_t dim, size_t dimsup)
{
  size_t intdimsup,realdimsup;
  opt_pk_array_t* op;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_EXPAND);
  if(oa->is_bottom || !oa->acl){
	if(destructive){
		oa->maxcols = oa->maxcols + dimsup;
		return oa;
	}
	else{
		op = opt_pk_bottom(man,oa->maxcols + dimsup-opk->dec,0);
		return op;
	}
  }
  unsigned short int maxcols = oa->maxcols;
  unsigned short int var = dim+opk->dec;
  unsigned short int num_var = maxcols - opk->dec;
  opt_pk_internal_realloc_lazy(opk,num_var+dimsup);
  unsigned short int j,k,nmaxcols;
  man->result.flag_best = man->result.flag_exact = true;   

  //if (dim<oa->intdim){
  //  intdimsup = dimsup;
  //  realdimsup = 0;
  //} else {
  //  intdimsup = 0;
  //  realdimsup = dimsup;
  //}
  nmaxcols = oa->maxcols + dimsup;
  //nrealdim = oa->realdim + realdimsup;

  if (dimsup==0){
    return (destructive ? oa : opt_pk_copy(man,oa));
  }
  array_comp_list_t *acla = oa->acl;
  unsigned short int num_compa = acla->size;
  if (destructive){
    op = oa;
    op->maxcols+=dimsup;
  }
  else {
   op = (opt_pk_array_t *)malloc(sizeof(opt_pk_array_t));
   op->maxcols = nmaxcols;
   op->acl = NULL; 
  }
  opt_pk_t ** poly_a = oa->poly;
  /* Get the constraints system, and possibly minimize */
    for(k=0; k<num_compa; k++){
	opt_pk_t * oak = poly_a[k];
	if (opk->funopt->algorithm<0){
    		opt_poly_obtain_C(man,oak,"expand operation");
	}
  	else{
    		opt_poly_chernikova(man,oak,"expand operation");
	}
  	if (opk->exn){
    	    opk->exn = ELINA_EXC_NONE;
    	    if (!oak->C){
     		man->result.flag_best = man->result.flag_exact = false;   
      		opt_poly_set_top(opk,op);
      		return op;
    	    }
    	   /* We can still proceed, although it is likely 
       		that the problem is only delayed
    	   */
  	}
  	/* if empty, return empty */
  	if (!oak->C){
    	    opt_poly_set_bottom(opk,op);
            return op;
        }
     }
  array_comp_list_t *tcla = copy_array_comp_list(acla);
  array_comp_list_t *acl = destructive ? acla : copy_array_comp_list(tcla);
  free_array_comp_list(tcla);
  comp_list_t * ecl = find(acl,var);
  if(ecl==NULL){
	if(!destructive){
		free(op);
		op = opt_pk_copy(man,oa);
		op->maxcols = nmaxcols;
	}
	return op;
  }
  for(j=maxcols; j < nmaxcols; j++){
	insert_comp(ecl,j);
  }
  short int ind = find_index(acl,ecl);
  opt_pk_t ** poly = destructive ? poly_a : (opt_pk_t **)malloc(num_compa*sizeof(opt_pk_t *));
  comp_list_t * cl = acl->head;
  for(k=0; k < num_compa; k++){
	if(k==ind){
  		opt_pk_t * oak = poly_a[k]; 
  		/* Prepare resulting matrix */
  		if (destructive){
			 if (oak->F){ opt_matrix_free(oak->F); oak->F = NULL; }
   			 if (oak->satF){ opt_satmat_free(oak->satF); oak->satF = NULL; }
   			 if (oak->satC){ opt_satmat_free(oak->satC); oak->satC = NULL; }
   		    	oak->nbeq  = 0;
    		    	oak->status &= ~opt_pk_status_consgauss & ~opt_pk_status_gengauss & ~opt_pk_status_minimaleps;
  		}
		else{
			unsigned short int comp_size = cl->size;
			poly[k] = opt_poly_alloc(comp_size,0);
		}
		unsigned short int ndim=0,l;
		unsigned short int * ca = to_sorted_array(cl,nmaxcols);
		for(l=0; l < cl->size; l++){
			if(ca[l]==(dim+opk->dec)){
				ndim = l;
				break;
			}
		}
		free(ca);
  		poly[k]->C = opt_matrix_expand(opk, destructive, oak->C, 
					ndim,cl->size - dimsup ,dimsup);
  		/* Minimize the result */
  		if (opk->funopt->algorithm>0){
                  opt_poly_minimize(man, poly[k]);
                  if (opk->exn) {
                    opk->exn = ELINA_EXC_NONE;
                    if (!poly[k]->C) {
                      man->result.flag_best = man->result.flag_exact = false;
                      opt_poly_set_top(opk, op);
                      return op;
                    }
    			}
  		}
	}
	else if(!destructive){
		poly[k] = opt_poly_alloc(cl->size,0);
		opt_poly_copy(poly[k],poly_a[k]);
	}
	cl = cl->next;
  }
   
   if(!destructive){
	op->poly = poly;
	op->acl = acl;
	op->is_bottom = false;
   }
  
  //assert(poly_check(pk,po));
  return op;
}

/* ********************************************************************** */
/* II. Fold */
/* ********************************************************************** */

/* ---------------------------------------------------------------------- */
/* Matrix */
/* ---------------------------------------------------------------------- */

/* Fold the last dimsup dimensions with dimension dim (not in the last dimsup
   ones) in the matrix */


/* the array tdim is assumed to be sorted */

void opt_pk_fold_new_comp(opt_pk_internal_t *opk, bool destructive, 
			  opt_pk_t ** poly, array_comp_list_t *acl, 
			  elina_dim_t *tdim, size_t size){
	size_t i;
        elina_dim_t dim = tdim[0];
        comp_list_t *cl = create_comp_list();
	insert_comp(cl,tdim[0] + opk->dec);
	insert_comp_list_tail(acl,cl);
	unsigned short int comp_size = acl->size;
	poly = (opt_pk_t **)realloc(poly,(comp_size+1)*sizeof(opt_pk_t*));
	opt_pk_t * factor = poly[comp_size];
	size_t fold_size = 0;
	unsigned short int k;
	for(k=0; k < comp_size; k++){
	        opt_pk_t * oak = poly[k];
		fold_size = fold_size + oak->F->nbrows;
	}
	factor->F = opt_matrix_alloc(fold_size,3,false);
        cl = acl->head;
}

static
opt_matrix_t* opt_matrix_fold_same_comp(opt_pk_internal_t* opk,
		      bool destructive,
		      opt_matrix_t* F, elina_dim_t dim,
		      elina_dim_t* tdim, size_t size)
{
  opt_matrix_t* nF;
  size_t i,j,row,col;
  size_t nbrows, nbcols, dimsup;
  elina_dimchange_t* dimchange;
  //printf("F\n");
  //matrix_fprint(stdout,F);
  dimsup = size;
  if (dimsup==0){
    return destructive ? F : opt_matrix_copy(F);
  }
  nbrows = F->nbrows;
  nbcols = F->nbcolumns;
  col = opk->dec + dim;

  nF = opt_matrix_alloc(size*nbrows+nbrows,
		     nbcols - dimsup,
		     false );
  dimchange = elina_dimchange_alloc(dimsup,0);
  for (i=0;i<dimsup;i++){
    dimchange->dim[i]=tdim[i];
  }
  row = 0;
  for(i=0; i<nbrows; i++){
    opt_vector_remove_dimensions(opk,nF->p[row],F->p[i],nbcols,
			     dimchange);
    opt_vector_normalize(opk,nF->p[row],nbcols-dimsup);
    row++;
    for (j=0;j<dimsup;j++){
      if (F->p[i][col] != F->p[i][opk->dec+tdim[j]]){
	  opt_vector_remove_dimensions(opk,
				 nF->p[row],F->p[i],nbcols,
				 dimchange);
	  nF->p[row][col]= F->p[i][opk->dec+tdim[j]];
	  opt_vector_normalize(opk,nF->p[row],nbcols-dimsup);
	  row++;
      }
    }
  }
  nF->nbrows = row;
  nF->_sorted = false;
  if (destructive){
    opt_matrix_free(F);
  }
  elina_dimchange_free(dimchange);
   //printf("nF\n");
  //matrix_fprint(stdout,nF);
  return nF;
}


static
opt_matrix_t* opt_matrix_fold_diff_comp(opt_pk_internal_t* opk,
		      bool destructive,
		      opt_matrix_t* F, opt_matrix_t ** fold_val,
		      elina_dim_t dim, size_t size)
{
  size_t i,i1,row;
  unsigned short int j,j1,col, dimsup, nbcols;
  size_t nbrows;
  dimsup = size-1;
  if (dimsup==0){
    return destructive ? F : opt_matrix_copy(F);
  }
 // printf("input\n");
  //opt_matrix_fprint(stdout,F);
  opt_matrix_t* nF;

  //size_t * num_vertex = (size_t *)calloc(size, sizeof(size_t));
  //num_vertex[0] = opt_generator_rearrange(F,NULL);
  //size_t num_vertex1 = num_vertex[0];
  size_t fold_size = 0;
  for(j=0; j < dimsup; j++){
	//num_vertex[j+1] = opt_generator_rearrange(fold_val[j],NULL);
	//num_vertex1 *= num_vertex[j+1];
	if(fold_val[j]){
		fold_size += fold_val[j]->nbrows;
	}
  }
  nbrows = F->nbrows;
  nbcols = F->nbcolumns;
  
  nF = opt_matrix_alloc(2*nbrows*fold_size + nbrows,
		     nbcols,
		     false);
  //size_t count = cartesian_product_vertices_fold(F,fold_val,num_vertex,size);
  

 
  col = opk->dec + dim;
  row = 0;
  for(i=0; i<nbrows; i++){
    opt_numint_t * si = F->p[i];
    opt_vector_copy(nF->p[row],si,nbcols);
    row++;
    for (j=0;j<dimsup;j++){
      opt_matrix_t *fv = fold_val[j];
	if(!fold_val[j]){
		continue;
	}
      for(i1=0; i1 < fv->nbrows; i1++){
	opt_numint_t * pi = fv->p[i1];
	opt_numint_t * di = nF->p[row];
      	//if (si[col]!= pi[2]){
		bool cond1 = is_vertex(pi) && is_vertex(si);
		bool cond2 = (is_line(pi) && is_line(si)) || (is_ray(pi) && is_ray(si)); 
		bool cond3 = is_line(pi) || is_ray(pi);
		if(cond1){
			opt_numint_t lcm = opt_numint_lcm(pi[1],si[1]);
			opt_numint_t lcm1 = lcm/pi[1];
			opt_numint_t lcm2 = lcm/si[1];
			opt_vector_copy(di,si,nbcols);
			di[1] = lcm;
			if(lcm2!=1){
				for(j1=2; j1 < nbcols; j1++){
					di[j1] = lcm2*si[j1];
				}
			}
			di[col] = lcm1*pi[2];
			opt_vector_normalize(opk,di,nbcols);
			row++;
		}
		
		else if(cond2){
			//opt_vector_copy(di,si,nbcols);
			di[0] = pi[0];
			di[col] = pi[2];
			opt_vector_normalize(opk,di,nbcols);
			row++;
			opt_vector_copy(nF->p[row],si,nbcols);
			nF->p[row][col] = 0;
			opt_vector_normalize(opk,nF->p[row],nbcols);			
    			row++;
		}
		else if(cond3){
			di[0] = pi[0];
			di[col] = pi[2];
			opt_vector_normalize(opk,di,nbcols);
			row++;
		}
        //}
      }
    }
  }
 
  nF->nbrows = row;
  nF->_sorted = false;
  remove_common_gen(opk,nF,0);
  if (destructive){
    opt_matrix_free(F);
  }
 // printf("output\n");
  //opt_matrix_fprint(stdout,nF);
  return nF;
}



opt_pk_array_t* opt_pk_fold(elina_manager_t* man,
	      bool destructive, opt_pk_array_t* oa,
	      elina_dim_t* tdim, size_t size)
{
 
  short int dimsup;
  opt_pk_array_t* op;
  opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_FOLD);
  man->result.flag_best = man->result.flag_exact = true;   
  unsigned short int j,k,k1,k2,l,maxcols, num_compa;
  dimsup = size - 1;
  
  maxcols = oa->maxcols;
  size_t i;
  array_comp_list_t *acla = oa->acl;
  if(oa->is_bottom || !oa->acl){
	man->result.flag_best = man->result.flag_exact = true;   
	if(destructive){
		oa->maxcols -= dimsup;
		return oa;
	}
	else{
          return opt_pk_bottom(man, oa->maxcols - dimsup, 0);
        }
  }
  num_compa = acla->size;
  opt_pk_t ** poly_a = oa->poly;
  if (destructive){
    op = oa;
    op->maxcols -= dimsup;
  }
  else {
    op = (opt_pk_array_t *)malloc(sizeof(opt_pk_array_t));
    op->maxcols = maxcols - dimsup;
    op->acl = NULL;
  }
  if(oa->acl->size==0){
	if(!destructive){
		opt_poly_set_top(opk,op);
	}
	return op;
  }
  
  
  for(k=0; k < num_compa; k++){
	  opt_pk_t * oak = poly_a[k];
	  if (opk->funopt->algorithm<0){
	    opt_poly_obtain_F(man,oak,"fold operation");
	  }
	  else{
	    opt_poly_chernikova(man,oak,"fold operation");
	  }
  		
	  if (opk->exn){
	    opk->exn = ELINA_EXC_NONE;
	    if (!oak->F){
	      man->result.flag_best = man->result.flag_exact = false;   
	      opt_poly_set_top(opk,op);
	      return op;
	    }
	  }
	  /* if empty, return empty */
	  if (!oak->F){
	    man->result.flag_best = man->result.flag_exact = true;   
	    opt_poly_set_bottom(opk,op);
	    return op;
	  }
  }
  array_comp_list_t *tacl = copy_array_comp_list(acla);
  array_comp_list_t *acl = destructive ? acla : copy_array_comp_list(tacl);
  free_array_comp_list(tacl);
  opt_pk_t ** poly = destructive? poly_a : (opt_pk_t **)malloc(num_compa*sizeof(opt_pk_t*));
 
  /*******************
	take cartesian product of factors corresponding to blocks containing
	tdim[1..size-1] and the factor corresponding to block containing
	tdim[0]
  ********************/
  /*size_t * num_vertex_a = (size_t *)calloc(num_compa,sizeof(size_t));
  char * map = (char *)calloc(num_compa, sizeof(char));
  comp_list_t * fl = create_comp_list();
  for(l=0; l < dimsup; l++){
	insert_comp(cl,tdim[l+1]+opk->dec);
  }
  comp_list_t * cl = acl->head;
  for(k=0; k < num_compa; k++){
	opt_pk_t * oak = poly[k];
	num_vertex_a[k] = opt_generator_rearrange(oak->F,oak->satF);
	if(is_disjoint(cl,fl,maxcols)){
		map[k] = 1;
	}
  }
  comp_list_t * res = find(acl,tdim[0]+opk->dec);
  opt_pk_t * src = opt_poly_alloc(res->size,0);
  cartesian_product_vertices_fold();*/
  unsigned short int var = tdim[0] + opk->dec;
  comp_list_t *fcl = find(acl,var);
  bool null_flag = false;
  if(fcl==NULL){
	// if there is no factor containing tdim[0];
	fcl = create_comp_list();
	insert_comp(fcl,var);
	insert_comp_list_tail(acl,fcl);
	//num_compa = acl->size;
	//poly = (opt_pk_t **)realloc(poly,acl->size*sizeof(opt_pk_t*));
	//poly[num_compa] = opt_poly_alloc(1,0);
	//poly[num_compa]->F = opt_matrix_alloc(1,3,true);
	//poly[num_compa]->F->p[0][0] = 1;
	//poly[num_compa]->F->p[0][1] = 1;
	//opt_pk_fold_new_comp(opk,true,poly,acl,tdim,size);
	//op->poly = poly;
	//op->acl = acl;
	//return op;
	null_flag = true;
  }
  
  unsigned short int ndim=tdim[0]+opk->dec;
  unsigned short int *fca = to_sorted_array(fcl,maxcols);
  
  for(l=0; l < fcl->size; l++){
	if(fca[l]==ndim){
		ndim = l;
		break;
	}
  }
  free(fca);
  opt_matrix_t ** fold_val = (opt_matrix_t **)malloc(dimsup * sizeof(opt_matrix_t *));
  for(k=0; k < dimsup; k++){
	fold_val[k]= NULL;
  }
  //opt_numint_t ** fold_val = (opt_numint_t **)malloc(dimsup*sizeof(opt_numint_t *));
  //opt_numint_t ** coeff_val = (opt_numint_t **)malloc(dimsup*sizeof(opt_numint_t *));
  k=0, k2 =0;
  comp_list_t * cl = acl->head;
  bool flag1 = false, flag2 = false;
  opt_matrix_t * tmp;
  while(k < num_compa){
	 unsigned short int * ca = to_sorted_array(cl,maxcols);
	 unsigned short int comp_size = cl->size;
	 elina_dim_t * tdimk = (elina_dim_t *)calloc(comp_size, sizeof(elina_dim_t));
	 unsigned short int dim_size = 0;
	 bool flag = true;
	 j=0, l=0;
	 opt_pk_t * oak = poly_a[k];
	 opt_matrix_t * F = oak->F;
	 size_t nbrows = F->nbrows;
 	 while((j < comp_size) && (l < dimsup) && flag){
		unsigned short int num = ca[j];
		unsigned short int var = tdim[l+1] + opk->dec;
		if(num==var){
			fold_val[l] = opt_matrix_alloc(nbrows,3,false);
			for(i=0; i < nbrows; i++){
				opt_numint_t * pi = F->p[i];
				opt_numint_t * di = fold_val[l]->p[i];
				di[0] = pi[0];
				di[1] = pi[1];
				di[2] = pi[j+opk->dec];
			}
			remove_comp(cl,var);
			tdimk[dim_size] = j;
			if(cl->size==0){
				flag = false;
			}
			j++;
			l++;
			dim_size++;
		}
		else if(num< var){
			j++;
		}
		else{
			l++;
		}
  	}
	if(!flag){
		comp_list_t * tcl = cl->next;
		remove_comp_list(acl,cl);
		cl = tcl;
		free(tdimk);
		if(destructive){
			opt_pk_t * tpoly = poly_a[k];
			for(k1=k; k1 < num_compa - 1; k1++){
				poly_a[k1] = poly_a[k1+1];
			}
			opt_poly_clear(tpoly);
			free(tpoly);
			num_compa--;	
		}
		else{
			k++;
		}
		free(ca);
		continue;
	}
	else if(dim_size){
		if(cl==fcl){
			if(destructive){
				poly_a[k]->F = opt_matrix_fold_same_comp(opk,true, poly_a[k]->F, ndim, tdimk, dim_size);
			}
			else{
				poly[k2] = opt_poly_alloc(oak->intdim - dim_size, oak->realdim);
				poly[k2]->F = opt_matrix_fold_same_comp(opk,false, poly_a[k]->F, ndim, tdimk, dim_size);
				// need to free poly[k2]->F afterwards
                                flag2 = true;
                                tmp = poly[k2]->F;
                        }
			if(dim_size==dimsup){
				// if all the folded dimensions are in the same block as tdim[0], 
				// then no need to apply fold_diff_comp
				flag1 = true;
			}
		}
		else{
			elina_dimchange_t dimchange1;
			dimchange1.dim = tdimk;
			dimchange1.intdim = dim_size;
			dimchange1.realdim = 0;
			if(destructive){
				if(oak->C){
				   opt_matrix_free(oak->C);
				   oak->C = NULL;
				}
				if(oak->satC){
				   opt_satmat_free(oak->satC);
				   oak->satC = NULL;
				}
				if(oak->satF){
				   opt_satmat_free(oak->satF);
				   oak->satF = NULL;
				}
				poly_a[k]->F = opt_matrix_remove_dimensions(opk,false,oak->F,&dimchange1);
				poly_a[k]->intdim = comp_size - dim_size;
				opt_matrix_free(oak->F);
				oak->F = NULL;
				poly_a[k]->nbline = oak->nbline;
				poly_a[k]->nbeq = 0;
				if (opk->funopt->algorithm>0){
				    opt_poly_chernikova(man,poly_a[k],"of the result");
				    if (opk->exn){
					opk->exn = ELINA_EXC_NONE;
				    }
				}
			}
			else{
				poly[k2] = opt_poly_alloc(oak->intdim - dim_size,
						              oak->realdim);
				poly[k2]->F = opt_matrix_remove_dimensions(opk,false,oak->F,&dimchange1);
				poly[k2]->nbeq = 0;
				poly[k2]->nbline = oak->nbline;
				if (opk->funopt->algorithm>0){
				    opt_poly_chernikova(man,poly[k2],"of the result");
				    if (opk->exn){
					opk->exn = ELINA_EXC_NONE;
				    }
				}
			}
		}			
	}
	else if(!dim_size){
		if(!destructive){
			poly[k2] = opt_poly_alloc(oak->intdim,oak->realdim);
			if(cl!=fcl){
				opt_poly_copy(poly[k2],poly_a[k]);
				poly[k2]->is_minimized = poly_a[k]->is_minimized;
			}
			else{
				poly[k2]->F = opt_matrix_copy(poly_a[k]->F);
			}
		}
	}
	k2++;
	k++;
	free(ca);
	free(tdimk);
	cl = cl->next;
  }
  
  if(!destructive){
	unsigned short int k3 = null_flag ? k2+1: k2;
	poly = (opt_pk_t **)realloc(poly,k3*sizeof(opt_pk_t*));
  }
  short int ind = find_index(acl,fcl);
 
  if(null_flag){
	poly[ind] = opt_poly_alloc(1,0);
	poly[ind]->F = opt_matrix_alloc(1,3,true);
	poly[ind]->F->p[0][0] = 1;
	poly[ind]->F->p[0][1] = 1; 
  }
   opt_pk_t * src = poly[ind];
  if(!flag1){
	  
	
	  /* Prepare resulting matrix */
	  if(destructive){
		if(src->C){
			opt_matrix_free(src->C);
			src->C = NULL;
	 	}
		if(src->satC){
			opt_satmat_free(src->satC);
			src->satC = NULL;
		}
		if(src->satF){
			opt_satmat_free(src->satF);
			src->satF = NULL;
		}
		src->nbeq = 0;
		src->status &= ~opt_pk_status_consgauss & ~opt_pk_status_gengauss & ~opt_pk_status_minimaleps;
	  }
	  opt_matrix_t * tmp2 = poly[ind]->F;
	 
	  poly[ind]->F = opt_matrix_fold_diff_comp(opk, destructive, src->F, fold_val,
			      ndim, size);
	  if(!destructive){
		opt_matrix_free(tmp2);
	  }
  }
	
  /* Minimize the result */
  if (opk->funopt->algorithm>0){
    opt_poly_chernikova(man,poly[ind],"of the result");
    if (opk->exn){
      opk->exn = ELINA_EXC_NONE;
      if (!poly[ind]->C){
	man->result.flag_best = man->result.flag_exact = false;   
	opt_poly_set_top(opk,op);
	return op;
      }
    }
  }


  man->result.flag_best = (dimsup==0);
  man->result.flag_exact = (dimsup==0);
  //assert(poly_check(pk,po));
  for(l=0; l < dimsup; l++){
	if(fold_val[l])
	opt_matrix_free(fold_val[l]);
  }
  free(fold_val);
  if(!destructive){
	op->poly = poly;
	op->acl = acl;
	op->is_bottom = false;
	//if(flag2 && !flag1){
	//	opt_matrix_free(tmp);
	//}
  }
  return op;
}

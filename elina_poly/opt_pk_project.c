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


#include "opt_pk_config.h"
#include "opt_pk.h"
#include "opt_pk_vector.h"
#include "opt_pk_project.h"
#include "opt_pk_representation.h"
#include "opt_pk_constructor.h"
#include "opt_pk_meetjoin.h"

int select_variable_gauss(opt_pk_internal_t *opk, opt_matrix_t *oc,  opt_numint_t *cons, size_t nbeq, elina_dim_t *tdim, size_t size, char *rmap){
	int res = -1;
	size_t nbcons = oc->nbrows;
	size_t nbcolumns = oc->nbcolumns;
	size_t i,j;
	size_t * p = (size_t *)calloc(size,sizeof(size_t));
	size_t * m = (size_t *)calloc(size,sizeof(size_t));
	for(i=nbeq; i < nbcons; i++){
		opt_numint_t *pi = oc->p[i];
		for(j=0; j < size; j++){
			elina_dim_t dim = opk->dec + tdim[j];
			if(!cons[dim] || rmap[j]){
				continue;
			}
			if(pi[dim]<0){
				m[j]++;
			}
			else if(pi[dim]>0){
				p[j]++;
			}
		}	
	}
	int maxadd = -nbcons*nbcons;
	//printf("maxadd: %d");
	int ind = -1;
	for(j=0; j < size; j++){
		elina_dim_t dim = opk->dec + tdim[j];
		if(!rmap[j] && cons[dim]){
			int nbadd = p[j]*m[j] - p[j] -m[j];
			if(nbadd >= maxadd){
				maxadd = nbadd;
				res = dim; 
				ind = j;
			}
		}
	}

	if(res!=-1){
		rmap[ind] = 1;
	}
	free(p);
	free(m);
	//printf("res: %d\n",res);	
	return res;
}

bool gauss_project(opt_pk_internal_t *opk, opt_pk_t *o, elina_dim_t *tdim, size_t *tsize, size_t * proj ){
	size_t nbeq = o->nbeq;
	opt_matrix_t * oc = o->C;
	size_t nbrows = oc->nbrows;
	size_t nbcolumns = oc->nbcolumns; 
	size_t i, j, l,  c2 = 0;
	int s;
	size_t size = *tsize;
	char *rmap = (char *)calloc(size,sizeof(char));
	opt_numint_t **p = oc->p;
	for(i=0; i < nbeq; i++){
		bool exch = false;
		opt_numint_t *pi = p[i];
		//s = -1;
		opt_numint_t * npi;
		//for(j = 0; j < size; j++){
		//	elina_dim_t dim = opk->dec + tdim[j];
		//	if(pi[dim]){
		//		s = dim;
		//		//printf("s: %d i: %d pi[dim]: %lld\n",s,i,pi[dim]);
		//		rmap[j] = 1;
		//		npi = opt_vector_neg(opk,pi,nbcolumns);
		//		if(pi[s] < 0){
		//			opt_numint_t *tmp = npi;
		//			npi = pi;
		//			pi = tmp; 
		//			exch  = true;
		//		}
		//		break;
		//	}
		//}
		s = select_variable_gauss(opk,oc,pi,nbeq,tdim,size,rmap); 
		if(s==-1){	
			opt_matrix_exch_rows(oc,i,c2);
			c2++;
		}
		else{
			npi = opt_vector_neg(opk,pi,nbcolumns);
			if(pi[s] < 0){
				opt_numint_t *tmp = npi;
				npi = pi;
				pi = tmp; 
				exch  = true;	
			}
			l=0;
			while(l < nbrows){
				if(l==i){
					l++;
					continue;
				}
				opt_numint_t * pl = p[l];
				if(pl[s]==0){
					l++;
					continue;
				}
				else{
					if(pl[s] > 0){
						opt_vector_combine(opk, pl, npi, pl, s, nbcolumns, true);
					}
					else{
						opt_vector_combine(opk, pl, pi, pl, s, nbcolumns,true);
					}
					opt_vector_simplify(opk, pl, nbcolumns);
					if(opt_vector_is_null_expr(opk, pl, nbcolumns)){
						/*** Null Inequality ***/
						if(pl[0]){
							if(pl[1] < 0){
								if(exch){
									free(pi);
								}
								else{
									free(npi);
								}
								return false;
							}
							else{
								nbrows--;
								opt_matrix_exch_rows(oc,l,nbrows);
							}
						}
						/**** Null Equality ***/
						else{
							if(pl[1]!= 0){
								if(exch){
									free(pi);
								}
								else{
									free(npi);
								}
								return false;
							}else{
								nbrows--;
								opt_matrix_exch_rows(oc,l,nbrows);
							}
						}
					}
					else{
						l++;
					}
					
				}
			
			}
			if(exch){
				free(pi);
			}
			else{
				free(npi);
			}
		}
	}
	
	for(i = c2; i < nbeq; i++){
		nbrows--;
		opt_matrix_exch_rows(oc,i,nbrows);
	}
	oc->nbrows = nbrows;
	
	size_t k = 0;
	for(i = 0; i < size; i++){
		if(!rmap[i]){
			tdim[k] = tdim[i]; 
			k++;
		}
	}
	*tsize = k;
	*proj = c2;
	o->nbeq = c2;
	free(rmap);
	return true;
}



void select_variable(opt_pk_internal_t * opk, elina_dim_t *tdim, size_t size, opt_matrix_t * oc,  int *nbadd, int *maxadd){
	//printf("Start %d %d\n",oc->nbrows,oc->nbcolumns);
	//opt_matrix_fprint(stdout,oc);
	//fflush(stdout);
	size_t i, j, k;
	size_t * p = (size_t *)calloc(size,sizeof(size_t));
	size_t * m = (size_t *)calloc(size,sizeof(size_t));
	size_t nbcolumns = oc->nbcolumns;
	size_t nbcons = oc->nbrows;
	//printf("start %d\n",nbcons);
	//fflush(stdout);
	for(i = 0; i < nbcons; i++){
		opt_numint_t *pi = oc->p[i];
		//printf("pi: %p %d %d\n",pi,size,tdim[0]);
		//fflush(stdout);
		for(j = 0; j < size; j++){
			//printf("pj: %d\n",tdim[j]);
			//fflush(stdout);
			opt_numint_t pj = pi[opk->dec + tdim[j]];
			
			if(pj > 0){
				p[j]++;
			}
			else if(pj < 0){
				m[j]++;
			}
		}
		//for(j=0; j < size; j++)
		//printf("%d %d %d\n",j,p[j],m[j]);
	}
	//printf("AA\n");
	//fflush(stdout);
	size_t ind = 0;
	*nbadd = nbcons*nbcons;
	*maxadd = nbcons*nbcons;
	for(i = 0; i < size; i++){
		int tmp = p[i]*m[i] - p[i] - m[i];
		if(tmp < *nbadd){
			*maxadd = p[i]*m[i];
			*nbadd = tmp;
			ind = i;
		}
	} 
	
	/**********
		Always put chosen variable at last position
	***********/
	size_t tmp = tdim[size-1];
	tdim[size-1] = tdim[ind];
	tdim[ind] = tmp;
	free(p);
	free(m);
	
}

bool fourier_motzkin(opt_pk_internal_t *opk, elina_dim_t dim, opt_matrix_t * oc){
	size_t nbcons = oc->nbrows;
	size_t nbcolumns = oc->nbcolumns;	
	size_t *ocm = (size_t*)calloc(nbcons,sizeof(size_t));	
	size_t *ocp = (size_t *)calloc(nbcons, sizeof(size_t));
	size_t *ocn = (size_t *)calloc(nbcons, sizeof(size_t));
	size_t nbcons1 = 0, nbcons2 = 0, nbcons3 = 0;
	/****
		Result matrix will add at most nbadd inequalities
	*****/
	//opt_matrix_t *ocn = opt_matrix_alloc(nbcons + nbadd, nbcolumns, false);
	size_t i, j;
	for(i= 0; i < nbcons; i++){
		opt_numint_t * pi = oc->p[i];
		if(pi[dim] > 0){
			ocp[nbcons2] = i;
			nbcons2++;
		}
		else if(pi[dim] < 0){
			ocm[nbcons1] = i;
			nbcons1++;
		}
		else{
			ocn[nbcons3] = i;
			nbcons3++;
		}
	}
	size_t s = 0;
	bool flag = false;
	for(i = 0; i < nbcons1; i++){
	   size_t mind = ocm[i];
	   opt_numint_t * pi = oc->p[mind];
	   opt_numint_t fi = -pi[dim];
	   for(j = 0; j < nbcons2; j++){
		size_t pind = ocp[j];
		opt_numint_t * pj = oc->p[pind];
		opt_numint_t fj = pj[dim];
		opt_numint_t * tmp1, *tmp2;
		//printf("start\n");
		//opt_vector_print(pi,nbcolumns);
		//opt_vector_print(pj,nbcolumns);
		if(opt_numint_abs(fi)!=opt_numint_abs(fj)){
		   tmp1 = opt_vector_scalar_product(opk, pi,fj, nbcolumns);
		   tmp2 = opt_vector_scalar_product(opk, pj,fi, nbcolumns);
		}else{
		   tmp1 = opt_vector_alloc(nbcolumns);
		   tmp2 = opt_vector_alloc(nbcolumns);
		   opt_vector_copy(tmp1,pi,nbcolumns);
		   opt_vector_copy(tmp2,pj,nbcolumns);
		}
		opt_vector_sum(opk, tmp1, tmp2, nbcolumns);
		//opt_vector_print(tmp1,nbcolumns);
		opt_vector_simplify(opk, tmp1, nbcolumns);
		//opt_vector_print(tmp1,nbcolumns);
		if(!opt_vector_is_null_expr(opk,tmp1,nbcolumns)){
			opt_numint_t * dst = oc->p[s+nbcons];
			opt_vector_copy(dst,tmp1,nbcolumns);
			s++;
		}
		else if(!flag && opt_vector_is_positivity_constraint(opk,tmp1,nbcolumns)){
			opt_numint_t * dst = oc->p[s+nbcons];
			opt_vector_copy(dst,tmp1,nbcolumns);
			flag = true;
			s++;
		}
		else if(tmp1[1] < 0){
			//printf("bottom here\n");
			//opt_vector_print(tmp1,nbcolumns);
			free(tmp1);
			free(tmp2);
			free(ocm);
			free(ocp);
			free(ocn);
			return false;
		}
		free(tmp1);
		free(tmp2);
	   }
	}
	for(i=0; i <nbcons3; i++){
		size_t nind = ocn[i];
		opt_matrix_exch_rows(oc,i,nind);
	}
	for(i=0; i < s; i++ ){
		opt_matrix_exch_rows(oc,nbcons3+i,nbcons+i);
	}
	oc->nbrows = s + nbcons3;
	free(ocm);
	free(ocp);
	free(ocn);
	return true;	
}


void opt_poly_projectforget_array(bool project,
			      elina_manager_t *man,	
			      opt_pk_t* op, opt_pk_t* oa, 
			      elina_dim_t* tdim, size_t size, bool destructive){
	
	opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_FORGET_ARRAY);
	opt_matrix_t *ocp = destructive ? oa->C : opt_matrix_copy(oa->C);
	if(!destructive){
		op->nbeq = oa->nbeq;
	}
	op->C = ocp;
	size_t psize = size;
	size_t i,j;
	elina_dim_t *pdim = (elina_dim_t *)calloc(psize,sizeof(elina_dim_t));
	for(i=0; i < psize; i++){
		pdim[i] = tdim[i];
	}
	size_t nbcons = ocp->nbrows;
	size_t nbcolumns = ocp->nbcolumns;
	opt_matrix_rearrange(ocp,op->nbeq);
	size_t proj=0;
	bool res = gauss_project(opk,op,tdim,&size,&proj);
	//printf("After Gauss\n");
	
	if(!res){
      		man->result.flag_best = man->result.flag_exact = true;
      		free(op->C);
		op->C=NULL;
      		return;
	}
	//quasi_removal(opk,op);
	nbcons = ocp->nbrows;
	int limit = nbcons;	
	int  nbadd;
	int maxadd;
	while(size > 0){

		select_variable(opk,tdim,size,ocp, &nbadd, &maxadd);
		size_t dim = tdim[size-1] + opk->dec;
		nbcons = ocp->nbrows; 
		//printf("nbcons: %d nbadd: %d dim: %d size: %d maxrows: %d maxadd: %d\n",nbcons,nbadd,dim,size,ocp->_maxrows,maxadd);
		opt_matrix_resize_rows(ocp,ocp->nbrows + maxadd);
		ocp->nbrows = nbcons;
		res = fourier_motzkin(opk,dim,ocp);
		if(!res){
			opt_matrix_fprint(stdout,op->C);
			man->result.flag_best = man->result.flag_exact = true;
    			free(op->C);
			op->C = NULL;
			return;
		}
		quasi_removal(opk,op);
		nbcons = ocp->nbrows; 
		size--;
		if(size>0){
			select_variable(opk,tdim,size,ocp, &nbadd, &maxadd);
			if(nbcons+nbadd > limit){
				//opt_matrix_fprint(stdout,op->C);
				//remove_redundancy(opk,op);
				if(op->C==NULL){
					return;
				}
				limit = ocp->nbrows;
				//printf("select variable %d %d %d\n",ocp->_maxrows,ocp->nbrows,nbadd);
				//fflush(stdout);
				//select_variable(opk,tdim,size,ocp, &nbadd, &maxadd);
			}
		}
	}
	
	if(project){
		ocp = op->C;
		nbcons = ocp->nbrows;
		opt_matrix_resize_rows(ocp,nbcons + psize);
		/**************************
			for each projected variable xi, add xi=0 to the projection
		***************************/
		for(i=0;i < psize;i++){
			size_t dim = pdim[i] + opk->dec;
			opt_numint_t * dst = ocp->p[nbcons];
			dst[0] = 0;
			for(j=1; j < nbcolumns; j++){
				dst[j] = 0;
			}
			dst[dim] = 1;
			nbcons++;
		}
		ocp->nbrows = nbcons;
		op->nbeq = op->nbeq + psize;
	}
	free(pdim);
	
	return;
}

opt_pk_array_t* opt_pk_forget_array(elina_manager_t* man, bool destructive, opt_pk_array_t* oa,
		      elina_dim_t* tdim, size_t size,
		      bool project){
	opt_pk_internal_t* opk = opt_pk_init_from_manager(man,ELINA_FUNID_FORGET_ARRAY);
        opt_pk_internal_realloc_lazy(opk,oa->maxcols - 2);
	array_comp_list_t * acla = oa->acl;
	unsigned short int maxcols = oa->maxcols;
        size_t i;
	size_t gg;
	#if defined(TIMING)
 	 		start_timing();
   	#endif 
	if(oa->is_bottom || !acla){
		man->result.flag_best = man->result.flag_exact = true;
		#if defined(TIMING)
 	 		record_timing(forget_array_time);
   		#endif 
		return opt_pk_bottom(man,maxcols -2 ,0);
	}
	
	
	/*******************************
		Minimize the input
	*********************************/
	unsigned short int num_compa = acla->size;
	unsigned short int j, k;
	
	opt_pk_t ** poly_a = oa->poly;
	for(k=0; k < num_compa; k++){
	    if(opk->funopt->algorithm>=0){
	       opt_poly_chernikova(man,poly_a[k],"cons to gen");
	    }
	    else{
		opt_poly_obtain_F(man,poly_a[k],"cons to gen");
	    }
	    if (opk->exn){
		opk->exn = ELINA_EXC_NONE;
		if (!poly_a[k]->F){
		     man->result.flag_best = man->result.flag_exact = false;
		     if(destructive){
			opt_poly_set_top(opk,oa);
			return oa;  
		     }
		     else{
			  return opt_pk_top(man,oa->maxcols-2,0);
		     }    
		}
  	    }
	   /* if empty, return empty */
	   if(!poly_a[k]->F){
	       man->result.flag_best = man->result.flag_exact = true;
	       if(destructive){
		  opt_poly_set_bottom(opk,oa);
		  return oa;
	       }
	       else{
		    return opt_pk_bottom(man,oa->maxcols-2,0);
	       }
	   }
	}			
	comp_list_t * cla = acla->head;
	/********************************
		Handle the destructive case
	*********************************/
	//char * exc_map = (char *)calloc(num_compa,sizeof(char));
	if(destructive){
		k = 0;
		while(k < num_compa){
			unsigned short int * ca = to_sorted_array(cla,maxcols);
			unsigned short int comp_size = cla->size;
			unsigned short int dim_size = 0;
			bool flag = true;
			elina_dim_t * tdimk = (elina_dim_t *)calloc(comp_size, sizeof(elina_dim_t));
			opt_pk_t * src = poly_a[k];
			for(i = 0; (i < size) && flag; i++){
				unsigned short int var = tdim[i] + opk->dec;
				unsigned short int l;
				//printf("remove var: %d\n",var);
				//fflush(stdout);
				for(l = 0; l < comp_size; l++){
					unsigned short int numa = ca[l];
					if(var==numa){
						tdimk[dim_size] = l;
						
						/**************************
							remove the component
						***************************/
						remove_comp(cla,var);
						if(cla->size==0){
							flag = false;
							break;	
						}
						dim_size++;
					}
				}
			}
			if(!flag){
				comp_list_t * tcla = cla->next;
				remove_comp_list(acla,cla);
				cla = tcla;
				free(ca);
				free(tdimk);
				opt_pk_t * tpoly = poly_a[k];
				unsigned short int k1;
				for(k1=k; k1 < num_compa - 1; k1++){
					poly_a[k1] = poly_a[k1+1];
				}
				opt_poly_clear(tpoly);
				free(tpoly);
				num_compa--;
				continue;
			}
			else if(dim_size){
				/**************************
					Project the blocks
				**************************/
				bool res = false;
				if(project){
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
				   opt_matrix_t * F = src->F;
				   for(i=0; i < F->nbrows; i++){
				       for(j=0; j < dim_size; j++){
					   F->p[i][opk->dec + tdimk[j]] = 0;
				       }
				       opt_matrix_normalize_row(opk,F,i);
				   }
				   src->status = 0;
				   if(opk->funopt->algorithm>=0){
				      opt_poly_chernikova(man,src,"of the result");
				   }
				}
				else {
   				     /* Forget */	
   				     opt_matrix_t *F = opt_matrix_alloc(size,src->F->nbcolumns,false);
				     for (i=0; i<dim_size; i++){
					  F->p[i][opk->dec+tdimk[i]] = 1;
				     }
				     opt_matrix_sort_rows(opk,F);
				     //opt_poly_dual(src);
				     //if (opk->funopt->algorithm>=0) {
				     //opt_poly_obtain_satC(src);
				     //}
				     opt_matrix_append_with(src->F,F);
				     //res = opt_poly_meet_matrix(false,opk->funopt->algorithm<0,man,src,src,F);
				     //opt_poly_dual(src);
				     opt_matrix_free(F);
				}
				/*if (res || opk->exn){
				      opk->exn = ELINA_EXC_NONE;
				      if (!src->F){
					   man->result.flag_best = man->result.flag_exact = false;
					   opt_poly_set_top(opk,oa);
					   printf("forget end9\n");
  	 				   fflush(stdout);
					   return oa;
				      }
				}*/
				if (opk->funopt->flag_best_wanted || opk->funopt->flag_exact_wanted){
				    bool real = true;
				    if (src->intdim>0){
					for (i=0; i<dim_size; i++){
					     if (tdimk[i]<src->intdim){
						  real = false;
						  break;
						}
					}
				    }
				    man->result.flag_best = man->result.flag_exact = 
					      real;
				}
				else {
				     man->result.flag_best = man->result.flag_exact =
					      src->intdim==0;
				}
				unsigned short int ncomp_size = comp_size - dim_size;
				unsigned short int * nca = to_sorted_array(cla,maxcols);
				poly_a[k] = opt_poly_alloc(ncomp_size,0);
				opt_matrix_t *src_mat, *dst_mat;
				src_mat = src->F;
				//dst_mat = poly_a[k]->F;
				poly_a[k]->nbeq = 0;
				poly_a[k]->nbline = src->nbline;
				size_t nbcons = src_mat->nbrows;
				dst_mat = opt_matrix_alloc(nbcons,ncomp_size+opk->dec,src_mat->_sorted);
				opt_numint_t ** tp = src_mat->p;
				opt_numint_t ** dp = dst_mat->p;
				for(i=0; i < nbcons; i++){
				    opt_numint_t * tpi = tp[i];
				    opt_numint_t * dpi = dp[i];
				    dpi[0] = tpi[0];
				    dpi[1] = tpi[1];
				    unsigned short int l = 0;
				    for(j=0; j < ncomp_size;j++){
					unsigned short int var = nca[j];
					while(ca[l]!=var){
					      l++;
					}
					dpi[j+2] = tpi[l+2];
					l++;
				    }
				}	
				poly_a[k]->is_minimized = false;
				opt_matrix_free(src_mat);
				free(src);
				free(nca);
				poly_a[k]->F = dst_mat;
				opt_poly_chernikova(man,poly_a[k],"forget destructive");
				if(opk->exn){
					opk->exn = ELINA_EXC_NONE;
					man->result.flag_best = man->result.flag_exact = false; 
					opt_poly_set_top(opk,oa);
					return oa; 
						//exc_map[k] = 1;
				}
			}
			cla = cla->next;
			k++;
		}
		#if defined(TIMING)
 	 		record_timing(forget_array_time);
   		#endif
		/*unsigned short int k1=0;
	  	unsigned short int bound = oa->acl->size;
	  	cl = oa->acl->head;
	  	for(k=0; k < num_comp;k++){
	  		opt_pk_t *oak = oa->poly[k1];
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
	      	}*/
		return oa;
	}
	/********************************
		Handle non destructive case
	**********************************/
	else{
		array_comp_list_t * acl = create_array_comp_list();
		/*****************************
			Initialize the blocks
		*****************************/
		opt_pk_t ** poly = (opt_pk_t **)malloc(num_compa*sizeof(opt_pk_t*));
		cla = acla->head;
		unsigned short int k1=0;
		for(k=0; k < num_compa; k++){
			unsigned short int * ca = to_sorted_array(cla,maxcols);
			unsigned short int comp_size = cla->size;
			unsigned short int dim_size = 0;
			elina_dim_t * tdimk = (elina_dim_t *)calloc(comp_size, sizeof(elina_dim_t));
			comp_list_t *cl = create_comp_list();
			unsigned short int l;
			opt_pk_t * src = poly_a[k];
			for(l = 0; l < comp_size; l++){
				bool flag = true;
				unsigned short int numa = ca[l];				
				for(i = 0; i < size; i++){
					unsigned short int var = tdim[i] + opk->dec;
					if(var==ca[l]){
						tdimk[dim_size] = l;
						dim_size++;
						flag = false;
						break;	
					}
				}
				if(flag){
					insert_comp(cl,numa);
				}
				
			}
			if(cl->size){
				insert_comp_list(acl,cl);
				if(dim_size){
					/**************************
						Project the blocks
					**************************/
					
					opt_pk_t * tmp = opt_poly_alloc(comp_size,0);
					bool res = false;
					if(project){
					   tmp->F = opt_matrix_copy(src->F);
					   opt_matrix_t * F = tmp->F;
					   for(i=0; i < F->nbrows; i++){
					       for(j=0; j < dim_size; j++){
						   F->p[i][opk->dec + tdimk[j]] = 0;
					       }
					       opt_matrix_normalize_row(opk,F,i);
					   }
					   tmp->status = 0;
					   if(opk->funopt->algorithm>=0){
					      opt_poly_chernikova(man,src,"of the result");
					   }
					  
					}
					else {
	   				     /* Forget */
	   				     opt_matrix_t *F = opt_matrix_alloc(size,src->F->nbcolumns,false);
					     for (i=0; i<dim_size; i++){
						  F->p[i][opk->dec+tdimk[i]] = 1;
					     }
					     opt_matrix_sort_rows(opk,F);
					     //opt_poly_dual(src);
					     //opt_poly_dual(tmp);
					     //if (opk->funopt->algorithm>=0){
					     //opt_poly_obtain_satC(src);
					     //}
					     tmp->F = opt_matrix_append(src->F,F);
					     //res = opt_poly_meet_matrix(false,opk->funopt->algorithm<0,man,src,src,F);
					     //opt_poly_dual(src);
					     //opt_poly_dual(tmp);	
					     opt_matrix_free(F);
					}
					/*if (res || opk->exn){
					      opk->exn = ELINA_EXC_NONE;
					      if (!src->F){
						   man->result.flag_best = man->result.flag_exact = false;
						   opt_poly_clear(tmp);
						   printf("forget end7\n");
  	 					   fflush(stdout);
						   return opt_pk_top(man,oa->maxcols-2,0);
					      }
					}*/
					if (opk->funopt->flag_best_wanted || opk->funopt->flag_exact_wanted){
					    bool real = true;
					    if (src->intdim>0){
						for (i=0; i<dim_size; i++){
						     if (tdimk[i]<src->intdim){
							  real = false;
							  break;
						     }
						}
					    }
					    man->result.flag_best = man->result.flag_exact = 
						      real;
					}
				        else {
					     man->result.flag_best = man->result.flag_exact =
						      src->intdim==0;
				        }
					unsigned short int ncomp_size = comp_size - dim_size;
					unsigned short int * nca = to_sorted_array(cl,maxcols);
					poly[k1] = opt_poly_alloc(ncomp_size,0);
					opt_matrix_t *src_mat, *dst_mat;
					src_mat = tmp->F;
					//dst_mat = poly[k1]->F;
					poly[k1]->nbline = tmp->nbline;
					poly[k1]->nbeq = 0;
					size_t nbcons = src_mat->nbrows;
					dst_mat = opt_matrix_alloc(nbcons,ncomp_size+opk->dec,src_mat->_sorted);
					
					opt_numint_t ** tp = src_mat->p;
					opt_numint_t ** dp = dst_mat->p;
					for(i=0; i < nbcons; i++){
						opt_numint_t * tpi = tp[i];
						opt_numint_t * dpi = dp[i];
						dpi[0] = tpi[0];
						dpi[1] = tpi[1];
						unsigned short int l = 0;
						for(j=0; j < ncomp_size;j++){
							unsigned short int var = nca[j];
							while(ca[l]!=var){
								l++;
							}
							dpi[j+2] = tpi[l+2];
							l++;
						}
					}
					poly[k1]->is_minimized = false;
					opt_matrix_free(src_mat);
					free(tmp);
					free(nca);
					poly[k1]->F = dst_mat;
					opt_poly_chernikova(man,poly[k1],"forget non destructive");
					if(opk->exn){
						opk->exn = ELINA_EXC_NONE;
						man->result.flag_best = man->result.flag_exact = false; 
						return opt_pk_top(man,oa->maxcols-2,0);
					}
				}
				else{
					poly[k1] = opt_poly_alloc(src->intdim,src->realdim);
					opt_poly_copy(poly[k1],src);
					poly[k1]->is_minimized = src->is_minimized;
				}
				k1++;
			}
			else{
				free_comp_list(cl);
			}
			free(ca);
			free(tdimk);
			cla = cla->next;
		}
		array_comp_list_t * res = copy_array_comp_list(acl);
		free_array_comp_list(acl);
		poly = (opt_pk_t **)realloc(poly,k1*sizeof(opt_pk_t*));
		opt_pk_array_t * op = opt_pk_array_alloc(poly,res,maxcols);
		 #if defined(TIMING)
 	 		record_timing(forget_array_time);
   		#endif
		return op;
	}
   
}


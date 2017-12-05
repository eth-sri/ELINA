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
/* opt_pk_cherni.c: Conversion from one representation to the dual one.  */
/* ********************************************************************** */

#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_bit.h"
#include "opt_pk_satmat.h"
#include "opt_pk_matrix.h"
#include "opt_pk_cherni.h"

/* ********************************************************************** */
/*  Conversion algorithm */
/* ********************************************************************** */
void opt_cherni_resize(opt_matrix_t* mat, opt_satmat_t* osc)
{
  assert(mat->nbrows==osc->nbrows);
  size_t nbrows = mat->nbrows;
  size_t currentsize = mat->_maxrows >= osc->_maxrows ? mat->_maxrows : osc->_maxrows;
  size_t addsize = currentsize < 20 ? 10 : currentsize / 2;
  opt_matrix_resize_rows(mat, currentsize+addsize);
  opt_satmat_resize_rows(osc, currentsize+addsize);
  mat->nbrows = osc->nbrows = nbrows;
  return;
}


size_t opt_cherni_conversion(opt_pk_internal_t* opk,
			 opt_matrix_t* con, size_t begin,
			 opt_matrix_t* ray, opt_satmat_t* osc, size_t nbline)
{
//  printf("%d %d\n",con->nbrows,con->nbcolumns);
 // fflush(stdout);
  size_t i,j,l,w;
  int is_inequality;
  size_t index_non_zero;
  size_t equal_bound,sup_bound,inf_bound,bound;
  int nbcommonconstraints;
  opt_bitindex_t k;
  opt_bitstring_t m,aux;
  bool redundant;
  opt_bitstring_t* bitstringp;
  const size_t nbcols = con->nbcolumns;
  const size_t satnbcols = opt_bitindex_size(con->nbrows);
  size_t nbrows = ray->nbrows;

  bitstringp = opt_bitstring_alloc(satnbcols);
 
  /* ================= Code ================== */
  k = opt_bitindex_init(begin);
  
  while (k.index < con->nbrows){
    /* Iteration for the constraints */
    is_inequality = opt_numint_sgn(con->p[k.index][0]);
	
    /* Scalar product and index: */
    /*
      We compute for the new considered constraint its scalar products
      with each frame and put them into the first coefficient of the
      frames. Moreover we set index_non_zero to the index of the
      first frame that does not saturate the constraint.
    */

    index_non_zero = nbrows;
//    unsigned short int * nz = build_index_vector(con->p[k.index],nbcols); 
    for (i=0; i<nbrows; i++){
      //ray->p[i][0] = opt_vector_product_with_index(opk, ray->p[i],
	//     con->p[k.index],nz);
       ray->p[i][0] = opt_vector_product(opk,ray->p[i],con->p[k.index],nbcols);
	if (opk->exn) {goto cherni_conversion_exit0;}
      if (index_non_zero == nbrows && opt_numint_sgn(ray->p[i][0])){
	index_non_zero = i;
      }
    }
  // free(nz);
    if (index_non_zero < nbline){
      /* If a line doesn't satisfy the constraint */
      /* Line becomes ray:
	 We transform the first line that does not saturate the constraint
	 into a ray and computes a new pointed cone. */

      /* remove it of lines and put it at index nbline */
      nbline--;
      if (index_non_zero != nbline)
	opt_matrix_exch_rows(ray,index_non_zero,nbline);
      /* compute new lineality space */
      opt_numint_t snbline = opt_numint_sgn(ray->p[nbline][0]); 
      for (i=index_non_zero; i<nbline; i++){
	if (opt_numint_sgn(ray->p[i][0])){
          bool flag;
          opt_numint_t si = opt_numint_sgn(ray->p[i][0]);
	  if(si==0){
		continue;
	  }
	  if(snbline * si < 0){
		flag = true;
	  }
	  else{
		flag = false;
	  }
	  opt_matrix_combine_rows(opk,ray,i,nbline,i,0,flag);
	  
	  if (opk->exn) {goto cherni_conversion_exit0;}
	}
       }
        
      /* orient the new ray */
      if (opt_numint_sgn(ray->p[nbline][0]) < 0){
	for (j=0; j<nbcols; j++){
	  ray->p[nbline][j] = -ray->p[nbline][j];
	}
      }
      /* compute the new pointed cone */
      snbline = opt_numint_sgn(ray->p[nbline][0]); 
      for (i=nbline+1; i<nbrows; i++){
	if (opt_numint_sgn(ray->p[i][0])){
	  bool flag;
	  opt_numint_t si = opt_numint_sgn(ray->p[i][0]);
	  if(si==0){
		continue;
	  }
	  if(snbline * si < 0){
		flag = true;
	  }
	  else{
		flag = false;
	  }
	  opt_matrix_combine_rows(opk,ray,i,nbline,i,0,flag);
	  if (opk->exn){ goto cherni_conversion_exit0;}
	}
      }
     
      /* For the saturation matrix, we only add a column, */
      /* so new bits are initialized to zero (see above) */
      if (is_inequality){
	/* rays saturates all apart the last inequality */
	opt_satmat_set(osc,nbline,k);
      } else {
	/* one remove the ray */
	nbrows --; ray->nbrows --; osc->nbrows--;
	opt_matrix_exch_rows(ray, nbline, nbrows);
	opt_satmat_exch_rows(osc, nbline, nbrows);
      }
      
      opt_bitindex_inc(&k);
	
    }

    else {
     
      /* if all lines saturates the constraint */
      /* Sort the rays: */
      /*
	Rays are sorted as follows:
	- nbline<= i < equal_bound:
	saturate the constraint;
	- equal_bound <= i < sup_bound:
	verify it;
	- sup_bound <= i < nbrows:
	do not verify it.
      */
      equal_bound=nbline;
      sup_bound=nbline;
      inf_bound=nbrows;
      while (inf_bound>sup_bound) {
	int s = opt_numint_sgn(ray->p[sup_bound][0]);
	if (s==0){
	  opt_matrix_exch_rows(ray, sup_bound, equal_bound);
	  opt_satmat_exch_rows(osc, sup_bound, equal_bound);
	  equal_bound++;
	  sup_bound++;
	} else if (s<0) {
	  inf_bound--;
	  opt_matrix_exch_rows(ray, sup_bound, inf_bound);
	  opt_satmat_exch_rows(osc, sup_bound, inf_bound);
	} else {
	  sup_bound++;
	}
      }
      if (is_inequality && sup_bound == nbrows){
	/* all rays satisfy the constraint:redundancy */
	con->nbrows--;
	opt_matrix_exch_rows(con, k.index, con->nbrows);
      }
      else {
	if (sup_bound==nbline){ /* no ray satisfies the constraint */
	  nbrows = ray->nbrows = osc->nbrows = nbline;
	}
	else { /* some rays do not satisfy the constraint */
	  /* Compute the new cones by combining adjacent constraints: */
	  bound = nbrows;
	  for (i=equal_bound; i<sup_bound; i++){
	    for(j=sup_bound; j<bound; j++){
	      /* For each pair R+,R-, */
	      /* compute the set of constraints saturated by both of them,
		 including equalities */
	      nbcommonconstraints=0;
	      for (w=0; w<k.word; w++) {
		aux = osc->p[i][w] | osc->p[j][w];
		bitstringp[w] = aux;
		for (m=opt_bitstring_msb; m!=0; m>>=1)
		  if ((aux & m)==0) nbcommonconstraints++;
	      }
	      aux = osc->p[i][k.word] | osc->p[j][k.word];
	      bitstringp[k.word] = aux;
	      for (m=opt_bitstring_msb; m!=k.bit; m>>=1){
		if ((aux & m)==0) nbcommonconstraints++;
	      }
	      if (nbcommonconstraints+nbline>=nbcols-3){ /* possibly adjacent */
		/* Does exist another ray saturating the same constraints ? */
		redundant=false;
		for (l=nbline; l<bound; l++){
		  if ((l!=i)&&(l!=j)){
		    for (w=0; w<=k.word; w++){
		      if (osc->p[l][w] & ~(bitstringp[w]))
			break;
		    }
		    if (w>k.word){
		      redundant=true; break;
		    }
		  }
		}
		if (!redundant){ /* if not */
		  if (opk->funopt->max_object_size && nbrows * (nbcols - opk->dec) >  opk->funopt->max_object_size){
		    /* out of space overflow */
		    opk->exn = ELINA_EXC_OUT_OF_SPACE;
		    goto cherni_conversion_exit0;
		  }
		  if (nbrows>=opt_matrix_get_maxrows(ray) || nbrows>=osc->_maxrows){
		    /* resize output matrices */
		    opt_cherni_resize(ray,osc);
		  }
		  /* Compute the new ray and put it at end */
		  bool flag;
		  if(opt_numint_sgn(ray->p[i][0])*opt_numint_sgn(ray->p[j][0]) < 0){
			flag = true;
		  }
		  else{
			flag = false;
		  }
		  //start_timing();
		  opt_matrix_combine_rows(opk,ray,j,i,nbrows,0,flag);
		  //record_timing(opt_conversion_time);
		  if (opk->exn){goto cherni_conversion_exit0;}
		  /* New row in saturation matrix */
		  for (w=0; w<=k.word; w++){
		    osc->p[nbrows][w] = bitstringp[w];
		  }
		  for (w=k.word+1; w<satnbcols; w++){
		    osc->p[nbrows][w] = 0;
		  }
		  nbrows ++; ray->nbrows ++; osc->nbrows ++;
		}
	      }
	    }
	  }
	  /* Remove non extremal rays by exchanging added rays */
	  /* with those that don't verify the constraint */
	  {
	    size_t l;
	    if (is_inequality){
	      j = sup_bound;
	      for (l=equal_bound; l<sup_bound; l++) opt_satmat_set(osc,l,k);
	    }
	    else {
	      j = equal_bound;
	    }
	    i = nbrows;
	    while ((j<bound)&&(i>bound)) {
	      i--;
	      opt_matrix_exch_rows(ray,i,j);
	      opt_satmat_exch_rows(osc,i,j);
	      j++;
	    }
	    nbrows = (j==bound) ? i : j;
	    ray->nbrows = osc->nbrows = nbrows;
	  }
	}
	opt_bitindex_inc(&k);
      }
	
    }
    
  }
  
  /* status coefficient */
  for (i=0; i<nbline; i++){
    ray->p[i][0] = 0;
  }
  for (i = nbline; i<nbrows; i++){
    ray->p[i][0] = 1;
  }
  ray->nbrows = osc->nbrows = nbrows;
  opt_bitstring_free(bitstringp);
  
  return nbline;

 cherni_conversion_exit0:
  opt_bitstring_free(bitstringp);
  
  return 0;
}


/* ====================================================================== */
/* III.3 Simplifying constraints (main function) */
/* ====================================================================== */

/*
We suppose that we just obtained ray and satc from con
with the chernikova algorithms. As a consequence the system of rays is
minimal. satf is the transposed matrix of satc, i.e. rows are
indexed by constraints. con is supposed to be non empty.

We have still to simplify con by detecting new equality constraints,
removing redundant inequalities, and obtaining a minimal system of
equalities. This is performed by gauss elimination.

Throw exception.
*/

int opt_cherni_simplify(opt_pk_internal_t* opk,
		    opt_matrix_t* con, opt_matrix_t* ray, opt_satmat_t* satf, size_t nbline)
{
  size_t i,j;
  long int nb,nbj;
  size_t nbeq,rank;
  size_t w;
  opt_bitstring_t m;
  bool redundant, is_equality;
  const size_t nbcols = con->nbcolumns;
  opt_bitindex_t nbrays = opt_bitindex_init(ray->nbrows);
  size_t nbcons = con->nbrows;
  con->_sorted = false;
  
  /* find the first inequality */
  for (nbeq=0; nbeq < nbcons; nbeq ++){
    if (opt_numint_sgn(con->p[nbeq][0])) break;
  }

  /* For each constraint,
     - put it with equalities if it detected as one
     - or set the status word to the number of rays that saturates the constraint,
  */
  for (i=nbeq; i < nbcons; i++){
    is_equality = (opt_numint_sgn(con->p[i][0])==0);
    if (!is_equality){
      is_equality = true;
	
      for (w=0; w<satf->nbcolumns; w++){
	if (satf->p[i][w]){
	  is_equality = false; break;
	}
      }
    }
    if (is_equality){
      /* we have an equality */
	
      con->p[i][0] = 0;
      opt_matrix_exch_rows(con, i,nbeq);
      opt_satmat_exch_rows(satf,i,nbeq);
      nbeq++;
    }
    else {
      /* we count the number of zero bits */
      nb = 0;
      for (w=0; w<nbrays.word; w++){
	for (m=opt_bitstring_msb; m!=0; m>>=1)
	  if ((satf->p[i][w] & m)==0) nb++;
      }
      for (m=opt_bitstring_msb; m!=nbrays.bit; m>>=1)
	if ((satf->p[i][nbrays.word] & m)==0) nb++;
        con->p[i][0] = (int)nb;
    }
  }
  /* remove redundant equalities and update nbeq */

  rank = opt_matrix_gauss_elimination(opk,con, nbeq); /* gauss pivot to simplify equalities */
  opk->exn = ELINA_EXC_NONE;

  /* remove redundants equations, located between rank and nbeq */
  if (rank<nbeq) {
    i = nbcons;
    j = rank;
    while( j < nbeq && i > nbeq ) {
      i--;
      opt_matrix_exch_rows(con, j,i);
      opt_satmat_exch_rows(satf,j,i);
      j++;
    }
    nbcons += rank - nbeq;
    nbeq = rank;
  }
  
  /* remove trivially redundants inequalities (nb < nbcols-nbeq-2) */
  i = nbeq;
  while (i < nbcons){
    //int_set_numint(&nb, con->p[i][0]);
    nb = con->p[i][0];
    if (nb < (long int)(nbcols-nbeq-2)){ /* redundant constraint */
      nbcons--;
      opt_matrix_exch_rows(con, i,nbcons);
      opt_satmat_exch_rows(satf,i,nbcons);
    }
    else
      i++;
  }
  /* remove others redundants inequalities */
  i=nbeq;
  while (i < nbcons){
    //int_set_numint(&nb,con->p[i][0]);
    nb = con->p[i][0];
    redundant = false;
    j = nbeq;
    while (j < nbcons){
      //int_set_numint(&nbj,con->p[j][0]);
     nbj = con->p[j][0];
      if (nbj > nb){
	/* does j saturates a strictly overset ? */
	redundant = true;
	for (w=0; w<satf->nbcolumns; w++){
	  if( ~(satf->p[i][w]) & satf->p[j][w] ){
	    redundant = false; break;
	  }
	}
	if (redundant)
	  break;
	else
	  j++;
      }
      else if (nbj == nb && j != i){
	  /* is j mutually redundant with i ? */
	  is_equality = true;
	  for (w=0; w<satf->nbcolumns; w++){
	    if( satf->p[i][w] != satf->p[j][w] ){
	      is_equality = false; break;
	    }
	  }
	  if (is_equality){
	    /* yes: we can remove j */
	    nbcons--;
	    opt_matrix_exch_rows(con, j,nbcons);
	    opt_satmat_exch_rows(satf,j,nbcons);
	  }
	  else
	    j++;
      }
      else
	j++;
    }
    if (redundant){
      nbcons--;
      opt_matrix_exch_rows(con, i,nbcons);
      opt_satmat_exch_rows(satf,i,nbcons);
    }
    else
      i++;
  }
  /* setting status coefficient */
  for (i=nbeq; i<nbcons; i++){
    con->p[i][0] = 1;
  }
  con->nbrows = satf->nbrows = nbcons;
  /* back substitution of remaining constraints */
  gauss_backsubstitute(opk, con, nbeq);
  opk->exn = ELINA_EXC_NONE;

  return nbeq;
}


/* ====================================================================== */
/*  Standard minimization */
/* ====================================================================== */

/* op is a polyhedron with a non-empty constraint system, and empty frame and
   saturation matrices. */

#define OPT_CHERNI_FACTOR 2

static inline size_t opt_uint_max(size_t a, size_t b)
{ return (a>b ? a : b); }

void opt_cherni_minimize(opt_pk_internal_t* opk,
		     bool con_to_ray,
		     opt_pk_t* op)
{
  
  size_t i;
  bool special;
  opt_matrix_t* C;
  opt_matrix_t* F;
  opt_satmat_t* satC;
  
  C = op->C;
  assert(opt_matrix_is_sorted(C) && 
	 op->F==NULL && op->satC==NULL && op->satF==NULL);
 // printf("cherni start %d %d\n",C->nbrows, C->nbcolumns);
  //opt_matrix_fprint(stdout,C);
  //fflush(stdout);
  /* initialization of F and sat */
  F = opt_matrix_alloc(OPT_CHERNI_FACTOR*opt_uint_max(C->nbrows, C->nbcolumns-1),
		   C->nbcolumns,false);
  satC = opt_satmat_alloc(F->nbrows, opt_bitindex_size(C->nbrows));
  
  F->nbrows = satC->nbrows = C->nbcolumns-1;
  /*bool flag = true;
  for(i=0; i < C->nbrows; i++){
	if(C->p[i][1]){
		flag = false;
		break;
	}
  }*/
  /*if(C->nbrows==1 && C->nbcolumns==3 && C->p[0][2]!=0 && con_to_ray){
	
	bool is_inequality = (C->p[0][0] == 1);
	F->p[0][0] = 1;
	F->p[0][1] = 1;
	F->p[0][2] = -C->p[0][1]*opt_numint_sgn(C->p[0][2]);
	if(is_inequality){
		F->p[1][0] = 1;
		F->p[1][1] = 0;
		F->p[1][2] = C->p[0][2];
		opt_satmat_set(satC,1,opt_bitindex_init(0));
	}
	else{
		F->nbrows--;
		satC->nbrows--;
	}
	op->nbline = 0;
  }*/
  //border case if Ax >= 0, then no need to embedd in higher dimensional space, APRON has a bug for this
  /*else if(flag && con_to_ray){
	for (i=0; i<C->nbcolumns-2; i++){
    		F->p[i][i+1] = 1;
  	}
	F->nbrows--;
	satC->nbrows--;
	opk->exn = ELINA_EXC_NONE;
  	op->nbline = opt_cherni_conversion(opk, C, 0, F, satC, C->nbcolumns-1);
  	if (opk->exn){
    	
    		opt_matrix_free(F);
    		opt_satmat_free(satC);
    		op->nbeq = op->nbline = 0;
  	}
  }*/
  //else{
  	for (i=0; i<C->nbcolumns-1; i++){
    		F->p[i][i+1] = 1;
  	}
  	/* conversion */
  	opk->exn = ELINA_EXC_NONE;
  	op->nbline = opt_cherni_conversion(opk, C, 0, F, satC, C->nbcolumns-1);
	
	/*if(C->nbrows>1){
		opt_matrix_exch_rows(C,0,C->nbrows-1);
 	}
  	C->nbrows--;*/
  	if (opk->exn){
    	/* out of space, overflow */
    		opt_matrix_free(F);
    		opt_satmat_free(satC);
    		op->nbeq = op->nbline = 0;
  	}
  	else {  
    	/* If con_to_ray, special case ? */
    	/* We have to test
       		- In non-strict mode, that $\xi$ can be strictly positive.
       		- In strict mode, that both $\xi$ and $\epsilon$ can be strictly
       		positive. Because $\xi\geq\epsilon$, we just need to check that
       		$\epslion$ can be strictly positive. */
		
    		if (con_to_ray){
      			special = true;
      			for (i = op->nbline; i< F->nbrows; i++){
				if (F->p[i][opk->dec-1]>0){
	  				special = false;
	  				break;
				}
      			}
      			if (special){
				/* this means we have an empty polyhedron */
				//printf("special\n");
				//fflush(stdout);
				opt_matrix_free(C);
				opt_matrix_free(F);
				opt_satmat_free(satC);
				op->C = 0;
				op->nbeq = op->nbline = 0;
				return;
      			}
    	     }
    	//}
    //}
    	op->F = F;
    
    	op->satF = opt_satmat_transpose(satC,C->nbrows);
    	opt_satmat_free(satC);
    	op->nbeq = opt_cherni_simplify(opk,C,F,op->satF,op->nbline);
	
    	if (F->_maxrows > 3*F->nbrows/2){
      		opt_matrix_resize_rows(F,F->nbrows);
      		opt_satmat_resize_cols(op->satF,opt_bitindex_size(F->nbrows));
    	}
   }
	
}


void opt_cherni_add_and_minimize(opt_pk_internal_t* opk, 
			     bool con_to_ray,
			     opt_pk_t* op,
			     size_t start)
{
  size_t i;
  bool special;
  opt_matrix_t* C;
  opt_matrix_t* F;
  opt_satmat_t* satC;
  size_t nbrows,nbcols;
 
  assert(opt_bitindex_size(op->C->nbrows)==op->satC->nbcolumns);

  C = op->C;
  F = op->F;
  satC = op->satC;
  nbrows = C->nbrows;
  nbcols = C->nbcolumns;

  assert(C!=NULL && F!=NULL && satC!=NULL);
  //printf("add and minimize start %d %d %d %d %d \n",start,op->nbeq,op->nbline,nbrows,nbcols);
  //opt_matrix_fprint(stdout,C);
  //opt_matrix_fprint(stdout,F);
  
  //fflush(stdout);
  if (op->satF!=NULL){ 
    opt_satmat_free(op->satF); 
    op->satF=NULL; 
  }

  /* saturation matrix */
  opt_satmat_resize_cols(satC, opt_bitindex_size(C->nbrows));
	
  /* conversion */
  F->_sorted = false;
  opk->exn = ELINA_EXC_NONE;
  
  op->nbline = opt_cherni_conversion(opk,
				 C,start,F,satC,op->nbline);
  
  //printf("output\n");
  //opt_matrix_fprint(stdout,C);
  //opt_matrix_fprint(stdout,F);
  //opt_satmat_fprint(stdout,satC);
  //fflush(stdout);
   
  if (opk->exn){
    /* out of space, overflow */
    opt_matrix_free(F);
    opt_satmat_free(satC);
    op->F = 0;
    op->satC = op->satF = 0;
    op->nbeq = op->nbline = 0;
  }
  else {  
    /* If con_to_ray, special case ? */
    /* We have to test
       - In non-strict mode, that $\xi$ can be strictly positive.
       - In strict mode, that both $\xi$ and $\epsilon$ can be strictly
       positive. Because $\xi\geq\epsilon$, we just need to check that
       $\epsilon$ can be strictly positive. */
    if (con_to_ray){ 
      special = true;
      for (i = op->nbline; i<F->nbrows; i++){
	if (opt_numint_sgn(F->p[i][opk->dec-1])>0){
	  special = false;
	  break;
	}
      }
      if (special){ /* this means we have an empty polyhedron */
	opt_matrix_free(C);
	opt_matrix_free(F);
	opt_satmat_free(satC);
	op->C = op->F = 0;
	op->satC = 0;
	op->nbeq = op->nbline = 0;
	return;
      }
    }
    op->satF = opt_satmat_transpose(satC,C->nbrows);
    opt_satmat_free(satC);
    op->satC = NULL;
    //printf("simplify input\n");
    //opt_matrix_fprint(stdout,C);
    //fflush(stdout);
    op->nbeq = opt_cherni_simplify(opk,C,F,op->satF,op->nbline);
     //printf("simplify output\n");
    //opt_matrix_fprint(stdout,C);
    //fflush(stdout);
    if (F->nbrows && (F->_maxrows > 3*F->nbrows/2)){
      opt_matrix_resize_rows(F,F->nbrows);
      opt_satmat_resize_cols(op->satF,opt_bitindex_size(F->nbrows));
    } 
  }
  //printf("add and minimize end\n");
  //opt_matrix_fprint(stdout,F);
  //opt_matrix_fprint(stdout,C);
  //opt_satmat_fprint(stdout,op->satF);
  //fflush(stdout);
}




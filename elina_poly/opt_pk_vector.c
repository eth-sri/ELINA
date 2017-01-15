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
/* opt_pk_vector.c: operations on vectors */
/* ********************************************************************** */


#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_internal.h"

/* ********************************************************************** */
/* I. Basic operations: creation, destruction, copying and printing */
/* ********************************************************************** */

/* Internal allocation function: the elements are not initialized. */
opt_numint_t* _opt_vector_alloc_int(unsigned short int size){
  assert(size>0);
  return (opt_numint_t*)malloc(size*sizeof(opt_numint_t));
}

/* Standard allocation function, with initialization of the elements. */
opt_numint_t* opt_vector_alloc(unsigned short int size)
{
  unsigned short int i;
  opt_numint_t* q;
  q = (opt_numint_t*)(malloc(size*sizeof(opt_numint_t)));
  for(i=0; i<size; i++){
    q[i] = 0;
  }
  return q;
}

/* Reallocation function, to change the dimension */
void opt_vector_realloc(opt_numint_t** pq, unsigned short int size, unsigned short int nsize)
{
  unsigned short int i;
  opt_numint_t* q;
  opt_numint_t* nq;

  q = *pq;
  for (i=nsize; i<size; i++){
    q[i] = 0;
  }
  nq = realloc(q, nsize*sizeof(opt_numint_t));
  for (i=size; i<nsize; i++){
    nq[i] = 0;
  }
  *pq = nq;
}

/* Copy/Assign function */
void opt_vector_copy(opt_numint_t * dst, opt_numint_t * src, unsigned short int size)
{
  unsigned short int i;
  #if defined(VECTOR)
	for(i=0; i < size/v_length; i++){
		v_opt_numint_t_type val = v_load_opt_numint_t(src + i*v_length);
		v_store_opt_numint_t(dst + i*v_length, val);
	}
  #else
  for (i=0; i<size/v_length; i++){
    dst[i] = src[i];
  }
  #endif
  for (i=size/v_length; i < size; i++){
	dst[i] = src[i];
  }
}

opt_numint_t * opt_vector_scalar_product(opt_pk_internal_t * opk, opt_numint_t * ov, opt_numint_t s, unsigned short int size){
	unsigned short int i;
	opt_numint_t * res = opt_vector_alloc(size);
	res[0] = ov[0];
	res[1] = ov[1]*s;
	#if defined(VECTOR)
		v_opt_numint_t_type vs = v_set1_opt_numint_t(s);
		for(i = opk->dec; i < size/v_length; i++){
			v_opt_numint_t_type val = v_load_opt_numint_t(ov + i*v_length);
			val = v_mul_opt_numint_t(val, vs);
			v_store_opt_numint_t(res + i*v_length, val);
		}
	#else
		for(i = opk->dec; i < size/v_length; i++){
			res[i] = s*ov[i];
		}
	#endif
	for(i = size/v_length; i < size; i++){
		res[i] = s*ov[i];
	}
	return res;
}

/********************************
	Assume both ov1 and ov2 are inequalities
*********************************/
void opt_vector_sum(opt_pk_internal_t * opk, opt_numint_t * ov1, opt_numint_t *ov2, unsigned short int size){
	unsigned short int i;
	ov1[1] = ov1[1] + ov2[1];
	#if defined(VECTOR)
		for(i = opk->dec; i < size/v_length; i++){
			v_opt_numint_t_type op1 = v_load_opt_numint_t(ov1 + i*v_length);
			v_opt_numint_t_type op2 = v_load_opt_numint_t(ov2 + i*v_length);
			v_opt_numint_t_type sum = v_add_opt_numint_t(op1,op2);
			v_store_opt_numint_t(ov1 + i*v_length, sum);
		}
	#else
		for(i = opk->dec; i < size/v_length; i++){
			ov1[i] = ov1[i] + ov2[i];
		}
	#endif
	for(i = size/v_length; i <size; i++){
		ov1[i] = ov1[i] + ov2[i];
	}
}
/* Deallocation function. */
void opt_vector_free(opt_numint_t* q, unsigned short int size)
{
  free(q);
}

/* Set all elements to zero. */
void opt_vector_clear(opt_numint_t * ov, unsigned short int size)
{
  unsigned short int i;
  #if defined(VECTOR)
	v_opt_numint_t_type zero = v_set1_opt_numint_t(0.0);
  	for(i = 0; i < size/v_length; i++){
		v_store_opt_numint_t(ov + i*v_length, zero);
	}
  #else
	for (i = 0; i<size/v_length; i++){
	 	ov[i] = 0;
  	}
  #endif
  for(i = size/v_length; i < size; i++){
	ov[i] = 0;
  }
}

/* Raw printing function. */
void opt_vector_print(opt_numint_t* ov, unsigned short int size)
{
  unsigned short int i;
  printf("vector %lld: ", (long)size);
  for (i=0; i<size; i++){
    printf("%lld ", ov[i]);
  }
  printf("\n");
}

/* ********************************************************************** */
/* II. Normalization */
/* ********************************************************************** */

/* ====================================================================== */
/* II.1 Pgcd computation */
/* ====================================================================== */

/* ---------------------------------------------------------------------- */
/* Searching the minimum non-zero coefficient */
/* ---------------------------------------------------------------------- */

/* The following functions search the index and the absolute value of the
   minimal non-zero coefficient of the vector opk->vector_numintp, of size size,
   supposed to contain positive values only.
   
   It returns its results with
   pointers index and min. If all coefficients are zero, then
   index is set to size and *min to 0.

   This function uses opk->vector_numintp and opk->vector_tmp[0]. */

static opt_numint_t
opt_vector_min_notzero(opt_pk_internal_t* opk,
		   unsigned short int size,
		   int* index)
{
  unsigned short int i;

  opt_numint_t * ov = opk->vector_numintp; 

  opt_numint_t min = 0;

  /* search the first non-zero coefficient
     and stores the index and the coeff in *index and *min */
  i = 0;
  while (i<size){
    if (ov[i] > 0){
      *index = i;
      min = ov[i];
      break;
    }
    i++;
  }
  i++;
  /* search now the minimum */
  while (i<size) {
    if (ov[i] > 0){
      if (min > ov[i]){
	*index = i;
	min = ov[i];
      }
    }
    i++;
  }
  return min;
}
/* This function computes the pgcd of a vector.

   This function uses opk->vector_numintp and opk->vector_tmp[0]. */

opt_numint_t opt_vector_gcd(opt_pk_internal_t* opk,
		opt_numint_t* ov, unsigned short int size)
{
  unsigned short int i;
  bool not_all_zero;
  opt_numint_t gcd;
  opt_numint_t * v = opk->vector_numintp; 

  for (i=0;i<size;i++){
    if(ov[i] < 0){
	v[i] = -ov[i];	
    }
    else{
	v[i] = ov[i];
    }
  }
  do {
    int index=0;
    gcd = opt_vector_min_notzero(opk,size,&index);
    if (gcd==0){
	 return gcd;
    }
    not_all_zero = false;
    for (i=0; i<size; i++)
      if ((int)i!=index){
	v[i] = v[i]%gcd;
	not_all_zero = not_all_zero || (v[i] != 0);
      }
  } while (not_all_zero);
  return gcd;
}

void opt_vector_simplify(opt_pk_internal_t *opk, opt_numint_t * ov, unsigned short int size){
	opt_numint_t gcd = opt_vector_gcd(opk, ov + 1, size-1);
	unsigned short int i;
	if(gcd){
		for(i = opk->dec - 1; i <size; i++){
			ov[i] = ov[i]/gcd;
		}
	}
	
}
/* ====================================================================== */
/* II.3 Main functions */
/* ====================================================================== */

/* The function vector_normalize normalizes the vector considered as
   a contraint or a generator. It does not modify q[0].

   This function use opk->vector_tmp[0..2] and opk->numintp. */

bool opt_vector_normalize(opt_pk_internal_t* opk,
		      opt_numint_t* ov, unsigned short int size)
{
  unsigned short int i;

  assert(size<=opk->maxcols);

  /*  computation of the pgcd */
  opk->vector_tmp[1] = opt_vector_gcd(opk,&ov[1],size-1);
  /* possible division */
  if (opk->vector_tmp[1]> 1){
    for (i=1; i<size; i++){
        ov[i] = ov[i]/opk->vector_tmp[1];
    }
    return true;
  }
  else
    return false;
}

/* The function vector_normalize normalizes the vector considered as
   an expression. It modifies q[0].

   This function use opk->vector_tmp[0..2] and opk->numintp. */

bool opt_vector_normalize_expr(opt_pk_internal_t* opk,
			   opt_numint_t* ov, unsigned short int size)
{
  unsigned short int i;

  assert(size<=opk->maxcols);

  /*  computation of the pgcd */
  opk->vector_tmp[1] = opt_vector_gcd(opk,&ov[0],size);
  /* possible division */
  if (opk->vector_tmp[1] > 1){
    for (i=0; i<size; i++){
      ov[i] = ov[i]/opk->vector_tmp[1];
    }
    return true;
  }
  else
    return false;
}

/* The function vector_normalize_constraint normalizes the vector considered as
   a constraint.

   - if strict mode, the epsilon coefficient is put to 0 or 1
 
   This function use opk->vector_tmp[0..1] and opk->numintp. */

bool opt_vector_normalize_constraint(opt_pk_internal_t* opk,
				 opt_numint_t* ov,
				 unsigned short int intdim, unsigned short int realdim)
{
  unsigned short int i;
  bool change = false;
  unsigned short int size = opk->dec+intdim+realdim;
  
  assert(opk->dec+intdim+realdim <= opk->maxcols);

    return opt_vector_normalize(opk,ov,size);
  //}
  return change;
}

/* The function vector_normalize_constraint_int normalizes the vector 
   considered as a constraint.

   - if it involves only integer dimensions, the constraint is tightened and
     renormalized.

   - it implies standard constraint normalization
    
   This function use opk->vector_tmp[0..1] and opk->numintp. */



bool opt_vector_normalize_constraint_int(opt_pk_internal_t* opk,
				       opt_numint_t* ov,
				       unsigned short int intdim, unsigned short int realdim)
{
  unsigned short int i;
  bool change = false;
  unsigned short int size = opk->dec+intdim+realdim;
  
  assert(opk->dec+intdim+realdim <= opk->maxcols);

  if (intdim>0 && 
      opt_vector_is_integer(opk,ov,intdim,realdim) &&
      !opt_vector_is_positivity_constraint(opk,ov,size)){
    if (opk->strict && (ov[opt_polka_eps] < 0)){
      change = true;
      ov[opt_polka_eps] = 0;
      ov[opt_polka_cst] = ov[opt_polka_cst] - 1;
    }
    /*  computation of the pgcd without constant (and epsilon, of course) */
   opk->vector_tmp[1] = opt_vector_gcd(opk, &ov[opk->dec], size-opk->dec);
    /* possible division */
    if (opk->vector_tmp[1] > 1){
      change = true;
      for (i=opk->dec; i<size; i++){
	ov[i] = ov[i]/opk->vector_tmp[1];
      }
      /* round the constant coefficient */
      if (ov[0]==0){
	ov[0] = ov[opt_polka_cst]%opk->vector_tmp[1];
	if (ov[0]!=0){
	  opt_vector_clear(ov,size);
	  ov[opt_polka_cst] = 1;
	} else {
	  ov[opt_polka_cst] = ov[opt_polka_cst]/opk->vector_tmp[1];
	}
      }
      else {
	ov[opt_polka_cst] = opt_numint_fdiv(ov[opt_polka_cst],opk->vector_tmp[1]);
      }
    }
  }
  else {
    change = opt_vector_normalize_constraint(opk,ov,intdim,realdim);
  }
  return change;
}

/* ********************************************************************** */
/* III. Comparison function */
/* ********************************************************************** */

/* Comparison function for vectors

The used order is the lexicographic order, with the exception that the
constant (and possibly epsilon) coefficient is considered last. As a
consequence, the equations or lines are classified before the
inequalities or rays when vectors are rows of a sorted matrix.

The meaning of the returned result res is:
- <0 : q1 is smaller than q2
- =0 : they are equal
- >0: q1 is greater than q2

This function uses opk->vector_tmp[0..3] and opk->vector_numintp.
*/

bool opt_vector_equal(opt_pk_internal_t* opk,
		   opt_numint_t* ov1, opt_numint_t* ov2,
		   unsigned short int size, unsigned short int* ind){
	if(ov1[*ind]!=ov2[*ind]){
		return false;
	}
	unsigned short int j;
	unsigned short int val = *ind;
	for(j=0; j < val; j++){
		if(ov1[j]!=ov2[j]){
			*ind = j;
			return false;
		}
	}
	for(j=val+1; j < size; j++){
		if(ov1[j]!=ov2[j]){
			*ind = j;
			return false;
		}
	}
	return true;
}

int opt_vector_compare_coeff(opt_pk_internal_t* opk,
		   opt_numint_t* ov1, opt_numint_t* ov2,
		   unsigned short int size){
	int res = 0;
	unsigned short int i;
	 
	for(i=opk->dec; i<size; i++){
		
    		res = (ov1[i] < ov2[i]? -1 : ov1[i] > ov2[i] ? 1 : 0);
    		if (res) {
			
			return res;
		}
  	}
	return res;
}

int opt_vector_compare(opt_pk_internal_t* opk,
		   opt_numint_t* ov1, opt_numint_t* ov2,
		   unsigned short int size)
{
  unsigned short int i;
  int res=1;
 
  assert(size<=opk->maxcols);

  /* bidirectional/unidirectional ? */
  res = (ov1[0] == ov2[0] ? 0 : (ov1[0]> ov2[0] ? 1 : -1));
  if (res) return res;
  /* comparison */
  
  res = opt_vector_compare_coeff(opk, ov1, ov2, size);
 
  if(res) return res;
  if (opt_polka_cst<size){
    res = (ov1[opt_polka_cst] < ov2[opt_polka_cst]? -1 : ov1[opt_polka_cst] > ov2[opt_polka_cst] ? 1 : 0);
    if (res) return res;
    //if (opk->strict && opt_polka_eps < size){
      //res = numint_cmp(q1[polka_eps],q2[polka_eps]);
    //}
  }
  return res;
}

/* ********************************************************************** */
/* IV. Combine function */
/* ********************************************************************** */

/* vector_combine computes a combination ov3 of ov1 and
   ov2 such that ov3[k]=0.  The first coefficient is never
   considered for computations, except when k==0.
   Assumes ov1[k] and ov2[k] have opposite signs.
   This function uses opk->vector_tmp[0..4] and opk->vector_numintp. */

void opt_vector_combine(opt_pk_internal_t* opk,
		    opt_numint_t* ov1, opt_numint_t * ov2,
		    opt_numint_t* ov3, size_t k, unsigned short int size, bool add)
{
  unsigned short int j;
  //ov3[0] = (ov1[0] | ov2[0]);
  //printf("input\n");
  //opt_vector_print(ov1,size);
  //opt_vector_print(ov2,size);
  //fflush(stdout);
  opk->vector_tmp[0] = opt_numint_gcd(ov1[k],ov2[k]);
 
  opk->vector_tmp[1] = opt_numint_abs(ov1[k]/opk->vector_tmp[0]);
  opk->vector_tmp[2] = opt_numint_abs(ov2[k]/opk->vector_tmp[0]);
  for (j=1;j<size;j++){
    if (j!=k){
      opk->vector_tmp[3] = opk->vector_tmp[2] * ov1[j];
      opk->vector_tmp[4] = opk->vector_tmp[1] * ov2[j];
      if(add)
      ov3[j] = opk->vector_tmp[3] + opk->vector_tmp[4];
      else 
      ov3[j] = opk->vector_tmp[3] - opk->vector_tmp[4];
    }
  }
  ov3[k] = 0;
  for (j=0; j<size; j++){
      if (opt_numint_abs(ov3[j]) > (ELINA_INT_MAX/2)){
	opt_vector_print(ov1,size);
	opt_vector_print(ov2,size);
	printf("exception %d %lld %lld\n",k,opt_numint_abs(ov3[j]),ELINA_INT_MAX/2);
	fflush(stdout);
	opk->exn = ELINA_EXC_OVERFLOW;
	return ;
      }
    }
  opt_vector_normalize(opk,ov3,size);
  //if (opk->max_coeff_size){
    
  //}
  //printf("output\n");
  //opt_vector_print(ov3,size);
  //fflush(stdout);
}


//************************************************************************ */
/* Resize Operations  */
/* *********************************************************************** */

/* Vector Permutation. */

void opt_vector_permute_dimensions(opt_pk_internal_t* opk,
			       opt_numint_t* nov, opt_numint_t* ov, unsigned short int size,
			       elina_dim_t* permut)
{
  bool destructive;
  unsigned short int j,newj;
  opt_numint_t* tmp;
  
  destructive = (nov==ov);
  
  /* Where to write ?
     If destructive, we write in a temporary vector
     otherwise we write in the destination
  */
  tmp = destructive ? opk->vector_numintp : nov;

  /* Fill the non-permuted fields */
  for (j=0; j<opk->dec && j<size; j++){
    	tmp[j] = ov[j];
  }
  /* Permutation itself */
  for (j=0; j<size-opk->dec; j++){
    newj = permut[j] + opk->dec;
    tmp[newj] = ov[j+opk->dec];
  }
  if (destructive){
    for(j=0; j<size; j++){
      nov[j]= tmp[j];
    }
  }
  return;
}

/* Vector Addition. */

void opt_vector_add_dimensions(opt_pk_internal_t* opk,
			   opt_numint_t* nov, 
			   opt_numint_t* ov, unsigned short int size,
			   elina_dimchange_t* dimchange)
{
  int i,k,dimsup;

  if (nov!=ov){ 
    for (i=0;i<(int)opk->dec && i<(int)size; i++) 
        nov[i] = ov[i];
  }
  dimsup = dimchange->intdim+dimchange->realdim;
  k = dimsup;
  for (i=(int)(size)-(int)opk->dec; i>=0; i--){
    if (i<(int)size-(int)opk->dec){
      nov[opk->dec+i+k] = ov[opk->dec+i];
    }
    while (k>=1 && dimchange->dim[k-1]==(elina_dim_t)i){
      k--;
      nov[opk->dec+i+k] = 0;
    }
  }
}



/* Vector Removal. */

void opt_vector_remove_dimensions(opt_pk_internal_t* opk,
			      opt_numint_t* nov, 
			      opt_numint_t* ov, unsigned short int size,
			      elina_dimchange_t* dimchange)
{
  size_t i,k,dimsup;
  
  if (nov!=ov){ 
    for (i=0;i<opk->dec && i<size; i++) 
      nov[i] = ov[i];
  }
  dimsup = dimchange->intdim+dimchange->realdim;
  k=0;
  for (i=0; i<size-dimsup-opk->dec; i++){
    while (k<dimsup && dimchange->dim[k]==i+k){
      k++;
    }
    nov[opk->dec+i] = ov[opk->dec+i+k];
  }
}

//* ********************************************************************** */
/* V. Algebraic operations */
/* ********************************************************************** */

/* Scalar product.

Compute the scalar product of q1 and q2 considered as vectors
of length size. The first coefficients are never considered.

This function uses opk->vector_tmp[0]. */

opt_numint_t * opt_vector_neg(opt_pk_internal_t *opk, opt_numint_t *ov, unsigned short int size){
	opt_numint_t * res = opt_vector_alloc(size);
	unsigned short int i;
	res[0] = ov[0];
	for(i = opk->dec - 1; i < size; i++){
		res[i] = -ov[i];
	}
	return res;
}


void opt_vector_mul_scalar(opt_numint_t *dst, opt_numint_t *src, opt_numint_t prod, unsigned short int size){
	unsigned short int i;
	for(i = 1; i < size; i++){
		dst[i] = src[i]*prod;
	}
}

void opt_vector_add(opt_numint_t *dst, opt_numint_t *op1, opt_numint_t *op2, unsigned short int size){
	unsigned short int i;
	for(i = 1; i < size; i++){
		dst[i] = op1[i] + op2[i]; 
	}
}

/*opt_numint_t opt_vector_product(opt_pk_internal_t* opk,
		    opt_numint_t* q1, opt_numint_t* q2, unsigned short int size)
{
  opt_numint_t prod;
  size_t j;
  prod = 0;
  for (j=1; j<size; j++){
    opk->vector_tmp[0] = q1[j]*q2[j];
    prod = prod + opk->vector_tmp[0];
  }
  return prod;
}*/

/* Same as previous function, but in case where opk->strict is
   true, the $\epsilon$ coefficients are not taken into account. */

opt_numint_t opt_vector_product_strict(opt_pk_internal_t* opk,
			   opt_numint_t* q1, opt_numint_t* q2, unsigned short int size)
{
  opt_numint_t prod;
  size_t j;
  if (opt_polka_cst<size){
    prod = q1[opt_polka_cst]*q2[opt_polka_cst];
  }
  else {
    prod = 0;
    return prod;
  }
  for (j=opk->dec; j<size; j++){
    opk->vector_tmp[0] = q1[j]*q2[j];
    prod = prod + opk->vector_tmp[0];
  }
  return prod;
}


opt_numint_t opt_vector_product_strict_comp(opt_pk_internal_t* opk,
			   opt_numint_t* q1, opt_numint_t* q2, unsigned short int * ind_map1,
			   unsigned short int *ind_map2, unsigned short int size1, 
			   unsigned short int size2, unsigned short int size)
{
  opt_numint_t prod;
  if (opt_polka_cst<size){
    prod = q1[opt_polka_cst]*q2[opt_polka_cst];
  }
  else {
    prod = 0;
    return prod;
  }
  opt_numint_t * tmp = opt_vector_alloc(size);
  unsigned short int j;
  //printf("check1\n");
  //fflush(stdout);
  for(j=opk->dec; j<size1; j++){
	unsigned short int ind = ind_map1[j-2];
	tmp[ind+2] = q1[j];
  }
  //printf("check2\n");
  //fflush(stdout);
  for(j=opk->dec; j<size2; j++){
	unsigned short int ind = ind_map2[j-2];
	tmp[ind+2] = tmp[ind+2]*q2[j];
  }
  //printf("check3\n");
  //fflush(stdout);
  for (j=opk->dec; j<size; j++){
   
    prod = prod + tmp[j];
  }
  //printf("check4\n");
  //fflush(stdout);
  free(tmp);
  return prod;
}
/* ********************************************************************** */
/* VI. Predicates */
/* ********************************************************************** */

/* The function tests if the given vector is null. */

bool opt_vector_is_null_expr(opt_pk_internal_t *opk,
			     opt_numint_t* ov, unsigned short int size){
	unsigned short int i;
	bool res = true;
	for(i = opk->dec; i < size; i++){
		if(ov[i] != 0){
			res = false;
			break;
		}
	}
	return res;
}

bool opt_vector_is_null(opt_pk_internal_t* opk,
		    opt_numint_t* ov, unsigned short int size)
{
  unsigned short int i;
  bool res = true;
 
  for (i=1; i<size; i++){
    if (ov[i]!=0){
      res = false;
      break;
    }
  }
  return res;
}



/* The function tests if the given vector represents a positivity
   constraint. */

bool opt_vector_is_positivity_constraint(opt_pk_internal_t* opk,
				     opt_numint_t* ov, unsigned short int size)
{
  if (size < opk->dec){
    return false;
  }
  else {
    unsigned short int i;
    bool res;

    res = (ov[0]>0);
    if (res){
      int s = ov[opt_polka_cst];
      if (s>0){
	/* Tests if it could be the positivity constraint */
	res = opk->strict ? ov[opt_polka_eps]<0 : true;
      }
      if (res){
	for (i=opk->dec; i<size; i++){
	  if (ov[i] != 0){
	    res = false;
	    break;
	  }
	}
      }
    }
    return res;
  }
}

/* The function tests if the given vector represents a positivity
   or a strictness constraint. */

bool opt_vector_is_dummy_constraint(opt_pk_internal_t* opk,
				opt_numint_t* ov, unsigned short int size)
{
  if (size < opk->dec){
    return false;
  }
  else {
    unsigned short int i;
    bool res;

    res = (ov[0]>0);
    if (res){
      int s = (ov[opt_polka_cst] > 0) ? 1 : (ov[opt_polka_cst] < 0) ? -1 : 0;
      if (s>0){
	/* Tests if it could be the positivity constraint */
	res = opk->strict ? (ov[opt_polka_eps]<0) : true;
      }
      else if (s==0){
	/* Tests if it could be the strictness constraint */
	res = opk->strict && (ov[opt_polka_eps]>0);
      }
      if (res){
	for (i=opk->dec; i<size; i++){
	  if (ov[i] != 0){
	    res = false;
	    break;
	  }
	}
      }
    }
    return res;
  }
}

/* Return true if all dimensions involved in the expression are integer
   dimensions */

bool opt_vector_is_integer(opt_pk_internal_t* opk, 
		       opt_numint_t* ov,
		       size_t intdim, size_t realdim)
{
  unsigned short int i;
  
  for (i=intdim; i<intdim+realdim; i++){
    if (ov[opk->dec+i] != 0){
      return false;
    }
  }
  return true;
}

// Assumes (size - 2) < comp_size
opt_numint_t * opt_map_vector(opt_numint_t *q, unsigned short int *map, unsigned short int comp_size, unsigned short int size){
	opt_numint_t * res = (opt_numint_t *)calloc(comp_size+2,sizeof(opt_numint_t));
	unsigned short int i;
	res[0] = q[0];
	res[1] = q[1];
	for(i=2; i < size; i++){
		unsigned short int ind = map[i-2];
		res[ind+2] = q[i];
	}
	return res;
}

opt_numint_t opt_vector_product(opt_pk_internal_t* opk,
		    opt_numint_t* q1, opt_numint_t* q2, size_t size)
{
  size_t j;
  opt_numint_t prod = 0;
  for (j=1; j<size; j++){
    prod = prod + q1[j]*q2[j];
  }
  return prod;
}

unsigned short int * build_index_vector(opt_numint_t *ov, unsigned short int size){
	unsigned short int i,s = 0;
	unsigned short int * res = (unsigned short int*)calloc(size+1,sizeof(unsigned short int));
	for(i=1; i < size; i++){
		if(ov[i]!=0){
			res[s+1] = i;
			s++;
		}
	}
	res[0] = s;
	return res;
}

opt_numint_t opt_vector_product_with_index(opt_pk_internal_t* opk,
		    opt_numint_t* q1, opt_numint_t* q2, unsigned short int *nz)
{
  size_t j;
  opt_numint_t prod = 0;
  unsigned short int size = nz[0];
  for (j=0; j<size; j++){
    unsigned short int j1 = nz[j+1];
    prod = prod + q1[j1]*q2[j1];
  }
  return prod;
}


bool is_vertex(opt_numint_t * v){
	return (v[0] && v[1]);
}

bool is_ray(opt_numint_t * v){
	return (v[0] && !v[1]);
}

bool is_line(opt_numint_t * v){
	return (!v[0] && !v[1]);
}


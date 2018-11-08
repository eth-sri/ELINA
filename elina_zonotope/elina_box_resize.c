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




#include "elina_box_internal.h"
#include "elina_box_resize.h"

elina_box_t* elina_box_forget_array(elina_manager_t* man,
                        bool destructive,
                        elina_box_t* a,
                        elina_dim_t* tdim,
                        size_t size,
                        bool project)
{
    elina_box_t* res;
    size_t i;
    
    man->result.flag_best = true;
    man->result.flag_exact = true;
    
    res = destructive ? a : elina_box_copy(man,a);
    if ((a->inf==NULL) && (a->sup==NULL)){
        return res;
    }
    if (project){
        for (i=0;i<size;i++){
            res->sup[tdim[i]]= 0;
            res->inf[tdim[i]] = 0;
        }
    }
    else {
        for (i=0;i<size;i++){
            res->inf[tdim[i]] = INFINITY;
            res->sup[tdim[i]] = INFINITY;
        }
    }
    return res;
}


elina_box_t* elina_box_add_dimensions(elina_manager_t* man,
			  bool destructive, elina_box_t* a,
			  elina_dimchange_t* dimchange,
			  bool project)
{
  elina_box_t* res;
  size_t size;
  size_t dimsup;
  int i,k;

  man->result.flag_best = true;  
  man->result.flag_exact = true;  
  res = destructive ? a : elina_box_copy(man,a);
  if ((a->inf==NULL) && (a->sup==NULL)){
    goto elina_box_add_dimensions_exit;
  }
  size = res->intdim+res->realdim;
  dimsup = dimchange->intdim+dimchange->realdim;
  res->inf = (double *)realloc(res->inf,(size+dimsup+1)*sizeof(double));
  res->sup = (double *)realloc(res->sup,(size+dimsup+1)*sizeof(double));
  k = dimsup;
  for (i=(int)size; i>=0; i--){
    if (i<(int)size){
      res->inf[i+k] = a->inf[i];
      res->sup[i+k] = a->sup[i];
    }
    while (k>=1 && dimchange->dim[k-1]==(elina_dim_t)i){
      k--;
      if (project){
	res->inf[i+k] = 0.0;
	res->sup[i+k] = 0.0;
      }
      else {
	res->inf[i+k] = INFINITY;
	res->sup[i+k] = INFINITY;
      }
    }
  }  
 elina_box_add_dimensions_exit:
  res->intdim = a->intdim+dimchange->intdim;
  res->realdim = a->realdim+dimchange->realdim;
  return res;
}

elina_box_t* elina_box_remove_dimensions(elina_manager_t* man,
			     bool destructive, elina_box_t* a,
			     elina_dimchange_t* dimchange)
{
  elina_box_t* res;
  size_t size;
  size_t dimsup;
  size_t i,k;
 
  man->result.flag_best = true;  
  man->result.flag_exact = true;  
  res = destructive ? a : elina_box_copy(man,a);
  if ((a->inf==NULL) && (a->sup==NULL)){
    goto elina_box_remove_dimensions_exit;
  }
  size = res->intdim+res->realdim;
  dimsup = dimchange->intdim+dimchange->realdim;
  k=0;
  for (i=0; i<size-dimsup; i++){
    while (k<dimsup && dimchange->dim[k]==i+k){
      k++;
    }
    res->inf[i] = a->inf[i+k];
    res->sup[i] = a->sup[i+k];
  }
  res->inf[size-dimsup] = 0.0;
  res->sup[size-dimsup] = 0.0;
  res->inf = (double*)realloc(res->inf,(size-dimsup+1)*sizeof(double));
  res->sup = (double*)realloc(res->sup,(size-dimsup+1)*sizeof(double));
 elina_box_remove_dimensions_exit:
  res->intdim = a->intdim-dimchange->intdim;
  res->realdim = a->realdim-dimchange->realdim;
  return res;
}

elina_box_t* elina_box_permute_dimensions(elina_manager_t* man,
                              bool destructive,
                              elina_box_t* a,
                              elina_dimperm_t* perm)
{
    elina_box_t* res;
    size_t size;
    size_t i;
    
    man->result.flag_best = true;
    man->result.flag_exact = true;
    if ((a->inf==NULL) && (a->sup==NULL)){
        return destructive ? a : elina_box_copy(man,a);
    }
    res = elina_box_copy(man,a);
    size = res->intdim+res->realdim;
    for (i=0;i<size;i++){
	res->inf[perm->dim[i]] = a->inf[i];
	res->sup[perm->dim[i]] = res->sup[i];
    }
    if (destructive) elina_box_free(man,a);
    return res;
}


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

#include "elina_box_internal.h"
#include "elina_box_resize.h"

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
  if (a->p == NULL) {
    goto elina_box_add_dimensions_exit;
  }
  size = res->intdim+res->realdim;
  dimsup = dimchange->intdim+dimchange->realdim;
  res->p = (elina_interval_t **)realloc(res->p, (size + dimsup + 1) *
                                                    sizeof(elina_interval_t *));
  for (i = (int)size + 1; i < (int)(size + dimsup + 1); i++) {
    res->p[i] = elina_interval_alloc();
  }
  k = dimsup;
  for (i=(int)size; i>=0; i--){
    if (i<(int)size){
      elina_interval_set(res->p[i + k], a->p[i]);
    }
    while (k>=1 && dimchange->dim[k-1]==(elina_dim_t)i){
      k--;
      if (project){
        elina_scalar_set_int(res->p[i + k]->inf, 0);
        elina_scalar_set_int(res->p[i + k]->sup, 0);
      }
      else {
        elina_interval_set_top(res->p[i + k]);
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
  if (a->p == NULL) {
    goto elina_box_remove_dimensions_exit;
  }
  size = res->intdim+res->realdim;
  dimsup = dimchange->intdim+dimchange->realdim;
  k=0;
  for (i=0; i<size-dimsup; i++){
    while (k<dimsup && dimchange->dim[k]==i+k){
      k++;
    }
    elina_interval_set(res->p[i], a->p[i + k]);
  }
  elina_interval_set_int(res->p[size - dimsup], 0, 0);
  for (i = size - dimsup + 1; i < size + 1; i++) {
    elina_interval_free(res->p[i]);
  }
  res->p = (elina_interval_t **)realloc(res->p, (size - dimsup + 1) *
                                                    sizeof(elina_interval_t *));
 elina_box_remove_dimensions_exit:
  res->intdim = a->intdim-dimchange->intdim;
  res->realdim = a->realdim-dimchange->realdim;
  return res;
}

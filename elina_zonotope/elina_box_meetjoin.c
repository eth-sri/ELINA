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

#include <string.h>
#include <stdio.h>

#include "elina_box_internal.h"
#include "elina_box_representation.h"
#include "elina_box_constructor.h"
#include "elina_box_meetjoin.h"
#include "elina_scalar_arith.h"

/* ============================================================ */
/* Meet and Join */
/* ============================================================ */



elina_box_t* elina_box_join(elina_manager_t* man, bool destructive, elina_box_t* a1, elina_box_t* a2)
{
  size_t i;
  size_t nbdims;
  elina_box_t* res;
  man->result.flag_best = true;
  man->result.flag_exact = false;
  res = destructive ? a1 : elina_box_alloc(a1->intdim,a1->realdim);
  if (a1->p == NULL) {
    if (a2->p != NULL) {
      man->result.flag_exact = true;
      elina_box_set(res,a2);
    }
    return res;
  } else if (a2->p == NULL) {
    man->result.flag_exact = true;
    if (!destructive) elina_box_set(res,a1);
    return res;
  }
  man->result.flag_exact = false;
  if (!destructive){
    elina_box_init(res);
  }
  
  nbdims = a1->intdim + a2->realdim;
  for (i=0; i<nbdims; i++){
    elina_scalar_min(res->p[i]->inf, a1->p[i]->inf, a2->p[i]->inf);
    elina_scalar_max(res->p[i]->sup, a1->p[i]->sup, a2->p[i]->sup);
  }
   
  return res;
}



/* ============================================================ */
/* Meet_lincons */
/* ============================================================ */

elina_box_t* elina_box_meet_lincons_array(elina_manager_t* man,
			      bool destructive,
			      elina_box_t* a,
			      elina_lincons0_array_t* array)
{
  elina_box_t* res;
  size_t kmax;
  elina_lincons0_array_t tlincons;
  elina_box_internal_t* intern = (elina_box_internal_t*)man->internal;

  res = destructive ? a : elina_box_copy(man,a);
  if (a->p == NULL) {
    man->result.flag_best = true;
    man->result.flag_exact = true;
  } else {

    man->result.flag_best = array->size==1;
    man->result.flag_exact = false;
    kmax = man->option.funopt[ELINA_FUNID_MEET_LINCONS_ARRAY].algorithm;
    if (kmax<1) kmax=2;
    tlincons = elina_lincons0_array_make(array->size);
    for(size_t i =0; i < array->size; i++){
	tlincons.p[i] = elina_lincons0_copy(&array->p[i]);
    }
    char tb = elina_lincons0_array_reduce_integer(&tlincons,a->intdim,ELINA_SCALAR_DOUBLE);
    if (tb==0){
      goto _elina_box_meet_lincons_array_bottom;
    }
    elina_boxize_lincons0_array(res->p, NULL, &tlincons, res->p, a->intdim,
                                kmax, false, ELINA_SCALAR_DOUBLE);
    if (elina_interval_is_bottom(res->p[0])) {
    _elina_box_meet_lincons_array_bottom:
      elina_box_set_bottom(res);
    }
    elina_lincons0_array_clear(&tlincons);
  }
  return res;
}




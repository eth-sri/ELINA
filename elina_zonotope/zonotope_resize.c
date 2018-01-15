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

#include "zonotope_resize.h"


/*********************/
/* Add dimensions */
/*********************/
zonotope_t* zonotope_add_dimensions(elina_manager_t* man, bool destructive, zonotope_t* z, elina_dimchange_t* dimchange, bool project)
{

    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ADD_DIMENSIONS);

    zonotope_t* res = destructive ? z : zonotope_copy(man, z);
    size_t intdim = z->intdim + dimchange->intdim;
    size_t dims = z->dims + (dimchange->intdim + dimchange->realdim);
    res->box = realloc(res->box, dims * sizeof(elina_interval_t *));
    res->paf = realloc(res->paf,dims*sizeof(zonotope_aff_t*));
    size_t i = 0;
    int j = 0;
    for (i=0; i<dimchange->intdim + dimchange->realdim; i++) {
      if (res->dims == dimchange->dim[i]) {
        /* add in the last place */
        res->box[res->dims] = elina_interval_alloc();
      } else {
        /* increment */
        for (j = (int)-1 + res->dims; j >= (int)dimchange->dim[i]; j--) {
          res->box[j + 1] = elina_interval_alloc();
          res->paf[j + 1] = res->paf[j];
          elina_interval_set(res->box[j + 1], res->box[j]);
        }
      }
      res->paf[dimchange->dim[i]] =
          project ? zonotope_aff_alloc_init(pr) : pr->top;
      res->paf[dimchange->dim[i]]->pby++;
      if (project)
        elina_interval_set_int(res->box[dimchange->dim[i]], 0, 0);
      else
        elina_interval_set_top(res->box[dimchange->dim[i]]);
      res->dims++;
    }
    res->intdim = intdim;

    return res;
}

/*********************/
/* Remove dimensions */
/*********************/
zonotope_t* zonotope_remove_dimensions(elina_manager_t* man, bool destructive, zonotope_t* z, elina_dimchange_t* dimchange)
{

    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_REMOVE_DIMENSIONS);

    zonotope_t *res = destructive ? z : zonotope_copy(man, z);
    size_t i = 0;
    size_t j = 0;
    for (i = 0; i < dimchange->intdim + dimchange->realdim; i++) {

      zonotope_aff_check_free(pr, res->paf[dimchange->dim[i]]);
      res->paf[dimchange->dim[i]] = NULL;
      for (j = dimchange->dim[i]; j < -1 + res->dims; j++) {
        res->paf[j] = res->paf[j + 1];
        elina_interval_set(res->box[j], res->box[j + 1]);
      }
    }
    res->intdim = z->intdim - dimchange->intdim;
    res->dims = z->dims - (dimchange->intdim + dimchange->realdim);
    res->box = realloc(res->box, res->dims * sizeof(elina_interval_t *));
    res->paf = realloc(res->paf, res->dims * sizeof(zonotope_aff_t *));

    return res;
}


zonotope_t* zonotope_permute_dimensions(elina_manager_t* man, bool destructive, zonotope_t* z, elina_dimperm_t* permutation)
{
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_PERMUTE_DIMENSIONS);
    zonotope_aff_t *tmp = NULL;
    zonotope_t *res = destructive ? z : zonotope_copy(man, z);
    size_t i = 0;
    for (i = 0; i < permutation->size; i++) {
      tmp = res->paf[i];
      res->paf[i] = res->paf[permutation->dim[i]];
      res->paf[permutation->dim[i]] = tmp;
      elina_interval_swap(res->box[i], res->box[permutation->dim[i]]);
    }
    return res;
}

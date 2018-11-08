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

#ifndef _ELINA_BOX_INTERNAL_H_
#define _ELINA_BOX_INTERNAL_H_

#include <string.h>
#include <stdio.h>


#include "elina_box.h"

#ifdef __cplusplus
extern "C" {
#endif

struct elina_box_t {
  double * inf;
  double * sup;
  size_t intdim;
  size_t realdim;
};

elina_manager_t* elina_box_manager_alloc(void);

typedef struct elina_box_internal_t {
  elina_interval_t * bound_linexpr_internal_itv;
  elina_interval_t * bound_linexpr_internal_itv2;
  elina_interval_t * bound_linexpr_itv;
  elina_interval_t * meet_lincons_internal_itv;
  elina_interval_t * meet_lincons_internal_itv2;
  elina_interval_t * meet_lincons_internal_itv3;
  elina_scalar_t* meet_lincons_internal_bound;
} elina_box_internal_t;

void elina_box_internal_init(elina_box_internal_t* intern);
void elina_box_internal_clear(elina_box_internal_t* intern);

elina_box_internal_t* elina_box_internal_alloc(void);
void elina_box_internal_free(elina_box_internal_t* intern);

/* Initializes some fields of pk from manager */
static inline elina_box_internal_t* elina_box_init_from_manager(elina_manager_t* man, elina_funid_t funid)
{
  elina_box_internal_t* itv = (elina_box_internal_t*)man->internal;
  return itv;
}

#ifdef __cplusplus
}
#endif

#endif

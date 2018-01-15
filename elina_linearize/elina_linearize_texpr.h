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

/* ************************************************************************* */
/* elina_linearize_aux.h: auxiliary functions for (quasi)linearisation */
/* ************************************************************************* */


/* Auxiliary module to elina_linearize, which contains functions depending of the
   number representation */

#ifndef _ELINA_LINEARIZE_TEXPR_H_
#define _ELINA_LINEARIZE_TEXPR_H_

#if defined (HAS_APRON)
#include "apron_wrapper.h"
#else
#include "elina_manager.h"
#include "elina_linexpr0.h"
#include "elina_abstract0.h"
#include "elina_linearize.h"
#endif

#include "elina_rat.h"
#include "elina_linearize.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  elina_interval_t *ulp;        /* [-1;1] * unit in the least place */
  elina_interval_t *min;        /* [-1;1] * minimum positive denormal */
  elina_interval_t *min_normal; /* [-1;1] * minimum positive normal */
  elina_interval_t *max;        /* [-1;1] * maximum non +oo  */
  elina_interval_t *max_exact;  /* [-1;1] * maximum exactly representable integer */
} elina_float_const;

/* ********************************************************************** */
/* I. Evaluation of tree expressions */
/* ********************************************************************** */

elina_interval_t* elina_eval_texpr0(elina_manager_t* man,
			      elina_abstract0_t* abs,
			      elina_texpr0_t* expr,
			      elina_scalar_discr_t discr,
			      bool* pexact);

/* ********************************************************************** */
/* II. Interval linearization of linear tree expressions */
/* ********************************************************************** */

/* Linearize a tree expression that is (syntaxically) interval linear with
   exact arithmetic.

   Compared to elina_intlinearize_texpr0() function below, this functions does
   not require a bounding box for dimensions.

   If the precondition is violated, returns NULL.
*/

bool elina_intlinearize_elina_texpr0_intlinear(elina_linexpr0_t** res, elina_texpr0_t* expr, elina_scalar_discr_t discr);

elina_linexpr0_t* elina_intlinearize_texpr0_intlinear(elina_manager_t* man,
						elina_texpr0_t* expr,
						elina_scalar_discr_t discr);

bool elina_boxize_lincons0_array(elina_interval_t** res,bool* tchange,
				      elina_lincons0_array_t* array,
				      elina_interval_t** env, size_t intdim,
				      size_t kmax,
				      bool intervalonly, elina_scalar_discr_t discr);

/* ********************************************************************** */
/* III. Interval linearization of tree expressions */
/* ********************************************************************** */

elina_texpr_rtype_t elina_interval_intlinearize_texpr0_rec(elina_texpr0_t* expr,
			    elina_interval_t** env, size_t intdim,
			    elina_linexpr0_t** lres /* out */, elina_interval_t *ires /* out */
			    ,elina_scalar_discr_t discr);

elina_linexpr0_t* elina_intlinearize_texpr0(elina_manager_t* man,
				      elina_abstract0_t* abs,
				      elina_texpr0_t* expr,
				      bool* pexact,
				      elina_scalar_discr_t discr,
				      bool quasilinearize);

elina_linexpr0_t** elina_intlinearize_texpr0_array(elina_manager_t* man,
					     elina_abstract0_t* abs,
					     elina_texpr0_t** texpr, size_t size,
					     bool* pexact,
					     elina_scalar_discr_t discr,
					     bool quasilinearize);

elina_lincons0_t elina_intlinearize_tcons0(elina_manager_t* man,
				     elina_abstract0_t* abs,
				     elina_tcons0_t* cons,
				     bool* pexact,
				     elina_scalar_discr_t discr,
				     bool quasilinearize, bool meet);

bool elina_intlinearize_elina_tcons0(elina_lincons0_t* res,
				   elina_tcons0_t* cons,
				   elina_interval_t** env, size_t intdim, elina_scalar_discr_t discr);

elina_lincons0_array_t elina_intlinearize_tcons0_array(elina_manager_t* man,
						 elina_abstract0_t* abs,
						 elina_tcons0_array_t* array,
						 bool* pexact,
						 elina_scalar_discr_t discr,
						 elina_linexpr_type_t type, bool meet,
						 bool boxize, size_t kmax, bool intervalonly);

#ifdef __cplusplus
}
#endif

#endif

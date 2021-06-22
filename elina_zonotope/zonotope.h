/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
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


#ifndef _ZONOTOPE_H_
#define _ZONOTOPE_H_

#include "elina_int.h"
#include "elina_rat.h"
#include "comp_list.h"


#ifdef __cplusplus
extern "C" {
#endif

#if defined (HAS_APRON)
#include "apron_wrapper.h"
//#include "num.h"
//#include "numint.h"
//#include "numrat.h"
//#include "bound.h"

//#include "itv.h"
//#include "itv_linexpr.h"
//#include "itv_linearize.h"

#else
#include "elina_coeff.h"
#include "elina_dimension.h"
#include "elina_linexpr0.h"
#include "elina_texpr0.h"
#include "elina_lincons0.h"
#include "elina_tcons0.h"
#include "elina_manager.h"
#include "elina_abstract0.h"
#endif

#include "elina_box.h"
#include "elina_box_meetjoin.h"
#include "elina_generic.h"
#include "elina_linearize_texpr.h"

//ap_manager_t* t1p_manager_alloc(ap_manager_t* manNS);
elina_manager_t* zonotope_manager_alloc(void);
void zonotope_relu(elina_manager_t *man, elina_abstract0_t * abs, elina_dim_t y, elina_dim_t x);

typedef struct _zonotope_t zonotope_t;

#include "zonotope_representation.h"

/****************/
/* Constructors */
/****************/
/* 1.Basic constructors */
/* 2.Accessors */
/* 3.Tests */
/* 4.Extraction of properties */
#include "zonotope_constructor.h"

/* Meet and Join */
/*****************/
/* 1.Meet */
/* 2.Join */
#include "zonotope_meetjoin.h"

/* Assign and Substitute */
/*************************/
#include "zonotope_assign.h"

/* Resize dimensions */
/*********************/
#include "zonotope_resize.h"

/* Other functions */
/*******************/
#include "zonotope_otherops.h"

#ifdef __cplusplus
}
#endif

#endif

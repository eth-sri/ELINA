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


#ifndef __OPT_ZONES_CLOSURE_H__
#define __OPT_ZONES_CLOSURE_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "opt_zones_internal.h"

bool check_negative_cycle(opt_zones_mat_t * oz, unsigned short int dim);

bool closure_dense(opt_zones_mat_t * oz, unsigned short int dim);

bool closure_dense_scalar(opt_zones_mat_t *oz, unsigned short int dim);

bool strengthening_zones(opt_zones_mat_t *oz, unsigned short int dim);

void compute_sparse_index(double *m, unsigned short int *ca, 
			  unsigned short int comp_size, unsigned short int *index1,
			  unsigned short int * index2, unsigned short int k, unsigned short int dim);

void closure_comp_sparse(opt_zones_mat_t *oz, unsigned short int dim);

int update_bounds_comp(double *m, unsigned short int * ca, 
			unsigned short int comp_size, unsigned short int n);

void strengthening_comp_zones(opt_zones_mat_t *oz,comp_list_t *cd, unsigned short int dim,int sgn);

bool strengthening_single_comp_zones(opt_zones_mat_t * oz, comp_list_t * cl, unsigned short int dim);

bool strengthening_intra_comp_zones(opt_zones_mat_t * oz, unsigned short int dim);

void strengthening_inter_comp_zones(opt_zones_mat_t *oz, array_comp_list_t *acl, char *map , unsigned short int dim);

void strengthening_incr_init_zones(opt_zones_mat_t *oz, array_comp_list_t *acl, char *map , unsigned short int dim);

#ifdef __cplusplus
}
#endif

#endif 

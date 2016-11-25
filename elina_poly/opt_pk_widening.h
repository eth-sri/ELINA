#ifndef __OPT_PK_WIDENING_H
#define __OPT_PK_WIDENING_H

#include "opt_pk_config.h"
#include "opt_pk_internal.h"
#include "opt_pk.h"

#ifdef __cplusplus
extern "C" {
#endif

bool is_vectors_equal_comp_list(opt_pk_internal_t *opk, opt_numint_t * v1, 
				opt_numint_t * v2, unsigned short int * ca1, 
				unsigned short int * ca2, unsigned short int comp_size1, 
				unsigned short int comp_size2);

void vector_copy_comp_list(opt_pk_internal_t *opk, opt_numint_t * dst, opt_numint_t * src, 
			   unsigned short int * ind_map, unsigned short int comp_size);

#ifdef __cplusplus
}
#endif

#endif

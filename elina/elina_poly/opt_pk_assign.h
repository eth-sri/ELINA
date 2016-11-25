/* ********************************************************************** */
/* opt_pk_assign.h: Assignements and Substitutions */
/* ********************************************************************** */



#ifndef _OPT_PK_ASSIGN_H_
#define _OPT_PK_ASSIGN_H_

#include "opt_pk_config.h"
#include "opt_pk.h"

#ifdef __cplusplus
extern "C" {
#endif

opt_pk_array_t* opt_poly_asssub_linexpr_det(bool assign, elina_manager_t* man,
			      bool destructive,
			      opt_pk_array_t* oa,
			      elina_dim_t dim, elina_linexpr0_t* linexpr);

opt_pk_array_t* opt_poly_asssub_linexpr_array_det(elina_manager_t* man,
				    bool destructive,
				    opt_pk_array_t* pa,
				    elina_dim_t* tdim, elina_linexpr0_t** texpr, 
				    size_t size);

#ifdef __cplusplus
}
#endif

#endif

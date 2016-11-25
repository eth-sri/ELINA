/* ********************************************************************** */
/* pk_constructor.h: constructors and accessors */
/* ********************************************************************** */



#ifndef _OPT_PK_CONSTRUCTOR_H_
#define _OPT_PK_CONSTRUCTOR_H_

#include "opt_pk_config.h"
#include "opt_pk_vector.h"

#include "opt_pk_matrix.h"
#include "opt_pk.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Fill the first (pk->dec-1) rows of the matrix with the constraints of the
   universe polyhedron */
void opt_matrix_fill_constraint_top(opt_pk_internal_t* opk, opt_matrix_t* oc, size_t start);

/* Assign with GMP semantics the given polyhedron with the empty
   (resp. universe) polyhedron, of same dimensions */
void opt_poly_set_bottom(opt_pk_internal_t* opk, opt_pk_array_t* op);
void opt_poly_set_top(opt_pk_internal_t* opk, opt_pk_array_t* op);

#ifdef __cplusplus
}
#endif

#endif

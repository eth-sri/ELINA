#ifndef __MAXPOOL_CONVEX_HULL_H__
#define __MAXPOOL_CONVEX_HULL_H__


#ifdef __cplusplus
extern "C" {
#endif


#include "setoper.h"
#include "cdd.h"


dd_MatrixPtr maxpool_deeppoly_approx(double *lb, double * ub, size_t pool_size);

#ifdef __cplusplus
}
#endif

#endif

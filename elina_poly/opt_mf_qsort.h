/* ********************************************************************** */
/* mf_qsort.h: quicksort */
/* ********************************************************************** */

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#define QSORT_EXTRA_CMP_ARGUMENT

typedef int (*opt_qsort2_cmp)(void*, const void *,const void *);

void opt_qsort2(void *base_ptr, size_t count, size_t size,
	    opt_qsort2_cmp cmp,
	    void *cmp_argument);

#ifdef __cplusplus
}
#endif


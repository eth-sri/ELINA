#include <stdio.h>
#include "itv.h"
#include "elina_interval.h"
#include "elina_interval_arith.h"
#include "zonotope_internal.h"

int main(){
	double a = 0.79844003347607328536;
	double b = 0.026801820391231023777;
	double c = 0.026801820391231023777;
	elina_interval_t *eitv1 = elina_interval_alloc();
	elina_interval_set_double(eitv1,a,a);
	elina_interval_t *eitv2 = elina_interval_alloc();
	elina_interval_set_double(eitv2,b,c);	
	elina_interval_t *eitv3 = elina_interval_alloc();
	elina_interval_mul(eitv3,eitv1,eitv2,ELINA_SCALAR_DOUBLE);
	elina_interval_fprint(stdout,eitv3);

	elina_interval_mul_double(eitv3,b,c,a,a);
	elina_interval_fprint(stdout,eitv3);

	ap_interval_t *aitv1 = ap_interval_alloc();
	ap_interval_set_double(aitv1,a,a);
	ap_interval_t *aitv2 = ap_interval_alloc();
	ap_interval_set_double(aitv2,b,c);
	
	itv_internal_t* intern = itv_internal_alloc();
	itv_t itv1, itv2, itv3;
	itv_set_ap_interval(intern,itv1,aitv1);
	itv_set_ap_interval(intern,itv2,aitv2);
	itv_mul(intern,itv3,itv1,itv2);
	itv_fprint(stdout,itv3);
	printf("\nsize: %d\n",sizeof(bound_t));
	//printf("%.*g %d %d",20,a*b+0.0,sizeof(double),sizeof(long double));

}

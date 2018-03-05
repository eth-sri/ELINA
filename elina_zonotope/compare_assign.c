#include "zonotope.h"
#include "t1p.h"
#include "t1p_representation.h" 

elina_manager_t *elina_man;
elina_linexpr0_t *elinexpr0;
elina_interval_t ** eitv;
elina_dimchange_t *edimchange;
elina_dim_t *edim;

ap_manager_t *ap_man;
ap_linexpr0_t *alinexpr0;
ap_interval_t ** aitv;
ap_dimchange_t *adimchange;
ap_dim_t *adim;

unsigned short int dim;


void generate_random_box(){
	//elina_interval_t ** interval = (elina_interval_t **)malloc(dim*sizeof(elina_interval_t *));
	unsigned short int i;
	
	for(i = 0; i < dim; i++){
		//interval[i] = elina_interval_alloc();
		double inf = (double)rand()/RAND_MAX*2.0-1.0;
		double sup = inf + ((double)rand()/RAND_MAX*2.0);
		elina_interval_set_double(eitv[i],inf,sup);
		ap_interval_set_double(aitv[i],inf,sup);
		
	}
	//return interval;
}


void generate_random_linexpr0(){
	elina_coeff_t *ecst, *ecoeff;
	ap_coeff_t *acst, *acoeff;
	unsigned short int j, k;
	//elina_linexpr0_t * linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,2);
	ecst = &elinexpr0->cst;
	acst = &alinexpr0->cst;
	elina_scalar_set_double(ecst->val.scalar,0);
	ap_scalar_set_double(acst->val.scalar,0);
	for(j=0; j<dim; j++){
		elina_linterm_t * elinterm = &elinexpr0->p.linterm[j];
		ap_linterm_t *alinterm = &alinexpr0->p.linterm[j];
		double val = (double)rand()/RAND_MAX*2.0-1.0;
		ecoeff = &elinterm->coeff;
		acoeff = &alinterm->coeff;
		elinterm->dim = j;
		alinterm->dim = j;
		elina_scalar_set_double(ecoeff->val.scalar,val);
		ap_scalar_set_double(acoeff->val.scalar,val);
	}
}


void generate_dimchange(){
	edimchange->dim[0] = dim;
	adimchange->dim[0] = dim;
	edim[0]= dim;
	adim[0] = dim;
}

void test_assign_elina(){
	//generate_random_box();
	//generate_random_linexpr0();
	zonotope_t *z1 = zonotope_of_box(elina_man,0,dim,eitv);
	zonotope_t *z2 = zonotope_add_dimensions(elina_man,false,z1,edimchange,false);
	printf("Zonotope assign inputs\n");
	zonotope_fprint(stdout,elina_man,z2,NULL);
	printf("x%d:= ",edim[0]);
	elina_linexpr0_fprint(stdout,elinexpr0,NULL);
	printf("\n");
	fflush(stdout);
	zonotope_t *z3 = zonotope_assign_linexpr_array(elina_man,false,z2,edim,&elinexpr0,1,NULL);
	zonotope_t *z4 = zonotope_join(elina_man,false,z3,zonotope_copy(elina_man,z3));
	printf("Zonotope assign output %d\n",zonotope_is_eq(elina_man,z3,z4));
	zonotope_fprint(stdout,elina_man,z3,NULL);
	zonotope_fprint(stdout,elina_man,z4,NULL);
	fflush(stdout);
	zonotope_free(elina_man,z1);
	zonotope_free(elina_man,z2);
	zonotope_free(elina_man,z3);
}


void test_assign_apron(){
	
	t1p_t *t1 = t1p_of_box(ap_man,0,dim,aitv);
	t1p_t *t2 = t1p_add_dimensions(ap_man,false,t1,adimchange,false);
	printf("T1p assign inputs\n");
	t1p_fprint(stdout,ap_man,t2,NULL);
	printf("x%d:= ",adim[0]);
	ap_linexpr0_fprint(stdout,alinexpr0,NULL);
	printf("\n");
	fflush(stdout);
	t1p_t *t3 = t1p_assign_linexpr_array(ap_man,false,t2,adim,&alinexpr0,1,NULL);
	t1p_t *t4 = t1p_join(ap_man,false,t3,t1p_copy(ap_man,t3));
	printf("T1p assign output %d\n",t1p_is_leq(ap_man,t3,t4));
	t1p_fprint(stdout,ap_man,t3,NULL);
	t1p_fprint(stdout,ap_man,t4,NULL);
	fflush(stdout);
	t1p_free(ap_man,t1);
	t1p_free(ap_man,t2);
	t1p_free(ap_man,t3);
}

int main(int argc, char **argv){
	//srand(time(NULL));
	dim = atoi(argv[1]);

	eitv = elina_interval_array_alloc(dim);
	aitv = ap_interval_array_alloc(dim);
	generate_random_box();

	elinexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,dim);
	alinexpr0 = ap_linexpr0_alloc(AP_LINEXPR_SPARSE,dim);
	generate_random_linexpr0();
	
	edimchange = elina_dimchange_alloc(0,1);
	adimchange = ap_dimchange_alloc(0,1);	
	
	edim = (elina_dim_t *)malloc(sizeof(elina_dim_t));
	adim = (ap_dim_t *)malloc(sizeof(ap_dim_t));

	generate_dimchange();
	elina_man = zonotope_manager_alloc();
	ap_man = t1p_manager_alloc();
    	test_assign_elina();
	test_assign_apron();
	elina_manager_free(elina_man);
	ap_manager_free(ap_man);
}

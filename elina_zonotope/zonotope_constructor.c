#include "zonotope.h"
#include "zonotope_internal.h"
#include "zonotope_constructor.h"

/****************/
/* Constructors */
/****************/
/* 1.Basic constructors */
zonotope_t* zonotope_bottom(elina_manager_t* man, size_t intdim, size_t realdim)
{
    start_timing();
    size_t i;
    size_t dims = intdim + realdim;
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_BOTTOM);
    zonotope_t* res = zonotope_alloc(man, intdim, realdim);
    for (i=0; i<dims; i++) {
	res->paf[i] = pr->bot;
	res->paf[i]->pby++;
	elina_interval_set_bottom(res->box[i]);
    }
    man->result.flag_best = true;
    man->result.flag_exact = true;
    record_timing(zonotope_bottom_time);
    return res;
}


zonotope_t* zonotope_top(elina_manager_t* man, size_t intdim, size_t realdim)
{
    start_timing();
    size_t i;
    size_t dims = intdim + realdim;
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_TOP);
    zonotope_t* res = zonotope_alloc(man, intdim, realdim);
    for (i=0; i<dims; i++) {
	res->paf[i] = pr->top;
	res->paf[i]->pby++;
	elina_interval_set_top(res->box[i]);
    }
    man->result.flag_best = true;
    man->result.flag_exact = true;
    record_timing(zonotope_top_time);
    return res;
}

/* Abstract an hypercube defined by the array of intervals of size intdim+realdim */
zonotope_t* zonotope_of_box(elina_manager_t* man, size_t intdim, size_t realdim, elina_interval_t** tinterval)
{
	//printf("input intervals\n");
	
    start_timing();
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_OF_BOX);
    zonotope_t* res = zonotope_alloc(man,intdim,realdim);
    size_t i = 0;
    for (i=0; i<intdim+realdim; i++) {
	//elina_interval_fprint(stdout,tinterval[i]);
	elina_interval_set(res->box[i], tinterval[i]);
	res->paf[i] = zonotope_aff_alloc_init(pr);
	if (elina_interval_is_bottom(tinterval[i])){
		 res->paf[i] = pr->bot;
	}
	else if (elina_interval_is_top(tinterval[i])){ 
		res->paf[i] = pr->top;
	}
	else if (!elina_scalar_infty(tinterval[i]->inf) && elina_scalar_equal(tinterval[i]->inf, tinterval[i]->sup)) {
		elina_interval_set(res->paf[i]->c, tinterval[i]);
	}
	else if (!(elina_scalar_infty(tinterval[i]->sup) || elina_scalar_infty(tinterval[i]->inf))){
		 zonotope_aff_add_itv(pr, res->paf[i], tinterval[i], IN);
	}
	else{
	 elina_interval_set(res->paf[i]->c, tinterval[i]);
	}
        //elina_interval_set(res->paf[i]->itv,tinterval[i]);
	res->paf[i]->pby++;
    }
    
    man->result.flag_best = true;
    man->result.flag_exact = true;
	//fflush(stdout);
    //printf("of box output\n");
    //zonotope_fprint(stdout,man,res,NULL);
    //fflush(stdout);
    record_timing(zonotope_of_box_time);
    return res;
}

elina_dimension_t zonotope_dimension(elina_manager_t* man, zonotope_t* z){
    elina_dimension_t res;
    res.intdim = z->intdim;
    res.realdim = z->dims - z->intdim;
    man->result.flag_best = true;
    man->result.flag_exact = true;
    return res;
}

/* 3.Tests */
bool zonotope_is_bottom(elina_manager_t* man, zonotope_t* z)
{
    start_timing();
    size_t i;
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_IS_BOTTOM);
    bool res = elina_interval_is_bottom(z->box[0]);
    if(res){
	return true;
    }
    for (i=1; i<z->dims; i++) {
	res = elina_interval_is_bottom(z->box[i]);
	if(res){
		return true;
	}
    }
    man->result.flag_best = true;
    man->result.flag_exact = true;
    record_timing(zonotope_is_bottom_time);
    return false;
}

bool zonotope_is_top(elina_manager_t* man, zonotope_t* z)
{
    start_timing();
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_IS_TOP);
    size_t i;
    bool res = elina_interval_is_top(z->box[0]);
    if(!res){
	return false;
    }
    for (i=1; i<z->dims; i++) {
	res = elina_interval_is_top(z->box[i]);
	if(!res){
		return false;
	}
    }
    man->result.flag_best = true;
    man->result.flag_exact = true;
    record_timing(zonotope_is_top_time);
    return true;
}

bool zonotope_is_eq(elina_manager_t* man, zonotope_t* z1, zonotope_t* z2)
{
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_IS_EQ);
    if(!z1 || !z2 || z1->dims != z2->dims){
	return false;
    }
    man->result.flag_best = true;
    man->result.flag_exact = true;
    elina_dimension_t dim1 = elina_abstract0_dimension(pr->manNS, z1->abs);
    elina_dimension_t dim2 = elina_abstract0_dimension(pr->manNS, z2->abs);
    if ( (dim1.intdim != dim2.intdim) || (dim1.realdim != dim2.realdim) ) return false;
    else if (!elina_abstract0_is_eq(pr->manNS, z1->abs, z2->abs)) return false;
    else if (z1 == z2) return true;
    else {
        start_timing();
	size_t i = 0;
	bool res = true;
	for (i=0; i<z1->dims; i++) {
	    if (z1->paf[i] != z2->paf[i] && elina_interval_equal(z1->box[i], z2->box[i])) {
		res = zonotope_aff_is_eq(pr, z1->paf[i], z2->paf[i]);
		if (!res) break;
	    } else {
		res = false;
		break;
	    }
	}
        record_timing(zonotope_is_equal_time);
	return res;
    }
}

elina_interval_t** zonotope_to_box(elina_manager_t* man, zonotope_t* z){
    start_timing();
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_TO_BOX);
    elina_interval_t **res = elina_interval_array_alloc(z->dims);
    size_t i = 0;
    for (i=0; i<z->dims; i++) {
	elina_interval_set(res[i], z->box[i]);
    }
    man->result.flag_best = true;
    man->result.flag_exact = true;
    record_timing(zonotope_to_box_time);
    return res;

}


elina_lincons0_array_t zonotope_to_lincons_array(elina_manager_t* man, zonotope_t* z)
    /* Same as APRON */
    /* TODO: use constraints in eps domain to deduce constraints on variable ? */
{
    size_t i;
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_TO_LINCONS_ARRAY);
    elina_lincons0_array_t array;

    size_t nbdims = z->dims;

    man->result.flag_best = true;
    man->result.flag_exact = true;
    if (nbdims==0){
	array = elina_lincons0_array_make(0);
    }
    else if (zonotope_is_bottom(man,z)){
	array = elina_lincons0_array_make(1);
	array.p[0] = elina_lincons0_make_unsat();
    }
    else {
	size_t size;
	elina_linexpr0_t* expr;
	elina_scalar_t* scalar;
	bool point;

	size = 0;
	for (i=0;i<nbdims;i++){
	    if (!elina_scalar_infty(z->box[i]->inf)) size++;
	    point = !elina_scalar_infty(z->box[i]->inf) && elina_scalar_equal(z->box[i]->inf,z->box[i]->sup);
	    if (!point && !elina_scalar_infty(z->box[i]->sup)) size++;
	}
	array = elina_lincons0_array_make(size);
	size = 0;
	for (i=0;i<nbdims;i++){
	    point = false;
	    
	    if (!elina_scalar_infty(z->box[i]->inf)){
		expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
		elina_coeff_set_scalar_int(&expr->p.linterm[0].coeff, 1);
		expr->p.linterm[0].dim = i;

		elina_coeff_reinit(&expr->cst,ELINA_COEFF_SCALAR,ELINA_SCALAR_DOUBLE);
		scalar = expr->cst.val.scalar;
		elina_scalar_set(scalar,z->box[i]->inf);

		point = !elina_scalar_infty(z->box[i]->inf) && elina_scalar_equal(z->box[i]->inf,z->box[i]->sup);
		array.p[size].constyp = point ? ELINA_CONS_EQ : ELINA_CONS_SUPEQ;
		array.p[size].linexpr0 = expr;
		size++;
	    }
	    if (!point && !elina_scalar_infty(z->box[i]->sup)){
		expr = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,1);
		elina_coeff_set_scalar_int(&expr->p.linterm[0].coeff, -1);
		expr->p.linterm[0].dim = i;

		elina_coeff_reinit(&expr->cst,ELINA_COEFF_SCALAR,ELINA_SCALAR_DOUBLE);
		elina_scalar_set(expr->cst.val.scalar,z->box[i]->sup);
		array.p[size].constyp = ELINA_CONS_SUPEQ;
		array.p[size].linexpr0 = expr;
		size++;
	    }
	}
     }
	return array;
}

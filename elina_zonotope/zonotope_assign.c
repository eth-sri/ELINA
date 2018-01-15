#include "zonotope_assign.h"

zonotope_t* zonotope_assign_linexpr_array(elina_manager_t* man,
		bool destructive,
		zonotope_t* z,
		elina_dim_t* tdim, elina_linexpr0_t** lexpr,
		size_t size,
		zonotope_t* dest)
{
    zonotope_internal_t* pr = zonotope_init_from_manager(man, ELINA_FUNID_ASSIGN_LINEXPR_ARRAY);
    zonotope_t* res = zonotope_copy(man, z);
    size_t i = 0;
    for (i=0; i<res->dims; i++) {
      elina_interval_set(res->paf[i]->itv, res->box[i]);
    }
    for (i=0; i<size; i++) {
	zonotope_aff_check_free(pr, res->paf[tdim[i]]);
	res->paf[tdim[i]] = zonotope_aff_from_linexpr0(pr, lexpr[i], z);
        zonotope_aff_reduce(pr, res->paf[tdim[i]]);
        if (zonotope_aff_is_top(pr, res->paf[tdim[i]])) {
	    zonotope_aff_check_free(pr, res->paf[tdim[i]]);
	    res->paf[tdim[i]] = pr->top;
	} else if (zonotope_aff_is_bottom(pr, res->paf[tdim[i]])) {
	    zonotope_aff_check_free(pr, res->paf[tdim[i]]);
	    res->paf[tdim[i]] = pr->bot;
	}
        elina_interval_set(res->box[tdim[i]], res->paf[tdim[i]]->itv);
        res->paf[tdim[i]]->pby++;
    }
    man->result.flag_best = false;
    man->result.flag_exact = false;
    return res;
}


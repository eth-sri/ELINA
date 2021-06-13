/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
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

#include "elina_box_assign.h"

elina_box_t* elina_box_assign_linexpr_array(elina_manager_t* man,
                                bool destructive,
                                elina_box_t* a,
                                elina_dim_t* tdim,
                                elina_linexpr0_t** linexpr,
                                size_t size,
                                elina_box_t* dest)
{
    bool exact;
    size_t i;
    elina_box_t* res;
    
    exact = true;
    if ((a->inf==NULL && a->sup==NULL) || (dest && dest->inf==NULL && dest->sup==NULL)){
        man->result.flag_best = true;
        man->result.flag_exact = true;
        return destructive ? a : elina_box_copy(man,a);
    }
    if (size==1){
        res = destructive ? a : elina_box_copy(man,a);
	exact = elina_double_interval_eval_elina_linexpr0(&res->inf[tdim[0]],&res->sup[tdim[0]],linexpr[0],a->inf, a->sup, ELINA_SCALAR_DOUBLE);
    }
    else {
        res = elina_box_copy(man,a);
        for (i=0;i<size;i++){
	     exact = elina_double_interval_eval_elina_linexpr0(&res->inf[tdim[i]],&res->sup[tdim[i]],linexpr[i],a->inf, a->sup, ELINA_SCALAR_DOUBLE) && exact;
        }
        if (destructive) elina_box_free(man,a);
    }
    if (dest)
        res = elina_box_meet(man,true,res,dest);
    man->result.flag_best = size==1 && exact;
    man->result.flag_exact = false;
    return res;
}






/*
	Copyright 2015 Department of Computer Science, ETH Zurich

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

		http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
*/


#ifndef __OPT_OCT_INCR_CLOSURE_COMP_SPARSE_H
#define __OPT_OCT_INCR_CLOSURE_COMP_SPARSE_H

#include "opt_oct_hmat.h"
#include "comp_list.h"
#include "opt_oct_closure_comp_sparse.h"

bool incremental_closure_comp_sparse(opt_oct_mat_t *oo, int dim, int v, bool is_int);

#endif

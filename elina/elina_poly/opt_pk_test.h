/*
	Copyright 2016 Software Reliability Lab, ETH Zurich

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

#ifndef __OPT_PK_TEST_H
#define __OPT_PK_TEST_H

#include "opt_pk_config.h"
#include "opt_pk_internal.h"
#include "opt_pk.h"

#ifdef __cplusplus
extern "C" {
#endif

bool opt_poly_leq(opt_pk_internal_t * opk, opt_matrix_t * C, opt_matrix_t * F);

#ifdef __cplusplus
}
#endif

#endif

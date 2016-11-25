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

#ifndef __APXX_OPT_PK_HH
#define __APXX_OPT_PK_HH

#include "apxx_manager.hh"

#include "opt_pk.h"

namespace apron {

//! Manager factory for the ELINA polyhedra library.
class opt_pk_manager : public manager {

public:

  /*! \brief Creates a new manager.
   *
   * \arg \c strict whether to allow strict inequalities as well as non-strict inequalities.
   * We only alloc non strict inequaities in the Polyhedra.
   */
  opt_pk_manager(bool strict = false);

  //! Copy operator.
  manager& operator=(const manager&);
};

#include "apxx_opt_pk_inline.hh"

}

#endif /* __APXX_OPT_PK_HH */

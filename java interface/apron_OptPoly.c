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

#include "japron.h"
#include "apron_OptPoly.h"
#include "opt_pk.h"

//////////////////////////////////////

/*
 * Class:     apron_optPoly
 * Method:    init
 * Signature: (Z)V
 */
JNIEXPORT void JNICALL Java_apron_OptPoly_init
  (JNIEnv *env, jobject o, jboolean strict)
{
  check_nonnull(o,);
  ap_manager_t* m = opt_pk_manager_alloc(strict);
  if (!m) { out_of_memory("cannot create manager"); return; }
  japron_manager_setup(m);
  set_manager(o, m);
}

/*
 * Class:     apron_OptPoly
 * Method:    class_init
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_apron_OptPoly_class_1init
  (JNIEnv *env, jclass cls)
{
  japron_cache(env);
}


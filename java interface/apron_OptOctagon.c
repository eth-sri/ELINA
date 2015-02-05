/*
	Copyright 2015 Software Reliability Lab, ETH Zurich

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
#include "apron_OptOctagon.h"
#include "opt_oct.h"

//////////////////////////////////////

/*
 * Class:     apron_OptOctagon
 * Method:    init
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_apron_OptOctagon_init
  (JNIEnv *env, jobject o)
{
  check_nonnull(o,);
  ap_manager_t* m = opt_oct_manager_alloc();
  if (!m) { out_of_memory("cannot create manager"); return; }
  japron_manager_setup(m);
  set_manager(o, m);
}

/*
 * Class:     apron_Octagon
 * Method:    class_init
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_apron_OptOctagon_class_1init
  (JNIEnv *env, jclass cls)
{
  japron_cache(env);
}


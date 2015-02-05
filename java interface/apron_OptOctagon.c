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


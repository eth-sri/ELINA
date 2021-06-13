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


/* ************************************************************************* */
/* elina_config.h */
/* ************************************************************************* */

#ifndef _ELINA_CONFIG_H_
#define _ELINA_CONFIG_H_

#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
#define HAS_BOOL
extern "C" {
#endif

#ifndef HAS_BOOL
#define HAS_BOOL
typedef char bool;
static const bool false = 0;
static const bool true  = 1;
#endif

#if !(defined __USE_SVID || defined __USE_BSD || defined __USE_XOPEN_EXTENDED || defined __APPLE__ || defined __CYGWIN__)

#endif

#ifdef __cplusplus
}
#endif

#endif

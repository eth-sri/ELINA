#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
#  This software is distributed under GNU Lesser General Public License Version 3.0.
#  For more information, see the ELINA project website at:
#  http://elina.ethz.ch
#
#  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
#  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
#  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
#  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
#  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY     
#  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
#  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
#  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
#  CONTRACT, TORT OR OTHERWISE).
#
#


from opt_pk_imports import *
from elina_manager_h import *

# ====================================================================== #
# Basics
# ====================================================================== #

def opt_pk_manager_alloc(strict):
    """
    Allocates an ElinaManager.
    
    Parameters
    ----------
    strict : c_bool
        if strict inequalities are required.    

    Returns
    -------
    man : ElinaManagerPtr
        Pointer to the newly allocated ElinaManager.

    """

    man = None
    try:
        opt_pk_manager_alloc_c = opt_pk_api.opt_pk_manager_alloc
        opt_pk_manager_alloc_c.restype = ElinaManagerPtr
        opt_pk_manager_alloc_c.argtypes = [c_bool]
        man = opt_pk_manager_alloc_c(strict)
    except:
        print('Problem with loading/calling "opt_pk_manager_alloc" from "liboptpoly.so"')

    return man

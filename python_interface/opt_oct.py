from opt_oct_imports import *
from elina_manager_h import *

# ====================================================================== #
# Basics
# ====================================================================== #

def opt_oct_manager_alloc():
    """
    Allocates an ElinaManager.

    Returns
    -------
    man : ElinaManagerPtr
        Pointer to the newly allocated ElinaManager.

    """

    man = None
    try:
        opt_oct_manager_alloc_c = opt_oct_api.opt_oct_manager_alloc
        opt_oct_manager_alloc_c.restype = ElinaManagerPtr
        opt_oct_manager_alloc_c.argtypes = None
        man = opt_oct_manager_alloc_c()
    except:
        print('Problem with loading/calling "opt_oct_manager_alloc" from "liboptoct.so"')

    return man

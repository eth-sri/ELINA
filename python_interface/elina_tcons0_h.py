# ************************************************************************* #
# elina_tcons0.h: tree expressions constraints and arrays
# ************************************************************************* #

from elina_coeff_h import *
from elina_texpr0_h import *
from elina_lincons0_h import *


# ======================================================================
# Datatypes */
# ======================================================================

class ElinaTcons0(Structure):
    """
    ElinaTcons0 ctype compatible with elina_tcons0_t from elina_tcons0.h.

    Fields
    ------
    texpr : ElinaTexpr0Ptr
        Pointer to the ElinaTexpr0.
    constyp : c_uint
        Enum defining the constraint type as described in ElinaConstyp.
    scalar : ElinaScalarPtr
        Pointer to the the ElinaScalar.

    """

    _fields_ = [('texpr0', ElinaTexpr0Ptr), ('constyp', c_uint), ('scalar', ElinaScalarPtr)]

ElinaTcons0Ptr = POINTER(ElinaTcons0)


class ElinaTcons0Array(Structure):
    """
    ElinaTcons0Array ctype compatible with elina_tcons0_array from elina_tcons_array.h.

    Fields
    ------
    p : ElinaTcons0Ptr
        Expression.
    size : c_size_t
        Size of the array.
        
    """

    _fields_ = [('p', ElinaTcons0Ptr), ('size', c_size_t)]

ElinaTcons0ArrayPtr = POINTER(ElinaTcons0Array)
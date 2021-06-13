#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
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

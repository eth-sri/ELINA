#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
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
# elina_abstract0.h: generic operations on numerical abstract values
# ************************************************************************* #

from elina_manager_h import *
from elina_texpr0_h import *
from elina_tcons0_h import *


class ElinaAbstract0(Structure):
    """
    ElinaAbstract0 ctype compatible with elina_abstract0_t from elina_manager.h
    
    Fields
    ------
    value : c_void_p
    man : ElinaManagerPtr
    
    """

    _fields_ = [('value', c_void_p), ('man', ElinaManagerPtr)]

ElinaAbstract0Ptr = POINTER(ElinaAbstract0)
ElinaAbstract0Array = POINTER(ElinaAbstract0Ptr)
